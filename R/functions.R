
check_required_packages <- function() {
  required <- c("VineCopula", "NSVineCopula", "evd", "fitdistrplus", "lubridate")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) stop("Missing required packages: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

clip01 <- function(u, eps = 1e-12) pmin(pmax(u, eps), 1 - eps)

extract_station_id <- function(file_path) {
  stem <- sub("\\.csv$", "", basename(file_path), ignore.case = TRUE)
  parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
  if (length(parts) < 3) stop("Could not parse station ID from: ", basename(file_path))
  parts[length(parts) - 2]
}

prepare_station_data <- function(data, ntr_col="POT_NTR", rf_col="Max_RF_in_5d", time_col="Time") {
  req <- c(ntr_col, rf_col, time_col)
  miss <- setdiff(req, names(data)); if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))
  data[[time_col]] <- as.Date(data[[time_col]])
  keep <- !is.na(data[[time_col]]) & is.finite(data[[ntr_col]]) & is.finite(data[[rf_col]]) & data[[rf_col]] != 0
  data <- data[keep, , drop=FALSE]
  data[order(data[[time_col]]), , drop=FALSE]
}

threshold_coverage_ok <- function(threshold_data, year_col="Year", threshold_col="Threshold", after_year=1950, min_fraction=0.80) {
  sub <- threshold_data[threshold_data[[year_col]] > after_year, , drop=FALSE]
  if (!nrow(sub)) return(FALSE)
  mean(!is.na(sub[[threshold_col]])) >= min_fraction
}

fit_threshold_model <- function(threshold_window, year_col="Year", threshold_col="Threshold") {
  d <- threshold_window[complete.cases(threshold_window[, c(year_col, threshold_col)]), , drop=FALSE]
  if (nrow(d) < 2) stop("Too few threshold observations in this window.")
  dd <- data.frame(Threshold=d[[threshold_col]], elapsed_years=d[[year_col]] - min(d[[year_col]]))
  lm(Threshold ~ elapsed_years, data=dd)
}

fit_time_varying_gpd <- function(x_obs, mu_obs, t_obs,
                                 guesses=list(c(0.2,0.01,0.1), c(0.2,0.001,0.1))) {
  y <- x_obs - mu_obs
  nll <- function(par) {
    a <- par[1]; b <- par[2]; xi <- par[3]
    sigma <- a + b*t_obs
    if (any(!is.finite(sigma)) || any(sigma <= 0)) return(1e12)
    ll <- tryCatch(evd::dgpd(y, loc=0, scale=sigma, shape=xi, log=TRUE), error=function(e) rep(NA_real_, length(y)))
    if (any(!is.finite(ll))) return(1e12)
    -sum(ll)
  }
  fits <- lapply(guesses, function(g) tryCatch(optim(g, nll, method="Nelder-Mead"), error=function(e) NULL))
  fits <- Filter(function(z) !is.null(z) && is.finite(z$value) && all(is.finite(z$par)), fits)
  if (!length(fits)) stop("Time-varying GPD fit failed.")
  fit <- fits[[which.min(vapply(fits, function(z) z$value, numeric(1)))]]
  a <- fit$par[1]; b <- fit$par[2]; xi <- fit$par[3]; sigma <- a + b*t_obs
  if (any(sigma <= 0)) stop("Non-positive fitted GPD scale.")
  list(fit=fit, a=a, b=b, xi=xi, sigma=sigma)
}

gpd_cdf_tv <- function(x, mu, sigma, xi) {
  if (abs(xi) < 1e-8) u <- 1 - exp(-(x-mu)/sigma) else {
    z <- 1 + xi*(x-mu)/sigma
    u <- ifelse(z <= 0, NA_real_, 1 - z^(-1/xi))
  }
  clip01(u)
}

gpd_quantile_tv <- function(p, mu, sigma, xi) {
  p <- clip01(p); alpha <- 1-p
  y <- if (abs(xi) < 1e-8) sigma*log(1/alpha) else (sigma/xi)*(alpha^(-xi)-1)
  mu + y
}

fit_best_rf_marginal <- function(x, distributions=c("gamma","lnorm","weibull","norm","exp")) {
  fits <- list()
  for (d in distributions) {
    f <- tryCatch(fitdistrplus::fitdist(x, d), error=function(e) NULL)
    if (!is.null(f)) fits[[d]] <- f
  }
  if (!length(fits)) stop("All RF marginal fits failed.")
  gof <- fitdistrplus::gofstat(fits)
  ia <- which.min(gof$aic); ib <- which.min(gof$bic)
  list(AIC_best=fits[[ia]], BIC_best=fits[[ib]], AIC_name=names(fits)[ia], BIC_name=names(fits)[ib], gof=gof)
}

rf_cdf <- function(x, fit) {
  d <- fit$AIC_best$distname; p <- fit$AIC_best$estimate
  u <- switch(d,
    exp=pexp(x, rate=as.numeric(p["rate"])),
    gamma=pgamma(x, shape=as.numeric(p["shape"]), rate=as.numeric(p["rate"])),
    lnorm=plnorm(x, meanlog=as.numeric(p["meanlog"]), sdlog=as.numeric(p["sdlog"])),
    weibull=pweibull(x, shape=as.numeric(p["shape"]), scale=as.numeric(p["scale"])),
    norm=pnorm(x, mean=as.numeric(p["mean"]), sd=as.numeric(p["sd"])),
    stop("Unsupported RF distribution: ", d))
  clip01(u)
}

rf_quantile <- function(u, fit) {
  u <- clip01(u); d <- fit$AIC_best$distname; p <- fit$AIC_best$estimate
  switch(d,
    exp=qexp(u, rate=as.numeric(p["rate"])),
    gamma=qgamma(u, shape=as.numeric(p["shape"]), rate=as.numeric(p["rate"])),
    lnorm=qlnorm(u, meanlog=as.numeric(p["meanlog"]), sdlog=as.numeric(p["sdlog"])),
    weibull=qweibull(u, shape=as.numeric(p["shape"]), scale=as.numeric(p["scale"])),
    norm=qnorm(u, mean=as.numeric(p["mean"]), sd=as.numeric(p["sd"])),
    stop("Unsupported RF distribution: ", d))
}

prepare_copula_matrix <- function(df) {
  u <- df[, c("U_NTR","U_RF"), drop=FALSE]
  u$U_NTR <- as.numeric(u$U_NTR); u$U_RF <- as.numeric(u$U_RF)
  u <- u[complete.cases(u) & u$U_NTR>0 & u$U_NTR<1 & u$U_RF>0 & u$U_RF<1, , drop=FALSE]
  u$U_NTR <- clip01(u$U_NTR); u$U_RF <- clip01(u$U_RF)
  m <- as.matrix(u); storage.mode(m) <- "double"; m
}

fit_nonstationary_copula <- function(u_data, primary_familyset=c(1,3), secondary_familyset=c(2,3,4,5)) {
  model <- tryCatch(NSVineCopula::NSBiCopEst(u_data, familyset=primary_familyset), error=function(e) NULL)
  if (is.null(model) || is.null(model$family) || model$family == 0) {
    model <- tryCatch(NSVineCopula::NSBiCopEst(u_data, familyset=secondary_familyset, indeptest=FALSE), error=function(e) NULL)
  }
  if (is.null(model) || is.null(model$family) || model$family == 0) stop("Non-stationary copula fit failed.")
  parfit <- NSVineCopula::NSBiCopPar(data=u_data, theta=model$beta, rhobar=model$TVTP[1,], family=model$family, Call="filtering")
  list(model=model, parameters=parfit)
}

copula_family_name <- function(fam) {
  map <- c(`1`="Gaussian", `2`="Student-t", `3`="Clayton", `4`="Gumbel", `5`="Frank",
           `13`="Survival Clayton", `14`="Survival Gumbel", `23`="Rotated Clayton 90",
           `24`="Rotated Gumbel 90", `33`="Rotated Clayton 270", `34`="Rotated Gumbel 270")
  out <- unname(map[as.character(as.integer(fam[1]))]); ifelse(is.na(out), paste0("Family_", fam[1]), out)
}

bound_copula_parameters <- function(par1, par2=0, fam) {
  fam <- as.integer(fam[1]); if (!is.finite(par1)) par1 <- 0; if (!is.finite(par2)) par2 <- 0
  if (fam == 1) par1 <- pmin(pmax(par1,-0.999),0.999)
  if (fam == 2) { par1 <- pmin(pmax(par1,-0.999),0.999); par2 <- pmax(par2,2.1) }
  if (fam %in% c(3,13)) par1 <- pmax(par1,1e-6)
  if (fam %in% c(23,33)) par1 <- pmin(pmax(-abs(par1),-28),-1e-6)
  if (fam %in% c(4,14)) par1 <- pmin(pmax(par1,1.000001),17)
  if (fam %in% c(24,34)) par1 <- pmin(pmax(-abs(par1),-17),-1.000001)
  if (fam == 5 && abs(par1)<1e-6) par1 <- ifelse(par1<0,-1e-6,1e-6)
  c(par1,par2)
}

generate_independence_contour <- function(target_return_period, events_per_year, n_points=100,
                                          rf_survival_min=0.005, rf_survival_max=0.5) {
  rhs <- 1/(events_per_year*target_return_period)
  s_rf <- seq(rf_survival_min, rf_survival_max, length.out=n_points)
  s_ntr <- rhs/s_rf
  ok <- is.finite(s_ntr) & s_ntr>0 & s_ntr<1
  data.frame(U_NTR=1-s_ntr[ok], U_RF=1-s_rf[ok])
}

nearest_parameter_index <- function(years, target, max_difference=2) {
  if (!length(years)) return(integer(0)); d <- abs(years-target); i <- which.min(d)
  if (!length(i) || !is.finite(d[i]) || d[i] > max_difference) integer(0) else i
}


evaluate_joint_return_periods <- function(contour, pars_df, fam, events_per_year, tau_test, center_year, ns_model,
                                          max_parameter_year_difference=2) {
  idx <- nearest_parameter_index(pars_df$Year, center_year, max_parameter_year_difference)
  if (!length(idx)) return(NULL)
  par1 <- pars_df$par1[idx]
  par2 <- if ("par2" %in% names(pars_df)) pars_df$par2[idx] else 0
  bp <- bound_copula_parameters(par1, par2, fam)
  Cval <- VineCopula::BiCopCDF(u1=contour$U_NTR, u2=contour$U_RF, family=fam, par=bp[1], par2=bp[2])
  pexc <- 1 - contour$U_NTR - contour$U_RF + Cval
  rp <- ifelse(is.finite(pexc) & pexc>0, 1/(events_per_year*pexc), NA_real_)
  data.frame(
    year=center_year, ntr=as.numeric(contour$NTR_Original), rf=as.numeric(contour$RF_Original),
    u_ntr_point=as.numeric(contour$U_NTR), u_rf_point=as.numeric(contour$U_RF), C_point=as.numeric(Cval),
    RP_Joint=as.numeric(rp), Kendall_Tau=as.numeric(tau_test$estimate), P_Value=as.numeric(tau_test$p.value),
    Copula=copula_family_name(fam), Copula_Family=as.integer(fam),
    Copula_AIC=as.numeric(ns_model$AIC), Copula_BIC=as.numeric(ns_model$BIC),
    Copula_LogLik=as.numeric(ns_model$Loglikelihood), stringsAsFactors=FALSE)
}

analyze_window <- function(data, threshold_data, center_year, window_half_width=20L,
                           target_return_period=100, events_per_year=5, contour_points=100,
                           ntr_col="POT_NTR", rf_col="Max_RF_in_5d", time_col="Time",
                           threshold_year_col="Year", threshold_col="Threshold",
                           primary_familyset=c(1,3), secondary_familyset=c(2,3,4,5)) {
  start_year <- center_year - window_half_width; end_year <- center_year + window_half_width
  obs_year <- as.integer(format(data[[time_col]], "%Y"))
  df <- data[obs_year > start_year & obs_year < end_year, , drop=FALSE]
  keep_t <- threshold_data[[threshold_year_col]] > start_year & threshold_data[[threshold_year_col]] < end_year
  keep_t[is.na(keep_t)] <- FALSE
  th <- threshold_data[keep_t, , drop=FALSE]
  if (nrow(df)<5 || nrow(th)<2) return(NULL)

  time_years <- as.numeric(difftime(df[[time_col]], min(df[[time_col]]), units="days"))/365.25
  lm_loc <- fit_threshold_model(th, threshold_year_col, threshold_col)
  mu_t <- predict(lm_loc, newdata=data.frame(elapsed_years=time_years))
  mask <- which(df[[ntr_col]] > mu_t)
  if (length(mask)<5) return(NULL)

  x_obs <- df[[ntr_col]][mask]; mu_obs <- mu_t[mask]; t_obs <- time_years[mask]
  gpd <- tryCatch(fit_time_varying_gpd(x_obs, mu_obs, t_obs), error=function(e) NULL)
  if (is.null(gpd)) return(NULL)
  u_ntr <- gpd_cdf_tv(x_obs, mu_obs, gpd$sigma, gpd$xi)

  rf_vals <- df[[rf_col]][mask]
  rf_fit <- tryCatch(fit_best_rf_marginal(rf_vals), error=function(e) NULL)
  if (is.null(rf_fit)) return(NULL)
  u_rf <- rf_cdf(rf_vals, rf_fit)

  copula_df <- data.frame(
    Datetime=df[[time_col]][mask], Year=lubridate::year(df[[time_col]][mask]),
    NTR_Original=x_obs, RF_Original=rf_vals, U_NTR=u_ntr, U_RF=u_rf)
  copula_df <- copula_df[complete.cases(copula_df), , drop=FALSE]
  if (nrow(copula_df)<5) return(NULL)
  tau_test <- suppressWarnings(cor.test(copula_df$U_NTR, copula_df$U_RF, method="kendall"))

  t_star <- median(time_years)
  mu_star <- as.numeric(predict(lm_loc, newdata=data.frame(elapsed_years=t_star)))
  sigma_star <- gpd$a + gpd$b*t_star
  if (!is.finite(sigma_star) || sigma_star<=0) return(NULL)

  contour <- generate_independence_contour(target_return_period, events_per_year, contour_points)
  contour$RF_Original <- rf_quantile(contour$U_RF, rf_fit)
  contour$NTR_Original <- gpd_quantile_tv(contour$U_NTR, mu_star, sigma_star, gpd$xi)

  u_data <- prepare_copula_matrix(copula_df)
  cop <- tryCatch(fit_nonstationary_copula(u_data, primary_familyset, secondary_familyset), error=function(e) NULL)
  if (is.null(cop)) return(NULL)
  pars_df <- as.data.frame(cop$parameters$TVTP)
  if (!ncol(pars_df)) return(NULL)
  names(pars_df) <- c("par1","par2")[seq_len(ncol(pars_df))]
  if (nrow(pars_df) != nrow(copula_df)) return(NULL)
  pars_df$Year <- copula_df$Year

  out <- evaluate_joint_return_periods(contour, pars_df, cop$model$family, events_per_year,
                                       tau_test, center_year, cop$model)
  if (is.null(out)) return(NULL)
  out$WindowStart <- start_year; out$WindowEnd <- end_year
  out$RF_Marginal <- rf_fit$AIC_name; out$GPD_a <- gpd$a; out$GPD_b <- gpd$b; out$GPD_xi <- gpd$xi
  out
}

analyze_station_file <- function(event_file, threshold_dir=".", output_dir=".",
                                 threshold_file_template="NTR_TimeVarying_Threshold_5_%s.csv",
                                 target_return_period=100, events_per_year=5, window_half_width=20L,
                                 contour_points=100, ntr_col="POT_NTR", rf_col="Max_RF_in_5d", time_col="Time",
                                 threshold_year_col="Year", threshold_col="Threshold",
                                 threshold_coverage_after_year=1950, minimum_threshold_coverage=0.80,
                                 primary_familyset=c(1,3), secondary_familyset=c(2,3,4,5), overwrite=FALSE) {
  station <- extract_station_id(event_file)
  threshold_file <- file.path(threshold_dir, sprintf(threshold_file_template, station))
  if (!file.exists(threshold_file)) { warning("Threshold file missing for station ", station); return(NULL) }
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  outfile <- file.path(output_dir, sprintf("KendallTau_RP_NTR_RF_%s_%syr_%sSamples_Windows%syears.csv",
                                          station, target_return_period, contour_points, 2*window_half_width))
  if (file.exists(outfile) && !overwrite) { message("Skipping station ", station, ": output exists."); return(invisible(outfile)) }

  data <- prepare_station_data(read.csv(event_file), ntr_col, rf_col, time_col)
  threshold_data <- read.csv(threshold_file)
  if (!threshold_coverage_ok(threshold_data, threshold_year_col, threshold_col,
                             threshold_coverage_after_year, minimum_threshold_coverage)) {
    message("Skipping station ", station, ": insufficient threshold coverage."); return(NULL)
  }

  years <- sort(unique(as.integer(format(data[[time_col]], "%Y"))))
  valid <- years[(years-window_half_width)>=min(years) & (years+window_half_width)<=max(years)]
  ans <- list()
  for (yr in valid) {
    message("Station ", station, " | center year ", yr)
    z <- tryCatch(analyze_window(data, threshold_data, yr, window_half_width, target_return_period,
                                 events_per_year, contour_points, ntr_col, rf_col, time_col,
                                 threshold_year_col, threshold_col, primary_familyset, secondary_familyset),
                  error=function(e) { message("  skipped: ", conditionMessage(e)); NULL })
    if (!is.null(z)) ans[[length(ans)+1]] <- z
  }
  if (!length(ans)) { warning("No valid windows for station ", station); return(NULL) }
  final <- do.call(rbind, ans)
  write.csv(final, outfile, row.names=FALSE)
  message("Saved: ", outfile)
  invisible(outfile)
}

run_all_stations <- function(input_dir=".", threshold_dir=input_dir, output_dir=file.path(input_dir,"results"),
                             event_file_pattern="POT_NTR_TimeVarying_RF_NOAA_GHCN_5_*.csv", ...) {
  check_required_packages()
  files <- Sys.glob(file.path(input_dir, event_file_pattern))
  if (!length(files)) stop("No event files matched the input pattern.")
  invisible(lapply(files, function(f) analyze_station_file(f, threshold_dir=threshold_dir, output_dir=output_dir, ...)))
}
