source(file.path("R", "functions.R"))

# ---------------- USER SETTINGS ----------------
INPUT_DIR <- "data"
THRESHOLD_DIR <- "data"
OUTPUT_DIR <- "results"

TARGET_RETURN_PERIOD <- 100
EVENTS_PER_YEAR <- 5
WINDOW_HALF_WIDTH <- 20L   # +/-20 years = nominal 40-year centered window
CONTOUR_POINTS <- 100

# ---------------- RUN ----------------
run_all_stations(
  input_dir = INPUT_DIR,
  threshold_dir = THRESHOLD_DIR,
  output_dir = OUTPUT_DIR,
  event_file_pattern = "POT_NTR_TimeVarying_RF_NOAA_GHCN_5_*.csv",
  threshold_file_template = "NTR_TimeVarying_Threshold_5_%s.csv",
  target_return_period = TARGET_RETURN_PERIOD,
  events_per_year = EVENTS_PER_YEAR,
  window_half_width = WINDOW_HALF_WIDTH,
  contour_points = CONTOUR_POINTS,
  ntr_col = "POT_NTR",
  rf_col = "Max_RF_in_5d",
  time_col = "Time",
  threshold_year_col = "Year",
  threshold_col = "Threshold",
  threshold_coverage_after_year = 1950,
  minimum_threshold_coverage = 0.80,
  primary_familyset = c(1, 3),
  secondary_familyset = c(2, 3, 4, 5),
  overwrite = FALSE
)
