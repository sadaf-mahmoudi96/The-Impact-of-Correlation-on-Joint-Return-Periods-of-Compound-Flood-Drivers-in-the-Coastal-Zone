# The-Impact-of-Correlation-on-Joint-Return-Periods-of-Compound-Flood-Drivers-in-the-Coastal-Zone

## Non-Stationary NTR–RF Copula Analysis

This repository contains reusable R and Python functions for finding the peak-over-threshold and analyzing the time-varying dependence between non-tidal residual water level (NTR) and rainfall (RF).

The workflow includes:

* time-varying NTR threshold estimation;
* non-stationary GPD fitting with

  $$
  \sigma(t)=a+bt
  $$
* RF marginal distribution fitting;
* transformation of NTR and RF to uniform marginals;
* non-stationary bivariate copula fitting;
* Kendall’s tau estimation; and
* calculation of joint return periods.

### Files

```text
R/functions.R
Marginals_Copulas.R
Install_NSVineCopula.md
POT_Function.ipynb
```


### Main parameters

Important settings include:

```r
TARGET_RETURN_PERIOD <- 100
EVENTS_PER_YEAR <- 5
WINDOW_HALF_WIDTH <- 20
CONTOUR_POINTS <- 100
```

The joint return period is calculated as:

$$
T =
\frac{1}
{\lambda\left[1-u_{NTR}-u_{RF}+C(u_{NTR},u_{RF})\right]}
$$

where \(\lambda\) is the average number of POT events per year.
