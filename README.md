# Extreme Partial Least-Squares (EPLS)


## Abstract

We propose a new approach, called Extreme-PLS, for dimension reduction in conditional extreme values settings. The objective is to find linear combinations of predictors that best explain the extreme values of the response variable in a non-linear inverse regression model. The asymptotic normality of the Extreme-PLS estimator is established in the single-index framework and under mild assumptions. The performance of the method is assessed on simulated data. A statistical analysis of French farm income data, considering extreme cereal yields, is provided as an illustration.

## Usage

- The following packages should be installed before running EPLS: 
  ```dplyr, data.table, stats, base, graphics, grDevices, copula, VineCopula.```
- The code should be downloaded in a directory named "Code" and the data in a directory named "Data".
- The first part of the script (lines 20-304) is used to reproduce the experiments on simulated data (Section 5),
  the default parametrization allows to draw the center-left panel of Figure 2 (Frank copula, theta=10, dimension p=3).
- The second part of the script (lines 318-995) is used to reproduce the experiments on real data (Section 6).

**Reference:** M. Bousebata, G. Enjolras & S. Girard, "Extreme Partial Least-Squares", https://hal.inria.fr/hal-03165399
