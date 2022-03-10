# EPLS = Extreme Partial Least-Squares

- The following packages should be installed before running EPLS: 
  devtools, dplyr, MASS, data.table, xtable, copula, VineCopula, fExtremes, RMThreshold, EnvStats, graphics, rgl, maxtrixcalc, ginv, VGAM, ggplot2, geometry, compositions.
- The code should be downloaded in a directory named "Code" and the data in a directory named "Data".
- The first part of the script (lines 16-298) is used to reproduce the experiments on simulated data (Section 5),
  the default parametrization allows to draw the center-left panel of Figure 2 (Frank copula, theta=10, dimension p=3).
- The second part of the script (lines 312-1087) is used to reproduce the experiments on real data (Section 6).

Reference: M. Bousebata, G. Enjolras & S. Girard, "Extreme Partial Least-Squares", https://hal.inria.fr/hal-03165399
