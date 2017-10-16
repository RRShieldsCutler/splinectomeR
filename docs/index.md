# splinectomeR
### R package of spline-based statistical analysis tools for longitudinal data
***
These functions are designed to provide statistical analyses in _real_ longitudinal data, which may be missing timepoints, possess limited data at some timepoints, have noisy biological variability, and variable numbers of observations per individual being measured. Comparisons can be made between two groups or within a single group for a non-zero change over the independent axis; both return a p-value based on a randomly permuted distribution of the real data. There is also a function for measuring significant differences at intervals across the entire x series (e.g. time) by interpolating splines from the original data.  
  
For information on installation of the package from GitHub, and a basic outline of the functions, [go here](installation_introduction.md).  
  
Happy splining, and have a _splinedid_ day!
  
#### Vignettes
There are two vignettes available for this package:
1. [Chick Weights](chickweights_web.html) - a simple exploration using the R dataset `ChickWeight`.
1. [Yassour et al. Antibiotics](yassour_antibiotics_web.html) - a more intensive application of the `splinectomeR` package to a large set of longitudinal microbiome data published by Yassour et al.[^1] following 38 babies over three years.
  
***  
  
[^1]: [Yassour et al. Natural history of the infant gut microbiome and impact of antibiotic treatment on bacterial strain diversity and stability. 8, 343RA81 (2016)](http://stm.sciencemag.org/content/8/343/343ra81.short)  
  
  
  
#### Disclaimer
Licensed under GPL-3. As such, the `splinectomer` package is freely offered, without any warrantee or guarantee. It may be used, redistributed, or modified for non-commercial purposes with appropriate citation of the original work. You, the user, are solely responsible for verifying the legitimacy of your results.  
  
***
#### References
