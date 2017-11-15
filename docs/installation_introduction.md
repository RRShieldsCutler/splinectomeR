---
title: "Installation and Introduction"
description: "A quick introduction to the splinectomeR package"
output: html_document
---
***
  
### Installation:
Recommended: Install v.0.9.0 from the pre-release, [here](https://github.com/RRShieldsCutler/splinectomeR/releases/tag/v0.9.0).

To install the latest development version from GitHub, you can use the `devtools` package. Beware, this version is actively worked on and may not be stable; therefore we highly recommend installing the pre-release. To install with devtools, run the following command:
```R
> devtools::install_github('RRShieldsCutler/splinectomeR')
> library(splinectomeR)
```
***
  
### Basic descriptions of the functions  
  
#### Permuspliner
This function tests for a greater-than-chance difference between two groups of interest over the x variable (e.g. time). Input data is a dataframe with the following columns at minimum:

* a column defining the individuals (`cases`, e.g. patient_id, user_name, etc)
* a column (`category =`) for the groups (two or more; if >2 select which with `groups = `)
* a column each for the independent (`x`) and dependent (`y`) variables (numeric, continuous)

An example using the `ChickWeight` dataset from the `datasets` package:
```R
> # Test for difference in weight change over time between Diet groups 1 and 2
> result <- permuspliner(data = ChickWeight, x = 'Time', y = 'weight', cases = 'Chick', category = 'Diet', groups = '1,2')
> result$pval
[1] 0.003
> # Test for difference in weight change over time between Diet groups 2 and 3
> result <- permuspliner(data = ChickWeight, x = 'Time', y = 'weight', cases = 'Chick', category = 'Diet', groups = '2,3')
> result$pval
[1] 0.159
```
##### Permuspliner plotting functions
There are two plotting functions built into the permuspliner test, `permuspliner.plot.permdistance` and `permuspliner.plot.permsplines`. The former plots the true measured distance between the two groups across the longitudinal axis, as well as all of the distances from the permuted data. These are the curves from which the p-value is derived, so this provides visual confidence of the p-value; if the real distance is greater than expected from random chance, its curve will be at or near the top (or bottom) of most of the permuted distances.  
  
The second plotting function, `permuspliner.plot.permsplines` shows the two groups' actual splines along with the individual permuted splines. Here, you can see how well the two groups are distinguished compared to all the random distributions, again, visualizing the strength of the difference. To run these, you use the result object from the main function.
```R
# Distance plot has a label argument for aesthetics - it just labels the axis
permuspliner.plot.permdistance(result, xlabel = 'Time')

# Splines plot requires the names of the longitudinal (x) and response (y) variables
# as were used in the original test function.
permuspliner.plot.permsplines(result, xvar = 'Time', yvar = 'weight')
```  

***
  
#### Trendyspliner
This function tests for a non-zero trend in the response over the x variable. Input data is a dataframe with the following columns at minimum:

* a column defining the individuals (`cases`, e.g. patient_id, user_name, etc)
* a column each for the independent (`x`) and dependent (`y`) variables (numeric, continuous)

If the dataframe contains multiple groups/populations, and a trend is sought for just one group, you can automatically subset to that group by defining the `category` (column name) and the `group` of interest, as in the example below.
```R
> # Test for non-zero trend in Chick weight in Diet group 1 over Time
> result <- trendyspliner(data = ChickWeight, x = 'Time', y = 'weight', cases = 'Chick', category = 'Diet', group = '1', perms = 999)
> result$pval
[1] 0.001
```
For further details and an example, check out the help docs for each function.
```R
> library(splinectomeR)
> ?permuspliner()
> ?trendyspliner()
```
##### trendyspliner.plot.perms
You can plot the actual trend spline with all the permuted splines using this function and the results object from the main test function. This allows you to visualize the confidence of your permutation test. It uses `ggplot2` and `reshape2` to plot the results, and you can run it simply as follows:
```R
trendyspliner.plot.perms(result, xlabel = 'time', ylabel = 'weight of chick')
```

***
  
#### Sliding_spliner
The objective of this function is to reveal whether two groups differ significantly at any point along the longitudinal axis. You might have a non-significant overall difference (determined by `permuspliner`), but perhaps the groups diverge dramatically for a short time in the middle of the study. This window may be of interest for followup hypotheses, for example. This function smooths each individual's datapoints and fills in gaps by interpolation of the splines. Then, it tests for significant difference between two groups at each step along the longitudinal axis (as long as there was enough user data at that x value in the original data; defaults to 3 in each group, can be modified). Input data is a dataframe with the following columns at minimum:

* a column defining the individuals (`cases`, e.g. patient_id, user_name, etc)
* a column each for the independent (`x`) and dependent (`y`) variables (numeric, continuous)

If the dataframe contains more than two groups/populations, you can automatically subset to those groups by defining the `groups` of interest, as a string, as in the example below.
```R
> # Test for significant differences in Chicks between Diets 1 and 2 at 100 Time intervals
> result <- sliding_spliner(data = ChickWeight, xvar = 'Time', yvar = 'weight', category = 'Diet', groups = '1,2', cases = 'Chick', ints = 100)
> str(result)
List of 3  # Result edited for clarity!
 $ pval_table     :'data.frame':	100 obs. of  3 variables:
 $ spline_table   :'data.frame':	100 obs. of  30 variables:
 $ spline_longform:'data.frame':	2900 obs. of  4 variables:
```
The output is a list with three items:
1. the pvalues and total sample size at each interval
2. the splines for each individual (as x,y)
3. the splines as a long-form dataframe, ready for plotting (see next...)

##### sliding_spliner.plot.pvals
This function uses the results object of the main test function to plot the trend in p-values over the x variable, with the points scaled by the number of observations (so points with greater _n_ in the significance test appear larger). This allows one to quickly visualize the trends.
```R
sliding_spliner.plot.pvals(result, cases = 'Chick')
```
