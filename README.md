# splinectomeR
### R package of spline-based statistical analysis tools for longitudinal data
***
These functions are designed to provide statistical analyses in _real_ longitudinal data, which may be missing timepoints or limited data at any given timepoint, have noisy biological variability, and variable numbers of samples per individual being measured. Comparisons can be made between two groups or within a single group for a non-zero change over the independent axis. Both return a p-value based on a randomly permuted distribution of the real data.

## Installation:
To install from GitHub, use the `devtools` package. Once you have that installed and loaded, run the following to install splinectomeR:
```R
> devtools::install_github('RRShieldsCutler/splinectomeR')
> library(splinectomeR)
```
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
This function uses the `result$pval_table` to plot the trend in p-values over the x variable, with the points scaled by the number of observations (so points with greater _n_ in the significance test appear larger). This allows one to quickly visualize the trends.
```R
sliding_spliner.plot.pvals(result$pval_table, cases = 'Chick')
```
##### sliding_spliner.plot.splines
To see what the spline smoothed data look like in comparing the two groups, this function uses the `spline_longform` element of the results. You can specify the category, xvar, and yvar to display on the plots (defaults are 'category', 'xvar', and 'yvar').
```R
sliding_spliner.plot.splines(result$spline_longform, category = 'Diet', xvar = 'Time', yvar = 'weight')
```
