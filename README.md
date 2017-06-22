# splinectomeR
### R package of spline-based statistical analysis tools for longitudinal data
***
These functions are designed to provide statistical analyses in _real_ longitudinal data, which may be missing timepoints or limited data at any given timepoint, have noisy biological variability, and variable numbers of samples per individual being measured. Comparisons can be made between two groups or within a single group for a non-zero change over the independent axis. Both return a p-value based on a randomly permuted distribution of the real data.

## Installation:
To install from GitHub Requires `devtools` package. Once you have that installed and loaded, run the following to install splinectomeR:
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
