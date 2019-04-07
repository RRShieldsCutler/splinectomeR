# splinectomeR
#### R package of spline-based statistical analysis tools for longitudinal data
[Explore on the splinectomeR website](https://rrshieldscutler.github.io/splinectomeR/)  
 
__To cite:__  
Shields-Cutler RR, Al-Ghalith GA, Yassour M, Knights D. (2018) SplinectomeR Enables Group Comparisons in Longitudinal Microbiome Studies. _Frontiers in Microbiology_ 9:785. doi: 10.3389/fmicb.2018.00785  
[![DOI](https://zenodo.org/badge/94937505.svg)](https://zenodo.org/badge/latestdoi/94937505)  
***
These functions are designed to provide statistical analyses in _real_ longitudinal data, which may be missing timepoints, possess limited data at some timepoints, have noisy biological variability, and variable numbers of observations per individual being measured. Comparisons can be made between two groups or within a single group for a non-zero change over the independent axis; both return a p-value based on a randomly permuted distribution of the real data. There is also a function for measuring significant differences at intervals across the entire x series (e.g. time) by interpolating splines from the original data.

### R package installation:
Recommended: Install v.0.9.3 from the releases, [here](https://github.com/RRShieldsCutler/splinectomeR/releases/tag/v0.9.3).

To install the latest development version from GitHub, you can use the `devtools` package. Beware, this version is actively worked on and may not be stable; therefore we highly recommend installing the pre-release. To install with devtools, run the following command:
```R
> devtools::install_github('RRShieldsCutler/splinectomeR')
> library(splinectomeR)
```
  
There is also a command line executable version (still requires R package installation); instructions located [below](#the-command-line-version). Note that the R package version is a little more robustly tested...
  
#### Permuspliner
This function tests for a greater-than-chance difference between two groups of interest over the x variable (e.g. time). Input data is a dataframe with the following columns at minimum:  

* a column defining the individuals (`cases`, e.g. patient_id, user_name, etc)
* a column (`category =`) for the groups (two or more; if >2 select which with `groups = `)
* a column each for the independent (`x`) and dependent (`y`) variables (numeric, continuous)
  
An example using the `ChickWeight` dataset from the `datasets` package:
```R
> # Test for difference in weight change over time between Diet groups 1 and 2
> result <- permuspliner(data = ChickWeight, x = 'Time', y = 'weight', perms=99,
                         cases = 'Chick', category = 'Diet', groups = c('1','3'))
> result$pval
[1] 0.02
> # Test for difference in weight change over time between Diet groups 2 and 3
> result <- permuspliner(data = ChickWeight, x = 'Time', y = 'weight', perms=99
                         cases = 'Chick', category = 'Diet', groups = c('2','3'))
> result$pval
[1] 0.19
```
#### Trendyspliner
This function tests for a non-zero trend in the response over the x variable. Input data is a dataframe with the following columns at minimum:

* a column defining the individuals (`cases`, e.g. patient_id, user_name, etc)
* a column each for the independent (`x`) and dependent (`y`) variables (numeric, continuous)

If the dataframe contains multiple groups/populations, and a trend is sought for just one group, you can automatically subset to that group by defining the `category` (column name) and the `group` of interest, as in the example below.
```R
> # Test for non-zero trend in Chick weight in Diet group 1 over Time
> result <- trendyspliner(data = ChickWeight, x = 'Time', y = 'weight',
                          cases = 'Chick', category = 'Diet', group = '1', perms = 999)
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

If the dataframe contains more than two groups/populations, you can automatically subset to those groups by defining the `groups` of interest as in the example below.
```R
> # Test for significant differences in Chicks between Diets 1 and 2 at 100 Time intervals
> result <- sliding_spliner(data = ChickWeight, xvar = 'Time',yvar = 'weight',
                            category = 'Diet', groups = c('1','2'), cases = 'Chick', ints = 100)
> str(result)
List of 3  # Result edited for clarity
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
sliding_spliner.plot.splines(result$spline_longform, category = 'Diet',
                             xvar = 'Time', yvar = 'weight')
```
  
***
  
## The command line version
This still requires splinectomeR installed in R and the same dependencies, but can be run from the command line. It will generate and save formatted plots if you do desire. These executable scripts are wrappers for the R package functions, so the tests being run are identical, just a different user interface.  

### Installation
You can clone the GitHub repository, then add the bin directory to your path for easiest execution.
```shell
# Change the path below to the absolute path to the bin directory
export PATH="$PATH:/path/to/bin"
# Restart the terminal
```
Or, you can open [the bin directory on GitHub](bin/), then right-click and save the linked file. You can now execute this file, so long as you have the R package installed. You may need to make the file executable:
```shell
chmod +x script_name.R
```
  
### Dependencies
These scripts have been developed and tested in R version 3.3.1 on Mac OSX, currently tested on 3.5.0 as well. The following R packages are required:

* The splinectomeR package
* Required for the permusplineR package:
  + devtools
  + ggplot2
  + reshape2
* For the command line versions:
  + optparse

***
  
### How to use
### permusplines.R
Given a longitudinal dataset, permusplinectomy reports a p-value for a binary categorical variable (e.g. Healthy vs Disease). It is robust to differing number and frequency of sampling between individuals ("units" in the code), and allows the user to alter key parameters controlling the stringency of the pre-test filtering. The input dataset should be a tab-delimited file in long format and have a column for the patient/individual ID, the independent and response variables (e.g. time and continuous measured variable), and the categorical variable. Before running, ensure there are no missing values in the measured variable or the categorical variable columns.  

The p-value is calculated by permutation of the categorical label across the individuals, and non-parametrically measuring the likelihood of the true area between the categories being a result of random chance.
```shell
# Sample run command
permusplines.R -i data_table.txt -x Years -y Blood_glucose -c Disease_status -s PATIENT_ID --perms 999

# To see usage and all commandline options
permusplines.R --help
```
  
### trendysplines.R
This tests whether the data you've selected to analyze follow a non-zero trend over the longitudinal axis. It does this by measuring the area between the group baseline and the group's spline over a set number of intervals. It then permutes (default 999 permutations) to find out whether the consistent change away from the baseline can be explained due to random chance, and reports the p-value. This can be a positive or negative trend; use the `--plot_path=my_plot.png` argument to visually assess the trend. Because in many cases, the individuals' starting points will vary greatly, there is a `--mean_center=T` option, which simply adjusts all the individuals' data to the same mean value in the y-axis before analysis; this can be useful for variable datasets, such as those from a human population.

```shell
# Sample run command with plot
trendysplines.R -i data_table.txt -x Years -y HbA1c -c Disease_status -s PATIENT_ID --group=diabetes --perms 9999 --mean_center=FALSE --plot_path=my_plot.pdf

# To see usage and all commandline options
trendyplines.R --help
```

### sliding_splines.R
Given the same dataset, the sliding_spline_test treats each individual separately, effectively converting their time series into a dense time series through extrapolation of a spline. Each interval (default = 100) is then tested for non-parametric significance between the two groups, provided there are enough data points (default = 3+ per group). The script produces three output files: a plot (png) of the splines for each individual, a plot of the negative log of p-values over the independent (x) variable where the size of the points is scaled by the number of data points contributing to that test (thicker line = greater n at that x), and a table of p-value at each interval. The user may provide a file prefix for each of these files, which are saved in the current working directory.  
```shell
# Sample run command
sliding_splines.R -i data_table.txt -x Years -y Blood_glucose -c Disease_status -s PATIENT_ID --spline_intervals=200 --prefix=Foo_

# To see usage and all command line options
sliding_splines.R --help
```
A good practice set is the `ChickWeight` dataset in R's datasets package. To export this data for use in the `splinectomy` scripts run the following from the R console:  
```R
> write.table(ChickWeight, file = 'ChickWeight.txt', sep='\t', quote =F, row.names = F)
```
Try it out! A few examples to get you started:
```shell
# Compare diets 1 and 2 overall (99 permutations for faster run):

permusplines.R -i ChickWeight.txt -s Chick -c Diet -x Time -y weight --perms=99 --plot_path=chicks_1v2.png --groups=1,2

# Should return a non-significant pvalue. But the separation isn't consistent across the time series... Try a sliding spline:

sliding_splines.R -i ChickWeight.txt -s Chick -c Diet -x Time -y weight --groups=1,2 --prefix=chicks_1vs2_

# Look at the pvalues plot...
# Cool! So, this suggests that the diets may lead to significantly different 
# chick weights at early timepoints, but this difference is not maintained 
# through the end of the experiment.
```


