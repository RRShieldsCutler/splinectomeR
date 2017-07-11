#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(cowplot))

usage = '\nPermutation test to determine whether two groups are significantly
different across longitudinal data that may have differing patterns
of data for each subject (time, n, etc). Applies a permutation test using
splines fit to the two groups. Categorical labels are shuffled among the
subjects to determine a random distribution for comparison.'

option_list = list(
  make_option(c('-i', '--input'),
              help='REQUIRED: Input is long-form tab-delimited table, prefiltered for missing values',
              default=NA, type = 'character'),
  make_option(c('-c', '--category'),
              help='REQUIRED: Categorical variable of interest; must be a column header in table',
              default=NA, type = 'character'),
  make_option(c('--groups'),
              help='If >2 groups in category, the two of interest: comma-delimited (e.g. groups=Healthy,Sick)',
              default=NA, type = 'character'),
  make_option(c('-x', '--x_variable'),
              help='REQUIRED: The independent variable (e.g. time); must be column header',
              default=NA, type = 'character'),
  make_option(c('-y', '--y_variable'),
              help='REQUIRED: The response variable; must be a column header',
              default=NA, type = 'character'),
  make_option(c('-p', '--unit_id'),
              help='REQUIRED: The column header for your grouping (e.g. Patient_ID, User_Name, etc)',
              default=NA, type = 'character'),
  make_option(c('--perms'),
              help='Number of permutation shuffles [default %default]',
              default=999, type = 'integer'),
  make_option(c('--cut'),
              help='Cut Unit IDs that occur less than this many times',
              default=NA, type = 'character'),
  make_option(c('--intervals'),
              help='Number of sampling intervals along spline [default %default]',
              default=10000, type = 'integer'),
  make_option(c('--spar'),
              help='The spar parameter when fitting splines (0 - 1); default is calculated from data',
              default=NULL),
  make_option(c('--plot'),
              help='Plot the data too! Provide a filename/path, with .png extension',
              default=NA, type = 'character')
)
opt = parse_args(OptionParser(usage=usage, option_list=option_list))

if (is.na(opt$input) | is.na(opt$category) | is.na(opt$x_variable) |
    is.na(opt$y_variable) | is.na(opt$unit_id)) {
  stop('Missing required parameters. See usage and options (--help)')
}

# Parse command line
infile = opt$input  # tsv file with data in long form
category = opt$category  # the column header label for the two groups
groups = opt$groups  # the two groups of interest, if column has >2
x.cat = opt$x_variable  # the time series label
y.cat = opt$y_variable  # the response variable label
unit.id = opt$unit_id  # the header label defining the individuals (patients, etc)
num.perm = as.numeric(opt$perms)  # default 999
cut.low = opt$cut
spar.param = opt$spar # default NULL
samp.intervals = opt$intervals 
plot.results = opt$plot  # name of the plot file.png
shuff.id = 'cat.shuff'


# TODO: Wrap the package function


