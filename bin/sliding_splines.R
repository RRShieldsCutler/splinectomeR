#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(cowplot))

# TODO:
# Optparsing
# Remove cowplot dependency
# Remove dplyr dependency

usage = '\nSliding spline test to compare a binary categorical variable across
longitudinal data that may be sparse or messy. Generates splines for each
individual and uses those points to compare groups across time. Result is
a set of non-parametric p-values at user-defined density across the time
series. Returns spline plot, p-value plot (png), and p-value table (txt) to
the current working directory.'

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
  make_option(c('--spar'),
              help='The spar parameter when fitting splines (0 - 1); default is calculated from data',
              default=NULL),
  make_option(c('--spline_intervals'),
              help='Number of intervals extrapolated across the splines [default %default]',
              default=100, type = 'integer'),
  make_option(c('--density'),
              help='Data points required for each sliding spline position to test [default %default]',
              default=3, type = 'integer'),
  make_option(c('--n_per_unit'),
              help='Data points required per unit (must be > 3) to keep in dataset [default %default]',
              default=4, type = 'integer'),
  make_option(c('--prefix'),
              help='Prefix for results files',
              default='', type = 'character')
)
opt = parse_args(OptionParser(usage=usage, option_list=option_list))

# TODO: Make sure this is complete:
if (is.na(opt$input) | is.na(opt$category) | is.na(opt$x_variable) |
    is.na(opt$y_variable) | is.na(opt$unit_id)) {
  stop('Missing required parameters. See usage and options (--help)')
}

# Parse command line
infile = opt$input  # tsv file with data in long form
category = as.character(opt$category)  # the column header label for the two groups
groups = opt$groups  # the two groups of interest, if column has >2
x.cat = opt$x_variable  # the time series label
y.cat = opt$y_variable  # the response variable label
unit.id = opt$unit_id  # the header label defining the individuals (patients, etc)
num.bits = as.numeric(opt$spline_intervals)  # default 100
sparsity = as.numeric(opt$density)  # cutoff for number of points at given time; default 3 (>= 3)
unit_number = as.numeric(opt$n_per_unit)  # min number of points needed per individual
spar.param = opt$spar  # default NULL
prefix = opt$prefix


# TODO: wrap the R function using optparse arguments



