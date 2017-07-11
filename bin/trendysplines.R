#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(cowplot))

usage = '\nPermutation to test whether there is a non-zero trend among a set
of individuals/samples over a continuous variable (such as time). So, there
does not need to be two groups in this test. The x variable datapoints are
permuated within each case/individual, thus maintaining the distribution in
the y component but shuffling the hypothesized trend.'

option_list = list(
  make_option(c('-i', '--input'),
              help='REQUIRED: Input is long-form tab-delimited table, prefiltered for missing values',
              default=NA, type = 'character'),
  make_option(c('-c', '--category'),
              help='REQUIRED: Categorical variable of interest; must be a column header in table',
              default=NA, type = 'character'),
  make_option(c('--group'),
              help='If >1 group in category, the group of interest (e.g. group=disease)',
              default=NA, type = 'character'),
  make_option(c('-x', '--x_variable'),
              help='REQUIRED: The independent variable (e.g. time); must be column header',
              default=NA, type = 'character'),
  make_option(c('-y', '--y_variable'),
              help='REQUIRED: The response variable; must be a column header',
              default=NA, type = 'character'),
  make_option(c('-p', '--unit_id'),
              help='REQUIRED: The column header for your case grouping (e.g. Patient_ID, User_Name, etc)',
              default=NA, type = 'character'),
  make_option(c('--perms'),
              help='Number of permutation shuffles [default %default]',
              default=999, type = 'integer'),
  make_option(c('--cut'),
              help='Cut Unit IDs that occur less than this many times',
              default=NA, type = 'character'),
  make_option(c('--mean_center'),
              help='Mean center the data by individual before processing [default %default]',
              default=FALSE),
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
category = opt$category  # the column header label for the group
test_grp = opt$group  # the group of interest, if column has >1
x.cat = opt$x_variable  # the time series label
y.cat = opt$y_variable  # the response variable label
unit.id = opt$unit_id  # the header label defining the individuals (patients, etc)
num.perm = as.numeric(opt$perms)  # default 999
cut.low = opt$cut
spar.param = opt$spar # default NULL
samp.intervals = opt$intervals 
plot.results = opt$plot  # name of the plot file.png
mean_adj = opt$mean_center
shuff.id = 'y_shuff'

# TODO: Wrap the package function here




