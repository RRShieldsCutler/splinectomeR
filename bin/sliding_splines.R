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
series. Returns p-value plot and p-value table (txt) to
the current working directory or to paths provided in the arguments.'

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
  make_option(c('-s', '--cases'),
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
  make_option(c('--cut_low'),
              help='Data points required per unit (must be > 3) to keep in dataset [default %default]',
              default=4, type = 'integer'),
  make_option(c('--result_fp'),
              help='File path for result table',
              default='sliding_splines.txt', type = 'character'),
  make_option(c('--plot_fp'),
              help='Plot the data too. Provide a filename/path,
              with extension explicitly in the name: pdf, png,
              tiff, jpg, bmp, eps, or ps.',
              default=NA, type = 'character')
)
opt = parse_args(OptionParser(usage=usage, option_list=option_list))

# Make sure this is complete:
if (is.na(opt$input) | is.na(opt$category) | is.na(opt$x_variable) |
    is.na(opt$y_variable) | is.na(opt$cases)) {
  stop('Missing required parameters. See usage and options (--help)')
}

require(splinectomeR)

# Parse command line
infile = opt$input  # tsv file with data in long form
category = as.character(opt$category)  # the column header label for the two groups
groups = opt$groups  # the two groups of interest, if column has >2
xvar = opt$x_variable  # the time series label
yvar = opt$y_variable  # the response variable label
cases = opt$cases  # the header label defining the individuals (patients, etc)
ints = as.numeric(opt$spline_intervals)  # default 100
sparsity = as.numeric(opt$density)  # cutoff for number of points at given time; default 3 (>= 3)
cut_low = as.numeric(opt$cut_low)  # min number of points needed per individual
spar_param = opt$spar  # default NULL
resulttab = as.character(opt$result_fp)
plot_path <- as.character(opt$plot_fp)  # name of the plot file.png

if (!is.na(groups)) {
  g1 <- strsplit(groups, ',')[[1]][1]
  g2 <- strsplit(groups, ',')[[1]][2]
  groups_in <- c(g1,g2)
}

# param data A dataframe object containing your data.
# param x The independent variable; is continuous, e.g. time.
# param y The dependent variable; is continuous, e.g. temperature.
# param category The column name of the category to be tested, if present.
# param cases The column name defining the individual cases, e.g. patients.
# param groups If more than two groups, the two groups to compare.
# param set_spar Set the spar parameter for splines
# param cut_sparse Remove sparse cases with fewer than __ data points (default 4)
# param test_density Minimum density of cases in each group to report p-value (default 3)
# param ints Number of x intervals over which to measure significance
# param quiet Silence all text outputs

# TODO: wrap the R function using optparse arguments

df <- read.delim(file = infile, header = 1, check.names = F, sep = '\t')

result <- sliding_spliner(data = df, xvar = xvar, yvar = yvar, category = category,
                          cases = cases, groups = groups_in, ints = ints,
                          set_spar = spar_param, cut_low = cut_low,
                          test_density = sparsity, quiet = F)

pval.p <- result['pval_table'][[1]]

write.table(x = pval.p, file = resulttab, quote = F,
            row.names = F, col.names = T, sep = '\t')
cat('Table written to the following file:\n')
cat(resulttab)

# Plot also if desired
if (!is.na(plot_fp)) {
  require(ggplot2)
  .norm_range <- function(x) {(x-min(x)) / (max(x)-min(x))}
  xvar <- names(pval.p)[1]
  pval.p$N.norm <- .norm_range(x = pval.p[,3])
  if (length(unique(pval.p[,3])) <= 1) {
    cat('Number of observations per interval is uniform; points will not be plotted.')
  }
  p <- ggplot(pval.p, aes(x=pval.p[,1], y=p_value)) + geom_line() +
    geom_point(shape = 20, size = (pval.p$N.norm * 2))  +
    geom_hline(aes(yintercept = 0.05), linetype='dashed') +
    xlab(xvar) + ylab('p-value (log)') +
    scale_y_continuous(trans = scales::log10_trans()) +
    theme_bw() + theme(legend.position='none',
                       plot.background = element_rect(color = 'white'),
                       panel.grid = element_blank(),
                       axis.text = element_text(color = 'black'))
  ggsave(p, filename = plot_fp, height = 3.5, width = 4, units = 'in', dpi = 600)
  cat('Plot saved\n')
}
  
  
  
  