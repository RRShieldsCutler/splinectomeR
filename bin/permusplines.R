#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(splinectomeR))

usage = '\nPermutation test to determine whether two groups are significantly
different across longitudinal data that may have differing patterns
of data for each subject (time, n, etc). Applies a permutation test using
splines fit to the two groups. Categorical labels are shuffled among the
subjects to determine a random distribution for comparison.'

option_list = list(
  make_option(c('-i', '--input'),
              help='REQUIRED: Input is long-form tab-delimited table,
                prefiltered for missing values',
              default=NA, type = 'character'),
  make_option(c('-c', '--category'),
              help='REQUIRED: Categorical variable of interest;
                must be a column header in table',
              default=NA, type = 'character'),
  make_option(c('--groups'),
              help='If >2 groups in category, the two of interest:
                comma-delimited (e.g. groups=Healthy,Sick)',
              default=NA, type = 'character'),
  make_option(c('-x', '--x_var'),
              help='REQUIRED: The independent variable (e.g. time);
                must be column header',
              default=NA, type = 'character'),
  make_option(c('-y', '--y_var'),
              help='REQUIRED: The response variable;
                must be acolumn header',
              default=NA, type = 'character'),
  make_option(c('-s', '--cases'),
              help='REQUIRED: The column header for your grouping
                (e.g. Patient_ID, User_Name, etc)',
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
              help='The spar parameter when fitting splines (0 - 1);
                default is calculated from data',
              default=NULL),
  make_option(c('--plot_path'),
              help='Plot the data too. Provide a filename/path,
                with extension explicitly in the name: pdf, png,
                tiff, jpg, bmp, eps, or ps.',
              default=NA, type = 'character')
)
opt <- parse_args(OptionParser(usage=usage, option_list=option_list))

if (is.na(opt$input) | is.na(opt$category) | is.na(opt$x_var) |
    is.na(opt$y_var) | is.na(opt$unit_id)) {
  stop('Missing required parameters. See usage and options (--help)')
}

# Parse command line
infile <- opt$input  # tsv file with data in long form
category <- opt$category  # the column header label for the two groups
groups <- opt$groups  # the two groups of interest, if column has >2
xvar <- opt$x_var  # the time series label
yvar <- opt$y_var  # the response variable label
cases <- opt$cases  # the header label defining the individuals (patients, etc)
perms <- as.numeric(opt$perms)  # default 999
cut_low <- opt$cut
spar_param <- opt$spar # default NULL
ints <- opt$intervals 
plot_path <- opt$plot  # name of the plot file.png

# TODO: Wrap the package function

df <- read.delim(file = infile, header = 1, check.names = F, sep = '\t')

result <- permuspliner(data = df, xvar = xvar, yvar = yvar, category = category,
                       cases = cases, groups = groups, ints = ints, perms = perms,
                       set_spar = spar_param, cut_low = cut_low)

if (!is.na(plot_path)) {
  if (is.na(groups)) {
    if (length(unique(df[, category])) > 2) {
      stop('More than two groups in category column. Define groups with "--groups=Name1,Name2"')
    }
    v1 <- unique(df[, category])[1]
    v2 <- unique(df[, category])[2]
  } else {
    v1 <- strsplit(groups, ',')[[1]][1]
    v2 <- strsplit(groups, ',')[[1]][2]
  }
  if (!is.na(cut_low)) {
    cut_low <- as.numeric(cut_low)
    keep.ids <- data.frame(table(df[, unit.id]))
    keep.ids <- as.character(keep.ids[keep.ids$Freq > cut_low, ]$Var1)
    df <- df[df[,unit.id] %in% keep.ids, ]
  }
  df.v1 <- df %>% filter(df[, category] == v1 & !is.na(df[, xvar]))
  df.v2 <- df %>% filter(df[, category] == v2 & !is.na(df[, xvar]))
  df.p <- rbind(df.v1, df.v2)
  df.pick <- c(xvar, category, yvar)
  plot.df <- df.p[, df.pick]
  plot.df <- plot.df[!is.na(plot.df[, xvar]), ]
  plot.df <- droplevels(plot.df)
  p <- ggplot(plot.df, aes(x=plot.df[,xvar], y=plot.df[,yvar],
                           color=as.character(plot.df[,category]))) +
    geom_point() + geom_smooth(span = spar_param, method = 'loess') + xlab(xvar) + ylab(yvar) +
    scale_color_manual(name=category, values = c("#0072B2","#D55E00")) +
    theme_bw() +
    theme(panel.background = element_blank(), panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_text(colour = 'black'), 
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'))
  
  ggsave(plot_path, height = 3.5, width = 4, units = 'in', dpi = 600)
  cat(paste('Plot saved\n'))
}

