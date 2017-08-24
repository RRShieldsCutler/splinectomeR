#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(splinectomeR))

usage = '\nPermutation to test whether there is a non-zero trend among a set
of individuals/samples over a continuous variable (such as time). So, there
does not need to be two groups in this test. The x variable datapoints are
permuated within each case/individual, thus maintaining the distribution in
the y component but shuffling the hypothesized trend.'

option_list = list(
  make_option(c('-i', '--input'),
              help='REQUIRED: Input is long-form tab-delimited table,
                prefiltered for missing values',
              default=NA, type = 'character'),
  make_option(c('-c', '--category'),
              help='REQUIRED: Categorical variable of interest;
                must be a column header in table',
              default=NA, type = 'character'),
  make_option(c('--group'),
              help='If >1 group in category, the group of interest
                (e.g. group=disease)',
              default=NA, type = 'character'),
  make_option(c('-x', '--x_var'),
              help='REQUIRED: The independent variable (e.g. time);
                must be column header',
              default=NA, type = 'character'),
  make_option(c('-y', '--y_var'),
              help='REQUIRED: The response variable; must be a column header',
              default=NA, type = 'character'),
  make_option(c('-s', '--cases'),
              help='REQUIRED: The column header for your case grouping
                (e.g. Patient_ID, User_Name, etc)',
              default=NA, type = 'character'),
  make_option(c('--perms'),
              help='Number of permutation shuffles [default %default]',
              default=999, type = 'integer'),
  make_option(c('--cut'),
              help='Cut Unit IDs that occur less than this many times',
              default=NA, type = 'character'),
  make_option(c('--mean_center'),
              help='Mean center the data by individual before
                processing [default %default]',
              default=TRUE),
  make_option(c('--intervals'),
              help='Number of sampling intervals along spline [default %default]',
              default=10000, type = 'integer'),
  make_option(c('--spar'),
              help='The spar parameter when fitting splines (0 - 1);
                default is derived from data',
              default=NULL),
  make_option(c('--plot_path'),
              help='Plot the data too. Provide a filename/path,
                with extension explicitly in the name: pdf, png,
                tiff, jpg, bmp, eps, or ps.',
              default=NA, type = 'character')
)
opt <- parse_args(OptionParser(usage=usage, option_list=option_list))

if (is.na(opt$input) | is.na(opt$category) | is.na(opt$x_var) |
    is.na(opt$y_var) | is.na(opt$cases)) {
  stop('Missing required parameters. See usage and options (--help)')
}

# Parse command line
infile <- opt$input  # tsv file with data in long form
category <- opt$category  # the column header label for the group
group <- opt$group  # the group of interest, if column has >1
xvar <- opt$x_var  # the time series label
yvar <- opt$y_var  # the response variable label
cases <- opt$cases  # the header label defining the individuals (patients, etc)
perms <- as.numeric(opt$perms)  # default 999
cut_low <- opt$cut
spar_param <- opt$spar # default NULL
ints <- opt$intervals 
mean_center <- opt$mean_center
plot_path <- opt$plot_path  # name of the plot file.png

# Read the data file
df <- read.delim(file = infile, header = 1, check.names = F, sep = '\t')
df_fun <- df

# Run the function and report p-value
result <- trendyspliner(data = df_fun, xvar = xvar, yvar = yvar, category = category,
              cases = cases, group = group, mean_center = mean_center, perms = perms,
              cut_low = cut_low, set_spar = spar_param, ints = ints)

# Plot results if called for
if (!is.na(plot_path)) {
  if (is.na(group)) {
    if (length(unique(df[, category])) > 1) {
      stop('More than one group in category column. Define group with "--group=Name1"')
    }
    v1 = unique(df[, category])[1]
  } else {
    v1 = as.character(group)
  }
    df.v1 = df %>% filter(df[, category] == v1 & !is.na(df[, xvar]))
    if (mean_center == TRUE) {
      all_ids <- as.character(unique(df.v1[, cases]))
      grp_mean <- mean(df.v1[, yvar])
      mean_center <- function(case_df, grp_mean) {
        offset <- (mean(case_df[, yvar]) - grp_mean)
        case_df$y_offset_adj <- (case_df[, yvar] - offset)
        case_df[, yvar] <- NULL
        names(case_df)[names(case_df) == 'y_offset_adj'] <- yvar
        return(case_df)
    }
    df.adj <- list()
    i = 1
    for (id in all_ids) {
      case_df <- df.v1[df.v1[, cases] == id, ]
      case_df <- mean_center(case_df, grp_mean)
      df.adj[[i]] <- case_df
      i = i + 1
    }
    df.v1 <- do.call(rbind, df.adj)
  }
  df.pick <- c(xvar, category, yvar)
  plot.df <- df.v1[, df.pick]
  plot.df <- plot.df[!is.na(plot.df[, xvar]), ]
  plot.df <- droplevels(plot.df)
  p <- ggplot(plot.df, aes(x=plot.df[, xvar], y=plot.df[, yvar], 
                           color=as.character(plot.df[, category]))) +
    geom_point() + geom_smooth(span = spar_param, method = 'loess') +
    xlab(xvar) + ylab(yvar) + theme_bw() +
    scale_color_manual(name=category, values = c("#0072B2")) +
    theme(panel.background = element_blank(), panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_text(colour = 'black'), 
          axis.line = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'))
  
  ggsave(plot_path, height = 3.5, width = 4, units = 'in', dpi = 600)
  cat(paste('Plot saved\n'))
}

