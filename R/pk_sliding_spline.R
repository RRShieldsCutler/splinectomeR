# The sliding spline function

#' @title Sliding_spliner
#' @usage Test usage
#' @details The data object needs to be organized with each observation as a
#' row, and have a column that identifies the case, patient, animal, etc, and
#' columns with the continuous x and y variables (row with x = NA will be
#' removed). If there are multiple groups in the data, you can filter to
#' the single group of interest with the category and group arguments. Otherwise
#' it assumes the entire dataset is the single population.
#' 
#' @rdname sliding_spliner
#' @description Test for a significant difference in two groups at imputed intervals
#' @param data A dataframe object containing your data.
#' @param x The independent variable; is continuous, e.g. time.
#' @param y The dependent variable; is continuous, e.g. temperature.
#' @param category The column name of the category to be tested, if present.
#' @param cases The column name defining the individual cases, e.g. patients.
#' @param groups If more than two groups, the two groups to compare.
#' @param set_spar Set the spar parameter for splines
#' @param cut_sparse Remove sparse cases with fewer than __ data points (default 4)
#' @param test_density Minimum density of cases in each group to report p-value (default 3)
#' @param ints Number of x intervals over which to measure significance
#' @param quiet Silence all text outputs
#' @export
#' @examples 
#' result <- trendyspliner(data = ChickWeight, xvar = 'Time',
#'              yvar = 'weight', category = 'Diet',
#'              cases = 'Chick', groups = '1,2')
#' result$pval
#' 

sliding_spliner <- function(data = NA, xvar = NA, yvar = NA, category = NA,
                            cases = NA, groups = NA, set_spar = NULL,
                            cut_sparse = 4, test_density = 3, ints = 100, quiet = FALSE) {
  
  require(dplyr)
  require(reshape2)
  
  reqs <- c(data, xvar, yvar, cases, category)
  if (any(is.na(reqs))) {
    stop('Missing required parameters. Run ?sliding_spliner to see help docs')
  }
  
  if (quiet == FALSE) {
    cat(paste('Running sliding spline test with', ints,
            'time points extrapolated from splines...\n'))
  }
  # Read infile and limit to ids matching sparsity parameter
  # df = read.delim(file = infile, header = 1, check.names = F, sep = '\t')
  # df = ChickWeight
  # category = 'Diet' # the column header label for the two groups
  # groups = '1,2'  # the two groups of interest, if column has >2
  # xvar = 'Time'  # the time series label
  # yvar = 'weight'  # the response variable label
  # cases = 'Chick'  # the header label defining the individuals (patients, etc)
  # ints = 100  # default 100
  # sparsity = 3  # cutoff for number of points at given time; default 3 (>= 3)
  # unit_number = 4  # min number of points needed per individual
  # set_spar = NULL  # default NULL
  
  df <- data
  
  if (is.na(groups)) {
    if (length(unique(df[, category])) > 2) {
      stop('More than two groups in category column. Define groups with "--groups=Name1,Name2"')
    }
    v1 <- unique(df[, category])[1]
    v2 <- unique(df[, category])[2]
  } else {
    groups <- as.character(groups)
    v1 <- strsplit(groups, ',')[[1]][1]
    v2 <- strsplit(groups, ',')[[1]][2]
  }
  
  # Make sure the df only contains data from categories tested
  cat.vars <- c(v1,v2)
  df <- df[df[, category] %in% cat.vars, ]
  cases.tab <- data.frame(table(df[, cases]))
  cases.notsparse <- cases.tab[cases.tab$Freq >= cut_sparse, ]
  cases.keep <- as.character(cases.notsparse$Var1)
  df <- df[df[, cases] %in% cases.keep, ]
  
  # Get the range of the independent variable (e.g. time) to set spline limits
  x.min <- min(as.numeric(df[, xvar]))
  x.max <- max(as.numeric(df[, xvar]))
  xrang <- seq(from = x.min, to = x.max, by = ((x.max - x.min) / (ints - 1)))
  
  # Save the group labels for each individual/unit
  df.groups <- df %>% 
              distinct_(cases, .keep_all = T) %>% 
              select_(cases, category)
  
  # Generate splines for each individual
  spl.table <- setNames(data.frame(xrang), c('x'))
  for (i in cases.keep) {
    unit.df <- subset(df, df[, cases]==i)
    unit.spl <- with(unit.df,
                    smooth.spline(x=unit.df[, xvar], y=unit.df[, yvar],
                                  spar = set_spar))
    xrang.i <- subset(xrang, xrang >= min(unit.spl$x) & xrang <= max(unit.spl$x))
    unit.spl.f <- data.frame(predict(unit.spl, xrang.i))
    colnames(unit.spl.f) <- c('x', i)
    spl.table <- merge(spl.table, unit.spl.f, by = 'x', all = T)
  }
  
  # Prepare the spline table for statistical testing
  spl.table.p <- tibble::column_to_rownames(df = spl.table, var = 'x')
  spl.table.p <- as.data.frame(t(spl.table.p))
  spl.table.p <- tibble::rownames_to_column(spl.table.p, var = cases)
  spl.table.p <- merge(spl.table.p, df.groups, by = cases, all = T)
  
  spline.table <- spl.table.p
  # Define the non-parametric test function
  mann_whitney_per_bit <- function(spline.table, n) {
    n <- as.character(n)
    x.pick <- c(cases, n, category)
    x.table <- spline.table[, x.pick]
    x.table <- x.table[!is.na(x.table[, n]), ]
    x.table <- droplevels(x.table)
    cat.freq <- data.frame(table(x.table[, category]))
    if (cat.freq$Freq[1] >= test_density & cat.freq$Freq[2] >= test_density) {
      mw <- wilcox.test(x.table[, n] ~ x.table[, category])
      pval <- mw$p.value
      num.pts <- cat.freq$Freq[1] + cat.freq$Freq[2]
      return(c(pval, n, num.pts))
    } else {}
  }
  
  # Run the stats test on each step of the spline on each data unit
  pvals.list <- list()
  pval.x <- list()
  pval.num.pts <- list()  # Save the number of total data points in each comparison
  for (n in xrang) {
    pval <- mann_whitney_per_bit(spl.table.p, n)
    pvals.list <- append(pvals.list, pval[1])
    pval.x <- append(pval.x, pval[2])
    pval.num.pts <- append(pval.num.pts, pval[3])
  }
  pval.df <- do.call(rbind, Map(data.frame,
                               x.series=as.numeric(pval.x),
                               p.value=as.numeric(pvals.list),
                               N=as.numeric(pval.num.pts)))
  colnames(pval.df) <- c(xvar, 'p_value', 'number_of_observations')
  
  plot.spline.data <- melt(data = spl.table, id.vars = 'x')
  colnames(plot.spline.data) <- c('x', cases, 'value')
  plot.spline.data <- merge(plot.spline.data, df.groups, by = cases, all = T)
  colnames(plot.spline.data) <- c('UNIT', 'x', 'value', 'category')
  
  result <- list('pval_table' = pval.df, 'spline_table' = spl.table,
                 'spline_longform' = plot.spline.data)
  return(result)
}

#' @title Plot group splines
#' @description Plot the individual splines grouped by two categories of interest
#' @rdname sliding_spliner.plot.splines
#' @param data The dataframe of imputed points generated from splines
#' @param xvar The independent variable; is continuous, e.g. time.
#' @param yvar The dependent variable; is continuous, e.g. temperature.
#' @param category The data category being compared
#' @export
#' @examples sliding_spliner.plot.splines(result$spline_longform,
#'  category = 'Diet', xvar = 'Time', yvar = 'weight')
#'

  # Plot the splines
sliding_spliner.plot.splines <- function(plot.spline.data, category = 'category',
                                         xvar = 'xvar', yvar = 'yvar') {
  require(ggplot2)
  ggplot(plot.spline.data, aes(x = x, y = value, group = UNIT,
    color = as.character(category))) + geom_line(na.rm = T) +
    scale_color_manual(name = category, values = c("#0072B2","#D55E00")) +
    xlab(xvar) + ylab(yvar) +
    theme_bw() + theme(legend.position='right',
                       plot.background = element_rect(color = 'white'),
                       panel.grid = element_blank())
}

#' @title Plot sliding spliner p-values
#' @description Plot the imputed p-values over the x variable range
#' @rdname sliding_spliner.plot.pvals
#' @param data The dataframe of p-values over the x variable
#' @param xvar The independent variable; is continuous, e.g. time.
#' @export
#' @examples sliding_spliner.plot.pvals(result$pval_table,
#'  xvar = 'Time')
#'


# Plot the p-values as function of the independent variable (e.g. time)
# Scale the size of the line and points according to the number of observations
sliding_spliner.plot.pvals <- function(data, xvar = 'xvar') {
  library(ggplot2)
  .norm_range <- function(x) {(x-min(x)) / (max(x)-min(x))}
  pval.p <- data
  xvar <- names(pval.p)[1]
  pval.p$N.norm <- .norm_range(x = pval.p[,3])
  ggplot(pval.p, aes(x=pval.p[,1], y=p_value)) + geom_line() +
    geom_point(shape = 20, size = (pval.p$N.norm * 2))  +
    geom_hline(aes(yintercept = 0.05), linetype='dashed') +
    xlab(xvar) + ylab('Mann-Whitney p-value') +
    scale_y_continuous(trans = scales::log10_trans()) +
    theme_bw() + theme(legend.position='none',
                       plot.background = element_rect(color = 'white'),
                       panel.grid = element_blank())
}

