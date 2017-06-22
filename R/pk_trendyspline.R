# The trendyspline function

#' @title Trendyspliner
#' @usage Permutation to test whether there is a non-zero trend among a set
#' of individuals/samples over a continuous variable -- such as time. So, there
#' does not need to be two groups in this test. The x variable datapoints are
#' permuated within each case/individual, thus maintaining the distribution in
#' the y component but shuffling the hypothesized trend.
#' @details The data object needs to be organized with each observation as a
#' row, and have a column that identifies the case, patient, animal, etc, and
#' columns with the continuous x and y variables (row with x = NA will be
#' removed). If there are multiple groups in the data, you can filter to
#' the single group of interest with the category and group arguments. Otherwise
#' it assumes the entire dataset is the single population.
#' 
#' @rdname trendyspliner
#' @description Test for a significant non-zero trend in a response over time.
#' @param data A dataframe object containing your data.
#' @param x The independent variable; is continuous, e.g. time.
#' @param y The dependent variable; is continuous, e.g. temperature.
#' @param category The column name of the category to be tested, if present.
#' @param cases The column name defining the individual cases, e.g. patients.
#' @param group If more than two groups, the two groups to compare.
#' @param mean_center Before processing, mean center data by individual/case
#' @param perms The number of permutations to generate
#' @param set_spar Set the spar parameter for splines
#' @param cut_low Remove data with fewer than __ points
#' @param ints Number of x intervals over which to measure area
#' @param quiet Silence all text outputs
#' @export
#' @examples 
#' trendyspliner(data = ChickWeight, xvar = 'Time',
#'              yvar = 'weight', category = 'Diet',
#'              cases = 'Chick', group = '1')

trendyspliner <- function(data = NA, xvar = NA, yvar = NA, category = NA,
                         cases = NA, group = NA, mean_center = TRUE, perms = 99, set_spar = NULL,
                         cut_low = NA, ints = 1000, quiet = FALSE) {
  
  require(dplyr)
  
  reqs = c(data, xvar, yvar, cases)
  if (any(is.na(reqs))) {
    stop('Missing required parameters. Run ?trendyspliner to see help docs')
  }
  
  perms <- as.numeric(perms)
  ints <- as.numeric(ints)
  cases <- as.character(cases)
  
  df <- data
  
  
  if (!is.na(category)) {
    if (is.na(group)) {
      group <- as.character(group)
      if (length(unique(df[, category])) > 1) {
        stop('More than one group in category column. Define group.')
      }
      v1 <- unique(df[, category])[1]
    } else {
      v1 <- group
    }
    df <- df %>% filter(df[, category] == v1)
  }
  
  df <- df %>% filter(!is.na(df[, xvar]))
  
  if (!is.na(cut_low)) {
    cut_low <- as.numeric(cut_low)
    keep.ids <- data.frame(table(df[, cases]))
    keep.ids <- as.character(keep.ids[keep.ids$Freq > cut_low, ]$Var1)
    df <- df[df[,cases] %in% keep.ids, ]
  }
  
  if (quiet == FALSE) {
    cat(paste('\nTesting for a non zero trend in', yvar, 'across', xvar, '...\n'))
    cat(paste('\nPerforming the trendyspline test with', perms, 'permutations...\n'))
  }
  
  
  if (mean_center == TRUE) {
    if (quiet == FALSE) {
      cat(paste('\nMean centering the data...\n'))
    }
    all_ids <- as.character(unique(df[, cases]))
    grp_mean <- mean(df[, yvar])
    mean_center <- function(case_df, grp_mean) {
      offset <- (mean(case_df[, yvar]) - grp_mean)
      case_df$y_offset_adj <- (case_df[, yvar] - offset)
      case_df[, yvar] <- NULL
      names(case_df)[names(case_df) == 'y_offset_adj'] <- yvar
      return(case_df)
    }
    df_adj <- list()
    i = 1
    for (id in all_ids) {
      case_df <- df[df[, cases] == id, ]
      case_df <- mean_center(case_df, grp_mean)
      df_adj[[i]] <- case_df
      i = i + 1
    }
    df <- do.call(rbind, df_adj)
  }
  
  
  ## First determine the group mean (null hypothesis for changing over time)
  # y.mean <- mean(df[, yvar])
  df.spl <- with(df,
                   smooth.spline(x=df[, xvar], y=df[, yvar],
                                 spar = set_spar))
  x0 <- min(df.spl$x)
  x1 <- max(df.spl$x)
  xby <- (x1 - x0) / (ints - 1)
  xx <- seq(x0, x1, by = xby)
  spl.fit <- data.frame(predict(df.spl, xx))
  colnames(spl.fit) <- c('x', 'var1')
  y_base = spl.fit$var1[1]
  real.spl.dist <- spl.fit
  spl.fit$y_base <- y_base
  spl.fit$distance <- (spl.fit$var1 - spl.fit$y_base)
  real.area <- sum(spl.fit$distance) / ints
  
  # Define the permutation function
  y_shuff <- 'y_shuff'
  .spline_permute <- function(randy, cases, xvar, yvar) {
    randy.meta <- randy %>% select_(cases, xvar)
    randy.meta$y_shuff <- sample(randy[, yvar])
    randy.spl <- with(randy.meta, smooth.spline(x = randy.meta[, xvar],
                                              y = randy.meta$y_shuff,
                                              spar = set_spar))
    x0 <- min(randy.spl$x)
    x1 <- max(randy.spl$x)
    xby <- (x1 - x0) / (ints - 1)
    xx <- seq(x0, x1, by = xby)
    randy.fit <- data.frame(predict(randy.spl, xx))
    colnames(randy.fit) <- c('x', 'var1')
    p_base = randy.fit$var1[1]
    randy.fit$p_base <- p_base
    randy.fit$distance <- (randy.fit$var1 - randy.fit$p_base)
    perm.area <- sum(randy.fit$distance) / ints
    permuted <- append(permuted, perm.area)
    return(permuted)
  }
  
  # Run the permutation over desired number of iterations
  permuted <- list()
  permuted <- replicate(perms, 
                       .spline_permute(randy = df, cases, xvar, yvar))
  pval <- (sum(lapply(permuted, abs) >= abs(real.area)) + 1) / (perms + 1)
  
  # Return the p-value
  if (quiet == FALSE) {
    cat(paste('\np-value =', round(pval, digits = 5), '\n\n'))
  }
  result <- list("pval" = pval)
  return(result)
}
