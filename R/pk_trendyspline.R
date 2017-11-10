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
#' @param group If more than one group in the data, the group to compare.
#' @param mean_center Before processing, mean center data by individual/case (default FALSE)
#' @param perms The number of permutations to generate
#' @param set_spar Set the spar parameter for splines
#' @param cut_low Remove data with fewer than __ points
#' @param ints Number of x intervals over which to measure area
#' @param quiet Silence all text outputs
#' @export
#' @examples 
#' result <- trendyspliner(data = ChickWeight, xvar = 'Time',
#'              yvar = 'weight', category = 'Diet',
#'              cases = 'Chick', group = '1')
#' result$pval

trendyspliner <- function(data = NULL, xvar = NULL, yvar = NULL, category = NULL,
                         cases = NULL, group = NULL, mean_center = FALSE, perms = 999, set_spar = NULL,
                         cut_low = NULL, ints = 1000, quiet = FALSE) {
  
  
  reqs = c(data, xvar, yvar, cases)
  if (any(is.null(reqs))) {
    stop('Missing required parameters. Run ?trendyspliner to see help docs')
  }
  
  perms <- as.numeric(perms)
  ints <- as.numeric(ints)
  cases <- as.character(cases)
  
  df <- data
  
  
  if (!is.null(category)) {
    if (is.null(group)) {
      group <- as.character(group)
      if (length(unique(df[, category])) > 1) {
        stop('More than one group in category column. Define group.')
      }
      v1 <- unique(df[, category])[1]
    } else {
      v1 <- group
    }
    df <- df[df[, category] %in% c(v1), ]
  } else if (is.null(category)) df <- df
  df <- df[!is.na(df[, xvar]), ]
  
  if (!is.null(cut_low)) {
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
  
  
  df.spl <- with(df,
                   smooth.spline(x=df[, xvar], y=df[, yvar],
                                 spar = set_spar))
  x0 <- min(df.spl$x)
  x1 <- max(df.spl$x)
  x0 <- x0 + ((x1 - x0) * 0.05)  # Trim the first and last 5% to avoid low-density artifacts
  x1 <- x1 - ((x1 - x0) * 0.05)
  xby <- (x1 - x0) / (ints - 1)
  xx <- seq(from = x0, to = x1, by = xby)
  spl.fit <- data.frame(predict(df.spl, xx))
  colnames(spl.fit) <- c('x', 'var1')
  y_base = spl.fit$var1[1]   # Null hypothesis: trend doesn't vary from zero over x axis
  real.spl.dist <- spl.fit
  spl.fit$y_base <- y_base
  spl.fit$distance <- (spl.fit$var1 - spl.fit$y_base)
  real.area <- sum(spl.fit$distance) / ints
  
  # Define the permutation function
  y_shuff <- 'y_shuff'
  .spline_permute <- function(randy) {
    randy.meta <- randy[, c(cases, xvar)]
    randy.meta$y_shuff <- sample(randy[, yvar])
    randy.spl <- with(randy.meta, smooth.spline(x = randy.meta[, xvar],
                                              y = randy.meta$y_shuff,
                                              spar = set_spar))
    randy.fit <- data.frame(predict(randy.spl, xx))
    colnames(randy.fit) <- c('x', 'var1')
    transfer.perms <- randy.fit
    colnames(transfer.perms)[2] <- c(paste0('perm_', ix))
    if (ix > 1) {
      perm_retainer <- perm_output$perm_retainer
      perm_retainer <- merge(perm_retainer, transfer.perms, by = 'x')
    } else perm_retainer <- transfer.perms
    perm_output$perm_retainer <- perm_retainer
    # if (ix > 1) colnames(perm_output)[ncol(perm_output)] <- paste0('perm_', ix)
    p_base = randy.fit$var1[1]
    randy.fit$p_base <- p_base
    randy.fit$distance <- (randy.fit$var1 - randy.fit$p_base)
    perm.area <- sum(randy.fit$distance) / ints
    permuted <- append(permuted, perm.area)
    perm_output$permuted <- permuted
    return(perm_output)
  }
  
  # Run the permutation over desired number of iterations
  permuted <- list()
  perm_output <- list()
  perm_retainer <- data.frame()
  for (ix in 1:perms) {
    perm_output <- .spline_permute(randy = df)
  }
  pval <- (sum(lapply(perm_output$permuted, abs) >= abs(real.area)) + 1) / (perms + 1)
  
  # Return the p-value
  if (quiet == FALSE) {
    cat(paste('\np-value =', round(pval, digits = 5), '\n\n'))
  }
  result <- list("pval" = pval, "imputed_curve" = spl.fit[, c(1,2,4)], "group_spline" = df.spl,
                 "permuted_splines" = perm_output$perm_retainer)
  return(result)
}



#' @title Plot permuted trends behind the real data
#' @description Compare how the permuted trend splines fit with the real data. Provides visual support for p values.
#' @rdname trendyspliner.plot.perms
#' @param data Required: The results object from the trendyspliner function
#' @param xlabel Optional: Title (as string) to print on the x axis
#' @param ylabel Optional: Title (as string) to print on the y axis
#' @export
#' @examples 
#' 
#' # Loads ggplot2 dependency
#' permsplot <- trendyspliner.plot.perms(results, xlabel='Time', ylabel='Weight')
#' 
#' # View the plot
#' permsplot
#' 
#' # Save the plot as PNG file
#' ggsave(permsplot, file = 'my_plot.png', dpi=300, height=4, width=4)
#' 


trendyspliner.plot.perms <- function(data = NULL, xlabel=NULL, ylabel=NULL) {
  require(ggplot2)
  require(reshape2)
  
  if (is.null(data)) stop('Missing required arugments.')
  if (is.null(xlabel)) xlabel <- 'longitudinal axis'
  if (is.null(ylabel)) ylabel <- 'response axis'
  
  permsplines <- data['permuted_splines'][[1]]
  num_perms <- ncol(permsplines)
  permsplines <- melt(permsplines, id.vars = 'x', variable.name = 'permutation',
                      value.name = 'y.par')
  true_df <- data['imputed_curve'][[1]][, 1:2]
  colnames(true_df) <- c('x','y')
  num_points <- length(true_df$x)
  if (num_perms > 1000) {
    alpha_level <- 0.002
  } else if (num_perms > 100) {
    alpha_level <- 0.02
  } else if (num_perms <= 100) {
    alpha_level <- 0.4
  }
  
  p <- ggplot() +
    geom_line(data=permsplines, aes(x=as.numeric(x), y=as.numeric(y.par),
                                    group=factor(permutation)),
              alpha=alpha_level, size=0.9) +
    geom_line(data=true_df, aes(x=as.numeric(x), y=as.numeric(y)),
              size=1.2, color = 'red') +
    theme_classic() + theme(axis.text = element_text(color='black')) +
    xlab(xlabel) + ylab(ylabel)
  return(p)
}


