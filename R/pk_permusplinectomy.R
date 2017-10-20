#' @title Permuspliner
#' @usage Permutation to test whether there is a non-zero trend among a set
#' of individuals/samples over a continuous variable -- such as time. So, there
#' does not need to be two groups in this test. The x variable datapoints are
#' permuated within each case/individual, thus maintaining the distribution in
#' the y component but shuffling the hypothesized trend.
#' @rdname permuspliner
#' @description Tests for a significant difference between two groups overall.
#' @param data A dataframe object containing your data.
#' @param xvar The independent variable; is continuous, e.g. time.
#' @param yvar The dependent variable; is continuous, e.g. temperature.
#' @param category The column name of the category to be tested.
#' @param cases The column name defining the individual cases, e.g. patients.
#' @param groups If more than two groups, the two groups to compare.
#' @param perms The number of permutations to generate
#' @param retain_perm Retain permuted spline data for permutation confidence interval plotting (default FALSE for speed)
#' @param test_direction Test whether the groups are significantly 'more' or significantly 'less' distinct than expected by random chance. Default is 'more'.
#' @param set_spar Set the spar parameter for splines
#' @param set_tol In rare cases, must manually set the tol parameter (default 1e-4)
#' @param cut_low Remove individual cases with fewer than _ observations
#' @param cut_sparse Minimum number of total observations necessary per group to fit a spline (default 4)
#' @param ints Number of x intervals over which to measure area
#' @param quiet Silence all text outputs
#' @export
#' @examples 
#' result <- permuspliner(data = ChickWeight, xvar = 'Time',
#'              yvar = 'weight', category = 'Diet',
#'              cases = 'Chick', groups = '1,2')
#' result$pval


permuspliner <- function(data, xvar = NULL, yvar = NULL, category = NULL,
                         cases = NULL, groups = NULL, perms = 999, retain_perm = FALSE,
                         test_direction = 'more', set_spar = NULL, cut_low = NA,
                         ints = 1000, quiet = FALSE, set_tol = 1e-4, cut_sparse = 4) {

  # reqs = c(data, category, xvar, yvar, cases)
  if (is.null(data) | is.null(category) | is.null(xvar) | is.null(yvar) | is.null(cases)) {
    stop('Missing required parameter(s). Run ?permuspliner to see help docs')
  }
  
  if ((test_direction == 'more' | test_direction == 'less') == FALSE) {
    stop('Error in test direction option: must be either "more" or "less"')
  }
  
  perms = as.numeric(perms)
  ints = as.numeric(ints)
  cases = as.character(cases)
  groups = as.character(groups)
  
  in_df <- data
  if (is.null(groups)) {
    if (length(unique(in_df[, category])) > 2) {
      stop('More than two groups in category column. Define groups with (groups = "Name1,Name2")')
    }
    v1 <- unique(in_df[, category])[1]
    v2 <- unique(in_df[, category])[2]
  } else if (!is.null(groups)) {
    v1 <- strsplit(groups, ',')[[1]][1]
    v2 <- strsplit(groups, ',')[[1]][2]
  }
  
  if (!is.na(cut_low)) {
    cut_low <- as.numeric(cut_low)
    keep_ids <- data.frame(table(in_df[, cases]))
    keep_ids <- as.character(keep_ids[keep_ids$Freq > cut_low, ]$Var1)
    in_df <- in_df[in_df[, cases] %in% keep_ids, ]
  }
  if (quiet == FALSE) {
    cat(paste('\nTesting between', v1, 'and', v2, 'for a difference in', yvar, '\n'))
    cat(paste('\nScalpel please: performing permusplinectomy with', perms, 'permutations...\n'))
  }
  
  # The experimentally reported response
  df_v1 <- in_df[in_df[, category] %in% c(v1) & !is.na(in_df[, xvar]), ]
  df_v2 <- in_df[in_df[, category] %in% c(v2) & !is.na(in_df[, xvar]), ]
  if (length(df_v1[, xvar]) < cut_sparse | length(df_v2[, xvar]) < cut_sparse) {
    stop('Not enough data in each group to fit spline')
  }
  df_v1_spl <- with(df_v1,
                   smooth.spline(x=df_v1[, xvar], y=df_v1[, yvar],
                                 spar = set_spar, tol = set_tol))
  df_v2_spl <- with(df_v2,
                   smooth.spline(x=df_v2[, xvar], y=df_v2[, yvar],
                                 spar = set_spar, tol = set_tol))
  x0 <- max(c(min(df_v1_spl$x)), min(df_v2_spl$x))
  x1 <- min(c(max(df_v1_spl$x)), max(df_v2_spl$x))
  xby <- (x1 - x0) / (ints - 1)
  xx <- seq(x0, x1, by = xby)
  v1_spl_f <- data.frame(predict(df_v1_spl, xx))
  colnames(v1_spl_f) <- c('x', 'var1')
  v2_spl_f <- data.frame(predict(df_v2_spl, xx))
  colnames(v2_spl_f) <- c('x', 'var2')
  real_spl_dist <- merge(v1_spl_f, v2_spl_f, by = 'x')
  real_spl_dist$abs.distance <- abs(real_spl_dist$var1 - real_spl_dist$var2)
  real_area <- sum(real_spl_dist$abs.distance) / ints
  
  # Define the permutation function
  case_shuff <- 'case_shuff'
  .spline_permute <- function(randy) {
    randy_meta <- randy[!duplicated(randy[, cases]), ]
    randy_meta$case_shuff <- sample(randy_meta[, category])
    randy_meta <- randy_meta[, c(cases, case_shuff)]
    randy <- merge(randy, randy_meta, by = cases, all = T)
    randy_v1 <- randy[randy[, case_shuff] %in% c(v1) & !is.na(randy[, xvar]), ]
    randy_v2 <- randy[randy[, case_shuff] %in% c(v2) & !is.na(randy[, xvar]), ]
    randy_v1_spl <- with(randy_v1,
                        smooth.spline(x=randy_v1[, xvar], y=randy_v1[, yvar],
                                      spar = set_spar, tol = set_tol))
    randy_v2_spl <- with(randy_v2,
                        smooth.spline(x=randy_v2[, xvar], y=randy_v2[, yvar],
                                      spar = set_spar, tol = set_tol))
    # x0 <- max(c(min(randy_v1_spl$x)), min(randy_v2_spl$x))
    # x1 <- min(c(max(randy_v1_spl$x)), max(randy_v2_spl$x))
    # xby <- (x1 - x0) / (ints - 1)
    # xx <- seq(x0, x1, by = xby)
    randy_v1_fit <- data.frame(predict(randy_v1_spl, xx))
    colnames(randy_v1_fit) <- c('x', 'var1')
    randy_v2_fit <- data.frame(predict(randy_v2_spl, xx))
    colnames(randy_v2_fit) <- c('x', 'var2')
    spl_dist <- merge(randy_v1_fit, randy_v2_fit, by = 'x')
    spl_dist$abs_distance <- abs(spl_dist$var1 - spl_dist$var2)
    if (retain_perm == TRUE) {
      transfer_perms <- spl_dist[, 2:4]
      colnames(transfer_perms) <- c(paste0('v1perm_',ix),
                                    paste0('v2perm_',ix),
                                    paste0('pdistance_',ix))
      if (ix > 1) perm_retainer <- perm_output$perm_retainer
      perm_retainer <- cbind(perm_retainer, transfer_perms)
      perm_output$perm_retainer <- perm_retainer
      # print(summary(xx))
      
      perm_area <- sum(spl_dist$abs_distance) / ints
      if (ix > 1) permuted <- perm_output$permuted
      permuted <- append(permuted, perm_area)
      perm_output$permuted <- permuted
      return(perm_output)
    } else if (retain_perm == FALSE) {
      # print(summary(xx))
      # spl_dist$abs_distance <- abs(spl_dist$var1 - spl_dist$var2)
      perm_area <- sum(spl_dist$abs_distance) / ints
      permuted <- append(permuted, perm_area)
      return(permuted)
    }
  }
  
  # Run the permutation over desired number of iterations
  in_rand <- rbind(df_v1, df_v2)
  permuted <- list()
  if (retain_perm == TRUE) {
    perm_output <- list()
    perm_retainer <- data.frame(row.names = xx)
    for (ix in 1:perms) {
      perm_output <- .spline_permute(randy = in_rand)
    }
    if (test_direction == 'more') {
      pval <- (sum(perm_output$permuted >= as.numeric(real_area)) + 1) / (perms + 1)
    } else if (test_direction == 'less') {
      pval <- (sum(perm_output$permuted <= as.numeric(real_area)) + 1) / (perms + 1)
    }
    # perm_output <- replicate(perms, 
    #                       .spline_permute(randy = in_rand))
  } else if (retain_perm == FALSE) {
    permuted <- replicate(perms, 
                       .spline_permute(randy = in_rand))
    if (test_direction == 'more') {
      pval <- (sum(permuted >= as.numeric(real_area)) + 1) / (perms + 1)
    } else if (test_direction == 'less') {
      pval <- (sum(permuted <= as.numeric(real_area)) + 1) / (perms + 1)
    }
  }
  
  # Return the p-value
  if (quiet == FALSE) {
    cat(paste('\np-value =', round(pval, digits = 5), '\n\n'))
  }
  v1_data <- df_v1; v2_data <- df_v2
  v1_data[, category] <- droplevels(v1_data[, category])
  v2_data[, category] <- droplevels(v2_data[, category])
  
  # Return the results list
  if (retain_perm == TRUE) {
    result <- list("pval" = pval,
                   "v1_interpolated" = v1_spl_f, "v2_interpolated" = v2_spl_f,
                   "v1_spline" = df_v1_spl, "v2_spline" = df_v2_spl,
                   "permuted_splines" = perm_output$perm_retainer,
                   "true_distance" = real_spl_dist,
                   "v1_data" = v1_data, "v2_data" = v2_data)
  } else if (retain_perm == FALSE) {
    result <- list("pval" = pval,
                   "v1_interpolated" = v1_spl_f, "v2_interpolated" = v2_spl_f,
                   "v1_spline" = df_v1_spl, "v2_spline" = df_v2_spl,
                   "v1_data" = v1_data, "v2_data" = v2_data)
  }
  return(result)
}

#' @title Plot permuted distance distribution
#' @description Compare the permuted distances to the true distance. Requires permuspliner run with the "retain_perm" option.
#' @rdname permuspliner.plot.permdistance
#' @param data The results object from the permuspliner function
#' @param xlabel Name for the x axis in the plot
#' @export
#' @examples 
#' distplot <- permuspliner.plot.permdistance(results, xlabel='Time')
#' 
#' # View the plot
#' distplot
#' 
#' # Save the plot as PNG file
#' ggsave(distplot, file = 'my_plot.png', dpi=300, height=4, width=4)
#' 

  # Plot the distributions
permuspliner.plot.permdistance <- function(data, xlabel=NULL) {
  require(ggplot2)
  require(reshape2)
  if (is.null(xlabel)) xlabel <- 'longitudinal parameter'
  dists <- data['permuted_splines'][[1]]
  dists <- dists[, grep('pdistance', colnames(dists))]
  dists$x.par <- rownames(dists); rownames(dists) <- NULL
  true_dist <- data['true_distance'][[1]]
  true_dist <- true_dist[, c(1,4)]
  dists <- melt(dists, id.vars = 'x.par', variable.name = 'permutation',
                value.name = 'permuted_distance')
  p <- ggplot() +
    geom_line(data=dists, aes(x=as.numeric(x.par), y=as.numeric(permuted_distance),
                              group=factor(permutation)), 
              color='black', alpha=0.2, size=1) +
    geom_line(data=true_dist, aes(x=as.numeric(x), y=as.numeric(abs.distance)),
              color='red', size=1.5) +
    theme_classic() + theme(axis.text = element_text(color='black')) +
    xlab(xlabel) + ylab('group spline distance')
  p
}


#' @title Plot permuted splines along the real data
#' @description Compare how the permuted splines fit with the real data. Provides visual support for p values.
#' @rdname permuspliner.plot.permsplines
#' @param data The results object from the permuspliner function
#' @param xvar Name (as string) of the longitudinal x variable in the data
#' @param yvar Name (as string) of the response/y variable in the data
#' @export
#' @examples 
#' permsplot <- permuspliner.plot.permsplines(results, xvar='Time', yvar='Weight')
#' 
#' # View the plot
#' permsplot
#' 
#' # Save the plot as PNG file
#' ggsave(permsplot, file = 'my_plot.png', dpi=300, height=4, width=4)
#' 

permuspliner.plot.permsplines <- function(data, xvar=NULL, yvar=NULL) {
  if (is.null(data) | is.null(xvar) | is.null(yvar)) {
    stop('Missing required arugments.')
  }
  require(ggplot2)
  require(reshape2)
  if (is.null(xlabel)) xlabel <- 'longitudinal parameter'
  permsplines <- data['permuted_splines'][[1]]
  permsplines <- permsplines[, grep('perm', colnames(permsplines))]
  permsplines$x.par <- rownames(permsplines); rownames(permsplines) <- NULL
  permsplines <- melt(permsplines, id.vars = 'x.par', variable.name = 'permutation',
                      value.name = 'y.par')
  # permsplines$group <- NA
  # permsplines$group[grep('v1', permsplines$permutation)] <- 'Group_1'
  # permsplines$group[grep('v2', permsplines$permutation)] <- 'Group_2'
  # permsplines_1 <- permsplines[permsplines$group=='Group_1', ]
  # permsplines_2 <- permsplines[permsplines$group=='Group_2', ]
  true_v1 <- data['v1_data'][[1]]
  spar_v1 <- data['v1_spline'][[1]]$spar
  true_v2 <- data['v2_data'][[1]]
  spar_v2 <- data['v2_spline'][[1]]$spar
  p <- ggplot() +
    geom_line(data=permsplines, aes(x=as.numeric(x.par), y=as.numeric(y.par),
                                      group=factor(permutation)), color='black', alpha=0.1, size=1.2) +
    # geom_line(data=permsplines_2, aes(x=as.numeric(x.par), y=as.numeric(y.par),
    #                                   group=factor(permutation)), color='light blue', alpha=0.4, size=0.5) +
    geom_smooth(aes(x=as.numeric(true_v1[, xvar]), y=as.numeric(true_v1[, yvar])),
                color='red',
                size=1.5, span = spar_v1, method='loess', show.legend = F, se=F) +
    geom_smooth(aes(x=as.numeric(true_v2[, xvar]), y=as.numeric(true_v2[, yvar])),
                color='blue', 
                size=1.5, span = spar_v2, method='loess', show.legend = F, se=F) +
        theme_classic() + theme(axis.text = element_text(color='black')) +
    xlab(xvar) + ylab(yvar)
  p
}






