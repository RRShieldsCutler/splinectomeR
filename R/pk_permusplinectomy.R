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
#' @param groups If more than two groups, the two groups to compare as character vector.
#' @param perms The number of permutations to generate
#' @param retain_perm Retain permuted spline data for permutation confidence interval plotting (set to FALSE for less memory)
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
#'              cases = 'Chick', groups = c(1,2))
#' result$pval


permuspliner <- function(data = NULL, xvar = NULL, yvar = NULL, category = NULL,
                         cases = NULL, groups = NA, perms = 999, retain_perm = TRUE,
                         test_direction = 'more', set_spar = NULL, cut_low = NA,
                         ints = 1000, quiet = FALSE, set_tol = 1e-4, cut_sparse = 4) {

  # reqs = c(data, category, xvar, yvar, cases)
  if (missing(data) | missing(category) | missing(xvar) | missing(yvar) | missing(cases)) {
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
  # Determine the two groups to compare
  if (is.na(groups[1])) {
    if (length(unique(in_df[, category])) > 2) {
      stop('More than two groups in category column. Define groups with (groups = c("Name1","Name2"))')
    }
    v1 <- as.character(unique(in_df[, category])[1])
    v2 <- as.character(unique(in_df[, category])[2])
  } else {
    v1 <- as.character(groups[1])
    v2 <- as.character(groups[2])
  }
  # Trim data if some cases have too few observations
  if (!is.na(cut_low)) {
    cut_low <- as.numeric(cut_low)
    keep_ids <- data.frame(table(in_df[, cases]))
    keep_ids <- as.character(keep_ids[keep_ids$Freq > cut_low, ]$Var1)
    in_df <- in_df[in_df[, cases] %in% keep_ids, ]
  }
  if (quiet == FALSE) {
    cat(paste('\nGroups detected:', v1, 'and', v2, '.\n'))
    cat(paste('\nNow testing between variables', v1, 'and', v2, 'for a difference in the response labeled', yvar, '\n'))
    cat(paste('\nScalpel please: performing permusplinectomy with', perms, 'permutations...\n'))
  }
  
  # The experimentally reported response
  df_v1 <- in_df[in_df[, category] %in% c(v1) & !is.na(in_df[, xvar]), ]
  df_v2 <- in_df[in_df[, category] %in% c(v2) & !is.na(in_df[, xvar]), ]
  if (length(df_v1[, xvar]) < cut_sparse | length(df_v2[, xvar]) < cut_sparse) {
    stop('Not enough data in each group to fit spline')
  }
  # Prevent issues arising from identical case labels across groups
  if (length(intersect(df_v1[, cases], df_v2[, cases])) > 0) {
    stop('\nIt appears there may be identically labeled cases in both groups.\n
         ...Please ensure that the cases are uniquely labeled between the two groups\n')
  }
  # Fit the splines for each group
  df_v1_spl <- with(df_v1,
                   smooth.spline(x=df_v1[, xvar], y=df_v1[, yvar],
                                 spar = set_spar, tol = set_tol))
  df_v2_spl <- with(df_v2,
                   smooth.spline(x=df_v2[, xvar], y=df_v2[, yvar],
                                 spar = set_spar, tol = set_tol))
  x0 <- max(c(min(df_v1_spl$x)), min(df_v2_spl$x))
  x1 <- min(c(max(df_v1_spl$x)), max(df_v2_spl$x))
  x0 <- x0 + ((x1 - x0) * 0.1)  # Trim the first and last 10% to avoid low-density artifacts
  x1 <- x1 - ((x1 - x0) * 0.1)
  xby <- (x1 - x0) / (ints - 1)
  xx <- seq(x0, x1, by = xby)  # Set the interval range
  v1_spl_f <- data.frame(predict(df_v1_spl, xx))  # Interpolate across the spline
  colnames(v1_spl_f) <- c('x', 'var1')
  v2_spl_f <- data.frame(predict(df_v2_spl, xx))
  colnames(v2_spl_f) <- c('x', 'var2')
  real_spl_dist <- merge(v1_spl_f, v2_spl_f, by = 'x')
  real_spl_dist$abs.distance <- abs(real_spl_dist$var1 - real_spl_dist$var2)  # Measure the real group distance
  real_area <- sum(real_spl_dist$abs.distance) / ints  # Calculate the area between the groups
  if (quiet == FALSE) {
    cat(paste('\nArea between groups successfully calculated, now spinning up permutations...\n'))
  }
  # Define the permutation function
  case_shuff <- 'case_shuff'  # Dummy label
  .spline_permute <- function(randy) {
    randy_meta <- randy[!duplicated(randy[, cases]), ]  # Pull out the individual IDs
    randy_meta$case_shuff <- sample(randy_meta[, category])  # Shuffle the labels
    randy_meta <- randy_meta[, c(cases, case_shuff)]
    randy <- merge(randy, randy_meta, by = cases, all = T)
    randy_v1 <- randy[randy[, case_shuff] %in% c(v1) & !is.na(randy[, xvar]), ]
    randy_v2 <- randy[randy[, case_shuff] %in% c(v2) & !is.na(randy[, xvar]), ]
    # Fit the splines for the permuted groups
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
    spl_dist$abs_distance <- abs(spl_dist$var1 - spl_dist$var2)  # Calculate the distance between permuted groups
    if (retain_perm == TRUE) {
      transfer_perms <- spl_dist[, 2:4]
      colnames(transfer_perms) <- c(paste0('v1perm_',ix),
                                    paste0('v2perm_',ix),
                                    paste0('pdistance_',ix))
      if (ix > 1) perm_retainer <- perm_output$perm_retainer
      perm_retainer <- cbind(perm_retainer, transfer_perms)
      perm_output$perm_retainer <- perm_retainer
      perm_area <- sum(spl_dist$abs_distance) / ints  # Calculate the area between permuted groups
      
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
    if (quiet == FALSE) {
      cat(paste('...permutations completed...\n'))
    }
    if (test_direction == 'more') {
      pval <- (sum(perm_output$permuted >= as.numeric(real_area)) + 1) / (perms + 1)
    } else if (test_direction == 'less') {
      pval <- (sum(perm_output$permuted <= as.numeric(real_area)) + 1) / (perms + 1)
    }
  } else if (retain_perm == FALSE) {
    permuted <- replicate(perms, 
                       .spline_permute(randy = in_rand))
    if (quiet == FALSE) {
      cat(paste('...permutations completed...\n'))
    }
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
  # Return the filtered data used for the splines and permutations
  v1_data <- df_v1; v2_data <- df_v2
  v1_data[, category] <- droplevels(factor(v1_data[, category]))
  v2_data[, category] <- droplevels(factor(v2_data[, category]))
  
  # Return the results list
  if (retain_perm == TRUE) {
    result <- list("pval" = pval, "category_1" = v1, "category_2" = v2,
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
  if (quiet == FALSE) {
    cat(paste('To visualize your results, try the following command, where "data" is your results object:'))
    cat(paste0('\npermuspliner.plot.permdistance(data, xlabel="', xvar,'")'))
    if (retain_perm == TRUE) {
      cat(paste0('\nor\npermuspliner.plot.permsplines(data, xvar="', xvar, '", yvar="', yvar, '")'))
    }
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
  num_perms <- ncol(dists)
  dists$x.par <- rownames(dists); rownames(dists) <- NULL
  true_dist <- data['true_distance'][[1]]
  true_dist <- true_dist[, c(1,4)]
  dists <- melt(dists, id.vars = 'x.par', variable.name = 'permutation',
                value.name = 'permuted_distance')
  if (num_perms > 1000) {
    alpha_level <- 0.002
  } else if (num_perms >= 100) {
    alpha_level <- 0.02
  } else if (num_perms < 100) {
    alpha_level <- 0.25
  }
  
  p <- ggplot() +
    geom_line(data=dists, aes(x=as.numeric(x.par), y=as.numeric(permuted_distance),
                              group=factor(permutation)), 
              color='black', alpha=alpha_level, size=1) +
    geom_line(data=true_dist, aes(x=as.numeric(x), y=as.numeric(abs.distance)),
              color='red', size=1.5) +
    theme_classic() + theme(axis.text = element_text(color='black')) +
    xlab(xlabel) + ylab('group spline distance')
  return(p)
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

permuspliner.plot.permsplines <- function(data = NULL, xvar=NULL, yvar=NULL) {
  if (is.null(data) | is.null(xvar) | is.null(yvar)) {
    stop('Missing required arugments.')
  }
  if (is.null(data['permuted_splines'][[1]])) {
    stop('Permuted data not in results. Did you set "retain_perms = TRUE" in the permuspliner() run?')
  }
  require(ggplot2)
  require(reshape2)
  permsplines <- data['permuted_splines'][[1]]
  permsplines <- permsplines[, grep('perm', colnames(permsplines))]
  num_perms <- (ncol(permsplines) / 2)
  permsplines$x.par <- rownames(permsplines); rownames(permsplines) <- NULL
  permsplines <- melt(permsplines, id.vars = 'x.par', variable.name = 'permutation',
                      value.name = 'y.par')
  var_1 <- as.character(data['category_1'][[1]])
  var_2 <- as.character(data['category_2'][[1]])
  true_v1 <- data['v1_interpolated'][[1]]
  true_v1$var <- var_1
  colnames(true_v1)[2] <- 'y'
  true_v2 <- data['v2_interpolated'][[1]]
  true_v2$var <- var_2
  colnames(true_v2)[2] <- 'y'
  true_data <- rbind(true_v1, true_v2)
  var_labels <- factor(true_data$var, ordered = T)
  true_data$var <- var_labels
  num_points <- length(true_v1$x)
  if (num_perms > 1000) {
    alpha_level <- 0.002
  } else if (num_perms > 100) {
    alpha_level <- 0.01
  } else if (num_perms <= 100) {
    alpha_level <- 0.31
  }
  
  p <- ggplot() +
    geom_line(data=permsplines, aes(x=as.numeric(x.par), y=as.numeric(y.par),
                                      group=factor(permutation)),
              alpha=alpha_level, size=1) +
    geom_line(data=true_data, aes(x=as.numeric(x), y=as.numeric(y), color=var_labels),
              size=1.2) +
    scale_color_manual(name='', values=c('red','blue')) +
        theme_classic() + theme(axis.text = element_text(color='black')) +
    xlab(xvar) + ylab(yvar)
  return(p)
}


