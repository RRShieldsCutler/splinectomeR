#' @title Permuspliner
#' @usage Permutation to test whether there is a non-zero trend among a set
#' of individuals/samples over a continuous variable -- such as time. So, there
#' does not need to be two groups in this test. The x variable datapoints are
#' permuated within each case/individual, thus maintaining the distribution in
#' the y component but shuffling the hypothesized trend.
#' @rdname permuspliner
#' @description Tests for a significant difference between two groups overall.
#' @param data A dataframe object containing your data.
#' @param x The independent variable; is continuous, e.g. time.
#' @param y The dependent variable; is continuous, e.g. temperature.
#' @param category The column name of the category to be tested.
#' @param cases The column name defining the individual cases, e.g. patients.
#' @param groups If more than two groups, the two groups to compare.
#' @param perms The number of permutations to generate
#' @param set_spar Set the spar parameter for splines
#' @param set_tol In rare cases, must manually set the tol parameter (default 1e-4)
#' @param cut_low Remove individual cases with fewer than _ observations
#' @param ints Number of x intervals over which to measure area
#' @param quiet Silence all text outputs
#' @export
#' @examples 
#' result <- permuspliner(data = ChickWeight, xvar = 'Time',
#'              yvar = 'weight', category = 'Diet',
#'              cases = 'Chick', groups = '1,2')
#' result$pval


permuspliner <- function(data = NA, xvar = NA, yvar = NA, category = NA,
                         cases = NA, groups = NA, perms = 999, set_spar = NULL,
                         cut_low = NA, ints = 1000, quiet = FALSE, set_tol = 1e-4) {
  
  require(dplyr)
  
  reqs = c(data, category, xvar, yvar, cases)
  if (any(is.na(reqs))) {
    stop('Missing required parameters. Run ?permusplinr to see help docs')
  }
  
  perms = as.numeric(perms)
  ints = as.numeric(ints)
  cases = as.character(cases)
  groups = as.character(groups)
  
  in_df <- data
  if (is.na(groups)) {
    if (length(unique(in_df[, category])) > 2) {
      stop('More than two groups in category column. Define groups with "--groups=Name1,Name2"')
    }
    v1 <- unique(in_df[, category])[1]
    v2 <- unique(in_df[, category])[2]
  } else {
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
  df_v1 <- in_df %>% filter(in_df[, category] == v1 & !is.na(in_df[, xvar]))
  df_v2 <- in_df %>% filter(in_df[, category] == v2 & !is.na(in_df[, xvar]))
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
  .spline_permute <- function(randy, cases, category, xvar, yvar) {
    randy_meta <- randy %>% distinct_(cases, .keep_all = T)
    randy_meta$case_shuff <- sample(randy_meta[, category])
    randy_meta <- randy_meta %>% select_(cases, case_shuff)
    randy <- merge(randy, randy_meta, by = cases, all = T)
    randy_v1 <- filter(randy, case_shuff == v1 & !is.na(randy[,xvar]))
    randy_v2 <- filter(randy, case_shuff == v2 & !is.na(randy[,xvar]))
    randy_v1_spl <- with(randy_v1,
                        smooth.spline(x=randy_v1[, xvar], y=randy_v1[, yvar],
                                      spar = set_spar, tol = set_tol))
    randy_v2_spl <- with(randy_v2,
                        smooth.spline(x=randy_v2[, xvar], y=randy_v2[, yvar],
                                      spar = set_spar, tol = set_tol))
    x0 <- max(c(min(randy_v1_spl$x)), min(randy_v2_spl$x))
    x1 <- min(c(max(randy_v1_spl$x)), max(randy_v2_spl$x))
    xby <- (x1 - x0) / (ints - 1)
    xx <- seq(x0, x1, by = xby)
    randy_v1_fit <- data.frame(predict(randy_v1_spl, xx))
    colnames(randy_v1_fit) <- c('x', 'var1')
    randy_v2_fit <- data.frame(predict(randy_v2_spl, xx))
    colnames(randy_v2_fit) <- c('x', 'var2')
    spl_dist <- merge(randy_v1_fit, randy_v2_fit, by = 'x')
    spl_dist$abs_distance <- abs(spl_dist$var1 - spl_dist$var2)
    perm_area <- sum(spl_dist$abs_distance) / ints
    permuted <- append(permuted, perm_area)
    return(permuted)
  }
  
  # Run the permutation over desired number of iterations
  permuted <- list()
  permuted <- replicate(perms, 
                       .spline_permute(randy = in_df, cases, category, xvar, yvar))
  pval <- (sum(permuted >= as.numeric(real_area)) + 1) / (perms + 1)
  
  # Return the p-value
  if (quiet == FALSE) {
    cat(paste('\np-value =', round(pval, digits = 5), '\n\n'))
  }
  result <- list("pval" = pval,
                 "v1_interpolated" = v1_spl_f, "v2_interpolated" = v2_spl_f,
                 "v1_spline" = df_v1_spl, "v2_spline" = df_v2_spl)
  return(result)
}
