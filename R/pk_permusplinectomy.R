# The permusplinectomy function
#
#' @title Permuspliner
#' @usage 
#' 
#' @description Tests for a significant difference between two groups overall
#' @param data A dataframe object containing your data
#' @param xvar The independent variable; is continuous, e.g. time
#' @param yvar The dependent variable; is continuous, e.g. temperature
#' @param category The column name of the category to be tested
#' @param unit The column name defining the individual cases, e.g. patients
#' @param groups If more than two groups, the two groups to compare
#' @param perms
#' @param spar
#' @param plot
#' @param cut
#' @param ints
#' 
#' @details 
#' 
#' 
#' @export
#' @examples
#' permuspliner(data = ChickWeight, perms = 99)
#' permuspliner(data = ChickWeight, perms = 99)


