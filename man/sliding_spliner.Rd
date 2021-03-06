% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_sliding_spline.R
\name{sliding_spliner}
\alias{sliding_spliner}
\title{Sliding_spliner}
\usage{
Test usage
}
\arguments{
\item{data}{A dataframe object containing your data.}

\item{category}{The column name of the category to be tested, if present.}

\item{cases}{The column name defining the individual cases, e.g. patients.}

\item{groups}{If more than two groups, the two groups to compare as character vector.}

\item{set_spar}{Set the spar parameter for splines}

\item{cut_low}{Remove low prevalence with fewer than __ data points (default 4)}

\item{test_density}{Minimum density of cases in each group to report p-value (default 3)}

\item{ints}{Number of x intervals over which to measure significance}

\item{quiet}{Silence all text outputs}

\item{x}{The independent variable; is continuous, e.g. time.}

\item{y}{The dependent variable; is continuous, e.g. temperature.}
}
\description{
Test for a significant difference in two groups at imputed intervals
}
\details{
The data object needs to be organized with each observation as a
row, and have a column that identifies the case, patient, animal, etc, and
columns with the continuous x and y variables (row with x = NA will be
removed). If there are multiple groups in the data, you can filter to
the single group of interest with the category and group arguments. Otherwise
it assumes the entire dataset is the single population.
}
\examples{
result <- trendyspliner(data = ChickWeight, xvar = 'Time',
             yvar = 'weight', category = 'Diet',
             cases = 'Chick', groups = c(1,2))
result$pval

}
