% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk_permusplinectomy.R
\name{permuspliner.plot.permdistance}
\alias{permuspliner.plot.permdistance}
\title{Plot permuted distance distribution}
\usage{
permuspliner.plot.permdistance(data, xlabel = NULL)
}
\arguments{
\item{data}{The results object from the permuspliner function}

\item{xlabel}{Name for the x axis in the plot}
}
\description{
Compare the permuted distances to the true distance. Requires permuspliner run with the "retain_perm" option.
}
\examples{
distplot <- permuspliner.plot.permdistance(results, xlabel='Time')

# View the plot
distplot

# Save the plot as PNG file
ggsave(distplot, file = 'my_plot.png', dpi=300, height=4, width=4)

}
