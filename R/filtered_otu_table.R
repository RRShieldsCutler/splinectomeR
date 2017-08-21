#' OTU table from Yassour et al. microbiome study
#'
#' The 16S OTU table from Yassour et al. (Science Trans Med, 2016), generated with
#' QIIME v1.8.0 and picked with Greengenes (October 2012).
#' 
#' 
#' @format A data frame with 1101 samples (columns) and 3091 taxa (rows):
#' \describe{
#'   \item{samples}{Each baby ID in format E999999, plus the timepoint in months}
#'   \item{OTU_ID}{Relative abundance of taxa, summarized at each level of taxonomy}
#'   
#' }
#' @source \url{https://pubs.broadinstitute.org/diabimmune/antibiotics-cohort}
"filtered_otu_table"