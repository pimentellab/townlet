#' Example dataset: testdata
#'
#' A village donor composition test data set.
#'
#' @format A data frame with 437 rows and 5 variables:
#' \describe{
#'   \item{donor}{character, donor name}
#'   \item{time}{numeric, day sampled}
#'   \item{replicate}{integer, village replicate number}
#'   \item{sex}{character, biological sex that can be included as optional covariate}
#'   \item{representation}{numeric, donor proportion where proportions grouped by replicate and sample time point will sum to 1}
#' }
#' @source Simulated
"testdata"
