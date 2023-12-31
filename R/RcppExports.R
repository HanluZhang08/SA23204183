# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title A probability matrix calculator
#' @description A probability matrix calculator
#' @param adjM the adjacency matrix
#' @param d the damping factor
#' @return a probability matrix
#' @examples
#' \dontrun{
#' A = matrix(c(0,1,0,1,0,0,0,1,1), ncol = 3)
#' P = PromC(A, 0.85)
#' }
#' @export
PromC <- function(adjM, d) {
    .Call('_SA23204183_PromC', PACKAGE = 'SA23204183', adjM, d)
}

