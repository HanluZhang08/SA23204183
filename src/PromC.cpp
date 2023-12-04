#include <Rcpp.h>
using namespace Rcpp;
//' @title A probability matrix calculator
//' @description A probability matrix calculator
//' @param adjM the adjacency matrix
//' @param d the damping factor
//' @return a probability matrix
//' @examples
//' \dontrun{
//' A = matrix(c(0,1,0,1,0,0,0,1,1), ncol = 3)
//' P = PromC(A, 0.85)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix PromC(const NumericMatrix& adjM, double d) {
  NumericVector cs = Rcpp::colSums(adjM);
  cs[cs == 0] = 1;
  int n = adjM.nrow();
  double delta = (1 - d) / n;
  NumericMatrix P(n, n);
  for (int i = 0; i < n; i++) {
    P(i, _) = delta + d * adjM(i, _) / cs;
  }
  return P;
}
