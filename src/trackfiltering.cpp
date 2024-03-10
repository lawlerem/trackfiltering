#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

/*
Get the "to" and "from" vertices from a graph node
*/

Eigen::VectorXi from_list(Rcpp::List edge_list) {
  return Rcpp::as<Eigen::VectorXi>(
    edge_list["from"]
  );
}
Eigen::VectorXi to_list(Rcpp::List edge_list) {
  return Rcpp::as<Eigen::VectorXi>(
    edge_list["to"]
  );
}
