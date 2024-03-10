template<class Type>
class spline {
  private:
    matrix<Type> x;
    vvint<Type> parents;
    matrix<Type> means;

    matrix<Type> nodes;
    vvint<Type> node_parents;
    vector<int> node_order;
    matrix<Type> node_means;
    matrix<Type> node_values;

    vector<covariance<Type> > covariance_functions;
    matrix<int> covariance_pairs;
    vector<matrix<Type> > covariance_matrices;

    int n_splines;
    int n_x;
    int n_nodes;
  public:
    matrix<Type> values;

    spline(
      matrix<Type> x,
      vvint<Type> parents,
      matrix<Type> means,
      matrix<Type> nodes,
      vvint<Type> node_parents,
      vector<int> node_order,
      matrix<Type> node_means,
      matrix<Type> node_values,
      vector<covariance<Type> > covariance_functions,
      matrix<int> covariance_pairs
    ) : x(x),
        parents(parents),
        means(means),
        nodes(nodes),
        node_parents(node_parents),
        node_order(node_order),
        node_means(node_means),
        node_values(node_values),
        covariance_functions(covariance_functions),
        covariance_pairs(covariance_pairs),
        n_splines(node_values.cols()),
        n_x(x.rows()),
        n_nodes(node_values.rows()) {
      // 1.) Compute covariance matrices, only filling in entries we need
      covariance_matrices.resize(node_values.cols());
      for(int s = 0; s < n_splines; s++) {
        matrix<Type> covmat(n_nodes, n_nodes);
        covmat.setZero();
        for(int pair = 0; pair < covariance_pairs.rows(); pair++) {
          int i = covariance_pairs(pair, 0);
          int j = covariance_pairs(pair, 1);
          covmat(i, j) = covariance_functions(s)(vector<Type>(nodes.row(i)), vector<Type>(nodes.row(j)));
          covmat(j, i) = covmat(i, j);
        }
        covariance_matrices(s) = covmat;
      }

      // 2.) Fill in x values
      means.resize(x.rows(), node_values.cols());
      values.resize(x.rows(), node_values.cols());
      for(int s = 0; s < n_splines; s++) {
        for(int i = 0; i < n_x; i++) {
          if( parents(i).size() == 1 ) {
            values(i, s) = node_values(parents(i)(0), s);
            continue;
          } else {}

          // Create covariance matrix needed for conditional mean
          matrix<Type> covmat(1 + parents(i).size(), 1 + parents(i).size());
          for(int row = 0; row < parents(i).size(); row++) {
            for(int col = 0; col < parents(i).size(); col++) {
              covmat(row + 1, col + 1) = covariance_matrices(s)(parents(i)(row), parents(i)(col));
            }
            covmat(row + 1, 0) = covariance_functions(s)(vector<Type>(nodes.row(parents(i)(row))), vector<Type>(x.row(i)));
            covmat(0, row + 1) = covmat(row + 1, 0);
          }
          covmat(0, 0) = covariance_functions(s)(vector<Type>(x.row(i)), vector<Type>(x.row(i)));

          vector<Type> these_values(1 + parents(i).size());
          vector<Type> these_means(1 + parents(i).size());
          these_values(0) = 0.0;
          these_means(0) = means(i, s);
          for(int p = 0; p < parents(i).size(); p++) {
            these_values(p + 1) = node_values(parents(i)(p), s);
            these_means(p + 1) = node_means(parents(i)(p), s);
          }

          conditional_normal<Type> cmvn(covmat, parents(i).size());
          values(i, s) = cmvn.conditional_mean(these_values, these_means)(0);
        }
      }
    }

  Type loglikelihood(Type robustness, int robust_code);
  spline<Type> simulate();
};


template<class Type>
Type spline<Type>::loglikelihood(Type robustness, int robust_code) {
  Type ll = 0.0;
  for(int s = 0; s < n_splines; s++) {
    for(int n = 0; n < node_order.size(); n++) {
      int i = node_order(n);
      // 1.) Get covariance matrix, node means, and node values
      // 2.) Compute conditional normal
      // 3.) Evaluate conditional normal density
      matrix<Type> covmat(1 + node_parents(i).size(), 1 + node_parents(i).size());
      for(int row = 0; row < node_parents(i).size(); row++) {
        for(int col = 0; col <  node_parents(i).size(); col++) {
          covmat(row + 1, col + 1) = covariance_matrices(s)(node_parents(i)(row), node_parents(i)(col));
        }
        covmat(row + 1, 0) = covariance_matrices(s)(node_parents(i)(row), i);
        covmat(0, row + 1) = covmat(row + 1, 0);
      }
      covmat(0, 0) = covariance_matrices(s)(i, i);

      vector<Type> these_values(1 + node_parents(i).size());
      vector<Type> these_means(1 + node_parents(i).size());
      these_values(0) = node_values(i, s);
      these_means(0) = node_means(i, s);
      for(int p = 0; p < node_parents(i).size(); p++) {
        these_values(p + 1) = node_values(node_parents(i)(p), s);
        these_means(p + 1) = node_means(node_parents(i)(p), s);
      }

      conditional_normal<Type> cmvn(covmat, node_parents(i).size());
      ll += robustify(cmvn.loglikelihood(these_values, these_means), robustness, robust_code);
    }
  }
  return ll;
}