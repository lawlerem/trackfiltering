template<class Type>
struct pings {
  matrix<Type> coordinate;
  vector<int> index;
  vector<int> quality_class;

  pings(SEXP r_list) :
    coordinate(asMatrix<Type>(VECTOR_ELT(r_list, 0))),
    index(asVector<int>(VECTOR_ELT(r_list, 1))),
    quality_class(asVector<int>(VECTOR_ELT(r_list, 2)))
    {};

  Type loglikelihood(
    const matrix<Type>& true_coordinate,
    vector<MVNORM_t<Type> >& quality_class_distribution,
    vector<Type>& quality_class_logdet,
    vector<Type>& loglikelihoods,
    vector<Type>& weights,
    Type robustness,
    int robust_code
  );
  matrix<Type> simulate(
    const matrix<Type>& true_coordinate,
    vector<MVNORM_t<Type> >& quality_class_distribution
  );
  int size() { return coordinate.rows(); };
};

template<class Type>
Type pings<Type>::loglikelihood(
    const matrix<Type>& true_coordinate,
    vector<MVNORM_t<Type> >& quality_class_distribution,
    vector<Type>& quality_class_logdet,
    vector<Type>& loglikelihoods,
    vector<Type>& weights,
    Type robustness,
    int robust_code
    ) {
  Type ans = 0.0;

  for(int i = 0; i < coordinate.rows(); i++) {
    loglikelihoods(i) = -1.0 * quality_class_distribution(quality_class(i))(
      vector<Type>(
        coordinate.row(i) - true_coordinate.row(index(i))
      )
    );
    weights(i) = robustify_weight_score(loglikelihoods(i) + quality_class_logdet(quality_class(i)), robustness, robust_code);
    ans += robustify(loglikelihoods(i) + quality_class_logdet(quality_class(i)), robustness, robust_code) - quality_class_logdet(quality_class(i));
  }

  return ans;
}

template<class Type>
matrix<Type> pings<Type>::simulate(
    const matrix<Type>& true_coordinate,
    vector<MVNORM_t<Type> >& quality_class_distribution
    ) {
  for(int i = 0; i < coordinate.rows(); i++) {
    coordinate.row(i) = quality_class_distribution(quality_class(i)).simulate() + vector<Type>(true_coordinate.row(index(i)));
  }

  return coordinate;
}
