template <class Type>
Type robustify(Type x, Type c, int code) {
	// c = tuning constant
	// code: 0 = ML, 1 = loglog, 2 = ssh
	switch(code){
		case 1: return log( Type(1.0) + exp(x + c) ) - log( Type(1.0) + exp(c) ); // log-logistic
		case 2: // rho = smoothed semi-Huber (ssh)
			return CppAD::CondExpGt(
          x,
          -c,
          x,
				  c * log(
            (x + c) / c + sqrt(
              1.0 + pow( (x + c) / c, 2)
            )
          ) - c
        );
		default:
      // ML; rho = identity
			return x;
	}
  return(0.0);
}

template<class Type>
Type robustify_weight_score(Type x, Type c, int code) {
  // c = tuning constant
  // code: 0 = ML, 1 = loglog, 2 = ssh
  switch(code) {
    case 1: return exp(x + c) / (1.0 + exp(x + c)); // log-logistic
    case 2: // rho = smoothed semi-Huber (ssh)
      return CppAD::CondExpGt(
        x,
        -c,
        Type(1.0),
        pow(
          1.0 + pow((x + c) / c, 2),
          -0.5
        )
      );
    default:
      // ML; rho = identity
      return Type(1.0);
  }
}

template<class Type>
vector<Type> robustify(vector<Type> x, Type c, int code) {
  vector<Type> ans(x.size());
  for(int i = 0; i < x.size(); i++) {
    ans(i) = robustify(x(i), c, code);
  }
  return ans;
}