// Class that computes covariances and covariance matrices
//
// Evaluate x(...) to compute covariances from distances
template<class Type>
class covariance {
  private:
    vector<Type> parameters;
    int covariance_code; // Which covariance function to use?
    template<typename T> T distance(T t1, T t2);
    template<typename T> T distance(vector<T> t1, vector<T> t2);

  public:
    // Constructor
    covariance(const vector<Type>& parameters, const int& covariance_code) :
      parameters{parameters}, covariance_code{covariance_code} {};
    covariance() : parameters{vector<Type>()}, covariance_code(0) {};

    // Compute covariances
    template<typename T> T operator() (T d);
    template<typename T> T operator() (T t1, T t2);
    template<typename T> matrix<T> operator() (vector<T> t);
    template<typename T> T operator() (vector<T> t1, vector<T> t2);
};


template<class Type>
template<typename T>
T covariance<Type>::distance(T t1, T t2) {
  return sqrt(pow(t1 - t2, 2));
}
template<class Type>
template<typename T>
T covariance<Type>::distance(vector<T> t1, vector<T> t2) {
  Type d2 = 0.0;
  for(int i = 0; i < t1.size(); i++) {
    d2 += pow(t1(i) - t2(i), 2);
  }
  return sqrt(d2);
}


template<class Type>
template<typename T>
T covariance<Type>::operator() (T d) {
  switch(covariance_code) {
    // Exponential [sd, range]
    case 0 : return (T) pow(parameters(0), 2) * parameters(1) * exp( -d / (T)parameters(1) );
    // Gaussian [marg_sd, range]
    case 1 : return (T) pow(parameters(0) * parameters(1), 2) * exp( -pow(d / (T)parameters(1), 2) );
    // Matern [sd, range, nu]
    case 2 : return (T) pow(parameters(0), 2) * pow(parameters(1), 2 * parameters(2)) * matern(d, parameters(1), parameters(2));
    // Matern32 [sd, range]
    case 3 : return (T) pow(parameters(0), 2) * pow(parameters(1), 3) * (1 + sqrt(3.0) * d / (T)parameters(1)) * exp( -sqrt(3.0) * d / (T)parameters(1) );
    // Exponential [sd, range]
    default : return (T) pow(parameters(0), 2) * parameters(1) * exp( -d / (T)parameters(1) );
  }
}

template<class Type>
template<typename T>
T covariance<Type>::operator() (T t1, T t2) {
  T d = distance(t1, t2);
  return operator()(d);
}

template<class Type>
template<typename T>
T covariance<Type>::operator() (vector<T> t1, vector<T> t2) {
  T d = distance(t1, t2);
  return operator()(d);
}

template<class Type>
template<typename T>
matrix<T> covariance<Type>::operator() (vector<T> t) {
  matrix<T> ans(t.size(), t.size());
  for(int i = 0; i < t.size(); i++) {
    for(int j = 0; j < t.size(); j++) {
      ans(i, j) = operator()(t(i), t(j));
      if( i != j ) {
        ans(i, j) *= 0.999;
      }
    }
  }
  return ans;
}
