template<class Type>
class track {
  private:
    matrix<Type> endpoint;
    matrix<Type> coordinate_mean;
    matrix<Type> coordinate;
    vector<Type> time;
    vvint<Type> parents;
    vector<vector<Type> > joint_time;
    vector<matrix<Type> > joint_mean;
    vector<covariance<Type> > covariance_function;


  public:
    track(
      const matrix<Type>& endpoint,
      const matrix<Type>& coordinate,
      const vector<Type>& time,
      vvint<Type> parents,
      const vector<covariance<Type> >& covariance_function
    ) : endpoint(endpoint),
        coordinate(coordinate),
        time(time),
        parents(parents),
        covariance_function(covariance_function) {
      coordinate_mean.resizeLike(coordinate);
      for(int t = 0; t < time.size(); t++) {
        Type percent_duration = (time(t) - time(0)) / (time(time.size() - 1) - time(0));
        for(int c = 0; c < coordinate.cols(); c++) {
          coordinate_mean(t, c) = (1.0 - percent_duration) * endpoint(0, c) + percent_duration * endpoint(1, c);
        }
      }

      joint_time.resize(time.size());
      for(int t = 0; t < time.size(); t++) {
        joint_time(t).resize(parents(t).size() + 1);
        joint_time(t)(0) = time(t);
        for(int p = 0; p < parents(t).size(); p++) {
          joint_time(t)(p + 1) = time(parents(t)(p));
        }
      }

      joint_mean.resize(time.size());
      for(int t = 0; t < time.size(); t++) {
        joint_mean(t).resize(parents(t).size() + 1, coordinate.cols());
        for(int c = 0; c < joint_mean(t).cols(); c++) {
          joint_mean(t)(0, c) = coordinate_mean(t, c);
        }
        for(int p = 0; p < parents(t).size(); p++) {
          for(int c = 0; c < joint_mean(t).cols(); c++) {
            joint_mean(t)(p + 1, c) = coordinate_mean(parents(t)(p), c);
          }
        }
      }
    };

    Type loglikelihood();
    matrix<Type> simulate();
    matrix<Type> simulate_posterior(matrix<Type> posterior_mean, matrix<Type> posterior_sd);
    int size() { return coordinate.rows(); }
};


template<class Type>
Type track<Type>::loglikelihood() {
  Type ll = 0.0;
  for(int t = 0; t < coordinate.rows(); t++) {
    for(int c = 0; c < coordinate.cols(); c++) {
      matrix<Type> joint_covariance = covariance_function(c)(joint_time(t));
      conditional_normal<Type> conditional_distribution(joint_covariance, parents(t).size());
      vector<Type> x(parents(t).size() + 1);
      x(0) = coordinate(t, c);
      for(int p = 0; p < parents(t).size(); p++) {
        x(p + 1) = coordinate(parents(t)(p), c);
      }
      vector<Type> mu = joint_mean(t).col(c);
      ll += conditional_distribution.loglikelihood(x, mu);
    }
  }
  return ll;
}

template<class Type>
matrix<Type> track<Type>::simulate() {
  for(int t = 0; t < coordinate.rows(); t++) {
    for(int c = 0; c < coordinate.cols(); c++) {
      matrix<Type> joint_covariance = covariance_function(c)(joint_time(t));
      conditional_normal<Type> conditional_distribution(joint_covariance, parents(t).size());
      vector<Type> x(parents(t).size() + 1);
      x(0) = 0.0;
      for(int p = 0; p < parents(t).size(); p++) {
        x(p + 1) = coordinate(parents(t)(p), c);
      }
      vector<Type> mu = joint_mean(t).col(c);
      coordinate(t, c) = conditional_distribution.simulate(x, mu)(0);
    }
  }
  return coordinate;
}

template<class Type>
matrix<Type> track<Type>::simulate_posterior(matrix<Type> posterior_mean, matrix<Type> posterior_sd) {
  for(int t = 0; t < coordinate.rows(); t++) {
    for(int c = 0; c < coordinate.cols(); c++) {
      matrix<Type> prior_covariance = covariance_function(c)(joint_time(t));
      matrix<Type> sd_deflation = prior_covariance;
      for(int i = 0; i < sd_deflation.rows(); i++) {
        for(int j = 0; j < sd_deflation.cols(); j++) {
          if( i == j ) {
            if( i == 0 ) {
              sd_deflation(i, j) = posterior_sd(t, c) / sqrt(prior_covariance(i, j));
            } else {
              sd_deflation(i, j) = posterior_sd(parents(t)(i - 1), c) / sqrt(prior_covariance(i, j));
            }
          } else {
            sd_deflation(i, j) = 0.0;
          }
        }
      }
      matrix<Type> posterior_covariance = sd_deflation * prior_covariance * sd_deflation;
      conditional_normal<Type> conditional_distribution(posterior_covariance, parents(t).size());
      vector<Type> x(parents(t).size() + 1);
      vector<Type> mu(parents(t).size() + 1);
      x(0) = 0.0;
      mu(0) = posterior_mean(t, c);
      for(int p = 0; p < parents(t).size(); p++) {
        x(p + 1) = coordinate(parents(t)(p), c);
        mu(p + 1) = posterior_mean(parents(t)(p), c);
      }
      coordinate(t, c) = conditional_distribution.simulate(x, mu)(0);
    }
  }
  return coordinate;
}