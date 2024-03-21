#define TMB_LIB_INIT R_init_trackfiltering
#include <TMB.hpp>
using namespace density;

#include "include/utilities.hpp"
#include "include/robust.hpp"
#include "include/covariance.hpp"
#include "include/conditional_normal.hpp"
#include "include/spline.hpp"
#include "include/pings.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  PARAMETER_MATRIX(end_points);
  DATA_VECTOR(end_times);
  DATA_MATRIX(ping_times);
  DATA_STRUCT(ping_parents, vvint);
  DATA_MATRIX(track_times);
  DATA_STRUCT(track_parents, vvint);
  DATA_IVECTOR(track_order);
  PARAMETER_MATRIX(track_coordinates);

  matrix<Type> ping_means(ping_times.rows(), track_coordinates.cols());
  for(int i = 0; i < ping_times.rows(); i++) {
    Type percent_t = (ping_times(i, 0) - end_times(0)) / (end_times(1) - end_times(0));
    ping_means(i, 0) = (1.0 - percent_t) * end_points(0, 0) + percent_t * end_points(1, 0);
    ping_means(i, 1) = (1.0 - percent_t) * end_points(0, 1) + percent_t * end_points(1, 1);
  }
  matrix<Type> track_means(track_times.rows(), track_coordinates.cols());
  for(int i = 0; i < track_times.rows(); i++) {
    Type percent_t = (track_times(i, 0) - end_times(0)) / (end_times(1) - end_times(0));
    track_means(i, 0) = (1.0 - percent_t) * end_points(0, 0) + percent_t * end_points(1, 0);
    track_means(i, 1) = (1.0 - percent_t) * end_points(0, 1) + percent_t * end_points(1, 1);
  }

  DATA_INTEGER(covariance_code);
  PARAMETER_MATRIX(working_covariance_parameters); // [par,var], columns may have trailing NA
  matrix<Type> covariance_parameters = working_covariance_parameters;
  covariance_parameters.col(0) = log(1 + exp(vector<Type>(working_covariance_parameters.col(0)))); // Variance Rate
  covariance_parameters.col(1) = log(1 + exp(vector<Type>(working_covariance_parameters.col(1)))); // Range
  covariance_parameters.col(2) = log(1 + exp(vector<Type>(working_covariance_parameters.col(2)))); // Smoothness
  REPORT(covariance_parameters);
  ADREPORT(covariance_parameters);

  vector<covariance<Type> > covariance_function(covariance_parameters.cols());
  for(int v = 0; v < covariance_function.size(); v++) {
    covariance_function(v) = covariance<Type> {
      vector<Type>(covariance_parameters.row(v)),
      covariance_code
    };
  }
  DATA_IMATRIX(covariance_pairs);
  DATA_SCALAR(track_robustness);
  DATA_INTEGER(track_robust_code);

  spline<Type> true_track(
    ping_times, // x
    ping_parents, // parents
    ping_means, // means
    track_times, // nodes
    track_parents, // node_parents
    track_order, // node_order
    track_means, // node_means
    track_coordinates, // node_values
    covariance_function, // covariance_functions
    covariance_pairs // covariance_pairs
  );
  Type track_ll = true_track.loglikelihood(track_robustness, track_robust_code);
  matrix<Type> ping_coordinates = true_track.values;
  REPORT(ping_coordinates);
  // ADREPORT(ping_coordinates);

  // Location Ping Likelihood
  PARAMETER(working_ping_correlation);
  Type ping_correlation = 2 * invlogit(working_ping_correlation) - 1.0;
  REPORT(ping_correlation);
  ADREPORT(ping_correlation);

  PARAMETER_VECTOR(working_ping_scaling);
  vector<Type> ping_scaling = exp(working_ping_scaling);
  REPORT(ping_scaling);
  ADREPORT(ping_scaling);

  matrix<Type> ping_covariance(2, 2);
  ping_covariance << 1.0, ping_correlation, ping_correlation, 1.0;

  vector<MVNORM_t<Type> > quality_class_distribution(ping_scaling.size());
  vector<Type> quality_class_logdet(ping_scaling.size());
  for(int quality_class = 0; quality_class < quality_class_distribution.size(); quality_class++) {
    matrix<Type> diagonal_adjustment(2, 2);
    diagonal_adjustment << ping_scaling(quality_class), 0.0, 0.0, ping_scaling(quality_class);
    matrix<Type> quality_class_covariance = diagonal_adjustment * ping_covariance * diagonal_adjustment;
    quality_class_distribution(quality_class) = MVNORM_t<Type>(quality_class_covariance);
    quality_class_logdet(quality_class) = atomic::logdet(diagonal_adjustment);
  }

  DATA_STRUCT(observed_pings, pings);
  DATA_SCALAR(ping_robustness);
  DATA_INTEGER(ping_robust_code);
  vector<Type> ping_weights(observed_pings.size());
  vector<Type> ping_loglikelihoods(observed_pings.size());

  Type ping_ll = observed_pings.loglikelihood(
    ping_coordinates,
    quality_class_distribution,
    quality_class_logdet,
    ping_loglikelihoods,
    ping_weights,
    ping_robustness,
    ping_robust_code
  );
  REPORT(ping_loglikelihoods);
  REPORT(ping_weights);
  Type ll = track_ll + ping_ll;

  DATA_INTEGER(posterior_simulation);
  DATA_MATRIX(posterior_coordinate_mean);
  DATA_MATRIX(posterior_coordinate_sd);

  // SIMULATE{
  //   if( posterior_simulation ) {
  //     true_coordinate = true_track.simulate_posterior(
  //       posterior_coordinate_mean,
  //       posterior_coordinate_sd
  //     );
  //   } else {
  //     true_coordinate = true_track.simulate();
  //   }
  //   matrix<Type> sim_pings = observed_pings.simulate(true_coordinate, quality_class_distribution);
  //   REPORT(true_coordinate);
  //   REPORT(sim_pings);
  // }

  return -1.0 * ll;
}
