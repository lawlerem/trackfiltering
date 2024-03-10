#' Fit a Gaussian process model to an sf data.frame of location pings.
#'
#' @param pings An sf data.frame with  point geometries, a "time" column, and a "quality_class" column
#' @param track_time Optional vector of times to use for modeling the true locations. Defaults to observed ping times.
#' @param time_units Units to use when computing time differences
#' @param covariance_code 0 = exponential, 1 = gaussian, 2 = matern, 3 = matern32
#' @param number_of_ping_parents Number of observation pings to use as parents
#' @param number_of_track_parents Number of directly prior track locations to use as parents
#' @param track_robustness Tuning parameters for robustifying process
#' @param track_robust_code Robustification type: 0 = ML, 1 = arctan, 2 = smooth semi-Huber
#' @param ping_robustness Tuning parameters for robustifying observations
#' @param ping_robust_code Robustification type: 0 = ML, 1 = arctan, 2 = smooth semi-Huber
#' @param starting_move_sd Starting values for std. dev. of movement velocity components
#' @param starting_move_range Starting values for movement range parameter
#' @param starting_ping_sd Starting values for observation error std. dev. (ARGOS class 3) (measured in m)
#' @param starting_ping_cor Starting values for observation error ellipse correlation
#' @param fix_move_sd Should these parameters be held fixed in the estimation?
#' @param fix_move_range Should these parameters be held fixed in the estimation?
#' @param fix_ping_sd Should these parameters be held fixed in the estimation?
#' @param fix_ping_cor Should these parameters be held fixed in the estimation?
#' @param k Length of sliding window for median filter on location estimates.
#' @param silent Should printing of tracing information be prevented
#' @param ... Options passed to nlminb
#'
#' @return A list with pings, estimated track, TMB object, optimization object, and sdreport object.
#'
#' @export
fit_track<- function(
  pings,
  track_time,
  time_units = "days",
  covariance_code = 3,
  number_of_ping_parents = 5,
  number_of_track_parents = 5,
  track_robustness = 6,
  track_robust_code = 0,
  ping_robustness = 6,
  ping_robust_code = 2,
  starting_move_sd = c(12, 12),
  starting_move_range = c(3, 3),
  starting_ping_sd = c(500, 500),
  starting_ping_cor = 0,
  fix_move_sd = FALSE,
  fix_move_range = TRUE,
  fix_ping_sd = FALSE,
  fix_ping_cor = TRUE,
  k = 5,
  silent = TRUE,
  ...
) {
  ping_time<- pings$time
  if( missing(track_time) ) {
    track_time<- ping_time
  } else {}
  track_time<- sort(unique(track_time))

  TMB_ping_time<- ping_time - track_time[[1]]
  units(TMB_ping_time)<- time_units
  TMB_ping_time<- as.numeric(TMB_ping_time)

  TMB_track_time<- track_time - track_time[[1]]
  units(TMB_track_time)<- time_units
  TMB_track_time<- as.numeric(TMB_track_time)

  end_times<- c(
    min(c(TMB_ping_time, TMB_track_time)),
    max(c(TMB_ping_time, TMB_track_time))
  )

  ping_parents<- lapply(
    TMB_ping_time,
    function(t) {
      d<- abs(t - TMB_track_time)
      if( any(d < .Machine$double.eps) ) {
        return( which(d < .Machine$double.eps)[1] )
      } else {}
      return(order(d)[seq(number_of_ping_parents)])
    }
  )
  track_parents<- lapply(
    seq_along(track_time),
    function(t) {
      track_parents<- seq(t - number_of_track_parents, t - 1)
      track_parents<- track_parents[track_parents > 0]
      return(sort(track_parents, decreasing = TRUE))
    }
  )
  track_order<- seq_along(track_time)

  track_pairs<- lapply(
    seq_along(track_parents),
    function(i) {
      parents<- track_parents[[i]]
      pairs<- matrix(i, ncol = 2)
      if( length(parents) > 0 ) {
        pairs<- rbind(
          pairs,
          t(combn(sort(c(i, parents)), 2))
        )
      } else {}
      pairs<- as.data.frame(pairs)
      colnames(pairs)<- c("i", "j")
      return(pairs)
    }
  )
  ping_pairs<- lapply(
    seq_along(ping_parents),
    function(i) {
      parents<- ping_parents[[i]]
      if( length(parents) < 2 ) {
        return(NULL)
      } else {}
      pairs<- t(combn(sort(parents), 2))
      pairs<- as.data.frame(pairs)
      colnames(pairs)<- c("i", "j")
      return(pairs)
    }
  )
  covariance_pairs<- do.call(
    rbind,
    c(track_pairs, ping_pairs)
  )
  covariance_pairs<- unique(covariance_pairs)


  default_class_starting_values<- units::set_units(starting_ping_sd, "m")
  default_class_starting_values<- units::set_units(
    default_class_starting_values,
    units(sf::st_distance(pings[1, ])),
    mode = "standard"
  )
  default_class_starting_values<- log(exp(units::drop_units(default_class_starting_values)) - 1)
  default_class_starting_values<- c(
    default_class_starting_values[[1]],
    starting_ping_cor,
    default_class_starting_values[[2]]
  )

  observed_pings<- list(
    coordinate = sf::st_coordinates(pings),
    index = seq_along(ping_time) - 1,
    quality_class = as.integer(pings$quality_class) - 1
  )
  data<- list(
    end_times = end_times,
    ping_times = matrix(TMB_ping_time, ncol = 1), # Times used to evaluate location spline
    ping_parents = lapply(ping_parents, `+`, -1), # List of parent nodes for each spline x
    track_times = matrix(TMB_track_time, ncol = 1), # Times used for spline knots
    track_parents = lapply(track_parents, `+`, -1), # List of parents nodes for each node
    track_order = track_order - 1, # Evaluation order for node likelihoods
    covariance_code = covariance_code,
    covariance_pairs = as.matrix(covariance_pairs) - 1, # List of all node pairs whose covariance needs to be computed
    track_robustness = track_robustness,
    track_robust_code = track_robust_code,
    quality_class_adjustment = as.matrix(loc_class_K[, 2:3]),
    observed_pings = observed_pings,
    ping_robustness = ping_robustness,
    ping_robust_code = 0,
    posterior_simulation = FALSE,
    posterior_coordinate_mean = 0.0 * observed_pings$coordinate,
    posterior_coordinate_sd = 0.0 * observed_pings$coordinate
  )
  para<- list(
    end_points = sf::st_coordinates(pings)[
      c(
        which.min(pings$time),
        which.max(pings$time)
      ),

    ],
    track_coordinates = matrix(0, nrow = length(track_time), ncol = 2),
    working_covariance_parameters = rbind(
      # c(-2.5, 4.2, 0.5),
      # c(-2.5, 4.5, 0.5)
      c(log(exp(starting_move_sd[[1]]) - 1), log(exp(starting_move_range[[1]]) - 1), log(exp(0.5) - 1)),
      c(log(exp(starting_move_sd[[2]]) - 1), log(exp(starting_move_range[[2]]) - 1), log(exp(0.5) - 1))
    ),
    working_ping_parameters = default_class_starting_values
  )

  map<- list()
  map$working_covariance_parameters<- factor(rbind(
    c(ifelse(fix_move_sd, NA, 1), ifelse(fix_move_range, NA, 3), NA),
    c(ifelse(fix_move_sd, NA, 2), ifelse(fix_move_range, NA, 4), NA)
  ))
  if( covariance_code == 2 ) {
    map$working_covariance_parameters<- factor(rbind(
      c(ifelse(fix_move_sd, NA, 1), ifelse(fix_move_range, NA, 3), 5),
      c(ifelse(fix_move_sd, NA, 2), ifelse(fix_move_range, NA, 4), 6)
    ))
  } else {}
  map$working_ping_parameters<- factor(c(
    ifelse(fix_ping_sd, NA, 1),
    ifelse(fix_ping_cor, NA, 2),
    ifelse(fix_ping_sd, NA, 3)
  ))

  # Fit non-robust ML to get starting values
  obj<- TMB::MakeADFun(
    data = data,
    para = para,
    map = map,
    random = "track_coordinates",
    DLL = "trackfiltering_TMB",
    silent = silent
  )
  opt<- nlminb(
    obj$par,
    obj$fn,
    obj$gr,
    ...
  )
  sdr<- sdreport(
    obj,
    opt$par,
    ignore.parm.uncertainty = TRUE,
    getReportCovariance = FALSE,
    skip.delta.method = TRUE
  )
  sdr_est<- as.list(sdr, "Est", report = FALSE)
  if( ping_robust_code != 0 ) {
    # Use ML Estimates as starting values and shrink towards robust estimate
    data$ping_robust_code<- ping_robust_code
    r_seq<- c(100, 50, 25, 10, 7, 5, 4, 3, 2, 1) * ping_robustness
    r_seq<- r_seq[r_seq < -min(obj$report()$ping_loglikelihoods, na.rm = TRUE)]
    good_data<- data
    good_obj<- obj
    good_opt<- opt
    for( r in r_seq ) {
      track_estimate<- sdr_est$track_coordinates
      track_estimate[-c(seq(floor(k/2)), seq(floor(k/2)) + nrow(track_estimate) - floor(k/2)), ]<- zoo::rollmedian(
        track_estimate,
        k = k,
        align = "center",
        fill = NULL
      )
      data$ping_robustness<- r
      para<- list(
        end_points = sdr_est$end_points,
        track_coordinates = track_estimate,
        working_covariance_parameters = sdr_est$working_covariance_parameters,
        working_ping_parameters = sdr_est$working_ping_parameters
      )
      obj<- TMB::MakeADFun(
        data = data,
        para = para,
        map = map,
        random = "track_coordinates",
        DLL = "trackfiltering_TMB",
        silent = silent
      )
      opt<- nlminb(
        obj$par,
        obj$fn,
        obj$gr,
        ...
      )
      if( opt$convergence == 0 ) {
        good_data<- data
        good_obj<- obj
        good_opt<- opt
      } else {
        data<- good_data
        obj<- good_obj
        opt<- good_opt
        break
      }
      sdr<- sdreport(
        obj,
        opt$par,
        ignore.parm.uncertainty = TRUE,
        getReportCovariance = FALSE,
        skip.delta.method = TRUE
      )
      sdr_est<- as.list(sdr, "Est", report = FALSE)
    }
  } else {}
  sdr<- sdreport(
    obj,
    opt$par,
    getJointPrecision = TRUE
  )
  sdr_est<- c(
    as.list(sdr, "Est", report = FALSE),
    as.list(sdr, "Est", report = TRUE)
  )
  sdr_se<- c(
    as.list(sdr, "Std", report = FALSE),
    as.list(sdr, "Std", report = TRUE)
  )

  track<- apply(
    sdr_est$track_coordinates,
    MARGIN = 1,
    sf::st_point,
    simplify = FALSE
  )
  track<- do.call(sf::st_sfc, track)
  track<- st_sf(
    time = track_time,
    x_se = sdr_se$track_coordinates[, 1],
    y_se = sdr_se$track_coordinates[, 2],
    geometry = track
  )
  sf::st_crs(track)<- sf::st_crs(pings)

  # Add location estimates and standard errors to pings
  ping_estimates<- apply(
    obj$report()$ping_coordinates,
    MARGIN = 1,
    sf::st_point,
    simplify = FALSE
  )
  ping_estimates<- do.call(sf::st_sfc, ping_estimates)
  ping_estimates<- st_sf(
    # x_se = sdr_se$ping_coordinates[, 1],
    # y_se = sdr_se$ping_coordinates[, 2],
    estimate = ping_estimates
  )
  pings<- cbind(pings, ping_estimates)

  return(
    list(
      pings = pings,
      track = track,
      obj = obj,
      opt = opt,
      sdr = sdr,
      TMB_input = list(
        data = data,
        para = para,
        map = map,
        random = "track_coordinates",
        DLL = "trackfiltering_TMB"
      )
    )
  )
}
