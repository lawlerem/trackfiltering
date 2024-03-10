#' Simulate a track with pings
#'
#' @param track_time A vector giving the timestamps of track locations
#' @param ping_time A vector giving the timestamps of location pings
#' @param endpoint The (expected) starting and ending points
#' @param covariance_code 0 = exponential, 1 = gaussian, 2 = matern, 3 = matern32
#' @param covariance_parameters A 2x3 matrix. Columns: sd, range, smoothness
#' @param ping_parameters A 3-vector giving the reference observation error
#' @param quality_class_probabilities A vector giving relative frequencies of quality classes
#' @param number_of_ping_parents Number of observation pings to use as parents
#' @param number_of_track_parents Number of directly prior track locations to use as parents
#' @param seed Random number seed
#'
#' @return A list
#'
#' @export
simulate_track<- function(
  track_time = seq(0, 1, length.out = 101),
  ping_time = seq(0, 1, length.out = 25),
  endpoint = rbind(
    c(0, 0),
    c(0, 0)
  ),
  covariance_code = 1,
  covariance_parameters = rbind(
    c(1, 0.2, Inf),
    c(1, 0.2, Inf)
  ),
  ping_parameters = c(0.1, 0.0, 0.1),
  quality_class_probabilities = c(
    0.1, 0.1, 0.2, 0.4, 0.2, 0.0, 0.0
  ),
  number_of_ping_parents = 5,
  number_of_track_parents = 5,
  seed
) {
  if( !missing(seed) ) set.seed(seed)
  track_time<- sort(unique(c(track_time, ping_time)))
  if( length(quality_class_probabilities) < 7 ) {
    quality_class_probabilities<- c(
      quality_class_probabilities,
      rep(0, 7 - length(quality_class_probabilities))
    )
  } else {}
  if( length(quality_class_probabilities) > 7 ) {
    warning("Using only the first 7 entries of quality_class_probabilities.")
    quality_class_probabilities<- head(quality_class_probabilities, 7)
  } else {}
  quality_class_probabilities<- quality_class_probabilities / sum(quality_class_probabilities)

  observed_pings<- list(
    coordinate = matrix(0, nrow = length(ping_time), ncol = 2),
    index = match(sort(ping_time), track_time) - 1,
    quality_class = sample(
      seq(1, 7) - 1,
      size = length(ping_time),
      replace = TRUE,
      prob = quality_class_probabilities
    )
  )
  track_parents<- lapply(
    seq_along(track_time),
    function(t) {
      ping_parents<- observed_pings$index + 1
      ping_parents<- unique(ping_parents[ping_parents < t])
      ping_parents<- tail(ping_parents, number_of_ping_parents)
      track_parents<- seq(t - number_of_track_parents, t - 1)
      track_parents<- track_parents[track_parents > 0]
      p<- unique(c(ping_parents, track_parents))
      return(sort(p, decreasing = TRUE))
    }
  )

  data<- list(
    track_time = track_time,
    track_parents = lapply(track_parents, `+`, -1),
    covariance_code = covariance_code,
    quality_class_adjustment = as.matrix(loc_class_K[, 2:3]),
    observed_pings = observed_pings,
    penalise_range = FALSE,
    range_hyperparameters = numeric(2)
  )
  para<- list(
    endpoint = endpoint,
    true_coordinate = matrix(0, nrow = length(track_time), ncol = 2),
    working_covariance_parameters = log(covariance_parameters),
    working_ping_parameters = c(
      log(ping_parameters[[1]]),
      qlogis(0.5 * (ping_parameters[[2]] + 1)),
      log(ping_parameters[[3]])
    )
  )
  obj<- TMB::MakeADFun(
    data = data,
    para = para,
    random = "true_coordinate",
    DLL = "trackfiltering_TMB",
    silent = TRUE
  )
  sim<- obj$simulate()

  track<- apply(
    sim$true_coordinate,
    MARGIN = 1,
    sf::st_point,
    simplify = FALSE
  )
  track<- do.call(sf::st_sfc, track)
  track<- st_sf(
    time = track_time,
    geometry = track
  )

  pings<- apply(
    sim$sim_pings,
    MARGIN = 1,
    sf::st_point,
    simplify = FALSE
  )
  pings<- do.call(sf::st_sfc, pings)
  pings<- st_sf(
    index = observed_pings$index + 1,
    quality_class = observed_pings$quality_class,
    time = track_time[observed_pings$index + 1],
    geometry = pings
  )
  pings$quality_class<- factor(
    c("G", 3, 2, 1, 0, "A", "B")[pings$quality_class + 1],
    levels = c("G", 3, 2, 1, 0, "A", "B"),
    ordered = TRUE
  )
  return(
    list(
      pings = pings,
      track = track
    )
  )
}
