#' Sample a track from the posterior of a fitted model
#'
#' @param track_fit The output of fit_track()
#' @param number_of_samples The number of posterior samples
#' @param method Either "joint_precision" or "posterior_nngp". If "joint_precision":
#'   posterior samples are obtained by inverting the joint precision matrix
#'   and sampling directly from the multivariate normal distribution with mean
#'   equal to the point estimates.
#' 
#'   If "posterior_nngp": correlation parameters are sampled from the point estimates
#'   and fixed effects precision matrix. The marginal mean and standard deviation 
#'   of each random effect is set to the point estimate and standard error from 
#'   the model fit. Then the above quantities are fed into the NNGP model for simulation.
#' @param n_cores Number of parallel cores to use.
#'
#' @return A list containing the posterior sampled tracks
#'
#' @export
sample_posterior_tracks<- function(
  track_fit,
  number_of_samples = 10,
  method = "joint_precision",
  n_cores = 1
) {
  if( method == "joint_precision" ) {
    coordinate_mean<- c(
      sf::st_coordinates(
        track_fit$track
      )
    )
    sdr<- track_fit$sdr
    coordinate_covariance<- solve(sdr$jointPrecision)
    coordinate_covariance<- coordinate_covariance[
      rownames(coordinate_covariance) == "true_coordinate",
      colnames(coordinate_covariance) == "true_coordinate"
    ]
    samples<- t(
      mvtnorm::rmvnorm(
        n = number_of_samples,
        mean = coordinate_mean,
        sigma = coordinate_covariance
      )
    )
    samples<- apply(
      samples,
      MARGIN = 2,
      function(x) {
        sample_track<- matrix(
          x,
          ncol = 2,
          byrow = FALSE
        )
        sample_track<- apply(
          sample_track,
          MARGIN = 1,
          sf::st_point,
          simplify = FALSE
        )
        sample_track<- do.call(sf::st_sfc, sample_track)
        sample_track<- sf::st_sf(
          time = track_fit$track$time,
          geometry = sample_track
        )
        sf::st_crs(sample_track)<- sf::st_crs(track_fit$track)
        return(sample_track)
      },
      simplify = FALSE
    )
  } else if( method == "posterior_nngp" ) {
    sdr<- track_fit$sdr
    par_samples<- t(
      mvtnorm::rmvnorm(
        n = number_of_samples,
        mean = sdr$par.fixed,
        sigma = sdr$cov.fixed
      )
    )
    coordinate_mean<- as.list(sdr, "Est")$true_coordinate
    coordinate_sd<- as.list(sdr, "Std")$true_coordinate
    data<- track_fit$TMB_input$data
    para<- as.list(sdr, "Est")
    map<- track_fit$TMB_input$map
    random<- track_fit$TMB_input$random
    DLL<- track_fit$TMB_input$DLL
    data$posterior_simulation<- TRUE
    data$posterior_coordinate_mean<- coordinate_mean
    data$posterior_coordinate_sd<- coordinate_sd
    obj<- TMB::MakeADFun(
      data = data,
      para = para,
      map = map,
      random = random,
      DLL = DLL,
      silent = TRUE
    )

    # samples<- parallel::mclapply(
    #   seq(ncol(par_samples)),
    #   function(i) {
    #     x<- par_samples[, i]
    #     obj$fn(x)
    #     sample_track<- obj$simulate()$true_coordinate
    #     sample_track<- apply(
    #       sample_track,
    #       MARGIN = 1,
    #       sf::st_point,
    #       simplify = FALSE
    #     )
    #     sample_track<- do.call(sf::st_sfc, sample_track)
    #     sample_track<- sf::st_sf(
    #       time = track_fit$track$time,
    #       geometry = sample_track
    #     )
    #     sf::st_crs(sample_track)<- sf::st_crs(track_fit$track)
    #     return(sample_track)
    #   }
    # )

    samples<- apply(
      par_samples,
      MARGIN = 2,
      function(x) {
        obj$fn(x)
        sample_track<- obj$simulate()$true_coordinate
        sample_track<- apply(
          sample_track,
          MARGIN = 1,
          sf::st_point,
          simplify = FALSE
        )
        sample_track<- do.call(sf::st_sfc, sample_track)
        sample_track<- sf::st_sf(
          time = track_fit$track$time,
          geometry = sample_track
        )
        sf::st_crs(sample_track)<- sf::st_crs(track_fit$track)
        return(sample_track)
      },
      simplify = FALSE
    )
  } else {
    stop("Method must be one of 'joint_precision' or 'posterior_nngp'.")
  }
  return(samples)
}

