#' Fit a model to a movement track
#'
#' @param pings
#'     An sf data frame with a POSIX or numeric column named "date",
#'     a factor column named "class", and point geometries. If class is an
#'     ordered factor, then the estimated observation error variances will
#'     increase with class.
#' @param time_mesh
#'     A vector of times used to discretize movement.
#' @param time_units
#'     Which time units should be used? Mostly matters for numerical stability.
#' @param independent_coordinates
#'     If true, then observation ellipses are assumed to have independent axes.
#' @param common_coordinate_correlation
#'     If true, then observation ellipses for different quality classes are
#'     assumed to have the same correlation structure.
#' @param ...
#'     Additional arguments, particularly those to pass to
#'     nnspline::create_lcspline, robustifyRTMB::robustly_optimize,
#'     and mcreportRTMB::mcreport.
#'
#' @return 
#'     A list with the following:
#'   * pings
#'       A copy of the provided location pings with an additional weights
#'       column giving the robust weight of each observation, an "est"
#'       column(s) giving the estimated locations, and an "se" column(s)
#'       giving the standard errors around the estimated locations.
#'   * error_matrices
#'       A list of estimated observation error covariance matrices for each
#'       quality class.
#'   * template_spline
#'       An nnspline object used to hold the nodes and node graph used for
#'       the model.
#'   * start_time
#'       The date used for t = 0.
#'   * robust_optimization
#'       The output of robustifyRTMB::robustly_optimize
#'   * mcreport
#'       The output of mcreportRTMB::mcreport for the spline parameters,
#'       node values, and coordinate means.
#'   * call
#'       The function call used.
#'
#' @export
fit_track <- function(
        pings,
        time_mesh = pings$date |> unique(),
        time_units = "days",
        ping_common_orientation = TRUE,
        ping_standard_orientation = FALSE,
        ping_common_shape = TRUE,
        ping_equal_shape = FALSE,
        ...
    ) {
    call <- match.call(expand.dots = TRUE) |> as.list()
    if( !("class" %in% (pings |> colnames())) ) {
        pings$class <- 1 |> rep(pings |> nrow())
    }
    if( !(pings$class |> is.factor()) ) pings$class <- pings$class |> as.factor()
    if( any( 0 < table(pings$class) & table(pings$class) <= 10) ) {
        warning(
            paste0(
                "One or more of the quality classes has 10 or ",
                "fewer observations, which may cause estimation issues."
            )
        )
    }
    
    arg_names<- nnspline::create_lcspline |> formals() |> names()
    spline_args<- call[names(call) %in% arg_names]
    spline_args$x<- time_mesh |>
        difftime(time_mesh[1], units = time_units) |>
        tail(-1)
    spline<- nnspline::create_lcspline |> do.call(spline_args)
    pings$spline_idx<- pings |> 
        _$date |> 
        cut(time_mesh, labels = FALSE, include.lowest = TRUE)

    arg_names<- robustifyRTMB::robustly_optimize |> formals() |> names()
    robopt_args <- call[(call |> names()) %in% arg_names]
    environment(nll) <- environment()
    robopt_args$func <- nll
    n_coords<- pings |> sf::st_coordinates() |> ncol()
    robopt_args$parameters <- list(
        qstretch = numeric(2),
        log_height = numeric(2),
        coordinates = matrix(
            pings |> sf::st_coordinates() |> apply(2, median),
            nrow = time_mesh |> length(),
            ncol = 2,
            byrow = TRUE
        ),
        working_ping_orientation = matrix(
            0,
            nrow = n_coords * (n_coords - 1) / 2,
            ncol = pings$class |> levels() |> length()
        ),
        working_ping_shape = matrix(
            0,
            nrow = n_coords - 1,
            ncol = pings$class |> levels() |> length(),
        ),
        working_ping_scale = pings$class |> levels() |> length() |> numeric()
    )
    map <- list(
        working_ping_orientation = matrix(
            robopt_args$parameters$working_ping_orientation |> seq_along(),
            nrow = robopt_args$parameters$working_ping_orientation |> nrow()
        ),
        working_ping_shape = matrix(
            robopt_args$parameters$working_ping_shape |> seq_along(),
            nrow = robopt_args$parameters$working_ping_shape |> nrow()
        ),
        working_ping_scale = robopt_args$parameters$working_ping_scale |> 
            seq_along()
    )
    missing_classes <- table(pings$class) == 0
    map$working_ping_orientation[, missing_classes] <- NA
    map$working_ping_shape[, missing_classes] <- NA
    robopt_args$parameters$working_ping_scale[missing_classes] <- -10
    map$working_ping_scale[missing_classes]<- NA
    if( ping_common_orientation ) {
        map$working_ping_orientation[] <- map$working_ping_orientation |>
            nrow() |> 
            seq()
    }
    if( ping_common_shape ) {
        map$working_ping_shape[] <- map$working_ping_shape |>
            nrow() |>
            seq()
    }
    if( ping_standard_orientation ) {
        map$working_ping_orientation[]<- NA
    }
    if( ping_equal_shape ) {
        map$working_ping_shape[]<- NA
    }

    map <- map |> lapply(as.factor)
    robopt_args$map <- map
    robopt_args$random <- "coordinates"
    robopt_args$smooth <- "coordinates"
    robopt_args$nodes <- time_mesh
    fit <- robustifyRTMB::robustly_optimize |> do.call(robopt_args)

    if( ((fit$sdr$cov.fixed |> diag()) <= 0) |>any() || 
        fit$sdr$cov.fixed |> is.nan() |> any() ) {
        # If any observation error component is very small it can mess up sdreport
        #   So I'll fix them at their estimated value and re-run sdreport.

        robopt_args$map$working_ping_diagonal <- 
            (fit$par$working_ping_diagonal <= -5) |> 
            ifelse(
                NA,
                fit$par$working_ping_diagonal |> seq_along()
            ) |>
            as.factor()
        robopt_args$map$robustness <- as.factor(NA)
        robopt_args$parameters <- fit$par
        robopt_args$silent <- TRUE
        fit$obj <- RTMB::MakeADFun |> do.call(robopt_args)
        fit$par <- fit$obj$env$parList()
        fit$opt$par <- fit$obj$par
        fit$sdr <- with(
            fit,
            obj |> RTMB::sdreport(opt$par, getJointPrecision = TRUE)
        )
    }

    arg_names<- mcreportRTMB::mcreport |> formals() |> names()
    mcreport_args <- call[(call |> names()) %in% arg_names]
    mcreport_args$obj <- fit$obj
    mcreport_args$sdr <- fit$sdr
    mc <- mcreportRTMB::mcreport |> do.call(mcreport_args)

    report <- fit$obj$report()
    sdr_est <- fit$sdr |> as.list("Est")
    sdr_se <- fit$sdr |> as.list("Std")

    pings <- cbind(
        weights = report$weights,
        est = sdr_est$coordinates[pings$spline_idx, ],
        se = sdr_se$coordinates[pings$spline_idx, ],
        pings
    )
    pings$spline_idx <- NULL

    names(report$Sigma_q) <- pings$class |> levels()
    report$Sigma_q <- report$Sigma_q |> 
        names() |>
        lapply(
            function(class) {
                if( missing_classes[[class]] ) return( NULL )
                return( report$Sigma_q[[class]] )
            }
        )
    names(report$Sigma_q) <- class |> levels()

    track<- data.frame(
            date = time_mesh,
            sdr_est$coordinates,
            std.error = sdr_se$coordinates
        ) |>
        sf::st_as_sf(coords = 2:3)

    return(
        list(
            pings = pings,
            track = track,
            error_matrices = report$Sigma_q,
            spline = spline,
            robust_optimization = fit,
            mcreport = mc,
            call = call
        )
    )
}



#' Interpolate a fitted track at new dates.
#'
#' @param fit 
#'     A fitted model from fit_track.
#' @param date 
#'     The dates at which to predict.
#' @param quantiles 
#'     The quantiles to predict. Defaults to 0.05, 0.5, and 0.95.
#'
#' @return 
#'     A list with the predicted quantiles of the coordinates at each date and
#'     the posterior track samples. If length(quantiles) < 1 then only the
#'     posterior samples are returned.
#'
#' @export
predict_track <- function(
        fit,
        date,
        quantiles = c(0.05, 0.5, 0.95)
    ) {
    orig_date <- date

    # Convert time to timescale used for spline
    date <- date |> 
        difftime(fit$start_time, units = fit$time_units) |> 
        as.numeric()

    # Recreate spline with same node but x is updated to date
    spline <- nnspline::create_nnspline(
        x = date,
        nodes = fit$template_spline$nodes,
        n_parents = fit$template_spline$n_parents,
        parameters = fit$template_spline$parameters,
        covariance_function = fit$template_spline$covariance_function,
        node_graph = fit$template_spline$node_graph,
        LT = fit$template_spline$LT
    )

    # Predict coordinates for each replicate
    mc <- fit$mcreport
    arg_names<- mcreportRTMB::mcreport |> formals() |> names()
    mcreport_args <- fit$call[(fit$call |> names()) %in% arg_names]
    parallel <- if( "parallel" %in% (mcreport_args |> names()) ) {
        mcreport_args$parallel
    } else {
        1
    }
    silent <- if( "silent" %in% (mcreport_args |> names()) ) {
        mcreport_args$silent
    } else {
        TRUE
    }
    if( (parallel > 1) && requireNamespace("parallel", quietly = TRUE) ) {
        lapplyfn <- parallel::mclapply
    } else {
        lapplyfn <- lapply
    }
    interpolate_replicate<- function(i, ...) {
        if( !silent ) {
            cat(
                paste0(
                    "\rUpdating mcreplicate: (", i, " / ", 
                    mc[[1]] |> ncol(), ")"
                )
            )
            if( i == (mc[[1]] |> ncol()) ) cat("\n")
            flush.console()
        }
        x_spline <- spline |>
            nnspline::update_spline(
                parameters = mc$spline_parameters[, 1, i],
                node_values = mc$node_values[, 1, i]
            )
        y_spline <- spline |>
            nnspline::update_spline(
                parameters = mc$spline_parameters[, 2, i],
                node_values = mc$node_values[, 2, i]
            )
        m <- cbind(
            mc$center[1, i] + x_spline$values,
            mc$center[2, i] + y_spline$values
        )
        m <- data.frame(
            replicate = i,
            date = orig_date,
            geometry = m
        )
        return( m )
    }
    replicates <- seq(ncol(mc[[1]])) |> 
        lapplyfn(
            interpolate_replicate,
            mc.cores = parallel,
            mc.preschedule = FALSE
        )
    sample_tracks <- rbind |> do.call(replicates)
    sample_tracks <- sample_tracks |> 
        sf::st_as_sf(
            coords = 3:4,
            crs = fit$pings |> sf::st_crs()
        )
    rownames(sample_tracks) <- NULL
    if( (quantiles |> length()) < 1 ) return( sample_tracks )


    quantiles <- quantiles |> sort()
    quantile_tracks <- abind::abind |> 
        do.call(
            c(
                lapply(replicates, `[`, 3:4),
                rev.along = 0
            )
        )
    quantile_tracks <- quantile_tracks |> 
        apply(
            MARGIN = 1:2,
            quantile,
            probs = quantiles
        )
    quantile_tracks <- quantile_tracks |> aperm(c(2:3, 1))
    quantile_tracks <- quantiles |> 
        seq_along() |>
        lapply(
            function(i) {
                m <- data.frame(
                    quantile = quantiles[i],
                    date = orig_date,
                    geometry = quantile_tracks[, , i]
                )
                return(m)
            }
        )
    quantile_tracks <- rbind |> do.call(quantile_tracks)
    quantile_tracks <- quantile_tracks |>
        sf::st_as_sf(
            coords = 3:4,
            crs = fit$pings |> sf::st_crs()
        )
    rownames(quantile_tracks) <- NULL

    return(
        list(
            quantiles = quantile_tracks,
            samples = sample_tracks
        )
    )
}
