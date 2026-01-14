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
#' @param ping_common_orientation
#'     If true, all classes' error ellipses will have the same orientation.
#' @param ping_standard_orientation
#'     If true, all classes' error ellipses will be oriented along the
#'         coordinate axes, i.e. have uncorrelated components.
#' @param ping_common_shape
#'     If true, all classes' error ellipses will have the same shape, i.e. the
#'         ratio of the standard deviation of coordinates over the standard
#'         deviation of the first coordinate.
#' @param ping_equal_shape
#'     If true, all classes' error ellipses will have unit shape for all
#'         coordinates, i.e. have the same standard deviation for each 
#'         coordinate.
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
        utils::tail(-1)
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
            pings |> sf::st_coordinates() |> apply(2, stats::median),
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
        sf::st_as_sf(
            coords = 2:3,
            crs = sf::st_crs(pings)
        )

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
