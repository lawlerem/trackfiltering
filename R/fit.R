#' Fit a model to a movement track
#' 
#' @param pings An sf data frame with a POSIX or numeric column named "date", a factor column named "class", and point geometries. If class is an ordered factor, then the estimated observation error variances will increase with class.
#' @param interpolation_time A vector of times at which to interpolate locations.
#' @param nodes A vector of times used as nodes in the underlying spline.
#' @param robustness A numeric tuning parameter. 0 corresponds to maximum likelihood, and increasing values correspond to more aggressive outlier detection.
#' @param robust_bandwidth See robustifyRTMB::robustly_optimize. Should be specified as a numeric values measured in units given by the time_units arguments.
#' @param independent_coordinates If true, then observation ellipses are assumed to have
#'   independent axes.
#' @param common_coordinate_correlation If true, then observation ellipses for different
#'   quality classes are assumed to have the same correlation structure.
#' @param time_units Which time units should be used? Mostly matters for numerical stability.
#' @param fix_range If set to a value, the range parameter will be fixed to the given value.
#' @param covariance_function See nnspline::create_nnspline. The parameters should be a vector of length two where the first element of the variance and the second element is the range.
#' @param silent Should optimization tracing be suppressed?
#' 
#' @return A list with the following:
#'   - pings A copy of the provided location pings with an additional weights column giving the robust weight of each observation, an "est" geometry column giving the estimated location, and an "se" geometry column giving the standard errors around the estimate.
#'   - interpolation An sf data.frame with a column giving the timestamp for the interpolation time, an "est" geometry column giving the estimated location, and an "se" geometry column giving the standard errors around the estimate.
#'   - error_matrices A list of estimated observation error covariance matrices for each quality class.
#'   - obj The RTMB::MakeADFun object used for the negative log-likelihood computation.
#'   - par The estimated parameter list.
#'   - opt The optimizer output.
#'   - sdr See RTMB::sdreport
#'   - sdr_est A list giving the estimates of model parameters and derived quantities
#'   - sdr_se A list giving the standard errors of model parameters and derived quantities
#'
#' @export
fit_track<- function(
        pings,
        interpolation_time = numeric(0),
        nodes = unique(pings$date),
        robustness = 0.05,
        robust_bandwidth = 1,
        independent_coordinates = FALSE,
        common_coordinate_correlation = TRUE,
        time_units = "days",
        fix_range,
        covariance_function = function(x1, x2, p) {
            d<- sqrt(sum((x1 - x2)^2))
            variance<- p[[1]]
            range<- p[[2]]
            poly<- 1 + sqrt(3) * (d / range)
            return(variance * range^3 * poly * exp( - sqrt(3) * (d / range)))
        },
        n_parents = 5,
        silent = TRUE
    ) {
    coordinates<- sf::st_coordinates(pings)
    time<- pings$date
    if( "class" %in% colnames(pings) ) {
        class<- pings$class
    } else {
        class<- rep(1, nrow(pings))
    }
    if( !is.factor(class) ) {
        class<- as.factor(class)
    } else {}
    ordered_classes<- is.ordered(class)
    
    min_time<- min(time)
    original_interpolation_time<- interpolation_time
    time<- difftime(
        time,
        min_time,
        units = time_units
    ) |> as.numeric()
    nodes<- difftime(
        nodes,
        min_time,
        units = time_units
    ) |> as.numeric()
    interpolation_time<- difftime(
        interpolation_time,
        min_time,
        units = time_units
    ) |> as.numeric()

    spline<- nnspline::create_nnspline(
        x = c(time, interpolation_time),
        nodes = nodes,
        n_parents = n_parents,
        parameters = c(1, 1),
        covariance_function = covariance_function
    )
    n_coordinates<- ncol(coordinates)
    n_class<- length(levels(class))

    if( !missing(fix_range) ) {
        fix_range<- log(fix_range)
        map_range<- as.factor(c(NA, NA))
    } else {
        fix_range<- 0
        map_range<- as.factor(c(1, 2))
    }
    pars<- list(
        working_variance = numeric(n_coordinates),
        working_range = fix_range + numeric(n_coordinates),
        node_values = cbind(
            spline$node_values,
            spline$node_values
        ),
        center = apply(
            sf::st_coordinates(pings),
            MARGIN = 2,
            stats::median
        ),
        working_ping_diagonal = matrix(
            0,
            nrow = n_coordinates,
            ncol = n_class
        ),
        ping_off_diagonal = matrix(
            0,
            nrow = n_coordinates * (n_coordinates - 1) / 2,
            ncol = n_class
        )
    )

    map<- list(
        working_range = map_range,
        center = as.factor(c(NA, NA)),
        working_ping_diagonal = matrix(
            seq_along(pars$working_ping_diagonal),
            nrow = nrow(pars$working_ping_diagonal)
        ),
        ping_off_diagonal = matrix(
            seq_along(pars$ping_off_diagonal),
            nrow = nrow(pars$ping_off_diagonal)
        )
    )

    if( any(table(class) > 0 & table(class) <= 10) ) {
        warning("One or more of the quality classes has 10 or fewer observations, which may cause estimation issues.")
    }
    missing_classes<- table(class) == 0
    pars$working_ping_diagonal[, missing_classes]<- -10
    map$working_ping_diagonal[, missing_classes]<- NA
    map$ping_off_diagonal[, missing_classes]<- NA
    
    if( common_coordinate_correlation ) {
        map$ping_off_diagonal[]<- seq(nrow(map$ping_off_diagonal))
    }
    if( independent_coordinates ) {
        pars$ping_off_diagonal[]<- 0
        map$ping_off_diagonal[]<- NA
    }
    map<- lapply(map, as.factor)

    environment(nll)<- environment()

    ### Bespoke robust optimization
    fitpar<- c(
        pars,
        list(
            robustness = robustness
        )
    )
    fitpar<- utils::as.relistable(fitpar)
    robust_map<- list(
        robustness = as.factor(NA)
    )
    spline_map<- list(
        node_values = as.factor(
            rep(
                NA,
                length(pars$node_values)
            )
        )
    )

    # Get good starting values for true path
    weight_fun<- Vectorize(
        function(d, bandwidth) exp( -d / bandwidth ),
        "d"
    )
    bandwidth<- as.difftime(robust_bandwidth, units = time_units) |> as.numeric()
    bandwidth<- bandwidth
    d<- stats::dist(
        c(
            spline$nodes,
            environment(nll)$time
        ),
        diag = TRUE,
        upper = TRUE
    )
    d<- as.matrix(d)
    weights<- d
    weights[]<- weight_fun(d, bandwidth)
    node_weights<- weights[seq(nrow(spline$nodes)), seq(nrow(spline$nodes))]
    node_weights<- sweep(
        node_weights,
        MARGIN = 1,
        STATS = rowSums(node_weights),
        FUN = "/"
    )
    ping_weights<- weights[seq(nrow(spline$nodes)), seq_along(environment(nll)$time) + nrow(spline$nodes)]
    ping_weights<- sweep(
        ping_weights,
        MARGIN = 1,
        STATS = rowSums(ping_weights),
        FUN = "/"
    )

    x_ping_ranks<- order(sf::st_coordinates(pings)[, 1])
    x_ping_weights<- ping_weights[, x_ping_ranks]
    y_ping_ranks<- order(sf::st_coordinates(pings)[, 2])
    y_ping_weights<- ping_weights[, y_ping_ranks]
    sorted_x_coord<- sf::st_coordinates(pings)[x_ping_ranks, 1]
    sorted_y_coord<- sf::st_coordinates(pings)[y_ping_ranks, 2]
    fitpar$node_values<- t(sapply(
        seq(nrow(fitpar$node_values)),
        function(i) {
            # 1.) Find the weight where sum w_->i >= 0.5
            # 2.) Linearly interpolate coordinate x_i and x_i+1
            x_cumweights<- cumsum(x_ping_weights[i, ])
            x_idx<- c(
                rev(which(x_cumweights < 0.5))[1],
                which(x_cumweights >= 0.5)[1]
            )
            if( length(x_idx) == 1 ) {
                x_val<- sorted_x_coord[x_idx]
            } else {
                w<- (x_cumweights[x_idx[2]] - 0.5) / (x_cumweights[x_idx[2]] - x_cumweights[x_idx[1]])
                x_val<- w * sorted_x_coord[x_idx[1]] + (1 - w) * sorted_x_coord[x_idx[2]]
            }

            y_cumweights<- cumsum(y_ping_weights[i, ])
            y_idx<- c(
                rev(which(y_cumweights < 0.5))[1],
                which(y_cumweights >= 0.5)[1]
            )
            if( length(y_idx) == 1 ) {
                y_val<- sorted_y_coord[y_idx]
            } else {
                w<- (y_cumweights[y_idx[2]] - 0.5) / (y_cumweights[y_idx[2]] - y_cumweights[y_idx[1]])
                y_val<- w * sorted_y_coord[y_idx[1]] + (1 - w) * sorted_y_coord[y_idx[2]]
            }

            return(c(x_val, y_val))
        }
    ))

    fitpar$node_values<- node_weights %*% fitpar$node_values
    fitpar$node_values<- sweep(fitpar$node_values, MARGIN = 2, fitpar$center)

    # Get good starting values for track parameters
    ff<- function(pars) {
        RTMB::getAll(pars)
        splinex<- nnspline::update_spline_covariance(
            spline,
            exp(c(working_variance[[1]], working_range[[1]])),
            only_node_covariance = TRUE
        )
        spliney<- nnspline::update_spline_covariance(
            spline,
            exp(c(working_variance[[2]], working_range[[2]])),
            only_node_covariance = TRUE
        )
        ll<- 0
        ll<- ll + nnspline::dspline(
            fitpar$node_values[, 1],
            splinex,
            log = TRUE
        )
        ll<- ll + nnspline::dspline(
            fitpar$node_values[, 2],
            spliney,
            log = TRUE
        )
        return( -ll )
    }
    parnames<- c("working_variance", "working_range")
    if( !silent ) cat("Estimating initial spline parameters.\n")
    obj<- RTMB::MakeADFun(
        ff,
        fitpar[parnames],
        map = map[names(map) %in% parnames],
        silent = silent
    )
    opt_spline<- stats::nlminb(
        obj$par,
        obj$fn,
        obj$gr
    )
    fitpar[parnames]<- obj$env$parList()

    # Get good starting values for observation error
    if( !silent ) cat("Estimating initial observation error parameters.\n")
    xpars<- exp(
        c(
            fitpar$working_variance[[1]],
            fitpar$working_range[[1]]
        )
    )
    x_spline<- nnspline::update_spline(
        spline,
        parameters = xpars,
        node_values = fitpar$node_values[, 1]
    )
    ypars<- exp(
        c(
            fitpar$working_variance[[2]],
            fitpar$working_range[[2]]
        )
    )
    y_spline<- nnspline::update_spline(
        spline,
        parameters = ypars,
        node_values = fitpar$node_values[, 2]
    )
    ping_pred<- cbind(
        fitpar$center[[1]] + nnspline::nns(time, x_spline),
        fitpar$center[[2]] + nnspline::nns(time, y_spline)
    )

    ff<- function(pars) {
        RTMB::getAll(pars)
        ping_diagonal<- exp(working_ping_diagonal)
        if( ordered_classes ) {
            for( i in seq(nrow(ping_diagonal)) ) {
                ping_diagonal[i, ]<- cumsum(ping_diagonal[i, ])
            }
        }
        Sigma_q<- lapply(
            seq(n_class),
            function(q) {
                cor<- theta2cor(ping_off_diagonal[, q])
                Sigma<- RTMB::diag(ping_diagonal[, q])
                Sigma<- Sigma %*% cor %*% Sigma
                return(Sigma)
            }
        )

        ll<- 0
        weights<- numeric(nrow(coordinates))
        coordinates<- RTMB::OBS(coordinates)
        for( i in seq(nrow(coordinates)) ) {
            this_ll<- RTMB::dmvnorm(
                coordinates[i, ],
                ping_pred[i, ],
                Sigma_q[[class[i]]],
                log = TRUE
            )
            ll<- ll + robustifyRTMB::robustify(
                this_ll,
                robustness,
                "ll"
            )
            weights[i]<- robustifyRTMB::robust_weight(
                this_ll,
                robustness,
                "ll"
            )
        }
        return( -ll )
    }
    parnames<- c("working_ping_diagonal", "ping_off_diagonal")
    obj<- RTMB::MakeADFun(
        ff,
        fitpar[parnames],
        map = map[names(map) %in% parnames],
        silent = silent
    )
    opt_obserr<- stats::nlminb(
        obj$par,
        obj$fn,
        obj$gr
    )
    fitpar[parnames]<- obj$env$parList()
    
    if( ordered_classes ) {
        map$working_ping_diagonal[fitpar$working_ping_diagonal < -3]<- NA
        map$working_ping_diagonal<- droplevels.factor(map$working_ping_diagonal)
    } else {}

    # Do a final optimization
    if( !silent ) cat("\n Running final optimization.\n")
    obj<- RTMB::MakeADFun(
        nll,
        fitpar,
        map = c(
            map,
            robust_map
        ),
        random = c("node_values"),
        silent = silent
    )
    opt<- with(obj, stats::nlminb(par, fn, gr))
    fitpar<- obj$env$parList()
    
    if( !silent ) cat("\n Computing standard errors.\n")
    sdr = RTMB::sdreport(
        obj,
        opt$par,
        getJointPrecision = TRUE
    )
    fit<- list(
        par = fitpar,
        obj = obj,
        opt = opt,
        sdr = sdr
    )
    ###

    report<- fit$obj$report()
    sdr_est<- c(
        as.list(fit$sdr, "Est", report = TRUE),
        as.list(fit$sdr, "Est", report = FALSE)
    )
    sdr_se<- c(
        as.list(fit$sdr, "Std", report = TRUE),
        as.list(fit$sdr, "Std", report = FALSE)
    )

    ping_est<- sf::st_as_sf(
        data.frame(sdr_est$ping_pred),
        coords = seq(n_coordinates),
        sf_column_name = "est"
    )
    ping_se<- sf::st_as_sf(
        data.frame(sdr_se$ping_pred),
        coords = seq(n_coordinates),
        sf_column_name = "se",
        na.fail = FALSE
    )
    pings<- cbind(
        weights = report$weights,
        pings,
        ping_est,
        ping_se
    )
    
    interpolate_est<- sf::st_as_sf(
        data.frame(sdr_est$interpolated_coordinates),
        coords = seq(n_coordinates),
        sf_column_name = "est"
    )
    interpolate_se<- sf::st_as_sf(
        data.frame(sdr_se$interpolated_coordinates),
        coords = seq(n_coordinates),
        sf_column_name = "se",
        na.fail = FALSE
    )
    interpolation<- cbind(
        date = original_interpolation_time,
        est = interpolate_est,
        se = interpolate_se
    )
    sf::st_crs(interpolation)<- sf::st_crs(pings)
    names(report$Sigma_q)<- levels(class)
    report$Sigma_q<- lapply(
        names(report$Sigma_q),
        function(class) {
            if( missing_classes[[class]] ) {
                return(NULL)
            } else {
                return(report$Sigma_q[[class]])
            }
        }
    )
    names(report$Sigma_q)<- levels(class)

    return(
        list(
            pings = pings,
            interpolation = interpolation,
            error_matrices = report$Sigma_q,
            obj = fit$obj,
            par = fit$par,
            opt = fit$opt,
            sdr = fit$sdr,
            sdr_est = sdr_est,
            sdr_se = sdr_se
        )
    )
}






#' Fit a model to a movement track
#' 
#' @param pings An sf data frame with a POSIX or numeric column named "date", a column named x_se, a column names y_se, and point geometries.
#' @param interpolation_time A vector of times at which to interpolate locations.
#' @param nodes A vector of times used as nodes in the underlying spline.
#' @param robustness A numeric tuning parameter. 0 corresponds to maximum likelihood, and increasing values correspond to more aggressive outlier detection.
#' @param robust_bandwidth See robustifyRTMB::robustly_optimize
#' @param time_units Which time units should be used? Mostly matters for numerical stability.
#' @param fix_range If set to a value, the range parameter will be fixed to the given value.
#' @param covariance_function See nnspline::create_nnspline
#'
#' @return A list with the following:
#'   - pings A copy of the provided location pings with an additional weights column giving the robust weight of each observation, an "est" geometry column giving the estimated location, and an "se" geometry column giving the standard errors around the estimate.
#'   - interpolation An sf data.frame with a column giving the timestamp for the interpolation time, an "est" geometry column giving the estimated location, and an "se" geometry column giving the standard errors around the estimate.
#'   - obj The RTMB::MakeADFun object used for the negative log-likelihood computation.
#'   - par The estimated parameter list.
#'   - opt The optimizer output.
#'   - sdr See RTMB::sdreport
#'   - sdr_est A list giving the estimates of model parameters and derived quantities
#'   - sdr_se A list giving the standard errors of model parameters and derived quantities
#' 
#' @export
fit_acoustic_track<- function(
        pings,
        interpolation_time = numeric(0),
        nodes = time,
        robustness = 0.1,
        robust_bandwidth = 1,
        time_units = "days",
        fix_range,
        covariance_function = function(x1, x2, p) {
            d<- sqrt(sum((x1 - x2)^2))
            variance<- p[[1]]
            range<- p[[2]]
            poly<- 1 + sqrt(5) * (d / range) + (5 / 3) * (d / range)^2
            return(variance * poly * exp( - sqrt(5) * (d / range)))
        }
    ) {
    coordinates<- sf::st_coordinates(pings)
    coordinates_se<- as.matrix(sf::st_drop_geometry(pings[, c("x_se", "y_se")]))
    time<- pings$date
    if( "class" %in% colnames(pings) ) {
        class<- pings$class
    } else {
        class<- rep(1, nrow(pings))
    }
    if( !is.factor(class) ) {
        class<- as.factor(class)
    } else {}
    if( any(table(class) <= 5) ) {
        warning("One of the quality classes has 5 or fewer observations, which may cause convergence issues.")
    }

    min_time<- min(time)
    original_interpolation_time<- interpolation_time
    time<- difftime(
        time,
        min_time,
        units = time_units
    ) |> as.numeric()
    nodes<- difftime(
        nodes,
        min_time,
        units = time_units
    ) |> as.numeric()
    interpolation_time<- difftime(
        interpolation_time,
        min_time,
        units = time_units
    ) |> as.numeric()

    spline<- nnspline::create_nnspline(
        x = c(time, interpolation_time),
        nodes = nodes,
        n_parents = 6,
        parameters = c(1, 1),
        covariance_function = covariance_function
    )
    n_coordinates<- ncol(coordinates)
    n_class<- length(levels(class))

    if( !missing(fix_range) ) {
        fix_range<- log(fix_range)
        map_range<- as.factor(c(NA, NA))
    } else {
        fix_range<- 0
        map_range<- as.factor(c(1, 2))
    }
    pars<- list(
        working_variance = numeric(n_coordinates),
        working_range = fix_range + numeric(n_coordinates),
        node_values = cbind(
            spline$node_values,
            spline$node_values
        ),
        center = apply(
            sf::st_coordinates(pings),
            MARGIN = 2,
            stats::median
        )
    )

    map<- list(
        working_range = map_range
    )
    map<- lapply(map, as.factor)

    environment(acoustic_nll)<- environment()
    
    ### Bespoke robust optimization
    fitpar<- c(
        pars,
        list(
            robustness = robustness
        )
    )
    fitpar<- utils::as.relistable(fitpar)
    robust_map<- list(
        robustness = as.factor(NA)
    )
    spline_map<- list(
        node_values = as.factor(
            rep(
                NA,
                length(pars$node_values)
            )
        )
    )

    # Get good starting values for true path
    weight_fun<- Vectorize(
        function(d, bandwidth) exp( -d / bandwidth ),
        "d"
    )
    bandwidth<- as.difftime(robust_bandwidth, units = time_units) |> as.numeric()
    bandwidth<- bandwidth
    d<- stats::dist(
        c(
            spline$nodes,
            environment(nll)$time
        ),
        diag = TRUE,
        upper = TRUE
    )
    d<- as.matrix(d)
    weights<- d
    weights[]<- weight_fun(d, bandwidth)
    ping_weights<- weights[seq(nrow(spline$nodes)), seq_along(environment(nll)$time) + nrow(spline$nodes)]
    ping_weights<- sweep(
        ping_weights,
        MARGIN = 1,
        STATS = rowSums(ping_weights),
        FUN = "/"
    )
    fitpar$node_values<- ping_weights %*% environment(nll)$coordinates
    fitpar$center<- colMeans(fitpar$node_values)
    fitpar$node_values<- sweep(fitpar$node_values, MARGIN = 2, fitpar$center)

    # Get good starting values for track parameters
    ff<- function(pars) {
        RTMB::getAll(pars)
        splinex<- nnspline::update_spline_covariance(
            spline,
            exp(c(working_variance[[1]], working_range[[1]])),
            only_node_covariance = TRUE
        )
        spliney<- nnspline::update_spline_covariance(
            spline,
            exp(c(working_variance[[2]], working_range[[2]])),
            only_node_covariance = TRUE
        )
        ll<- 0
        ll<- ll + nnspline::dspline(
            fitpar$node_values[, 1],
            splinex,
            log = TRUE
        )
        ll<- ll + nnspline::dspline(
            fitpar$node_values[, 2],
            spliney,
            log = TRUE
        )
        return( -ll )
    }
    parnames<- c("working_variance", "working_range")
    obj<- RTMB::MakeADFun(
        ff,
        fitpar[parnames],
        map = map[names(map) %in% parnames],
        silent = TRUE
    )
    opt_spline<- stats::nlminb(
        obj$par,
        obj$fn,
        obj$gr
    )
    fitpar[parnames]<- obj$env$parList()

    # Do a final optimization
    obj<- RTMB::MakeADFun(
        nll,
        fitpar,
        map = c(
            map,
            robust_map
        ),
        random = c("node_values"),
        silent = TRUE
    )
    opt<- with(obj, stats::nlminb(par, fn, gr))
    fitpar<- obj$env$parList()
    sdr = RTMB::sdreport(
        obj,
        opt$par,
        getJointPrecision = TRUE
    )
    fit<- list(
        par = fitpar,
        obj = obj,
        opt = opt,
        sdr = sdr
    )
    ###

    report<- fit$obj$report()
    sdr_est<- c(
        as.list(fit$sdr, "Est", report = TRUE),
        as.list(fit$sdr, "Est", report = FALSE)
    )
    sdr_se<- c(
        as.list(fit$sdr, "Std", report = TRUE),
        as.list(fit$sdr, "Std", report = FALSE)
    )

    ping_est<- sf::st_as_sf(
        data.frame(sdr_est$ping_pred),
        coords = seq(n_coordinates),
        sf_column_name = "est"
    )
    ping_se<- sf::st_as_sf(
        data.frame(sdr_se$ping_pred),
        coords = seq(n_coordinates),
        sf_column_name = "se",
        na.fail = FALSE
    )
    pings<- cbind(
        weights = report$weights,
        pings,
        ping_est,
        ping_se
    )
    
    interpolate_est<- sf::st_as_sf(
        data.frame(sdr_est$interpolated_coordinates),
        coords = seq(n_coordinates),
        sf_column_name = "est"
    )
    interpolate_se<- sf::st_as_sf(
        data.frame(sdr_se$interpolated_coordinates),
        coords = seq(n_coordinates),
        sf_column_name = "se",
        na.fail = FALSE
    )
    interpolation<- cbind(
        date = original_interpolation_time,
        est = interpolate_est,
        se = interpolate_se
    )
    sf::st_crs(interpolation)<- sf::st_crs(pings)

    return(
        list(
            pings = pings,
            interpolation = interpolation,
            obj = fit$obj,
            opt = fit$opt,
            sdr = fit$sdr,
            sdr_est = sdr_est,
            sdr_se = sdr_se
        )
    )
}




