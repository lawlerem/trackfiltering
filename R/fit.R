#' Fit a model to a movement track
#' 
#' @param pings An sf data frame with a POSIX or numeric column named "date", an optional factor column named "class", and point geometries.
#' @param interpolation_time A vector of times at which to interpolate locations.
#' @param nodes A vector of times used as nodes in the underlying spline.
#' @param robust_schedule See robustifyRTMB::robustly_optimize
#' @param robust_function See robustifyRTMB::robustify
#' @param robust_bandwidth See robustifyRTMB::robustly_optimize. Should be specified as a numeric values measured in units given by the time_units arguments.
#' @param independent_coordinates If true, then observation ellipses are assumed to have
#'   independent axes.
#' @param common_coordinate_correlation If true, then observation ellipses for different
#'   quality classes are assumed to have the same correlation structure.
#' @param time_units Which time units should be used? Mostly matters for numerical stability.
#' @param fix_range If set to a value, the range parameter will be fixed to the given value.
#' @param correlation_function See nnspline::create_nnspline
#'
#' @export
fit_track<- function(
        pings,
        interpolation_time = numeric(0),
        nodes = unique(pings$date),
        robust_schedule = seq(0, 0.05, by = 0.01),
        robust_function = "ll",
        robust_bandwidth = 1,
        independent_coordinates = FALSE,
        common_coordinate_correlation = TRUE,
        time_units = "days",
        fix_range,
        correlation_function = function(x1, x2, p) {
            d<- sqrt(sum((x1 - x2)^2))
            range<- p[[1]]
            # (0.05 * (d == 0) + 0.95 ) * exp( -(d / range)^2 )
            exp( -(d / range)^2 )
        }
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
        n_parents = 4,
        correlation_function = correlation_function
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
            median
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
        working_ping_diagonal = matrix(
            seq_along(pars$working_ping_diagonal),
            nrow = nrow(pars$working_ping_diagonal)
        ),
        ping_off_diagonal = matrix(
            seq_along(pars$ping_off_diagonal),
            nrow = nrow(pars$ping_off_diagonal)
        )
    )
    missing_classes<- table(class) == 0
    pars$working_ping_diagonal[, missing_classes]<- NA
    map$working_ping_diagonal[, missing_classes]<- NA
    map$ping_off_diagonal[, missing_classes]<- NA
    if( common_coordinate_correlation ) {
        map$ping_off_diagonal[]<- seq(nrow(map$ping_off_diagonal))
    } else {}
    if( independent_coordinates ) {
        pars$ping_off_diagonal[]<- 0
        map$ping_off_diagonal[]<- NA
    } else {}
    map<- lapply(map, as.factor)

    environment(nll)<- environment()
    # fit<- robustifyRTMB::robustly_optimize(
    #     nll,
    #     pars,
    #     random = "node_values",
    #     map = map,
    #     smooth = "node_values",
    #     nodes = spline$nodes,
    #     robust_schedule = robust_schedule,
    #     bandwidth = as.difftime(robust_bandwidth, units = time_units) |> as.numeric()
    # )

    ### Bespoke robust optimization
    weight_fun<- Vectorize(
        function(d, bandwidth) exp( -d / bandwidth ),
        "d"
    )
    bandwidth<- as.difftime(robust_bandwidth, units = time_units) |> as.numeric()
    fitpar<- c(
        pars,
        list(
            robustness = max(robust_schedule)
        )
    )
    fitpar<- as.relistable(fitpar)
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
    d<- dist(
        c(
            spline$nodes,
            environment(nll)$time
        ),
        diag = TRUE,
        upper = TRUE
    )
    d<- as.matrix(d)
    within<- d <= 5 * bandwidth
    weights<- 0 * d
    weights[within]<- weight_fun(d[within], bandwidth)
    ping_weights<- weights[seq(nrow(spline$nodes)), seq_along(environment(nll)$time) + nrow(spline$nodes)]
    zero_rows<- rowSums(ping_weights) == 0
    ping_weights[!zero_rows, ]<- sweep(
        ping_weights[!zero_rows, ],
        MARGIN = 1,
        STATS = rowSums(ping_weights[!zero_rows, ]),
        FUN = "/"
    )
    fitpar$node_values[!zero_rows, ]<- ping_weights[!zero_rows, ] %*% environment(nll)$coordinates

    # node_weights<- weights[seq(nrow(spline$nodes)), seq(nrow(spline$nodes))]

    # Get good starting values for track parameters
    opt_spline<- nlminb(
        c(fitpar$working_variance),
        function(x) {
            splinex<- nnspline::update_spline_covariance(
                spline,
                exp(x[[1]]),
                exp(fitpar$working_range[[1]])
            )
            spliney<- nnspline::update_spline_covariance(
                spline,
                exp(x[[2]]),
                exp(fitpar$working_range[[2]])
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
    )
    fitpar$working_variance<- opt_spline$par

    # Get good starting values for observation error
    obj<- RTMB::MakeADFun(
        nll,
        fitpar,
        map = c(list(working_variance = as.factor(c(NA, NA))), spline_map, map, robust_map),
        random = "node_values",
        silent = TRUE
    )
    opt<- with(obj, nlminb(par, fn, gr))
    fitpar<- obj$env$parList()
    report<- obj$report()
    residuals<- report$ping_pred - environment(nll)$coordinates
    class_sd<- by(
        residuals,
        environment(nll)$class,
        function(x) {
            return(robust::covRob(x)$cov)
        }
    )

    fitpar$working_ping_diagonal<- sapply(
        class_sd,
        function(x) {
            if( is.null(x) ) return(c(0, 0)) else return(log(sqrt(c(x[1, 1], x[2, 2]))))
        }
    )
    fitpar$ping_off_diagonal[]<- 0

    # Do a final optimization
    obj<- RTMB::MakeADFun(
        nll,
        fitpar,
        map = c(map, robust_map),
        random = "node_values",
        silent = TRUE
    )
    opt<- with(obj, nlminb(par, fn, gr))
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
#' @param interpolate_time A vector of times at which to interpolate locations.
#' @param nodes A vector of times used as nodes in the underlying spline.
#' @param robust_schedule See robustifyRTMB::robustly_optimize
#' @param robust_function See robustifyRTMB::robustify
#' @param robust_bandwidth See robustifyRTMB::robustly_optimize
#' @param time_units Which time units should be used? Mostly matters for numerical stability.
#'
#' @export
fit_acoustic_track<- function(
        pings,
        interpolation_time = numeric(0),
        nodes = time,
        robust_schedule = seq(0, 0.05, by = 0.01),
        robust_function = "ll",
        robust_bandwidth = 0.001 * diff(range(time)),
        time_units = "days"
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
        n_parents = 5
    )
    n_coordinates<- ncol(coordinates)
    n_class<- length(levels(class))

    pars<- list(
        working_variance = numeric(n_coordinates),
        working_range = numeric(n_coordinates),
        node_values = cbind(
            spline$node_values,
            spline$node_values
        ),
        center = apply(
            sf::st_coordinates(pings),
            MARGIN = 2,
            median
        )
    )

    environment(acoustic_nll)<- environment()
    fit<- robustifyRTMB::robustly_optimize(
        acoustic_nll,
        pars,
        random = "node_values",
        smooth = "node_values",
        nodes = spline$nodes,
        robust_schedule = robust_schedule,
        bandwidth = robust_bandwidth
    )

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
            error_matrices = report$Sigma_q,
            obj = fit$obj,
            opt = fit$opt,
            sdr = fit$sdr,
            sdr_est = sdr_est,
            sdr_se = sdr_se
        )
    )
}




