theta2cor<- function(theta) {
    n<- floor(sqrt(2 * length(theta))) + 1
    m<- RTMB::matrix(0, nrow = n, ncol = n)
    m[lower.tri(m)]<- theta
    h<- (exp(m) - 1) / (exp(m) + 1)
    
    scalings<- sqrt(1 - h^2)
    scalings[upper.tri(scalings, diag = TRUE)]<- 0
    scalings<- t(RTMB::apply(scalings, MARGIN = 1, cumprod))
    scalings<- scalings[, c(ncol(scalings), 1:(ncol(scalings) - 1))]
    scalings[, 1]<- 1
    
    L<- h
    diag(L)[]<- 1
    L<- L * scalings
    cor<- L %*% t(L)
    return(cor)
}

nll<- function(pars) {
    "[<-"<- RTMB::ADoverload("[<-")
    RTMB::getAll(pars)
    variance<- exp(working_variance)
    range<- exp(working_range)

    RTMB::ADREPORT(variance)
    RTMB::ADREPORT(range)

    ll<- 0

    # Track likelihood and predictions
    xpars<- c(variance[[1]], range[[1]])
    x_spline<- nnspline::update_spline_covariance(
        spline,
        parameters = xpars
    )
    ll<- ll + nnspline::dspline(
        node_values[, 1],
        x_spline,
        log = TRUE
    )
    x_spline<- nnspline::update_spline_values(
        x_spline,
        node_values[, 1]
    )

    ypars<- c(variance[[2]], range[[2]])
    y_spline<- nnspline::update_spline_covariance(
        spline,
        parameters = ypars
    )
    ll<- ll + nnspline::dspline(
        node_values[, 2],
        y_spline,
        log = TRUE
    )
    y_spline<- nnspline::update_spline_values(
        y_spline,
        node_values[, 2]
    )
    
    interpolated_coordinates<- cbind(
        c(center[1]) + c(nnspline::nns(
            interpolation_time,
            x_spline
        )),
        center[2] + nnspline::nns(
            interpolation_time,
            y_spline
        )
    )
    ping_pred<- cbind(
        center[1] + nnspline::nns(
            time,
            x_spline
        ),
        center[2] + nnspline::nns(
            time,
            y_spline
        )
    )
    RTMB::REPORT(interpolated_coordinates)
    RTMB::ADREPORT(interpolated_coordinates)
    RTMB::REPORT(ping_pred)
    RTMB::ADREPORT(ping_pred)

    # Ping likelihood and robust weights
    ping_diagonal<- exp(working_ping_diagonal)
    ping_diagonal[]<- ping_diagonal[]
    Sigma_q<- lapply(
        seq(n_class),
        function(q) {
            cor<- theta2cor(ping_off_diagonal[, q])
            Sigma<- RTMB::diag(ping_diagonal[, q])
            Sigma<- Sigma %*% cor %*% Sigma
            return(Sigma)
        }
    )
    RTMB::REPORT(Sigma_q)

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
    RTMB::REPORT(weights)
    return(-1.0 * ll)
}



acoustic_nll<- function(pars) {
    "[<-"<- RTMB::ADoverload("[<-")
    RTMB::getAll(pars)
    variance<- exp(working_variance)
    range<- exp(working_range)

    RTMB::ADREPORT(variance)
    RTMB::ADREPORT(range)

    ll<- 0

    # Track likelihood and predictions
    xpars<- c(variance[[1]], range[[1]])
    x_spline<- nnspline::update_spline_covariance(
        spline,
        parameters = xpars
    )
    ll<- ll + nnspline::dspline(
        node_values[, 1],
        x_spline,
        log = TRUE
    )
    x_spline<- nnspline::update_spline_values(
        x_spline,
        node_values[, 1]
    )

    ypars<- c(variance[[2]], range[[2]])
    y_spline<- nnspline::update_spline_covariance(
        spline,
        parameters = ypars
    )
    ll<- ll + nnspline::dspline(
        node_values[, 2],
        y_spline,
        log = TRUE
    )
    y_spline<- nnspline::update_spline_values(
        y_spline,
        node_values[, 2]
    )
    
    interpolated_coordinates<- cbind(
        c(center[1]) + c(nnspline::nns(
            interpolation_time,
            x_spline
        )),
        center[2] + nnspline::nns(
            interpolation_time,
            y_spline
        )
    )
    ping_pred<- cbind(
        center[1] + nnspline::nns(
            time,
            x_spline
        ),
        center[2] + nnspline::nns(
            time,
            y_spline
        )
    )
    RTMB::REPORT(interpolated_coordinates)
    RTMB::ADREPORT(interpolated_coordinates)
    RTMB::REPORT(ping_pred)
    RTMB::ADREPORT(ping_pred)

    # Ping likelihood and robust weights
    weights<- numeric(nrow(coordinates))
    coordinates<- RTMB::OBS(coordinates)
    for( i in seq(nrow(coordinates)) ) {
        this_var<- (coordinates_se[i, ])^2
        Sigma<- RTMB::diag(this_var, nrow = 2, ncol = 2)
        this_ll<- RTMB::dmvnorm(
            coordinates[i, ],
            ping_pred[i, ],
            Sigma,
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
    RTMB::REPORT(weights)
    return(-1.0 * ll)
}




