theta2cor<- function(theta) {
    n<- floor(sqrt(2 * length(theta))) + 1
    m<- RTMB::AD(matrix(0, nrow = n, ncol = n))
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
    RTMB::getAll(pars)
    ll<- 0

    # Track likelihood and predictions
    spline_parameters = exp(working_spline_parameters)
    x_spline<- nnspline::update_spline(
        spline,
        parameters = spline_parameters[, 1],
        node_values = node_values[, 1]
    )
    ll<- ll + nnspline::dspline(
        node_values[, 1],
        x_spline,
        log = TRUE
    )
    y_spline<- nnspline::update_spline(
        spline,
        parameters = spline_parameters[, 2],
        node_values = node_values[, 2]
    )
    ll<- ll + nnspline::dspline(
        node_values[, 2],
        y_spline,
        log = TRUE
    )
    predicted_coordinates<- cbind(
        center[1] + x_spline$values,
        center[2] + y_spline$values
    )
    mcreportRTMB::MCREPORT(spline_parameters)
    mcreportRTMB::MCREPORT(node_values)
    mcreportRTMB::MCREPORT(center)
    RTMB::ADREPORT(predicted_coordinates)

    # Ping likelihood and robust weights
    ping_diagonal<- exp(working_ping_diagonal)
    if( is.ordered(pings$class) ) {
        for( i in seq(nrow(ping_diagonal)) ) {
            ping_diagonal[i, ]<- cumsum(ping_diagonal[i, ])
        }
    }
    Sigma_q<- lapply(
        seq(ncol(ping_diagonal)),
        function(q) {
            cor<- theta2cor(ping_off_diagonal[, q])
            Sigma<- RTMB::diag(ping_diagonal[, q])
            Sigma<- Sigma %*% cor %*% Sigma
            return(Sigma)
        }
    )
    RTMB::REPORT(Sigma_q)

    raw_ll<- RTMB::AD(numeric(nrow(pings)))
    for( i in seq(nrow(pings)) ) {
        raw_ll[i]<- RTMB::dmvnorm(
            sf::st_coordinates(pings)[i, ],
            predicted_coordinates[pings$spline_idx[i], ],
            Sigma_q[[pings$class[i]]],
            log = TRUE
        )
    }
    ll<- ll + sum(robustifyRTMB::robustify(raw_ll, robustness, "ll"))
    weights<- robustifyRTMB::robust_weight(raw_ll, robustness, "ll")
    RTMB::REPORT(weights)
    return( -ll )
}

