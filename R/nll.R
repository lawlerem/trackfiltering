theta2cor<- function(theta) {
    n<- floor(sqrt(2 * length(theta))) + 1
    m<- matrix(0, nrow = n, ncol = n) |> RTMB::AD()
    m[m |> lower.tri()]<- theta
    h<- (exp(m) - 1) / (exp(m) + 1)
    
    scalings<- sqrt(1 - h^2)
    scalings[scalings |> upper.tri(diag = TRUE)]<- 0
    scalings<- scalings |> RTMB::apply(MARGIN = 1, cumprod) |> t()
    scalings<- scalings[, c(scalings |> ncol(), 1:(ncol(scalings) - 1))]
    scalings[, 1]<- 1
    
    L<- h
    diag(L)[]<- 1
    L<- L * scalings
    cor<- L %*% t(L)
    return( cor )
}

nll<- function(pars) {
    pars |> RTMB::getAll()
    ll<- 0

    # Track likelihood and predictions
    spline_parameters = working_spline_parameters |> exp()
    x_spline<- spline |>
        nnspline::update_spline(
            parameters = spline_parameters[, 1],
            node_values = node_values[, 1]
        )
    ll<- ll + nnspline::dspline(
        node_values[, 1],
        x_spline,
        log = TRUE
    )
    y_spline<- spline |>
        nnspline::update_spline(
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
    spline_parameters |> mcreportRTMB::MCREPORT()
    node_values |> mcreportRTMB::MCREPORT()
    center |> mcreportRTMB::MCREPORT()
    predicted_coordinates |> RTMB::ADREPORT()

    # Ping likelihood and robust weights
    ping_diagonal<- working_ping_diagonal |> exp()
    if( pings$class |> is.ordered() ) {
        for( i in ping_diagonal |> nrow() |> seq() ) {
            ping_diagonal[i, ]<- ping_diagonal[i, ] |> cumsum()
        }
    }
    Sigma_q<- ping_diagonal |>
        ncol() |> 
        seq() |> 
        lapply(
            function(q) {
                cor<- ping_off_diagonal[, q] |> theta2cor()
                Sigma<- ping_diagonal[, q] |> RTMB::diag()
                Sigma<- Sigma %*% cor %*% Sigma
                return(Sigma)
            }
        )
    Sigma_q |> RTMB::REPORT()

    raw_ll<- pings |> nrow() |> numeric() |> RTMB::AD()
    for( i in pings |> nrow() |> seq() ) {
        raw_ll[i]<- RTMB::dmvnorm(
            sf::st_coordinates(pings)[i, ],
            predicted_coordinates[pings$spline_idx[i], ],
            Sigma_q[[pings$class[i]]],
            log = TRUE
        )
    }
    ll<- ll + sum(robustifyRTMB::robustify(raw_ll, robustness, "ll"))
    weights<- raw_ll |> robustifyRTMB::robust_weight(robustness, "ll")
    weights |> RTMB::REPORT()
    return( -ll )
}

