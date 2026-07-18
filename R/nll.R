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
    ll<- ll + sum(RTMB::dnorm(qstretch, 0, 1, log = TRUE))
    ll<- ll + RTMB::dflat(coordinates[1, 1], log = TRUE)
    ll<- ll + RTMB::dflat(coordinates[1, 2], log = TRUE)
    ll<- ll + nnspline::dlcspline(
        diff(coordinates[, 1]) / diff(num_time_mesh),
        spline,
        RTMB::plogis(qstretch[1]),
        exp(lheight[1]),
        log = TRUE
    )
    ll<- ll + nnspline::dlcspline(
        diff(coordinates[, 2]) / diff(num_time_mesh),
        spline,
        RTMB::plogis(qstretch[2]),
        exp(lheight[2]),
        log = TRUE
    )
    coordinates |> mcreportRTMB::MCREPORT()

    # Ping likelihood and robust weights
    ping_scale<- working_ping_scale |> exp()
    if( pings$class |> is.ordered() ) {
        for( i in ping_diagonal |> nrow() |> seq() ) {
            ping_diagonal[i, ]<- ping_diagonal[i, ] |> cumsum()
        }
    }
    Sigma_q<- ping_scale |>
        seq_along() |> 
        lapply(
            function(q) {
                orient<- working_ping_orientation[, q] |> theta2cor()
                shape<- working_ping_shape[, q] |>
                    exp() |>
                    (\(x) RTMB::diag(c(RTMB::AD(1), x)))()
                scale<- ping_scale[q]
                Sigma<- scale * (shape %*% orient %*% shape)
                return(Sigma)
            }
        )
    Sigma_q |> RTMB::REPORT()

    ping_ll<- pings |> nrow() |> numeric() |> RTMB::AD()
    ping_coordinates<- RTMB::OBS(ping_coordinates)
    for( i in pings |> nrow() |> seq() ) {
        ping_ll[[i]]<- RTMB::dmvnorm(
            # sf::st_coordinates(pings)[i, ],
            ping_coordinates[i, ],
            coordinates[pings$spline_idx[i], ],
            Sigma_q[[pings$class[i]]],
            log = TRUE
        )
    }
    ll<- ll + sum(robustifyRTMB::robustify(ping_ll, robustness, "ll"))
    weights<- ping_ll |> robustifyRTMB::robust_weight(robustness, "ll")
    weights |> RTMB::REPORT()
    return( -ll )
}

