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

    ll<- ll + sum(RTMB::dnorm(qstretch, 0.5, 1, TRUE))
    stretch<- RTMB::plogis(qstretch)
    height<- exp(log_height)

    # Track likelihood and predictions
    ll<- ll + nnspline::dlcspline(
        diff(coordinates[, 1]) / as.numeric(diff(time_mesh)),
        spline,
        stretch[1],
        height[1],
        log = TRUE
    )
    ll<- ll + nnspline::dlcspline(
        diff(coordinates[, 2]) / as.numeric(diff(time_mesh)),
        spline,
        stretch[2],
        height[2],
        log = TRUE
    )
    coordinates |> mcreportRTMB::MCREPORT()

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

    ping_ll<- pings |> nrow() |> numeric() |> RTMB::AD()
    for( i in pings |> nrow() |> seq() ) {
        ping_ll[[i]]<- RTMB::dmvnorm(
            sf::st_coordinates(pings)[i, ],
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

