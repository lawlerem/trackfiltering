library(sf)
# library(trackfiltering)
devtools::load_all()
# crs<- readRDS("track_crs.rds")
# tracks<- sf::st_read("tracks.gpkg") |>
#     sf::st_transform(crs)
# i<- 0

# # i<- i + 1
# deployid<- sample(unique(tracks$deployid), 1)
# print(deployid)
# # track<- subset(tracks, deployid == unique(deployid)[[i]])
# track<- tracks[tracks$deployid == deployid, ]
track<- sf::st_read("example_track.gpkg")
track$class<- factor(
    track$class,
    levels = c("G", "3", "2", "1", "0", "A", "B", "Z"),
    ordered = TRUE
)
small_classes<- names(which(table(track$class) <= 5))
track<- subset(
    track,
    !(class %in% small_classes) &
    !is.na(date)
)
print(track)
fit<- trackfiltering::fit_track(
    track,
    interpolation_time = seq(
        min(track$date),
        max(track$date),
        by = "1 hours"
    ),
    nodes = seq(
        min(track$date),
        max(track$date),
        by = "1 days"
    ),
    # fix_range = 15,
    robust_schedule = 0.05,
    robust_function = "ll",
    # independent_coordinates = FALSE,
    # common_coordinate_correlation = FALSE,
    robust_bandwidth = 3
    # correlation_function = function(x1, x2, p) {
    #     d<- sqrt(sum((x1 - x2)^2))
    #     range<- p[[1]]
    #     poly<- 1 + sqrt(5) * (d / range) + (5 / 3) * (d / range)^2
    #     return(poly * exp( -sqrt(5) * (d / range)))
    # }
)

fit$opt
fit$sdr
fit$error_matrices
by(
    fit$pings$weights,
    fit$pings$class,
    function(x) c(round(summary(x), 2), n = length(x))
)

est<- sf::st_coordinates(fit$interpolation$est)
se<- sf::st_coordinates(fit$interpolation$se)
lower<- est - 1.96 * se
upper<- est + 1.96 * se

par(mfrow = c(2, 1))
    plot(
        x = range(fit$interpolation$date),
        y = range(c(lower[, 1], upper[, 1])),
        type = "n",
        xlab = "Date",
        ylab = "X coordinate"
    )
    lines(
        x = fit$interpolation$date,
        y = est[, 1]
    )
    lines(
        x = fit$interpolation$date,
        y = lower[, 1],
        lty = 2
    )
    lines(
        x = fit$interpolation$date,
        y = upper[, 1],
        lty = 2
    )
    points(
        x = fit$pings$date,
        y = sf::st_coordinates(fit$pings)[, 1],
        col = factor(
            fit$pings$class,
            levels = c("G", "3", "2", "1", "0", "A", "B", "Z"),
            ordered = TRUE
        ),
        pch = 19
    )
    with(
        subset(fit$pings, weights < 0.2),
        text(
            x = date,
            y = sf::st_coordinates(geom)[, 1],
            labels = round(weights, 1),
            pos = 3
        )
    )

    plot(
        x = range(fit$interpolation$date),
        y = range(c(lower[, 2], upper[, 2])),
        type = "n",
        xlab = "Date",
        ylab = "Y coordinate"
    )
    lines(
        x = fit$interpolation$date,
        y = est[, 2]
    )
    lines(
        x = fit$interpolation$date,
        y = lower[, 2],
        lty = 2
    )
    lines(
        x = fit$interpolation$date,
        y = upper[, 2],
        lty = 2
    )
    points(
        x = fit$pings$date,
        y = sf::st_coordinates(fit$pings)[, 2],
        col = factor(
            fit$pings$class,
            levels = c("G", "3", "2", "1", "0", "A", "B", "Z"),
            ordered = TRUE
        ),
        pch = 19
    )
    with(
        subset(fit$pings, weights < 0.2),
        text(
            x = date,
            y = sf::st_coordinates(geom)[, 2],
            labels = round(weights, 1),
            pos = 3
        )
    )

