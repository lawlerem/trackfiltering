library(sf)
library(trackfiltering)

i<- 2
# i<- 14

# i<- 0
# i<- i + 1

track.dir<- "~/research/analysis/grey_seal_space_use/data/derived/tracks/"
track.file<- paste0(track.dir, list.files(track.dir)[[i]])
pings<- st_read(track.file)[, c("date", "class")]
pings<- subset(
    pings,
    !is.na(date) &
    !is.na(class)
)
pings$class<- factor(
    pings$class,
    levels = c("G", "3", "2", "1", "0", "A", "B", "Z"),
    ordered = TRUE
)
interpolate_time<- seq(
    min(pings$date),
    max(pings$date),
    by = "1 days"
)
full_pings<- pings

classes<- c()
classes<- c(classes, "G")
classes<- c(classes, "3")
classes<- c(classes, "2")
classes<- c(classes, "1")
classes<- c(classes, "0")
classes<- c(classes, "A")
classes<- c(classes, "B")
classes<- c(classes, "Z")
pings<- subset(full_pings, class %in% classes)
good_classes<- names(which(table(pings$class) >= 10))

fit<- fit_track(
    subset(pings, class %in% good_classes),
    interpolate_time,
    nodes = seq(
        min(pings$date),
        max(pings$date),
        by = "2 days"
    ),
    independent_coordinates = TRUE,
    robust_schedule = seq(0, 0.05, by = 0.01)
)



pred<- data.frame(
    date = fit$interpolation$date,
    x_est = sf::st_coordinates(fit$interpolation$est)[, 1],
    y_est = sf::st_coordinates(fit$interpolation$est)[, 2],
    x_se = sf::st_coordinates(fit$interpolation$se)[, 1],
    y_se = sf::st_coordinates(fit$interpolation$se)[, 2]
)
pred<- within(
    pred,
    {
        x_lower<- x_est - 1.96 * x_se
        x_upper<- x_est + 1.96 * x_se
        y_lower<- y_est - 1.96 * y_se
        y_upper<- y_est + 1.96 * y_se
    }
)

ylims<- apply(
    rbind(
        as.matrix(sf::st_coordinates(pings)),
        as.matrix(pred[, c("x_lower", "y_lower")]),
        as.matrix(pred[, c("x_upper", "y_upper")])
    ),
    MARGIN = 2,
    range
)

outlier_weight_cutoff<- 0.8
opar<- par(
    mfrow = c(2, 1),
    mar = c(4.2, 4.2, 0.4, 0.4)
)
plot(
    y = sf::st_coordinates(fit$pings)[, 1],
    x = fit$pings$date,
    pch = 19,
    col = fit$pings$class,
    xlab = "Time",
    ylab = "X Coordinate",
    ylim = ylims[, 1]
)
lines(
    y = pred$x_est,
    x = pred$date
)
lines(
    y = pred$x_lower,
    x = pred$date,
    lty = 2
)
lines(
    y = pred$x_upper,
    x = pred$date,
    lty = 2
)
legend(
    "topright",
    legend = levels(fit$pings$class),
    col = seq_along(levels(fit$pings$class)),
    pch = 19
)
if( any(fit$pings$weight < outlier_weight_cutoff) ) {
    text(
        y = sf::st_coordinates(fit$pings)[fit$pings$weight < outlier_weight_cutoff, 1],
        x = fit$pings$date[fit$pings$weight < outlier_weight_cutoff],
        labels = round(fit$pings$weight, 2)[fit$pings$weight < outlier_weight_cutoff],
        pos = 3
    )
}

plot(
    y = sf::st_coordinates(fit$pings)[, 2],
    x = fit$pings$date,
    pch = 19,
    col = fit$pings$class,
    xlab = "Time",
    ylab = "Y Coordinate",
    ylim = ylims[, 2]
)
lines(
    y = pred$y_est,
    x = pred$date
)
lines(
    y = pred$y_lower,
    x = pred$date,
    lty = 2
)
lines(
    y = pred$y_upper,
    x = pred$date,
    lty = 2
)
if( any(fit$pings$weight < outlier_weight_cutoff) ) {
    text(
        y = sf::st_coordinates(fit$pings)[fit$pings$weight < outlier_weight_cutoff, 2],
        x = fit$pings$date[fit$pings$weight < outlier_weight_cutoff],
        labels = round(fit$pings$weight, 2)[fit$pings$weight < outlier_weight_cutoff],
        pos = 3
    )
}
par(opar)



