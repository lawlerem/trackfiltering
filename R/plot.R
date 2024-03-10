#' Plot the coordinate profiles of a track
#'
#' @param track A list with a pings sf data.frame and a track sf data.frame.
#' @param x_se, y_se Optional vector of standard errors for the x and y coordinate.
#' @param ping_weights An optional vector of score-level weights for each ping
#' @param lower_weight_limit Pings with weight below this limit will not be plotted.
#' @param point_size_scale Scaling factor for point sizes
#'
#' @return A tmap object
#'
#' @export
plot_track_profile<- function(track, x_se, y_se, ping_weights, lower_weight_limit, point_size_scale = 1) {
  original_par<- par(
    mfrow  = c(2, 1),
    mar = c(1, 3, 2, 1)
  )
  pings<- track$pings
  track<- track$track

  if( missing(ping_weights) ) {
    ping_weights<- rep(1, nrow(pings))
  } else {}
  pings<- pings[ping_weights >= lower_weight_limit, ]
  ping_weights<- ping_weights[ping_weights >= lower_weight_limit]
  weight_labels<- which(ping_weights < 0.99)
  ping_weights<- round(ping_weights, digits = 2)

  ping_coordinates<- sf::st_coordinates(pings)
  track_coordinates<- sf::st_coordinates(track)
  all_coordinates<- rbind(ping_coordinates, track_coordinates)

  tlim<- c(min(track$time), max(track$time))
  xlim<- c(
    min(all_coordinates[, 1]),
    max(all_coordinates[, 1])
  )
  xlim<- xlim + c(-1, 1) * 0.05 * diff(xlim)
  ylim<- c(
    min(all_coordinates[, 2]),
    max(all_coordinates[, 2])
  )
  ylim<- ylim + c(-1, 1) * 0.05 * diff(ylim)
  point_size<- point_size_scale * 0.5 * log(loc_class_K[, 2] + 3)[pings$quality_class]

  plot(
    x = track$time,
    y = track_coordinates[, 1],
    xlim = tlim,
    ylim = xlim,
    main = "",
    ylab = "X coordinate",
    xlab = "Time",
    type = "n"
  )
  points(
    x = pings$time,
    y = ping_coordinates[, 1],
    pch = 20,
    cex = point_size,
    col = pings$quality_class
  )
  if( length(weight_labels) > 0 ) {
    text(
      x = pings$time[weight_labels],
      y = ping_coordinates[weight_labels, 1],
      labels = ping_weights[weight_labels],
      pos = 1
    )
  } else {}
  lines(
    x = track$time,
    y = track_coordinates[, 1]
  )
  if( !missing(x_se) ) {
    lines(
      x = track$time,
      y = track_coordinates[, 1] + 1.96 * x_se,
      lty = 2
    )
    lines(
      x = track$time,
      y = track_coordinates[, 1] - 1.96 * x_se,
      lty = 2
    )
  } else {}
  legend(
    x = "topright",
    legend = c("G", 3, 2, 1, 0, "A", "B"),
    col = seq(7),
    pt.cex = 2 * log(loc_class_K[, 2] + 1.5),
    pch = 20
  )

  plot(
    x = track$time,
    y = track_coordinates[, 2],
    xlim = tlim,
    ylim = ylim,
    main = "",
    ylab = "Y coordinate",
    xlab = "Time",
    type = "n"
  )
  points(
    x = pings$time,
    y = ping_coordinates[, 2],
    pch = 20,
    cex = point_size,
    col = pings$quality_class
  )
  if( length(weight_labels) > 0 ) {
    text(
      x = pings$time[weight_labels],
      y = ping_coordinates[weight_labels, 2],
      labels = ping_weights[weight_labels],
      pos = 1
    )
  } else {}
  lines(
    x = track$time,
    y = track_coordinates[, 2]
  )
  if( !missing(y_se) ) {
    lines(
      x = track$time,
      y = track_coordinates[, 2] + 1.96 * y_se,
      lty = 2
    )
    lines(
      x = track$time,
      y = track_coordinates[, 2] - 1.96 * y_se,
      lty = 2
    )
  } else {}

  par(original_par)
  return(invisible())
}

#' Plot a track
#'
#' @param track A list with a pings sf data.frame and a track sf data.frame.
#' @param ping_weights An optional vector of score-level weights for each ping
#' @param lower_weight_limit Pings with weight below this limit will not be plotted.
#' @param point_size_scale Scaling factor for point sizes
#' @param return_type Either "plot" to return a plot, or "list" to return a list of plot elements.
#'
#' @return A tmap object
#'
#' @export
plot_track<- function(track, ping_weights, lower_weight_limit = 0.01, point_size_scale = 1, return_type = "plot") {
  has_tmap<- requireNamespace("tmap", quietly = TRUE)
  if( !has_tmap ) {
    stop("Must have package tmap installed to use plot_track.")
  } else {}
  pings<- track$pings
  track<- track$track
  point_size<- 0.2 * log(loc_class_K[, 2] + 3)

  if( missing(ping_weights) ) {
    ping_weights<- rep(1, nrow(pings))
  } else {}
  pings<- pings[ping_weights >= lower_weight_limit, ]
  ping_weights<- ping_weights[ping_weights >= lower_weight_limit]
  weight_labels<- which(ping_weights < 0.99)
  ping_weights<- round(ping_weights, digits = 2)
  pings$weight<- ping_weights

  bounding_box<- st_as_sfc(st_bbox(
    c(
      st_as_sfc(st_bbox(pings)),
      st_as_sfc(st_bbox(track))
    )
  ))

  track_linestring<- sf::st_cast(
    sf::st_combine(
      track
    ),
    "LINESTRING"
  )
  pings_linestring<- lapply(
    seq(nrow(pings)),
    function(i) {
      line<- c(
        sf::st_geometry(pings)[[i]],
        pings$estimate[[i]]
      )
      line<- sf::st_cast(line, "LINESTRING")
      return(line)
    }
  )
  pings_linestring<- do.call(sf::st_sfc, pings_linestring)
  pings_linestring<- sf::st_sf(
    quality_class = pings$quality_class,
    geometry = pings_linestring
  )
  sf::st_crs(pings_linestring)<- sf::st_crs(pings)


  track_plot<- tmap::tm_shape(track_linestring) + tmap::tm_lines()
  pings_plot<- tmap::tm_shape(pings) + tmap::tm_bubbles(
    size = "quality_class",
    size.scale = tmap::tm_scale_ordinal(
      values = point_size,
      values.scale = point_size_scale
    ),
    size.legend = tmap::tm_legend_combine(
      "fill"
    ),
    fill = "quality_class",
    fill.scale = tmap::tm_scale_ordinal(
      values = palette()
    ),
    fill.legend = tmap::tm_legend(
      title = "Quality Class"
    )
  )
  if( length(weight_labels) > 0 ) {
    ping_weights_plot<- tmap::tm_shape(pings[weight_labels, ]) + tmap::tm_text(
      text = "weight"
    )
  } else {
    ping_weights_plot<- NULL
  }
  pings_linestring_plot<- tmap::tm_shape(pings_linestring) + tmap::tm_lines(
    col = "quality_class",
    col.scale = tmap::tm_scale_ordinal(
      values = palette()
    ),
    col.legend = tmap::tm_legend(
      show = FALSE
    ),
    col_alpha = 0.3
  )
  bbox_plot<- tmap::tm_shape(bounding_box, is.main = TRUE) + tmap::tm_borders(
    col_alpha = 0
  )
  if( return_type == "list" ) {
    return(
      list(
        pings_linestring = pings_linestring_plot,
        pings = pings_plot,
        ping_weights = ping_weights_plot,
        track = track_plot,
        bounding_box = bounding_box
      )
    )
  } else {}
  return(
    pings_linestring_plot + pings_plot + ping_weights_plot + track_plot + bbox_plot
  )
}

#' Plot the coordinate profiles of a track
#'
#' @param track A list with a pings sf data.frame and a track sf data.frame.
#' @param ping_weights An optional vector of score-level weights for each ping
#' @param lower_weight_limit Pings with weight below this limit will not be plotted.
#' @param posterior_tracks The output of sample_posterior_tracks()
#' @param alpha Transparency value for posterior samples. 0 is invisible, 1 is opaque.
#' @param point_size_scale Scaling factor for point sizes
#' @param which_sample Index vector of which samples to plot
#'
#' @return A tmap object
#'
#' @export
plot_posterior_track_profiles<- function(
  track,
  posterior_tracks,
  ping_weights,
  lower_weight_limit = 0.01,
  alpha = 0.2,
  point_size_scale = 1,
  which_samples = seq_along(posterior_tracks)
) {
  original_par<- par(
    mfrow  = c(2, 1),
    mar = c(1, 3, 2, 1)
  )
  pings<- track$pings
  track<- track$track

  if( missing(ping_weights) ) {
    ping_weights<- rep(1, nrow(pings))
  } else {}
  pings<- pings[ping_weights >= lower_weight_limit, ]
  ping_weights<- ping_weights[ping_weights >= lower_weight_limit]
  weight_labels<- which(ping_weights < 0.99)
  ping_weights<- round(ping_weights, digits = 2)

  ping_coordinates<- sf::st_coordinates(pings)
  track_coordinates<- sf::st_coordinates(track)
  all_coordinates<- rbind(ping_coordinates, track_coordinates)
  posterior_tracks<- posterior_tracks[which_samples]

  tlim<- c(min(track$time), max(track$time))
  xlim<- c(
    min(all_coordinates[, 1]),
    max(all_coordinates[, 1])
  )
  xlim<- xlim + c(-1, 1) * 0.05 * diff(xlim)
  ylim<- c(
    min(all_coordinates[, 2]),
    max(all_coordinates[, 2])
  )
  ylim<- ylim + c(-1, 1) * 0.05 * diff(ylim)
  point_size<- point_size_scale * 0.5 * log(loc_class_K[, 2] + 3)[pings$quality_class]

  plot(
    x = track$time,
    y = track_coordinates[, 1],
    xlim = tlim,
    ylim = xlim,
    main = "",
    ylab = "X coordinate",
    xlab = "Time",
    type = "n"
  )
  points(
    x = track$time[pings$index],
    y = ping_coordinates[, 1],
    pch = 20,
    cex = point_size,
    col = pings$quality_class
  )
  if( length(weight_labels) > 0 ) {
    text(
      x = track$time[pings$index[weight_labels]],
      y = ping_coordinates[weight_labels, 1],
      labels = ping_weights[weight_labels],
      pos = 1
    )
  } else {}
  lines(
    x = track$time,
    y = track_coordinates[, 1]
  )
  for( i in seq_along(posterior_tracks) ) {
    lines(
      x = track$time,
      y = sf::st_coordinates(posterior_tracks[[i]])[, 1],
      col = rgb(0, 0, 0, alpha)
    )
  }
  legend(
    x = "topright",
    legend = c("G", 3, 2, 1, 0, "A", "B"),
    col = seq(7),
    pt.cex = 2 * log(loc_class_K[, 2] + 1.5),
    pch = 20
  )

  plot(
    x = track$time,
    y = track_coordinates[, 2],
    xlim = tlim,
    ylim = ylim,
    main = "",
    ylab = "Y coordinate",
    xlab = "Time",
    type = "n"
  )
  points(
    x = track$time[pings$index],
    y = ping_coordinates[, 2],
    pch = 20,
    cex = point_size,
    col = pings$quality_class
  )
  if( length(weight_labels) > 0 ) {
    text(
      x = track$time[pings$index[weight_labels]],
      y = ping_coordinates[weight_labels, 2],
      labels = ping_weights[weight_labels],
      pos = 1
    )
  } else {}
  lines(
    x = track$time,
    y = track_coordinates[, 2]
  )
  for( i in seq_along(posterior_tracks) ) {
    lines(
      x = track$time,
      y = sf::st_coordinates(posterior_tracks[[i]])[, 2],
      col = rgb(0, 0, 0, alpha)
    )
  }

  par(original_par)
  return(invisible())
}

#' Plot an estimates track with posterior samples
#'
#' @param track A list with a pings sf data.frame and a track sf data.frame.
#' @param ping_weights An optional vector of score-level weights for each ping
#' @param lower_weight_limit Pings with weight below this limit will not be plotted.
#' @param posterior_tracks The output of sample_posterior_tracks()
#' @param alpha Transparency value for posterior samples. 0 is invisible, 1 is opaque.
#' @param point_size_scale Scaling factor for point sizes
#' @param which_sample Index vector of which samples to plot
#'
#' @return A tmap object
#'
#' @export
plot_posterior_tracks<- function(
  track,
  posterior_tracks,
  ping_weights,
  lower_weight_limit = 0.01,
  alpha = 0.2,
  point_size_scale = 1,
  which_samples = seq_along(posterior_tracks)
) {
  track_plot<- plot_track(
    track,
    ping_weights,
    lower_weight_limit,
    point_size_scale = point_size_scale,
    return_type = "list"
  )
  posterior_tracks<- posterior_tracks[which_samples]
  posterior_plots<- lapply(
    posterior_tracks,
    function(sample) {
      sample_linestring<- sf::st_cast(
        sf::st_combine(
          sample
        ),
        "LINESTRING"
      )
      sample_plot<- tmap::tm_shape(sample_linestring) + tmap::tm_lines(col_alpha = alpha)
      return(sample_plot)
    }
  )
  bboxes<- lapply(
    posterior_tracks,
    function(sample) {
      return(sf::st_as_sfc(sf::st_bbox(sample)))
    }
  )
  bounding_box<- st_as_sfc(st_bbox(
    c(
      track_plot$bounding_box,
      do.call(c, bboxes)
    )
  ))
  bbox_plot<- tmap::tm_shape(bounding_box, is.main = TRUE) + tmap::tm_borders(
    col_alpha = 0
  )
  return(
    bbox_plot + Reduce(`+`, c(posterior_plots, track_plot[1:3]))
  )
}
