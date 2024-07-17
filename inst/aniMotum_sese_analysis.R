library(trackfiltering)

# Load aniMotum to get an example dataset
# we also need to convert the data into
# an sf object. We'll also rename the columns
# and make  quality class a factor.
library(aniMotum)
pings<- as.data.frame(aniMotum::sese)
pings<- sf::st_as_sf(
    pings,
    coords = c("lat", "lon")   
)
sf::st_crs(pings)<- sf::st_crs("+proj=stere +lon_0=68 +units=km +datum=WGS84")
colnames(pings)[2:3]<- c("time", "quality_class")
pings<- subset(pings, id == "ct109-085-14" & quality_class != "Z")
pings$id<- NULL
pings$quality_class<- ordered(
  pings$quality_class,
  levels = c("G", "3", "2", "1", "0", "A", "B")
)

# The dataset provided by aniMotum was already cleaned by hand,
#   so we'll make some artificial outliers by randomly swapping
#   some of the time stamps around.
set.seed(140)
outlier_index<- sample(
    nrow(pings),
    size = 50,
    replace = FALSE
)
pings$time[sort(outlier_index)]<- pings$time[outlier_index]
pings<- pings[order(pings$time), ]

# Fit a model to the track using the fit_track function
#
# For this function to work, the column names of your data need to be
# the exact same as mine: "time" and "quality_class", although the name
# of the geometry column can be different.
#
# The track_time argument is optional, but the default values don't
# always work well. My usual advice is to pick a regular time step
# from the beginning of your track to the end of the track and pick
# as small a time step as you want. Smaller time steps give you higher
# resolution predictions at the cost of higher computational burden and
# possibly worse convergence issues.
#
# ping_robustness is a tuning parameter. Lower values will make the model
# more aggressive in labelling points as outliers. If you're seeing points
# that should be labelled as outliers but are not, make the value lower.
# If you're seeing too many points labelled as outliers, increase the value.
#
# ping_robust_code is a choice between different ways to find outliers.
# If you set it equal to 0 you get regular maximum likelihood without
#   removing any outliers
# If you set it equal to 1 or 2 you will remove outliers. I recommend 2.
track_fit<- trackfiltering::fit_track(
    pings = pings,
    track_time = seq(
      from = min(pings$time),
      to = max(pings$time),
      by = "2 hours"
    ),
    ping_robustness = 3,
    ping_robust_code = 2
)

# Check for convergence of the model
#
# If you see "relative convergence" or "X convergence" 
#   you're in good shape. If you see "false convergence"
#   you can try a larger time step for track_time, or by
#   trying different starting values with the starting_*
#   arguments or even fixing some parameters fixed at
#   their starting values with the fix_* arguments. See
#   help("fit_track") for all the arguments you can give.
track_fit$opt$message

# You can get a copy of the supplied pings back. It will have
#   an additional geometry column giving the model-based estimate
#   of the true location at each observed ping.
track_fit$pings

# You can get an interpolated track that gives the estimated
#   location at every time included in the track_time argument.
#   The interpolated track here also gives you estimates of the
#   standard errors / uncertainty in the estimate. If you want
#   standard errors for your estimates of the observed locations,
#   you should add the observed times to the track_time argument.
track_fit$track

# There are also some basic plotting functions. You need the "tmap"
#   package installed to use the plot_track function, but
#   the plot_track_profile function doesn't need any additional packages.
#   This particular track is very vertical, so you may get an error
#   about the margins being too large when using plot_track.
#
# The lower_weight_limit argument is a cutoff such that any points labelled
#   as an outlier with weight less than that won't be plotted.
plot_track_profile(
    track_fit,
    x_se = track_fit$track$x_se,
    y_se = track_fit$track$y_se,
    ping_weights = track_fit$obj$report()$ping_weights,
    lower_weight_limit = 0
)
plot_track(
    track_fit,
    ping_weights = track_fit$obj$report()$ping_weights,
    lower_weight_limit = 0
)



