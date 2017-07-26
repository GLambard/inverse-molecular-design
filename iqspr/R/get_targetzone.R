#'get a targeted zone
#' @description Get a temporary targeted space of properties (a.k.a physico-chemical properties here) to populate.
#' This space can be on the path to reach the final desired target, or the final targeted space by itself. This depends on
#' how the temporary targeted space fits the population of observables constrains delimited by the parameters min.points and
#' max.points described below.
#' @param targ.min is a d-dimensional vector of minimum values delimiting the targeted space, where d is the number of properties.
#' @param targ.max is a d-dimensional vector of maximum values delimiting the targeted space, where d is the number of properties.
#' @param feat is a n rows x d columns matrix, where d and n are the number of properties and datapoints respectively.
#' @param min.points minimum number of points required in the targeted space of properties.
#' @param max.points maximum number of points required in the targeted space of properties.
#'
#' @examples
#' \dontrun{data(qspr.data)
#' ty <- as.matrix(qspr.data[,c(2,5)])
#' targ.zone <- get_targetzone(targ.min=c(50,9),targ.max=c(200,11),
#'                              feat=ty,
#'                              min.points=100,max.points=150)}
#'
#' @return the boundaries of a targeted zone.
#'
#' @import stats
#' @importFrom stats density
#'
#' @export get_targetzone
get_targetzone <- function(targ.min=c(50,9),
                           targ.max=c(200,11),
                           feat=NULL,
                           min.points=100,
                           max.points=150){

  # number of properties
  D = dim(feat)[2]

  # maximum densities coordinates for each property, initialisation
  max.dens <- vector()

  # parameters of the global target zone
  # half length of the box
  targ.dim = (targ.max - targ.min) * 0.5

  # center coordinates
  targ.cent = targ.min + targ.dim
  featcut <- feat

  # calculate density distribution and
  # maxima's coordinates for each property
  # cut off unnecessary data off the feat to form featcut
  for (d in 1:D) {
    feat.dens <- density(feat[, d])
    max.dens[d] = feat.dens$x[which(feat.dens$y == max(feat.dens$y))]
  }

  # length of the segment between
  # the maximum density location and the final targeted space center
  segm = targ.cent - max.dens
  npoints = 0
  n = 1
  targ.zone = NULL

  repeat {
    # every next step is twice shorter than the previous
    n = n * 0.5

    # eternal loop prevention
    if (n < 1.e-6) {
      cat("The step is dangerously small")
      break
    }

    # if number of points is inside of the desired space, finish
    if (npoints <= max.points & npoints >= min.points) {
      targ.zone = cbind(temp.targ.min, temp.targ.max)
      break
    } else if (npoints > max.points) {
      # if more points than needed, move towards the final targeted space
      targ.cent = targ.cent + segm * n

      temp.targ.min = targ.cent - targ.dim
      temp.targ.max = targ.cent + targ.dim
    } else if (npoints < min.points) {
      # if less points than needed, move towards maximum density
      targ.cent = targ.cent - segm * n

      temp.targ.min = targ.cent - targ.dim
      temp.targ.max = targ.cent + targ.dim
    }

    # count the number of points inside of the temporary space
    temp.p <- t(t(featcut) < temp.targ.max & t(featcut) > temp.targ.min)
    temp.p <- apply(temp.p, 1, function(row) all(row != 0))
    npoints = sum(temp.p)

  }
  return(targ.zone)
}
