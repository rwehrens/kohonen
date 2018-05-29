object.distances <- function(kohobj, type = c("data", "codes"), whatmap) {
  if (missing(whatmap)) {
    whatmap <- kohobj$whatmap
  } else {
    whatmap <- check.whatmap(kohobj, whatmap)
  }
  
  weights <- kohobj$user.weights[whatmap] * kohobj$distance.weights[whatmap]
  maxNA.fraction <- kohobj$maxNA.fraction
  distanceFunctions <- kohobj$dist.fcts[whatmap]
  dist.ptrs <- getDistancePointers(distanceFunctions,
                                   maxNA.fraction = maxNA.fraction)
  
  type <- match.arg(type)
  data <- kohobj[[type]][whatmap]
  if (any(factor.idx <- sapply(data, is.factor)))
    data[factor.idx] <- lapply(data[factor.idx], classvec2classmat)

  nvars <- sapply(data, ncol)
  nobjects <- nrow(data[[1]])

  nNA <- getnNA(data, maxNA.fraction, nobjects)

  datamat <- matrix(unlist(data), ncol = nobjects, byrow = TRUE)
##  datamat <- matrix(unlist(data[whatmap]), ncol = nobjects, byrow = TRUE)
  res <- ObjectDistances(data = datamat,
                         numVars = nvars,
                         numNAs = nNA,
                         distanceFunctions = dist.ptrs,
                         weights = weights)

  attributes(res) <- list("Size" = nobjects, "Diag" = FALSE, "Upper" = FALSE,
                          "method" = "supersom", "call" = match.call(),
                          "class" = "dist")
  
  res
}
