"map" <- function(x, ...)
{
  UseMethod("map")
}

map.kohonen <- function(x,
                        newdata,
                        whatmap = NULL,
                        user.weights = NULL,
                        maxNA.fraction = x$maxNA.fraction, ...)
{
  ## ##########################################################################
  ## Get relevant info from kohonen object
  codes <- x$codes
  nlayers <- length(codes)
  
  if (missing(newdata) & !is.null(x$data)) {
    newdata <- x$data
  } else {
    newdata <- check.data(newdata)
  }
  
  if (is.null(user.weights)) {
    user.weights <- x$user.weights
    useTrainingWeights <- TRUE
  } else {
    useTrainingWeights <- FALSE
  }
  if (length(user.weights) == 1) user.weights <- rep(1, nlayers)
 
  dist.ptrs <- getDistancePointers(x$dist.fcts,
                                   maxNA.fraction = maxNA.fraction)

  ## ##########################################################################
  ## Check whatmap
  if (is.null(whatmap)) {
    whatmap <- whatmap.tr <- x$whatmap
  } else {
    whatmap.tr <- check.whatmap(x, whatmap)
    whatmap <- check.whatmap(newdata, whatmap)
  }

  ## ##########################################################################
  ## Apply whatmap.
  newdata <- newdata[whatmap]
  nachecks <- check.data.na(newdata, maxNA.fraction = maxNA.fraction)
  newdata <- remove.data.na(newdata, nachecks)
            
  if (useTrainingWeights & any(user.weights[whatmap.tr] < 1e-8))
    warning("Mapping new data using data layers not involved in training")

  ## ##########################################################################
  ## Apply whatmap to other objects
  dist.ptrs <- dist.ptrs[whatmap.tr]
  codes <- codes[whatmap.tr]
  user.weights.orig <- user.weights
  user.weights <- user.weights[whatmap.tr]
  if (length(whatmap.tr) == 1) {
    user.weights <- 1
  } else {
    if (sum(user.weights >= 1e-8) == 0)
      stop("Only user.weights of zero given")
  }

  ## ##########################################################################
  ## Final preparations
  nvars <- sapply(codes, ncol)
  ncodes <- nrow(codes[[1]])
  nobjects <- nrow(newdata[[1]])
  
  nNA <- getnNA(newdata, maxNA.fraction, nobjects)

  newdata <- matrix(unlist(newdata), ncol=nobjects, byrow=TRUE)
  codes <- matrix(unlist(codes), ncol=ncodes, byrow=TRUE)

  weights <- user.weights * x$distance.weights[whatmap.tr]
  weights <- weights / sum(weights)

  ## ##########################################################################
  ## Go!
  res <- RcppMap(data = newdata,
                 numVars = nvars,
                 numNAs = nNA,
                 codes = codes,
                 weights = weights,
                 distanceFunctions = dist.ptrs)

  if (length(nachecks) > 0) {
    unit.classif <- distances <- rep(NA, nobjects + length(nachecks))
    unit.classif[-nachecks] <-  res$winners + 1
    distances[-nachecks] <- res$unitdistances

    list(unit.classif = unit.classif,
         distances = distances,
         whatmap = whatmap,
         user.weights = user.weights.orig)
  } else {
    list(unit.classif = res$winners + 1,
         distances = res$unitdistances,
         whatmap = whatmap,
         user.weights = user.weights.orig)
  }
}
