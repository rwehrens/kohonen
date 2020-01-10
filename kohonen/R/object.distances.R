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

## layer.distances calculates for each unit the average distance of
## objects mapped to that unit and the corresponding codebook vector.

layer.distances <- function(kohobj, whatmap, data, classif = NULL) {
  if (is.null(classif)) {
    if (is.null(kohobj$unit.classif)) {
      stop("No classification information present")
    } else {
      classif <- kohobj$unit.classif
    }
  }
  
  d2wus <- dist2WU(kohobj, whatmap, data, classif = classif)
  aggdists <- aggregate(d2wus, list(classif), mean)
  
  unit.distances <- rep(NA, nunits(kohobj))
  unit.distances[aggdists[,1]] <- aggdists[,2]

  unit.distances
}

## dist2WU returns the distances of objects to the
## corresponding winning units. If whatmap equals the kohobj whatmap
## element this is already available in slot unit.distances (unless no
## data were saved withthe object) - this function will usually be
## called with a single layer as whatmap argument.

dist2WU <- function(kohobj, whatmap, data, classif = NULL) {
  if (is.null(classif)) {
    if (is.null(kohobj$unit.classif)) {
      stop("No classification information present")
    } else {
      classif <- kohobj$unit.classif
    }
  }

  if (missing(data)) {
    if (!is.null(kohobj$data)) {
      data <- kohobj$data
    } else {
      stop("No data present")
    }
  }

  if (missing(whatmap)) {
    whatmap <- kohobj$whatmap
  } else {
    whatmap <- check.whatmap(kohobj, whatmap)
  }
  if (all(whatmap == kohobj$whatmap) &
      !is.null(kohobj$distances))
    return(kohobj$distances)

  ## the following statement is not strictly necessary, but it is hard
  ## to see the use of a whatmap argument other than a single layer,
  ## or all layers used in training.
  if (length(whatmap) > 1 & !all(whatmap == kohobj$whatmap))
    stop("Incorrect whatmap argument")
  
  weights <- kohobj$user.weights[whatmap] * kohobj$distance.weights[whatmap]
  maxNA.fraction <- kohobj$maxNA.fraction
  distanceFunctions <- kohobj$dist.fcts[whatmap]
  dist.ptrs <- getDistancePointers(distanceFunctions,
                                   maxNA.fraction = maxNA.fraction)
  
  data <- data[whatmap]
  codes <- kohobj$codes[whatmap]
  if (any(factor.idx <- sapply(data, is.factor)))
    data[factor.idx] <- lapply(data[factor.idx], classvec2classmat)

  nvars <- sapply(data, ncol)
  nobjects <- nrow(data[[1]])

  nNA <- getnNA(data, maxNA.fraction, nobjects)

  datamat <- matrix(unlist(data), ncol = nobjects, byrow = TRUE)
  codemat <- matrix(unlist(codes), ncol = nunits(kohobj), byrow = TRUE)
  
  LayerDistances(data = datamat,
                 codes = codemat,
                 uclassif = kohobj$unit.classif - 1, # C numbering
                 numVars = nvars,
                 numNAs = nNA,
                 distanceFunctions = dist.ptrs,
                 weights = weights)
}
