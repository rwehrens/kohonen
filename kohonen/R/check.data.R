## several checks for input data of the supersom and map functions

check.data <- function(data) {
  if (!is.list(data) | is.data.frame(data))
    data <- list(data)

  ## Check whether data contains only vectors, matrices or factors
  vector.idx <- sapply(data, is.vector, mode = "numeric")
  matrix.idx <- sapply(data, is.matrix)
  factor.idx <- sapply(data, is.factor)
  
  if (!all(vector.idx | matrix.idx | factor.idx))
    stop("Argument 'data' should be a list of numeric vectors or matrices, or factors")
  
  ## Convert vectors to one-column matrices in layers
  if (any(vector.idx)) {
    vector.idx <- which(vector.idx)
    data[vector.idx] <- lapply(data[vector.idx],
                               function(x) matrix(x, ncol = 1))
  }
  if (any(factor.idx))
    data[factor.idx] <- lapply(data[factor.idx], classvec2classmat)

  ## Check whether data is numeric
  if (!all(sapply(data, is.numeric)))
    stop("Argument data should be numeric")

  nobjects <- unique(sapply(data, nrow))
  if (length(nobjects) > 1)
    stop("Unequal numbers of objects in data list")

  data
}

## Objective: identify rows with too many NA values in individual data
## layers. Data is a list of matrices. We cannot incorporate this in
## check.data because in the map.kohonen function we need to keep track of the
## records that have been removed, so narows is essential information.

check.data.na <- function(data, maxNA.fraction) {
  narows <-
    lapply(data,
           function(x)
             which(apply(x, 1,
                         function(y)
                         (sum(is.na(y)) / length(y)) > maxNA.fraction)))

  unique(unlist(narows))
}

remove.data.na <- function(data, narows) {
  for (i in seq(along = data)) {
    if (length(narows) > 0) 
      data[[i]] <- data[[i]][-narows, , drop=FALSE]
  }
  
  ## check to see if there are any empty data layers
  ## because of the maxNA.fraction
  if (0 %in% c(sapply(data, dim)))
    stop("Empty data layer - check maxNA.fraction argument")
  
  data
}

getnNA <- function(data, maxNA.fraction, nobjects) {
  if (maxNA.fraction > 0L) {
    t(sapply(data,
             function(x)
               apply(x, 1, function(y) sum(is.na(y)))))
  } else {
    matrix(0, length(data), nobjects)
  }
}


## we should allow for rows containing NA values
is.factor.matrix <- function(datamat, tolerance = 1e-8, completeThreshold = 5) {
  idx <- apply(datamat, 1, function(x) !any(is.na(x)))
  if (sum(idx) > completeThreshold) { 
    datamat <- datamat[idx,]
    
    is.matrix(datamat) &&
      min(datamat, na.rm = TRUE) >= 0 &&
      max(datamat, na.rm = TRUE) <= 1 &&
      all(abs(rowSums(datamat) - 1) < tolerance)
  } else {
    FALSE
  }
}
