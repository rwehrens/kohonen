## several checks for input data of the supersom and map functions

check.data <- function(data, maxNA.fraction) {
  ## Check whether data is a list of data matrices or factors
  if (!is.list(data) |
      !all(sapply(data, class) %in% c("numeric", "matrix", "factor")))
    stop("Argument data should be a list of numeric vectors or matrices, or factors")
  
  ## Convert vectors to one-column matrices in layers
  if (any (vector.idx <- sapply(data, is.vector, mode = "numeric"))) {
    vector.idx <- which(vector.idx)
    data[vector.idx] <- lapply(data[vector.idx],
                               function(x) matrix(x, ncol = 1))
  }
  if (any(factor.idx <- sapply(data, is.factor)))
    data[factor.idx] <- lapply(data[factor.idx], classvec2classmat)

  ## Check whether data is numeric
  if (!all(sapply(data, is.numeric)))
    stop("Argument data should be numeric")

  nobjects <- unique(sapply(data, nrow))
  if (length(nobjects) > 1)
    stop("Unequal numbers of objects in data list")
  
  data
}

## Objective: identify rows and columns with too many NA values. Data
## is a list of matrices.
check.data.na <- function(data, maxNA.fraction) {
  nacols <-
    lapply(data,
           function(x)
             which(apply(x, 2,
                         function(y)
                         (sum(is.na(y)) / length(y)) > maxNA.fraction)))
  narows <-
    lapply(data,
           function(x)
             which(apply(x, 1,
                         function(y)
                         (sum(is.na(y)) / length(y)) > maxNA.fraction)))
  narows <- unique(unlist(narows))

  list(narows, nacols)
}

remove.data.na <- function(data, nachecks) {
  narows <- nachecks[[1]]
  nacols <- nachecks[[2]]
  
  for (i in seq(along = data)) {
    if (length(narows) > 0) 
      data[[i]] <- data[[i]][-narows, , drop=FALSE]

    if (length(nacols[[i]]) > 0) 
      data[[i]] <- data[[i]][, -nacols[[i]], drop=FALSE]
  }
  
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
