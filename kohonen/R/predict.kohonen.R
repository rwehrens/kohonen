## workflow:
## if unit.predictions not given, calculate them:
##    if trainingdata given: check with codes, determine whatmap
##    else trainingdata <- data, no checks required
##    map trainingdata to obtain unit.predictions
## else check unit.predictions
## check newdata with unit.predictions
## map newdata

predict.kohonen <- function(object,
                            newdata = NULL,
                            unit.predictions = NULL,
                            trainingdata = NULL,
                            whatmap = NULL,
                            threshold = 0,
                            maxNA.fraction = object$maxNA.fraction,
                            ...)
{
  codeNcols <- sapply(object$codes, ncol)
  ncodes <- nrow(object$codes[[object$whatmap[1]]])

  if (is.null(whatmap)) {
    whatmap <- object$whatmap
  } else {
    whatmap <- check.whatmap(object, whatmap)
  }
  
  if (!is.null(unit.predictions)) {
    if (!is.list(unit.predictions))
      unit.predictions <- list(unit.predictions)
    
    if (any(isFactorUnPred <- sapply(unit.predictions, is.factor)))
      unit.predictions[isFactorUnPred] <-
        lapply(unit.predictions[isFactorPred], classvec2classmat)
  } else {
    ## calculate unit.predictions from mapping the trainingdata
    if (is.null(trainingdata)) {
      trainingdata <- object$data
    } else {
      if (is.null(trainingdata))
        stop("Missing trainingdata argument, no data available in kohonen object either")
      trainingdata <- check.data(trainingdata)
    }
    
    whatmap.tr <- check.whatmap(trainingdata, whatmap)
    if (!all(whatmap.tr == whatmap))
      stop("Data layer mismatch between map and trainingdata")

    ## remove only rows that contain too many NAs in the whatmap layers. 
    narows <- check.data.na(trainingdata[whatmap], maxNA.fraction)
    trainingdata <- remove.data.na(trainingdata, narows)
    
    if (any(is.null(object$codes[whatmap])))
      stop("Attempt to map training data on the basis of data layers not present in the SOM")
    
    if (!checkListVariables(trainingdata[whatmap], codeNcols[whatmap]))
      stop("Number of columns of trainingdata do not match ",
           "codebook vectors")

    object$data <- trainingdata  # saves a check in map.kohonen
    mappingX <- map(object,
                    whatmap = whatmap,
                    maxNA.fraction = maxNA.fraction,
                    ...)$unit.classif
    
    ## now calculate unit averages for ALL layers in the training data
    unit.predictions.tmp <-
      lapply(trainingdata,
             function(trd)
               aggregate(x = trd,
                         by = list(mappingX),
                         FUN = mean, drop = FALSE))
    ## if some units are empty we could interpolate from neighbouring
    ## units but for now we return NA
    unit.predictions <- lapply(unit.predictions.tmp,
                               function(x) {
                                 varnames <- colnames(x)[-1]
                                 x <- as.matrix(x)
                                 res <- matrix(NA, ncodes, ncol(x) - 1)
                                 colnames(res) <- varnames
                                 res[x[,1],] <- x[, -1, drop=FALSE]
                                 res
                               })
  }

  ## ###############################################################
  ## map newdata
  if (is.null(newdata)) {
    newdata <- object$data
    if (is.null(newdata))
      stop("Missing newdata argument, no data available in kohonen object either")
    
    newrownames <- rownames(newdata[[1]])
    if (is.null(newrownames)) newrownames <- 1:nrow(newdata[[1]])
    nnewrows <- length(newrownames)

    newmapping <- object$unit.classif

    if (any(factorNew <- sapply(newdata, is.factor)))
      newdata[factorNew] <- lapply(newdata[factorNew], classvec2classmat)
  } else {
    newdata <- check.data(newdata)

    ## Data layers of newdata may be a subset of whatmap, but only it
    ## the names agree.
    if (is.null(newnames <- names(newdata))) {
      whatmap.new <- whatmap
    } else {
      ## we keep whatmap.new as a vector of strings, now. We can use
      ## that to index codes and trainingdata as well.
      whatmap.new <- intersect(newnames, names(object$codes)[whatmap])
      dummy <- check.whatmap(newdata, whatmap.new)
    }
    
    ## agreement between map and newdata. If no names are given we
    ## assume the layers are in the same order.
    if (!checkListVariables(newdata[whatmap.new], codeNcols[whatmap.new]))
      stop("Number of columns of newdata do not match codebook vectors")

    ## check all layers but only remove records from whatmap layers

    newrownames <- rownames(newdata[[1]])
    if (is.null(newrownames)) newrownames <- 1:nrow(newdata[[1]])
    nnewrows <- length(newrownames)

    narows <- check.data.na(newdata[whatmap.new],
                            maxNA.fraction = maxNA.fraction)
    newdata <- remove.data.na(newdata, narows)
    
    ## finally: calculate mapping of new data
    newmapping <- map(object,
                      newdata = newdata,
                      whatmap = whatmap.new,
                      ...)$unit.classif
  }

  if (length(narows) > 0) {
    rowidx <- (1:nnewrows)[-narows]
  } else {
    rowidx <- 1:nnewrows
  }
  nonNA <- rowidx[which(!is.na(newmapping))]
  predictions <- lapply(unit.predictions,
                        function(x) {
                          pred <- matrix(NA, nnewrows, ncol(x))
                          pred[nonNA,] <-
                            x[newmapping[nonNA],,drop=FALSE]
                          pred
                          })
  for (i in seq(along = predictions))
    dimnames(predictions[[i]]) <- list(newrownames,
                                       colnames(unit.predictions[[i]]))

  ## codebooks are always matrices, also when they really signify
  ## factors, so we have to either store that info in the kohonen
  ## object or check for it. The danger of the current solution is
  ## that we may convert things into factors that should not be
  ## converted; the advantage is that we don't need to store useless
  ## stuff.
  if (any(isFactorPred <- sapply(predictions, is.factor.matrix)))
    predictions[isFactorPred] <- lapply(predictions[isFactorPred],
                                        classmat2classvec,
                                        threshold = threshold)
      
  ## return the predictions associated to the map units closest to
  ## the new data points
  list(predictions = predictions,
       unit.classif = newmapping,
       unit.predictions = unit.predictions,
       whatmap = whatmap)
}


## aux functions to check dimensions of lists of vectors and matrices
checkListObjects <- function(mylist, targetlength) {
  lengths <- sapply(mylist,
                    function(x)
                      if (is.vector(x) | is.factor(x)) {
                        length(x)
                      } else {
                        nrow(x)
                      })
  if (missing(targetlength)) targetlength <- lengths[1]
  all(targetlength == lengths)
}

checkListVariables <- function(mylist, targetlength) {
  lengths <-
    sapply(mylist,
           function(x)
             if (is.factor(x)) {
               nlevels(x)
             } else {
               if (is.vector(x)) {
                 1
               } else {
                 if (is.matrix(x)) {
                   ncol(x)
                 } else {
                   stop("Data type not allowed: should be a matrix or a factor")
                 }
               }
             })
  
  if (missing(targetlength)) targetlength <- lengths[1]
  all(targetlength == lengths)
}
