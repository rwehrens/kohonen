summary.kohonen <- function(object, ...)
{
  cat("SOM of size ", object$grid$xdim, "x", object$grid$ydim,
      " with a ", object$grid$topo,
      if (object$grid$toroidal) "toroidal", " topology and a ",
      as.character(object$grid$neighbourhood.fct),
      " neighbourhood function.", sep="")

  if (!is.null(object$data)) {
    cat("\nTraining data included of",
        nrow(object$data[[1]]), "objects")
    cat("\nThe number of layers is", length(object$data))
    if (length(object$data) > length(object$whatmap))
      cat(", of which", length(object$whatmap),
          ifelse(length(object$whatmap) > 1, "have", "has"),
          "been used in training.")
    cat("\nMean distance to the closest unit in the map:",
        mean(object$distances, na.rm = TRUE))
  } else {
    cat("\nNo training data included in the object.")
  }
  
  cat("\n")
  
  invisible()
}

print.kohonen <- function(x, ...)
{
  cat("SOM of size ", x$grid$xdim, "x", x$grid$ydim,
      " with a ", x$grid$topo, if (x$grid$toroidal) " toroidal",
      " topology.", sep="")
  if (!is.null(x$data))
    cat("\nTraining data included.")
  cat("\n")
  
  invisible()
}
