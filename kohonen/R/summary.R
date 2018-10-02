summary.kohonen <- function(object, ...)
{
  cat("SOM of size ", object$grid$xdim, "x", object$grid$ydim,
      " with a ", object$grid$topo,
      if (object$grid$toroidal) "toroidal", " topology and a ",
      as.character(object$grid$neighbourhood.fct),
      " neighbourhood function.", sep="")

  cat("\nThe number of data layers is ", length(object$codes), ".", sep = "")
  cat("\nDistance measure(s) used: ",
      paste(object$dist.fcts, collapse = ", "), ".", sep = "")
  
  if (!is.null(object$data)) {
    cat("\nTraining data included:",
        nrow(object$data[[1]]), "objects.")
    if (length(object$data) > length(object$whatmap))
      cat(", of which", length(object$whatmap),
          ifelse(length(object$whatmap) > 1, "have", "has"),
          "been used in training.")
    cat("\nMean distance to the closest unit in the map: ",
        round(mean(object$distances, na.rm = TRUE), 3), ".", sep = "")
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
