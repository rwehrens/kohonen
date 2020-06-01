"som" <- function (X, ...)
{
  supersom(list(X), ...)
}

"nunits" <- function(kohobj) {
  nrow(kohobj$grid$pts)
}
