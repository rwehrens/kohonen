check.whatmap <- function(x, whatmap)
{
  whatmap <- unique(whatmap)

  checkpar <- NULL
  if (inherits(x, "kohonen")) {
    checkpar <- getCodes(x)
  } else {
    if (is.list(x)) # not foolproof!
      checkpar <- x
  }
  if (is.null(checkpar))
    stop("no possibility to check argument 'whatmap'!")
  
  if (is.null(whatmap)) {
    if (is.null(x$whatmap)) {
      return(1:length(checkpar))  # no selection, return all layers
    } else {
      return(x$whatmap)
    }
  }

  if (is.numeric(whatmap) && all(whatmap %in% 1:length(checkpar)))
    return(sort(whatmap))

  if (is.character(whatmap)) {
    idx <- match(whatmap, names(checkpar)) ## works also when comparing to NULL
    if (any(!is.na(idx))) 
      return(sort(idx))
  }

  stop("incorrect whatmap argument") # invalid selection
}
