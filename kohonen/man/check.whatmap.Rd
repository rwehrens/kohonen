\name{check.whatmap}
\alias{check.whatmap}
\title{Check the validity of a whatmap argument}
\description{Not meant to be called directly by the user.}
\usage{
check.whatmap(x, whatmap)
}
\arguments{
  \item{x}{A \code{kohonen} object, or a
    list of data matrices that can be used as input data for SOM functions.}
  \item{whatmap}{An indication of a subset of the data; either by naming
    the elements, or giving indices. If \code{whatmap} equals NULL, the
    selection of \code{x} is used if \code{x} is a \code{kohonen}
    object, or else no selection is performed.}
}
\value{
  Returns a numerical vector with the indices of the selected layers. An
  invalid selection leads to an error.
}
\author{Ron Wehrens}
\keyword{classif}
