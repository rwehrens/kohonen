\name{check.whatmap}
\alias{check.whatmap}
\title{Check the validity of a whatmap argument}
\description{Not meant to be called directly by the user.}
\usage{
check.whatmap(x, whatmap)
}
\arguments{
  \item{x}{Either a \code{kohonen} object from \code{supersom}, or a
    list of data matrices that can be used as input data for \code{supersom}.}
  \item{whatmap}{An indication of a subset of the data; either by naming
    the elements, or giving indices. If \code{whatmap} equals NULL, no
    selection is performed.}
}
\value{
  Returns a numerical vector with the indices of the selected layers.
}
\author{Ron Wehrens}
\keyword{classif}
