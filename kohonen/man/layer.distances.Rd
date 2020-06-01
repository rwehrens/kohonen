\name{layer.distances}
\alias{layer.distances}
\alias{dist2WU}
\title{Assessing distances to winning units}
\description{
  Given a trained SOM, distances of individual objects to their closest
  units may be calculated with function \code{dist2WU}. Aggregation on
  the unit level is obtained through the function
  \code{layer.distances}. The latter function is the workhorse for the
  "quality" plots in function \code{plot.kohonen}.
}
\usage{
layer.distances(kohobj, whatmap, data, classif = NULL)
dist2WU(kohobj, whatmap, data, classif = NULL)
}
\arguments{
  \item{kohobj}{A trained \code{kohonen} object. Data and mapping
    results should be included.}
  \item{whatmap}{What layers to take into account - default is to consider
    all layers used in training. Also single layers may be chosen. Note
    that although the underlying C code can also calculate results for
    any subset, currently subsets larger than one are forbidden.}
  \item{data}{Data to use - default is to use the data from the trained SOM.}
  \item{classif}{Classification vector, corresponding to the
    \code{unit.classif} element of a \code{kohonen} object. It can also
    be provided explicitly.}
}
\details{The results will be weighted using both the user weights and
  distance weights. Summing all the results for individual layers
  therefore would lead to the \code{unit.classif} vector of the
  \code{kohonen} object.}
\value{Function \code{dist2WU} returns a vector, representing for each
  object the distance to its winning unit. Function
  \code{layer.distances} returns (as a vector) for each unit the average
  distance of objects for which it is the winning unit.
}
\author{Ron Wehrens}
\seealso{
Quality plots from \code{\link{plot.kohonen}}.
}
\examples{
library(kohonen)
data(wines)
wines.sc <- scale(wines)
set.seed(7)
xyf.wines <- xyf(wines.sc, vintages, grid = somgrid(5, 5, "hexagonal"))
dist2WU(xyf.wines, whatmap = 1)
plot(xyf.wines, "quality", whatmap = 1)
plot(xyf.wines, "property",
     property = layer.distances(xyf.wines, whatmap = 1))
}
\keyword{manip}
\keyword{cluster}
