\name{object.distances}
\alias{object.distances}
\title{Calculate distances between object vectors in a SOM}
\description{
This function calculates the distance between objects using the distance
functions, weights and other attributes of a trained SOM. This function
is used in the calculation of the U matrix in function
\code{plot.kohonen} using the \code{type = "dist.neighbours" argument.}
}
\usage{
object.distances(kohobj, type = c("data", "codes"), whatmap)
}
\arguments{
  \item{kohobj}{An object of class \code{kohonen}.}
  \item{type}{Whether to calculate distances between the data objects,
    or the codebook vectors.}
  \item{whatmap}{What data layers to use. If unspecified the data layers
    defined in the kohonen object are used.}
}
\value{An object of class \code{dist}, which can be directly fed into
  (e.g.) a hierarchical clustering.}
\references{R. Wehrens and J. Kruisselbrink, submitted, 2017.}
\author{Ron Wehrens}
\seealso{
  \code{\link{unit.distances}}, \code{\link{supersom}}
}
\examples{
data(wines)
set.seed(7)
sommap <- supersom(list(measurements = scale(wines),
                        vintages = vintages),
                   grid = somgrid(6, 4, "hexagonal"))
obj.dists <- object.distances(sommap, type = "data")
code.dists <- object.distances(sommap, type = "codes")
}
\keyword{classif}
