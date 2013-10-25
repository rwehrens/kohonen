\name{topo.error}
\alias{topo.error}
\title{Measures of topographical error}
\description{
The function implements two forms of what is known as 'topographical
error', basically a measure of continuity of the mapping. This is done
by assessing the distances between units containing similar
patterns. See 'Details'.}
\usage{
topo.error(somobj, type = c("nodedist", "bmu"), data) 
}
\arguments{
  \item{somobj}{an object of class \code{"kohonen".}}
  \item{type}{which kind of topographical error to return. See 'Details'.}
  \item{data}{used in case \code{type = "bmu"} - the objects for which
    the best and second best matching units are calculated. Usually this
    is done for training data. If not explicitly given, the function
    will use any data present in \code{somobj}; if no data are available
    there either, an error message is given.}
}
\value{A number representing the topological error.}
\details{
  Two types of topographical error are implemented, depending on the
  value of the \code{type} argument:
  \enumerate{
  \item \code{"nodedist"}: the average distance, in terms of (x, y)
    coordinates in the map, between all pairs of most similar codebook
    vectors.
  \item \code{"bmu"}: the average distance, in terms of (x, y)
    coordinates in the map, between the best matching unit and the
    second best matching unit, for all data points.}
  In both cases, low values indicate a smooth mapping, with similar
  objects mapped to close-by units.
}
\references{"Topology preservation in self-organising maps",
  K. Kiviluoto, IEEE Int. Conf. on Neural Networks, pp. 294-299 (1996)

  "Self-organizing maps", 3rd Ed., T. Kohonen, New York: Springer (2001) 
}
\examples{
data(yeast)
## take only complete cases
X <- yeast[[3]][apply(yeast[[3]], 1, function(x) sum(is.na(x))) == 0,]
yeast.som <- som(X, somgrid(5, 8, "hexagonal"))

## quantization error:
mean(yeast.som$distances)
## topographical error measures:
topo.error(yeast.som, "nodedist")
topo.error(yeast.som, "bmu")
}
\keyword{classif}
