\name{supersom}
\alias{supersom}
\title{Super-organising maps}
\description{
  An extension of \code{xyf} maps to multiple data layers, possibly with
  different numbers of variables (though equal numbers of objects). NAs
  are allowed (see below). A weighted distance over all layers is
  calculated to determine the winning units during training.
}

\usage{
supersom(data, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
    radius = quantile(nhbrdist, 0.67) * c(1, -1), 
    contin, toroidal = FALSE, n.hood, whatmap = NULL, weights = 1,
    maxNA.fraction = .5, keep.data = TRUE)
}

\arguments{
  \item{data}{list of data matrices.}
  \item{grid}{a grid for the representatives: see \code{\link[class]{somgrid}}.}
  \item{rlen}{the number of times the complete data set will be
    presented to the network.}
  \item{alpha}{learning rate, a vector of two numbers indicating the
    amount of change. Default is to decline linearly from 0.05 to 0.01
    over \code{rlen} updates.}
  \item{radius}{the radius of the neighbourhood, either given as a
    single number or a vector (start, stop). If it is given as a single
    number the radius will run from the given number to the negative
    value of that number; as soon as the neighbourhood gets smaller than
    one only the winning unit will be updated. The default is to start with a
    value that covers 2/3 of all unit-to-unit distances.}
  \item{contin}{parameter indicating whether data are continuous or
    categorical, i.e. a logical vector. If only \code{TRUE} or
    \code{FALSE} is given, this is taken to hold for all elements in
    \code{data}. The default is to check whether row sums in the data
    matrices are equal to 1: in that case the corresponding
    \code{contin} element is \code{FALSE}.}
  \item{toroidal}{if TRUE, the edges of the map are joined. Note
    that in a hexagonal toroidal map, the number of rows must be even.}
  \item{n.hood}{the shape of the neighbourhood, either "circular" or
    "square". The latter is the default for rectangular maps, the former
    for hexagonal maps.}
  \item{whatmap}{For \code{supersom} maps: what layers to use in the
    mapping.}
  \item{weights}{the weights given to individual layers. Default is 1/n,
    with n the number of layers.}
  \item{maxNA.fraction}{the maximal fraction of values that may be NA to
    prevent the row or column to be removed.}
  \item{keep.data}{save data in return value.}
}

\value{
  an object of class "kohonen" with components
  \item{data}{data matrix, only returned if \code{keep.data == TRUE}.}
  \item{contin}{parameter indicating whether elements of \code{data}
    are continuous or categorical.}
  \item{na.rows}{indices of objects (rows) that are removed because at
    least one of the layers has to many NAs for these objects.}
  \item{unit.classif}{winning units for all data objects,
    only returned if \code{keep.data == TRUE}.}
  \item{distances}{distances of objects to their corresponding winning
    unit, only returned if \code{keep.data == TRUE}.}
  \item{grid}{the grid, an object of class \code{\link[class]{somgrid}}.}
  \item{codes}{a list of matrices containing codebook vectors.}
  \item{changes}{matrix of mean average deviations from code vectors;
    every map corresponds with one column.}
  \item{toroidal}{whether a toroidal map is used.}
  \item{n.hood}{the shape of the neighbourhood, either "circular" or
    "square". The latter is the default for rectangular maps, the former
    for hexagonal maps.}
  \item{weights}{For \code{supersom} maps: weights of layers uses in the
    mapping.}
  \item{whatmap}{For \code{supersom} maps: what layers to use in the
    mapping.}
  \item{method}{type of map, here "supersom".}
}

\references{R. Wehrens and L.M.C. Buydens, J. Stat. Softw. 21 (5), 2007}
\author{Ron Wehrens}

\seealso{\code{\link[class]{somgrid}}, \code{\link{plot.kohonen}}}
\examples{
data(yeast)
yeast.supersom <- supersom(yeast, somgrid(6, 6, "hexagonal"), whatmap = 3:6)
obj.classes <- as.integer(yeast$class)
colors <- c("yellow", "green", "blue", "red", "orange")
plot(yeast.supersom, type = "mapping", col = colors[obj.classes],
     pch = obj.classes, main = "yeast data", keepMargins = TRUE)
}

\keyword{classif}

