\name{xyf}
\alias{xyf}
\title{Supervised version of Kohonen's self-organising maps}
\description{
  Supervised version of self-organising maps for mapping
  high-dimensional spectra or patterns to 2D. The name stands for X-Y
  fused SOMs. One vector for each object is created by concatenating X
  and Y, and a SOM is trained in the usual way, with one exception: the
  distance of an object to a unit is the sum of separate distances for X
  and Y spaces. Prediction is done only using the X-space. For
  continuous Y, the Euclidean distance is used; for categorical Y the
  Tanimoto distance.
}

\usage{
xyf(data, Y, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
    radius = quantile(nhbrdist, 0.67) * c(1, -1), xweight = 0.5,
    contin, toroidal = FALSE, n.hood, keep.data = TRUE)
}

\arguments{
  \item{data}{a matrix, with each row representing an object.}
  \item{Y}{property that is to be modelled. In case of classification, Y
    is a matrix of zeros, with exactly one '1' in each row indicating
    the class. For prediction of continuous properties, Y is a vector. A
    combination is possible, too, but one then should take care of
    appropriate scaling.}
  \item{grid}{a grid for the representatives: see \code{somgrid}.}
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
  \item{xweight}{the weight given to the X map in the
    calculation of distances for updating Y. Default is 0.5.}
  \item{contin}{parameter indicating whether Y is continuous or
    categorical. The default is to check whether all row sums of Y equal
    1: in that case \code{contin} is \code{FALSE}.}
  \item{toroidal}{if TRUE, the edges of the map are joined. Note
    that in a hexagonal toroidal map, the number of rows must be even.}
  \item{n.hood}{the shape of the neighbourhood, either "circular" or
    "square". The latter is the default for rectangular maps, the former
    for hexagonal maps.}
  \item{keep.data}{save data in return value.}
}

\value{
  an object of class "kohonen" with components
 \item{data}{data matrix, only returned if \code{keep.data == TRUE}.}
  \item{Y}{Y, only returned if \code{keep.data == TRUE}.}
  \item{contin}{parameter indicating whether Y is continuous or
    categorical.}
  \item{grid}{the grid, an object of class "somgrid".}
  \item{codes}{list of two matrices, containing codebook vectors for X
    and Y, respectively.}
  \item{changes}{matrix containing two columns of mean average
    deviations from code vectors. Column 1 contains deviations used for
    updating Y; column 2 for updating X.}
  \item{toroidal}{whether a toroidal map is used.}
  \item{unit.classif}{winning units for all data objects,
    only returned if \code{keep.data == TRUE}.}
  \item{distances}{distances of objects to their corresponding winning
    unit, only returned if \code{keep.data == TRUE}.}
   \item{method}{the type of som, here "xyf".}
}

\references{W.J. Melssen, R. Wehrens, and L.M.C. Buydens.
Chemom. Intell. Lab. Syst., 83, 99-113 (2006).}
\author{Ron Wehrens}

\seealso{\code{\link{som}}, \code{\link{bdk}},
  \code{\link{plot.kohonen}}, \code{\link{predict.kohonen}}}
\examples{
### Wine example
data(wines)
set.seed(7)

training <- sample(nrow(wines), 120)
Xtraining <- scale(wines[training,])
Xtest <- scale(wines[-training,],
               center = attr(Xtraining, "scaled:center"),
               scale = attr(Xtraining, "scaled:scale"))

xyf.wines <- xyf(Xtraining,
                 factor(wine.classes[training]),
                 grid = somgrid(5, 5, "hexagonal"))

xyf.prediction <- predict(xyf.wines, newdata=Xtest)
table(wine.classes[-training], xyf.prediction$prediction)
}

\keyword{classif}
