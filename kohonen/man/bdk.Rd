\name{bdk}
\alias{bdk}
\title{Supervised version of Kohonen's self-organising maps}
\description{
  Supervised version of self-organising maps for mapping
  high-dimensional spectra or patterns to 2D: the Bi-Directional Kohonen
  map. This is an alternating training of the X-space and the Y-space
  of the map, where for updating the X-space more weight is given to the
  features in Y, and vice versa. Weights start by default with values of
  (0.75, 0.25) and during training go to (0.5, 0.5). Prediction is done
  only using the X-space. For continuous Y, the Euclidean distance is
  used; for categorical Y the Tanimoto distance. 
}

\usage{
bdk(data, Y, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
    radius = quantile(nhbrdist, 0.67) * c(1, -1), xweight = 0.75,
    contin, toroidal = FALSE, n.hood, keep.data = TRUE)
}

\arguments{
  \item{data}{a matrix, with each row representing an object.}
  \item{Y}{property that is to be modelled. In case of classification, Y
    is a matrix with exactly one '1' in each row indicating
    the class, and zeros elsewhere. For prediction of continuous
    properties, Y is a vector. A combination is possible, too, but one
    then should take care of appropriate scaling.}
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
  \item{xweight}{the initial weight given to the X map in the
    calculation of distances for updating Y, and to the Y map for
    updating X. This will linearly go to 0.5 during training. Defaults
    to 0.75.}
  \item{contin}{parameter indicating whether Y is continuous or
    categorical. The default is to check whether all row sums of Y equal
    1: in that case \code{contin} is \code{FALSE}.}
  \item{toroidal}{if \code{TRUE}, the edges of the map are joined. Note
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
  \item{method}{the type of som, here "bdk"}
}

\references{W.J. Melssen, R. Wehrens, and L.M.C. Buydens.
Chemom. Intell. Lab. Syst., 83, 99-113 (2006).}
\author{Ron Wehrens}

\seealso{\code{\link{som}}, \code{\link{xyf}},
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

bdk.wines <- bdk(Xtraining,
                 factor(wine.classes[training]),
                 grid = somgrid(5, 5, "hexagonal"))

bdk.prediction <- predict(bdk.wines, newdata=Xtest)
table(wine.classes[-training], bdk.prediction$prediction)
}
\keyword{classif}
