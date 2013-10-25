\name{som}
\alias{som}
\title{Kohonen's self-organising maps}
\description{
  Self-organising maps for mapping high-dimensional spectra or patterns
  to 2D; Euclidean distance is used. Modelled after the SOM function in
  package \code{class}.
}

\usage{
som(data, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
    radius = quantile(nhbrdist, 0.67) * c(1, -1), init,
    toroidal = FALSE, n.hood, keep.data = TRUE)
}

\arguments{
  \item{data}{a matrix, with each row representing an object.}
  \item{grid}{a grid for the representatives: see 'somgrid'.}
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
  \item{init}{the initial representatives, represented as a matrix.  If
    missing, chosen (without replacement) randomly from 'data'.}
  \item{toroidal}{if TRUE, the edges of the map are joined. Note
    that in a hexagonal toroidal map, the number of rows must be even.}
  \item{n.hood}{the shape of the neighbourhood, either "circular" or
    "square". The latter is the default for rectangular maps, the former
    for hexagonal maps.}
  \item{keep.data}{save data in return object.}
}

\value{
  an object of class "kohonen" with components
  \item{data}{data matrix, only returned if \code{keep.data == TRUE}.}
  \item{grid}{the grid, an object of class "somgrid".}
  \item{codes}{a matrix of code vectors.}
  \item{changes}{vector of mean average deviations from code vectors.}
  \item{unit.classif}{winning units for all data objects,
    only returned if \code{keep.data == TRUE}.}
  \item{distances}{distances of objects to their corresponding winning
    unit, only returned if \code{keep.data == TRUE}.}
  \item{toroidal}{whether a toroidal map is used.}
  \item{method}{the type of som, here "som".}
}

\references{
  "Self-organizing maps", 3rd Ed., T. Kohonen, New York: Springer (2001) 
}
\author{Ron Wehrens}

\seealso{\code{\link{xyf}}, \code{\link{bdk}}, \code{\link{plot.kohonen}}}
\examples{
data(wines)
set.seed(7)

training <- sample(nrow(wines), 120)
Xtraining <- scale(wines[training, ])
Xtest <- scale(wines[-training, ],
               center = attr(Xtraining, "scaled:center"),
               scale = attr(Xtraining, "scaled:scale"))

som.wines <- som(Xtraining, grid = somgrid(5, 5, "hexagonal"))

som.prediction <- predict(som.wines, newdata = Xtest,
          trainX = Xtraining,
          trainY = factor(wine.classes[training]))
table(wine.classes[-training], som.prediction$prediction)
}
\keyword{classif}
