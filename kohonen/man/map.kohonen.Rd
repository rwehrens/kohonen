\name{map.kohonen}
\alias{map.kohonen}
\alias{map}
\title{Map data to a supervised or unsupervised SOM}
\description{
  Map a data matrix onto a trained SOM.
}
\usage{
\method{map}{kohonen}(x, newdata, whatmap = NULL, user.weights = NULL,
                      maxNA.fraction = x$maxNA.fraction, ...)
}
\arguments{
  \item{x}{An object of class \code{kohonen}.}
  \item{newdata}{list of data matrices (numerical) of factors, equal to
    the \code{data} argument of the \code{supersom} function. No
    \code{data.frame} objects are allowed.}
  \item{whatmap, user.weights, maxNA.fraction}{parameters that usually will
    be taken from the \code{x} object, but can be supplied by the user
    as well. Note that it is not possible to change distance functions
    from the ones used in training the map.
    See \code{\link{supersom}} for more information.}
  \item{\dots}{Currently ignored.}
}
\value{
  A list with elements
  \item{unit.classif}{a vector of units that are closest to the objects
    in the data matrix.}
  \item{distances}{distances of the objects to the closest
    units. Distance measures are the same ones used in training the map.}
  \item{whatmap, user.weights}{Values used for these arguments.}
}
\author{Ron Wehrens}
\seealso{\code{\link{predict.kohonen}}, \code{\link{supersom}}}
\examples{
data(wines)
set.seed(7)

training <- sample(nrow(wines), 150)
Xtraining <- scale(wines[training, ])
somnet <- som(Xtraining, somgrid(5, 5, "hexagonal"))

map(somnet,
    scale(wines[-training, ],
          center=attr(Xtraining, "scaled:center"),
          scale=attr(Xtraining, "scaled:scale")))
}
\keyword{classif}
