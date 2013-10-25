\name{map.kohonen}
\alias{map.kohonen}
\alias{map}
\title{Map data to a supervised or unsupervised SOM}
\description{
  Map a data matrix onto a trained SOM.
}
\usage{
\method{map}{kohonen}(x, newdata, whatmap = NULL, weights,
            scale.distances = (nmaps > 1), ...)
}
\arguments{
  \item{x}{A trained supervised or unsupervised SOM obtained from
    functions \code{som}, \code{xyf}, or \code{bdk}.}
  \item{newdata}{Data matrix, with rows corresponding to objects.}
  \item{whatmap}{For \code{supersom} maps: the layers to take into
    account.}
  \item{weights}{For \code{supersom} maps: weights of the layers that
    are used for mapping.}
  \item{scale.distances}{whether to rescale distances per layer in the
    case of \code{supersom} maps (default): if \code{TRUE} the maximal
    distance of each layer equals one. If the absolute values of the
    distances per layer should be used, this argument should be set to
    \code{FALSE}. Note that in that case, when mapping the training
    data, the result returned by \code{map.kohonen} will differ from the
    mapping present in the map.}
  \item{\dots}{Currently ignored.}
}
\value{
  A list with elements
  \item{unit.classif}{a vector of units that are closest to the objects
    in the data matrix.}
  \item{dists}{distances (currently only Euclidean distances) of the
    objects to the units.}
  \item{whatmap,weights,scale.distances}{Values used for these arguments.}
}
\author{Ron Wehrens}
\seealso{\code{\link{predict.kohonen}}}
\examples{
data(wines)
set.seed(7)

training <- sample(nrow(wines), 120)
Xtraining <- scale(wines[training, ])
somnet <- som(Xtraining, somgrid(5, 5, "hexagonal"))

mapping <- map(somnet,
               scale(wines[-training, ],
                     center=attr(Xtraining, "scaled:center"),
                     scale=attr(Xtraining, "scaled:scale")))
}
\keyword{classif}
