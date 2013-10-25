\name{predict.kohonen}
\alias{predict.kohonen}
\title{Predict properties using a trained Kohonen map}
\description{
  Map objects to a trained Kohonen map, and return for each object the
  property associated with the corresponding winning unit. For
  \code{som} and \code{supersom} maps, the unit properties are
  calculated using explicit arguments \code{trainX} and \code{trainY};
  for \code{xyf} and \code{bdk} maps, the predicted properties are the
  Y-codebookvectors. Note that in the latter case only the X-space is
  used for prediction.
}
\usage{
\method{predict}{kohonen}(object, newdata, trainX, trainY, unit.predictions,
                 threshold = 0, whatmap = NULL, weights = 1, ...)
}
\arguments{
  \item{object}{Trained network.}
  \item{newdata}{Data matrix for which predictions are to be made. If
    not given, defaults to the training data (when available).}
  \item{trainX}{Training data for obtaining predictions for unsupervised
    maps; necessary for \code{som} maps trained with the
    \code{keep.data = FALSE} option.}
  \item{trainY}{Values for the dependent variable for the training data;
    necessary for \code{som} and \code{supersom} maps.}
  \item{unit.predictions}{Possible override of the predictions for each
    unit.}
  \item{threshold}{Used in class predictions; see
    \code{\link{classmat2classvec}}.}
  \item{whatmap}{For \code{supersom} maps: what layers to use in the
    mapping.}
  \item{weights}{For \code{supersom} maps: weights of layers uses in the
    mapping.}
  \item{\dots}{Currently not used.} 
}
\value{
  Returns a list with components
  \item{prediction}{predicted values for the properties of
    interest. When multiple values are predicted, this element is a
    list, otherwise a vector or a matrix.}
  \item{unit.classif}{unit numbers to which objects in the data matrix are
    mapped.}
  \item{unit.predictions}{mean values associated with map units. Again,
    when multiple properties are predicted, this is a list.}
}
\author{Ron Wehrens}
\seealso{\code{\link{som}},\code{\link{xyf}},\code{\link{bdk}},
  \code{\link{supersom}}, \code{\link{map}}}
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

### more complicated examples
data(yeast)

### only consider complete cases
missings <- (apply(cbind(yeast$alpha, yeast$cdc15), 1,
                   function(x) any(is.na(x))))

yeast2 <- list(alpha = yeast$alpha[!missings,],
               cdc15 = yeast$cdc15[!missings,],
               class = yeast$class[!missings])

set.seed(7)
training.indices <- sample(nrow(yeast2$alpha), 300)
training <- rep(FALSE, nrow(yeast2$alpha))
training[training.indices] <- TRUE

## unsupervised mapping
yeast2.som <- som(yeast2$alpha[training,], somgrid(4, 6, "hexagonal"))
yeast2.som.prediction <- predict(yeast2.som,
                                 newdata = yeast2$alpha[!training,],
                                 trainY = yeast2$class[training])
table(yeast2$class[!training], yeast2.som.prediction$prediction)

## supervised mapping (XYF) - trainY is no longer necessary
yeast2.xyf <- xyf(yeast2$alpha[training,], yeast2$class[training],
                  somgrid(4, 6, "hexagonal"))
yeast2.xyf.prediction <- predict(yeast2.xyf,
                                 newdata = yeast2$alpha[!training,])
table(yeast2$class[!training], yeast2.xyf.prediction$prediction)

## supervised mapping (BDK)
yeast2.bdk <- bdk(yeast2$alpha[training,], yeast2$class[training],
                  somgrid(4, 6, "hexagonal"))
yeast2.bdk.prediction <- predict(yeast2.bdk,
                                 newdata = yeast2$alpha[!training,])
table(yeast2$class[!training], yeast2.bdk.prediction$prediction)

## unsupervised mapping (supersom): prediction of data layer not used
## in training
yeast2.ssom <- supersom(lapply(yeast2, function(x) subset(x, training)),
                        grid = somgrid(4, 6, "hexagonal"),
                        whatmap = 1)
yeast2.ssom.prediction <- predict(yeast2.ssom,
                                  newdata = lapply(yeast2,
                                    function(x) subset(x, !training)),
                                  trainY = list(class = yeast2$class[training]))
table(yeast2$class[!training], yeast2.ssom.prediction$prediction)

## supervised mapping (supersom): prediction of a data layer that has
## been used in training - trainY is not necessary
yeast2.ssom2 <- supersom(lapply(yeast2, function(x) subset(x, training)),
                         grid = somgrid(4, 6, "hexagonal"),
                         whatmap = c(1,3))
yeast2.ssom2.prediction <- predict(yeast2.ssom2, 
                                   newdata = lapply(yeast2,
                                     function(x) subset(x, !training)),
                                   whatmap = 1)
table(yeast2$class[!training], yeast2.ssom2.prediction$prediction)
}
\keyword{classif}
