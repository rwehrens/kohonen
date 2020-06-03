\name{predict.kohonen}
\alias{predict.kohonen}
\title{Predict properties using a trained Kohonen map}
\description{
  Map objects to a trained Kohonen map, and return for each object the
  desired property associated with the corresponding winning
  unit. These properties may be provided explicitly (argument
  \code{unit.predictions}) or implicitly (by providing
  \code{trainingdata}, that will be mapped to the SOM - the averages of
  the winning units for the trainingdata then will be used as
  unit.predictions). If not given at all, the codebook vectors of the
  map will be used.
}
\usage{
\method{predict}{kohonen}(object,
                          newdata = NULL,
                          unit.predictions = NULL,
                          trainingdata = NULL,
                          whatmap = NULL,
                          threshold = 0,
                          maxNA.fraction = object$maxNA.fraction,
                          ...)
}
\arguments{
  \item{object}{Trained network, containing one or more information layers.}
  \item{newdata}{List of data matrices, or one single data matrix, for
    which predictions are to be made. The data layers should match those
    in the trained map. If not presented, the training data in the map
    will be used. No \code{data.frame} objects are allowed.}
  \item{unit.predictions}{Explicit definition of the predictions for each
    unit. Should be a list of matrices, vectors or factors, of the same
    length as \code{object$codes}.}
  \item{trainingdata}{List of data matrices, or one single data matrix,
    determining the mapping of the training data. Normally, data stored
    in the \code{kohonen} object will be used for this, but one can also
    specify this argument explicitly. Layers should match the trained map.}
  \item{whatmap, maxNA.fraction}{parameters that usually will
    be taken from the \code{x} object, but can be supplied by the user
    as well. See \code{\link{supersom}} for more information.}
  \item{threshold}{Used in converting class predictions back into
    factors; see \code{\link{classmat2classvec}}.}
  \item{\dots}{Further arguments to be passed to \code{map.kohonen}, in
    particular \code{user.weights}. If not provided will be taken from
    \code{object}.}
}
\value{
  Returns a list with components
  \item{prediction}{predicted values for the properties of
    interest. When multiple values are predicted, this element is a
    list, otherwise a vector or a matrix.}
  \item{unit.classif}{vector of unit numbers to which objects in the newdata
    object are mapped.}
  \item{unit.predictions}{prediction values associated with map
    units. Again, when multiple properties are predicted, this is a
    list.}
  \item{whatmap}{the numbers of the data layers in the kohonen object
    used in the mapping on which the predictions are based.}
}
\details{
  The new data are mapped to the trained SOM using
  the layers indicated by the \code{whatmap} argument. The predictions
  correspond to the \code{unit.predictions}, normally corresponding to
  the averages of the training data mapping to individual units. If no
  \code{unit.predictions} are provided, the \code{trainingdata} will be
  used to calculate them - if \code{trainingdata} is not provided by the
  user and the \code{kohonen} object contains data, these will be used.
  If no objects of the training data are mapping to a particular unit,
  the prediction for that unit will be NA.
}
\author{Ron Wehrens}
\seealso{\code{\link{som}},\code{\link{xyf}},
  \code{\link{supersom}}, \code{\link{map}}}
\examples{
data(wines)

training <- sample(nrow(wines), 120)
Xtraining <- scale(wines[training, ])
Xtest <- scale(wines[-training, ],
               center = attr(Xtraining, "scaled:center"),
               scale = attr(Xtraining, "scaled:scale"))
trainingdata <- list(measurements = Xtraining,
                     vintages = vintages[training])
testdata <- list(measurements = Xtest, vintages = vintages[-training])

mygrid = somgrid(5, 5, "hexagonal")
som.wines <- supersom(trainingdata, grid = mygrid)

## ################################################################
## Situation 0: obtain expected values for training data (all layers,
## also if not used in training) on the basis of the position in the map
som.prediction <- predict(som.wines)

## ################################################################
## Situation 1: obtain predictions for all layers used in training

som.prediction <- predict(som.wines, newdata = testdata)
table(vintages[-training], som.prediction$predictions[["vintages"]])


## ################################################################
## Situation 2: obtain predictions for the vintage based on the mapping
## of the sample characteristics only. There are several ways of doing this:

som.prediction <- predict(som.wines, newdata = testdata,
                          whatmap = "measurements")
table(vintages[-training], som.prediction$predictions[["vintages"]])

## same, but now indicated implicitly
som.prediction <- predict(som.wines, newdata = testdata[1])
table(vintages[-training], som.prediction$predictions[["vintages"]])

## if no names are present in the list elements whatmap needs to be
## given explicitly; note that the order of the data layers needs to be
## consistent with the kohonen object
som.prediction <- predict(som.wines, newdata = list(Xtest), whatmap = 1)
table(vintages[-training], som.prediction$predictions[["vintages"]])

## for xyf: explicitly indicate which layer is to be used for the mapping
xyf.wines <- xyf(Xtraining, vintages[training], grid = mygrid)
xyf.prediction <- predict(xyf.wines, Xtest, whatmap = 1)
table(vintages[-training], xyf.prediction$predictions[[2]])

## ###############################################################
## Situation 3: predictions for layers not present in the original
## data. Training data need to be provided for those layers.
som.wines <- supersom(Xtraining, grid = mygrid)
som.prediction <- predict(som.wines, newdata = testdata,
                          trainingdata = trainingdata)
table(vintages[-training], som.prediction$predictions[["vintages"]])

## ################################################################
## yeast examples, including NA values

data(yeast)
training.indices <- sample(nrow(yeast$alpha), 300)
training <- rep(FALSE, nrow(yeast$alpha))
training[training.indices] <- TRUE

## unsupervised mapping, based on the alpha layer only. Prediction
## for all layers including alpha
yeast.som <- supersom(lapply(yeast, function(x) subset(x, training)),
                      somgrid(4, 6, "hexagonal"),
                      whatmap = "alpha", maxNA.fraction = .5)
yeast.som.prediction <-
  predict(yeast.som,
          newdata = lapply(yeast, function(x) subset(x, !training)))

table(yeast$class[!training], yeast.som.prediction$prediction[["class"]])

## ################################################################
## supervised mapping - creating the map is now based on both
## alpha and class, prediction for class based on the mapping of alpha.
yeast.som2 <- supersom(lapply(yeast, function(x) subset(x, training)),
                       grid = somgrid(4, 6, "hexagonal"),
                       whatmap = c("alpha", "class"), maxNA.fraction = .5)
yeast.som2.prediction <-
  predict(yeast.som2,
          newdata = lapply(yeast, function(x) subset(x, !training)),
          whatmap = "alpha")
table(yeast$class[!training], yeast.som2.prediction$prediction[["class"]])
}
\keyword{classif}
