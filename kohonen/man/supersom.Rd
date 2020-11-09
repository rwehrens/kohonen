\name{supersom}
\alias{som}
\alias{xyf}
\alias{supersom}
\alias{nunits}
\title{Self- and super-organising maps}
\description{
  A supersom is an extension of self-organising maps (SOMs) to multiple
  data layers, possibly with different numbers and different types of
  variables (though equal numbers of objects). NAs are allowed. A
  weighted distance over all layers is calculated to determine the
  winning units during training.
  Functions \code{som} and \code{xyf} are simply wrappers for supersoms
  with one and two layers, respectively. Function \code{nunits} is a
  utility function returning the number of units in the map.
}

\usage{
som(X, ...)
xyf(X, Y, ...)
supersom(data, grid=somgrid(), rlen = 100, alpha = c(0.05, 0.01),
         radius = quantile(nhbrdist, 2/3), 
         whatmap = NULL, user.weights = 1, maxNA.fraction = 0L,
         keep.data = TRUE, dist.fcts = NULL,
         mode = c("online", "batch", "pbatch"), cores = -1, init,
         normalizeDataLayers = TRUE)
nunits(kohobj)
}

\arguments{
  \item{X, Y}{numerical data matrices, or factors. No \code{data.frame}
    objects are allowed - convert them to matrices first.}
  \item{data}{list of data matrices (numerical) of factors. If a vector
    is entered, it will be converted to a one-column matrix. No
    \code{data.frame} objectss are allowed.}
  \item{grid}{a grid for the codebook vectors:
    see \code{somgrid}.}
  \item{rlen}{the number of times the complete data set will be
    presented to the network.}
  \item{alpha}{learning rate, a vector of two numbers indicating the
    amount of change. Default is to decline linearly from 0.05 to 0.01
    over \code{rlen} updates. Not used for the batch algorithm.}
  \item{radius}{the radius of the neighbourhood, either given as a
    single number or a vector (start, stop). If it is given as a single
    number the radius will change linearly from \code{radius} to zero; as
    soon as the neighbourhood gets smaller than one only the winning unit
    will be updated. Note that the default before version 3.0 was to run
    from \code{radius} to \code{-radius}. If nothing is supplied, the
    default is to start with a value that covers 2/3 of all unit-to-unit
    distances.} 
  \item{whatmap}{What data layers to use. If unspecified all layers are used.}
  \item{user.weights}{the weights given to individual layers. This can
    be a single number (all layers have the same weight, the default), a
    vector of the same length as the \code{whatmap} argument, or a vector
    of the same length as the \code{data} argument. In xyf maps, this
    argument provides the same functionality as the now-deprecated
    \code{xweight} argument that was used prior to version 3.0.}
  \item{maxNA.fraction}{the maximal fraction of values that may be NA to
    prevent the row to be removed.}
  \item{keep.data}{if TRUE, return original data and mapping
    information. If FALSE, only return the trained map (in essence the
    codebook vectors).}
  \item{dist.fcts}{vector of distance functions to be used for the
    individual data layers, of the same length as the \code{data}
    argument, or the same length of the \code{whatmap} argument. If the
    length of this vector is one, the  
    same distance will be used for all layers. Admissable values
    currently are "sumofsquares", "euclidean", "manhattan", and
    "tanimoto". Default is to use "sumofsquares" for continuous data,
    and "tanimoto" for factors.}
  \item{mode}{type of learning algorithm.}
  \item{cores}{number of cores to use in the "pbatch" learning mode. The
    default, -1, corresponds to using all available cores.}
  \item{init}{list of matrices, initial values for the
    codebook vectors. The list should have the same length as the data
    list, and corresponding numbers of variables (columns). Each list
    element should have a number of rows corresponding to the number of
    units in the map.}
  \item{normalizeDataLayers}{boolean, indicating whether
    \code{distance.weights} should be calculated (see details section).
    If \code{normalizeDataLayers == FALSE} the user weights 
    are applied to the data immediately.}
  \item{kohobj}{an object of class \code{kohonen}.}
  \item{\dots}{Further arguments for the \code{supersom} function
    presented to the \code{som} or \code{xyf} wrappers.}
}
\value{
  An object of class "kohonen" with components
  \item{data}{data matrix, only returned if \code{keep.data == TRUE}.}
  \item{unit.classif}{winning units for all data objects,
    only returned if \code{keep.data == TRUE}.}
  \item{distances}{distances of objects to their corresponding winning
    unit, only returned if \code{keep.data == TRUE}.}
  \item{grid}{the grid, an object of class \code{somgrid}.}
  \item{codes}{a list of matrices containing codebook vectors.}
  \item{changes}{matrix of mean average deviations from code vectors;
    every map corresponds with one column.}
  \item{na.rows}{vector of row numbers with too many NA values
    (according to argument \code{maxNA.fraction}.}
  \item{alpha, radius, user.weights, whatmap, maxNA.fraction}{input arguments
    presented to the function.}
  \item{distance.weights}{if \code{normalizeDataLayers} weights to
    equalize the influence of the individual data layers, else a vector
    of ones.}
  \item{dist.fcts}{distance functions corresponding to all layers of the
    data, not just the ones indicated by the whatmap argument.}
}
\details{In order to avoid some layers to overwhelm others, simply
  because of the scale of the data points, the \code{supersom} function
  by default applies internal weights to balance this. The \code{user.weights}
  argument is applied on top of that: the result is that when a user
  specifies equal weights for all layers (the default), all layers
  contribute equally to the global distance measure. For large data
  sets (defined as containing more than 500 records), a sample of size
  500 is used to calculate the mean distances in each data layer. If
  \code{normalizeDataLayers == FALSE} the user weights are applied
  directly to the data (\code{distance.weights} are set to 1).

  Various definitions of the Tanimoto distance exist in the
  literature. The implementation here returns (for two binary vectors of
  length n) the fraction of cases in which the two vectors
  disagree. This is basically the Hamming distance divided by n - the
  incorrect naming is retained (for the moment) to guarantee backwards
  compatibility. If the vectors are not binary, they will be
  converted to binary strings (with 0.5 as the class boundary). This
  measure should not be used when variables are outside the range [0-1];
  a check is done to make sure this is the case.
}
\references{R. Wehrens and L.M.C. Buydens, J. Stat. Softw. 21 (5), 2007;
R. Wehrens and J. Kruisselbrink, submitted, 2017.}
\author{Ron Wehrens and Johannes Kruisselbrink}

\seealso{\code{somgrid}, \code{\link{plot.kohonen}},
  \code{\link{predict.kohonen}}, \code{\link{map.kohonen}}}
\examples{
data(wines)

## som
som.wines <- som(scale(wines), grid = somgrid(5, 5, "hexagonal"))
summary(som.wines)
nunits(som.wines)

## xyf
xyf.wines <- xyf(scale(wines), vintages, grid = somgrid(5, 5, "hexagonal"))
summary(xyf.wines)

## supersom example
data(yeast)
yeast.supersom <- supersom(yeast, somgrid(6, 6, "hexagonal"),
                           whatmap = c("alpha", "cdc15", "cdc28", "elu"),
                           maxNA.fraction = .5)

plot(yeast.supersom, "changes")

obj.classes <- as.integer(yeast$class)
colors <- c("yellow", "green", "blue", "red", "orange")
plot(yeast.supersom, type = "mapping", col = colors[obj.classes],
     pch = obj.classes, main = "yeast data")
}

\keyword{classif}

