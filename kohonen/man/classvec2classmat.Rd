\name{classvec2classmat}
\alias{classvec2classmat}
\alias{classmat2classvec}
\title{Convert a classification vector into a matrix or the other way around.}
\description{
  Functions toggle between a matrix representation, where class
  membership is indicated with one '1' and for the rest zeros at each
  row, and a factor. The classification matrix contains one column per
  class. Conversion from a class matrix to a class vector assigns each
  row to the column with the highest value. An optional argument can be
  used to assign only those objects that have a probability higher than
  a certain threshold (default is 0).
}
\usage{
classvec2classmat(yvec)
classmat2classvec(ymat, threshold=0)
}

\arguments{
  \item{yvec}{class vector. Usually a factor; if it is a vector of
    integer values, it will be converted to a factor.}
  \item{ymat}{class matrix: every column corresponds to a class.}
  \item{threshold}{only classify into a class if the probability is
    larger than this threshold.}
}
\value{
  \code{classvec2classmat} returns the classification matrix, where each
  column consists of zeros and ones; \code{classmat2classvec} returns a
  factor.
}
\author{Ron Wehrens}

\seealso{\code{\link{som}},\code{\link{xyf}},\code{\link{supersom}}}
\examples{
classes <- c(rep(1, 5), rep(2, 7), rep(3, 9))
classmat <- classvec2classmat(classes)
classmat
classmat2classvec(classmat)
}

\keyword{classif}
