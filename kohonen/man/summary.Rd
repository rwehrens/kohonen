\name{summary.kohonen}
\alias{summary.kohonen}
\alias{print.kohonen}
\title{Summary and print methods for kohonen objects}
\description{
  Summary and print methods for \code{kohonen} objects. The \code{print}
  method shows the dimensions and the topology of the map; if
  information on the training data is included, the \code{summary}
  method additionally prints information on the size of the data, the
  distance functions used, and the
  mean distance of an object to its closest codebookvector, which is an
  indication of the quality of the mapping.}
\usage{
\method{summary}{kohonen}(object, \dots)
\method{print}{kohonen}(x, \dots)
}
\arguments{
  \item{x, object}{a \code{kohonen} object}
  \item{\dots}{Not used.}
}
\author{Ron Wehrens}
\seealso{\code{\link{som}}, \code{\link{xyf}}, \code{\link{supersom}}}
\examples{
data(wines)
xyf.wines <- xyf(scale(wines), classvec2classmat(vintages),
                 grid = somgrid(5, 5, "hexagonal"))
xyf.wines
summary(xyf.wines)
}
\keyword{classif}
