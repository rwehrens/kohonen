\name{getCodes}
\alias{getCodes}
\title{Extract codebook vectors from a kohonen object}
\description{Utility function for extracting codebook vectors. These
  are present as a list element in a \code{kohonen} object, and
  themselves are a list as well, with one entry for each data
  layer. This function returns either a list of codebook matrices (if
  more layers are selected), or just one matrix (if one layer is
  selected).
}
\usage{
getCodes(x, idx = 1:length(codes))
}
\arguments{
  \item{x}{An object of class \code{kohonen}.}
  \item{idx}{Indices of the layer(s) for which codebook vectors are
    returned.}
}
\value{If \code{idx} is a single number, a matrix of codebook vectors;
  if it is a vector of numbers, a list of codebook matrices.}
\author{Ron Wehrens}
\seealso{\code{\link{supersom}}}
\examples{
data(wines)
set.seed(7)
som.wines <- som(scale(wines), grid = somgrid(5, 5, "hexagonal"))
dim(getCodes(som.wines))
}
\keyword{manip}
