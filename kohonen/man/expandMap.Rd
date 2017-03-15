\name{expandMap}
\alias{expandMap}
\title{Expand a self-organising map}
\description{
Double the size of a map, imputing the codebookvectors of the new units
by averiging their immediate neighbours.
}
\usage{
expandMap(kohobj)
}
\arguments{
  \item{kohobj}{Object of class \code{"kohonen"})}
}
\value{A new kohonen object, with a double size.}
\author{Ron Wehrens}

\examples{
data(yeast)
yeast.supersom <- supersom(yeast, somgrid(4, 4, "hexagonal"),
                           whatmap = 3:6, maxNA.fraction = .5)
yeast.supersom2 <- expandMap(yeast.supersom)
yeast.supersom3 <- supersom(yeast, yeast.supersom2$grid,
                            whatmap = 3:6, maxNA.fraction = .5,
                            init = yeast.supersom2$codes[3:6])
}
\keyword{classif}
