\name{unit.distances}
\alias{unit.distances}
\alias{somgrid}
\title{SOM-grid related functions}
\description{
  Function \code{somgrid} (modified from the version in the class
  package) sets up a grid of units, of a specified size 
  and topology. Distances between grid units are calculated by function
  \code{unit.distances}.
}
\usage{
somgrid(xdim = 8, ydim = 6, topo = c("rectangular", "hexagonal"),
        neighbourhood.fct = c("bubble", "gaussian"), toroidal = FALSE)
unit.distances(grid, toroidal)
}
\arguments{
  \item{xdim, ydim}{dimensions of the grid.}
  \item{topo}{choose between a hexagonal or rectangular topology.}
  \item{neighbourhood.fct}{choose between bubble and gaussian
    neighbourhoods when training a SOM.}
  \item{toroidal}{logical, whether the grid is toroidal or not. If not
    provided to the \code{unit.distances} function, the information in
    the \code{grid} object will be used.}
  \item{grid}{an object of class \code{somgrid}.}
}
\value{
  Function \code{somgrid} returns an object of class "somgrid", with
  elements \code{pts}, and the input arguments to the function.
  
  Function \code{unit.distances} returns a (symmetrical) matrix
  containing distances. When \code{grid$n.hood} equals "circular",
  Euclidean distances are used; for \code{grid$n.hood} is "square"
  maximum distances. For toroidal maps (joined at the edges) distances
  are calculated for the shortest path.
}
\author{Ron Wehrens}
\examples{
mygrid <- somgrid(5, 5, "hexagonal")
fakesom <- list(grid = mygrid)
class(fakesom) <- "kohonen"

par(mfrow = c(2,1))
dists <- unit.distances(mygrid)
plot(fakesom, type="property", property = dists[1,],
     main="Distances to unit 1", zlim=c(0,6),
     palette = rainbow, ncolors = 7)

dists <- unit.distances(mygrid, toroidal=TRUE)
plot(fakesom, type="property", property = dists[1,],
     main="Distances to unit 1 (toroidal)", zlim=c(0,6),
     palette = rainbow, ncolors = 7)
}
\keyword{classif}
