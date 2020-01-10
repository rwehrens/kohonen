\name{plot.kohonen}
\alias{plot.kohonen}
\alias{identify.kohonen}
\alias{add.cluster.boundaries}
\title{Plot kohonen object}
\description{
  Plot objects of class \code{kohonen}. Several types
  of plots are supported.
}
\usage{
\method{plot}{kohonen}(x, type = c("codes", "changes", "counts",
                       "dist.neighbours", "mapping", "property", "quality"),
                       whatmap = NULL, classif = NULL, labels = NULL,
                       pchs = NULL, main = NULL, palette.name = NULL,
                       ncolors, bgcol = NULL, zlim = NULL,
                       heatkey = TRUE, property, codeRendering = NULL,
                       keepMargins = FALSE, heatkeywidth = .2,
                       shape = c("round", "straight"), border = "black",
                       na.color = "gray", \dots)
\method{identify}{kohonen}(x, \dots)
add.cluster.boundaries(x, clustering, lwd = 5, ...)
}
\arguments{
  \item{x}{kohonen object.}
  \item{type}{type of plot. (Wow!)}
  \item{whatmap}{For a "codes" plot: what maps
    to show; for the "dist.neighbours" plot: what maps to take into
    account when calculating distances to neighbouring units.}
  \item{classif}{classification object, as returned by
    \code{predict.kohonen}, or vector of unit numbers.
    Only needed if \code{type} equals "mapping" and
    "counts".}
  \item{labels}{labels to plot when \code{type} equals
    "mapping".}
  \item{pchs}{symbols to plot when \code{type} equals
    "mapping".}
  \item{main}{title of the plot.}
  \item{palette.name}{colors to use as unit background for "codes",
    "counts", "prediction", "property", and "quality" plotting types.}
  \item{ncolors}{number of colors to use for the unit
    backgrounds. Default is 20 for continuous data, and the number of
    distinct values (if less than 20) for categorical data.}
  \item{bgcol}{optional argument to colour the unit backgrounds for the
    "mapping" and "codes" plotting type. Defaults to "gray" and
    "transparent" in both types, respectively.}
  \item{zlim}{optional range for color coding of unit backgrounds.}
  \item{heatkey}{whether or not to generate a heatkey at the left side
    of the plot in the "property" and "counts" plotting types.}
  \item{property}{values to use with the "property" plotting
    type.}
  \item{codeRendering}{How to show the codes. Possible choices:
    "segments", "stars" and "lines".}
  \item{keepMargins}{if \code{FALSE} (the default), restore the original
    graphical parameters after plotting the kohonen map. If \code{TRUE},
    one retains the map coordinate system so that one can add symbols to
    the plot, or map unit numbers using the \code{identify} function.}
  \item{heatkeywidth}{width of the colour key; the default of 0.2 should
    work in most cases but in some cases, e.g. when plotting multiple
    figures, it may need to be adjusted.}
  \item{shape}{kind shape to be drawn: "round" (circle) or "straight".
    Choosing "straight" produces a map of squares when the grid is
    "rectangular", and produces a map of hexagons when the grid is
    "hexagonal".}
  \item{border}{color of the shape's border.}
  \item{na.color}{background color matching NA - default "gray".}
  \item{lwd, \dots}{other graphical parameters.}
  \item{clustering}{cluster labels of the map units.}
}
\details{Several different types of plots are supported:
  \describe{
    \item{"changes"}{shows the mean distance to the closest codebook vector
      during training.}
    \item{"codes"}{shows the codebook vectors.}
    \item{"counts"}{shows the number of objects mapped to the
      individual units. Empty units are depicted in gray.}
    \item{"dist.neighbours"}{shows the sum of the distances to all
      immediate neighbours. This kind of visualisation is also known as
      a U-matrix plot. Units near a class boundary can be expected to
      have higher average distances to their neighbours. Only available
      for the "som" and "supersom" maps, for the moment.}
    \item{"mapping"}{shows where objects are mapped. It
      needs the "classif" argument, and a "labels" or "pchs"
      argument.}
    \item{"property"}{properties of each unit can be calculated and
      shown in colour code. It can be used to visualise the similarity
      of one particular object to all units in the map, to show the mean
      similarity of all units and the objects mapped to them,
      etcetera. The parameter \code{property} contains the numerical
      values. See examples below.}
    \item{"quality"}{shows the mean distance of objects mapped to a
      unit to the codebook vector of that unit. The smaller the
      distances, the better the objects are represented by the codebook
      vectors. It is possible to visualize this for the complete set of
      layers used in training, or for individual layers only (using the
      \code{whatmap} argument).}
  }
  Function \code{identify.kohonen} shows the number of a unit that is
  clicked on with the mouse. The tolerance is calculated from the ratio
  of the plotting region and the user coordinates, so clicking at any
  place within a unit should work.

  Function \code{add.cluster.boundaries} will add to an existing plot of
  a map thick lines, visualizing which units would be clustered
  together. In toroidal maps, boundaries at the edges will only be shown
  on the top and right sides to avoid double boundaries.
}
\value{Several types of plots return useful values (invisibly): the
  \code{"counts"}, \code{"dist.neighbours"}, and \code{"quality"} return
  vectors corresponding to the information visualized in the plot (unit
  background colours and heatkey).
}
\author{Ron Wehrens}
\seealso{\code{\link{som}}, \code{\link{supersom}}, \code{\link{xyf}},
  \code{\link{predict.kohonen}}}

\examples{
data(wines)
set.seed(7)

kohmap <- xyf(scale(wines), vintages,
              grid = somgrid(5, 5, "hexagonal"), rlen=100)
plot(kohmap, type="changes")
counts <- plot(kohmap, type="counts", shape = "straight")

## show both sets of codebook vectors in the map
par(mfrow = c(1,2))
plot(kohmap, type="codes", main = c("Codes X", "Codes Y"))

par(mfrow = c(1,1))
similarities <- plot(kohmap, type="quality", palette.name = terrain.colors)
plot(kohmap, type="mapping",
     labels = as.integer(vintages), col = as.integer(vintages),
     main = "mapping plot")

## add background colors to units according to their predicted class labels
xyfpredictions <- classmat2classvec(getCodes(kohmap, 2))
bgcols <- c("gray", "pink", "lightgreen")
plot(kohmap, type="mapping", col = as.integer(vintages),
     pchs = as.integer(vintages), bgcol = bgcols[as.integer(xyfpredictions)],
     main = "another mapping plot", shape = "straight", border = NA)

## Show 'component planes'
set.seed(7)
sommap <- som(scale(wines), grid = somgrid(6, 4, "hexagonal"))
plot(sommap, type = "property", property = getCodes(sommap, 1)[,1],
     main = colnames(getCodes(sommap, 1))[1])

## Show the U matrix
Umat <- plot(sommap, type="dist.neighbours", main = "SOM neighbour distances")
## use hierarchical clustering to cluster the codebook vectors
som.hc <- cutree(hclust(object.distances(sommap, "codes")), 5)
add.cluster.boundaries(sommap, som.hc)

## and the same for rectangular maps
set.seed(7)
sommap <- som(scale(wines),grid = somgrid(6, 4, "rectangular"))
plot(sommap, type="dist.neighbours", main = "SOM neighbour distances")
## use hierarchical clustering to cluster the codebook vectors
som.hc <- cutree(hclust(object.distances(sommap, "codes")), 5)
add.cluster.boundaries(sommap, som.hc)
}
\keyword{classif}
