somgrid <- function (xdim = 8, ydim = 6,
                     topo = c("rectangular", "hexagonal"),
                     neighbourhood.fct = c("bubble", "gaussian"),
                     toroidal = FALSE) 
{
  topo <- match.arg(topo)
  x <- 1L:xdim
  y <- 1L:ydim
  pts <- as.matrix(expand.grid(x = x, y = y))

  if (topo == "hexagonal") {
    pts[, 1L] <- pts[, 1L] + 0.5 * (pts[, 2L]%%2)
    pts[, 2L] <- sqrt(3)/2 * pts[, 2L]
  }

  ## Check neighbourhood function
  neighbourhood.fct <- match.arg(neighbourhood.fct)
  neighbourhood.fct <- factor(neighbourhood.fct,
                              levels = c("bubble", "gaussian"))

  res <- list(pts = pts, xdim = xdim, ydim = ydim, topo = topo,
              neighbourhood.fct = neighbourhood.fct,
              toroidal = toroidal)
  class(res) <- "somgrid"

  res
}
