library(kohonen)
library(grid)

## Initialize a matrix of rgb colors with each rgb level scaled between 0 and 1
colors <- rbind(c(1,0,0), c(1,1,0), c(0,1,0), c(0,1,1), c(0,0,1), c(1,0,1))

## Initialize a somgrid grid and create an initial matrix of codebook vectors
## filled with colors randomly sampled from the colors matrix
cells_x <- 50
cells_y <- 50
somgrid <- somgrid(cells_x,cells_y, "rectangular")
init_codes <- replicate(3, runif(cells_x * cells_y))

## Start training the rgb somgrid
som_rgb <- som(
  X = colors,
  grid = somgrid,
  init = init_codes,
  rlen = 100,
  keep.data = FALSE)

## Render the rgb values of the codebook vectors of the trained som
map <- rgb(som_rgb$codes[[1]][,1], som_rgb$codes[[1]][,2], som_rgb$codes[[1]][,3])
col <- matrix(map,nrow = som_rgb$grid$xdim,ncol = som_rgb$grid$ydim)
grid.raster(col, interpolate=FALSE)

