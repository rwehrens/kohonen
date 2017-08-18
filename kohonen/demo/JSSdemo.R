### Demo presenting the application results from the 2017 JSS
### publication on the kohonen package v. 3.0. Training the SOM for
### the pepper image takes some time, on my laptop approx 40 seconds.

require("kohonen")
require("RColorBrewer")

## function to plot the pepper picture
plot.picture <- function(dmat, what = c("col", "truth"),
                         image.dims = c(600, 800), mycols, ...) {
  what <- match.arg(what)

  if (what == "truth") { # dmat is a factor
    if (missing(mycols)) {
      require(lattice)    
      mycols <- brewer.pal(nlevels(dmat), "Set1")
    }
    myat <- 0:nlevels(dmat) + .5
    immat <- matrix(as.numeric(dmat), image.dims[1], image.dims[2])
    levelplot(t(immat)[,nrow(immat):1], at = myat,
              colorkey = list(at = myat,
                              space = "top",
                              col = mycols,
                              labels = list(at = myat[-1] - .5,
                                            labels = levels(dmat))),
              col.regions = mycols, xlab = "", ylab = "", ...)
   } else { # dmat is a matrix with columns "red", "green" and "blue"
    require("grid")
    immat <- rgb(dmat[,"red"], dmat[,"green"], dmat[,"blue"],
                 maxColorValue = 255)
    dim(immat) <- image.dims
    grid.raster(immat, interpolate = FALSE, ...)
  }
}

### ####################################################################
### Image segmentation example

## Prepare input data
data(peppaPic)
X <- peppaPic[,c("red", "green", "blue")]
Y <- factor(peppaPic[,"truth"])
classes <- c("background", "leaves", "peppers", "peduncles", "stems",
             "sideshoots", "wires", "cuts")
levels(Y) <- classes
coords <- cbind(rep(1:800, each = 600), rep(1:600, 800))
colnames(coords) <- c("x", "y")

imdata <- list("rgb" = X, "truth" = Y, "coords" = coords)

## do segmentation
mygrid <- somgrid(4, 8, "hexagonal")
set.seed(15)
system.time(somsegm1 <- supersom(imdata, whatmap = c("rgb", "coords"),
                                 user.weights = c(1, 9),
                                 grid = mygrid))
## reconstruct the original image
reconstrIm <- predict(somsegm1, newdata = imdata)
X32 <- reconstrIm$predictions[["rgb"]]

plot.picture(X, x = unit(.025, "npc"),
             just = "left", width = unit(.45, "npc"))
plot.picture(X32, x = unit(.525, "npc"), 
             just = "left", width = unit(.45, "npc"))

## find class predictions for map units, concentrate only on peppers
## and leaves
truthpredictions <-
  classmat2classvec(predictions$unit.predictions$truth)
trclasses <- levels(truthpredictions)
trclasses[!(trclasses %in% c("leaves", "peppers"))] <- "other"
levels(truthpredictions) <- trclasses
pepperunits <- which(truthpredictions == "peppers")

## visualize the rgb part of the map, using unit class labels as
## background colors
rgbpalette <- function(n) {
  if (n != 3) stop("Hey! Count again!")
  c("red", "green", "blue")
}
cluscolors <- brewer.pal(nlevels(truthpredictions), "Set3")
unitcols <- cluscolors[as.integer(truthpredictions)]
plot(somsegm1, "codes", whatmap = "rgb", bgcol = unitcols, 
     palette.name = rgbpalette, main = "")
## add unit numbers and a legend on top
idx <- (0:7)*4 + 1
text(somsegm1$grid$pts[idx,1] - .5, somsegm1$grid$pts[idx,2], 
     labels = idx, pos = 2)
add.cluster.boundaries(somsegm1, truthpredictions)
legend("top", legend = rev(levels(truthpredictions)), bty="n",
       pch = rep(21, 3), pt.bg = rev(cluscolors), ncol = 3)

## visualize the pixel location part of the map
dev.new()
plot(somsegm1$codes[["coords"]][,1], 
     max(somsegm1$codes[["coords"]][,2]) - somsegm1$codes[["coords"]][,2],
     pch = 21, cex = 3, axes = FALSE, xlab = "", ylab = "",
     bg = rgb(somsegm1$codes[["rgb"]][,1],
                          somsegm1$codes[["rgb"]][,2],
                          somsegm1$codes[["rgb"]][,3], maxColorValue = 255))
## numbers of pepper units are printed in red
text(somsegm1$codes[["coords"]][,1],
     max(somsegm1$codes[["coords"]][,2]) - 
     somsegm1$codes[["coords"]][,2],
     col = (1:nrow(somsegm1$codes[["coords"]]) %in% pepperunits) + 1,
     labels = 1:nrow(somsegm1$codes[["coords"]]), pos = 1, offset = .8)
box()

## visualize the pixels that are projected to the pepper units
dev.new()
oneclass <- rep("background", 480000)
for (ii in seq(along = pepperunits))
  oneclass[predictions$unit.classif == pepperunits[ii]] <-
    paste("unit", pepperunits[ii])
classcols <- c("lightgray", brewer.pal(length(pepperunits), "Set1"))
plot.picture(factor(oneclass), what = "truth", mycols = classcols,
             scales = list(draw = FALSE))

### Ecology example using Bray-Curtis distance
## definition of the distance measure in C++
BCcode <- '
  #include <Rcpp.h>
  typedef double (*DistanceFunctionPtr)(double *, double *, int, int);
  
  double brayCurtisDissim(double *data, double *codes, int n, int nNA) {
    if (nNA > 0) return NA_REAL;

    double num = 0.0, denom = 0.0;
    for (int i = 0; i < n; i++) {
        num += std::abs(data[i] - codes[i]);
        denom += data[i] + codes[i];
    }
 
    return num/denom;
  }

  // [[Rcpp::export]]
  Rcpp::XPtr<DistanceFunctionPtr> BrayCurtis() {
    return Rcpp::XPtr<DistanceFunctionPtr>(new DistanceFunctionPtr(&brayCurtisDissim));
  } '
library("Rcpp")
sourceCpp(code = BCcode)

## alternatively, use the BC.cpp file in the inst/Distances
## subdirectory of the package:
## sourceCpp(file = paste(path.package("kohonen"), "Distances/BC.cpp", sep = "/"))

## the example from the paper - load the data, and split in training
## and test sets. For the soil data standardization is performed based
## on the means and standard deviations from the training set.
data("varespec", package = "vegan")
data("varechem", package = "vegan")
varechem <- as.matrix(varechem)
varespec <- as.matrix(varespec)
n <- nrow(varechem)
tr.idx <- 1:(n-5)
tst.idx <- (n-4):n
chemmat.tr <- scale(varechem[tr.idx,])
chemmat.tst <- scale(varechem[tst.idx,], 
                     center = colMeans(varechem[tr.idx,]),
                     scale = apply(varechem[tr.idx,], 2, sd))
specmat.tr <- varespec[tr.idx,]
specmat.tst <- varespec[tst.idx,]

## train the map using the training data
set.seed(101)
varesom <- supersom(list(species = specmat.tr, soil = chemmat.tr), 
                    somgrid(4, 3, "hexagonal"),
                    dist.fcts = c("BrayCurtis", "sumofsquares"))

## now get predictions for the test data. Predictions may be based on
## either the soil data, the species data, or both.
soilBasedPred <- predict(varesom, 
                         newdata = list(soil = chemmat.tst))
soilBasedPred$predictions$species[,1:5]
speciesBasedPred <- predict(varesom,
                            newdata = list(species = specmat.tst))
speciesBasedPred$predictions$soil[,1:5]
bothPred <- predict(varesom,
                    newdata = list(species = specmat.tst,
                                   soil = chemmat.tst))
cbind(bothPred$predictions$species[,1:3],bothPred$predictions$soil[,1:3])

### ########################################################################
### Example for mapping X-ray powder patterns to a SOM
## get the data
data("degelder")
mydata <- list(patterns = degelder$patterns,
               CellVol = log(degelder$properties[,"cell.vol"]))
## compile the function definition of the WCCd dissimilarity function
require(Rcpp)
sourceCpp(paste(path.package("kohonen"), "Distances/wcc.cpp", sep = "/"))

## train the map
set.seed(7)
powsom <- supersom(data = mydata,
                   grid = somgrid(6, 4, "hexagonal"),
                   dist.fcts = c("WCCd", "sumofsquares"),
                   keep.data = TRUE)

## show codebook vectors
par(mfrow = c(1,2))
plot(powsom, type = "codes", bgcol = "lightblue", 
     main = c("Diffraction patterns", "Cell volume"))

##  show predictions
cellPreds <- predict(powsom, newdata = mydata, whatmap = "patterns")
names(cellPreds)
cellPreds$predictions$CellVol[1:5,]

