### Demo presenting the code to replicate the second example from the
### 2017 JSS publication on the kohonen package v. 3.0, mapping X-ray
### powder patterns to a SOM using a specialized distance function.

require("kohonen")
require("Rcpp")

## get the data
data("degelder")
mydata <- list(patterns = degelder$patterns,
               CellVol = log(degelder$properties[,"cell.vol"]))

par(mfrow = c(2,1))
matplot(degelder$thetas, t(degelder$patterns[c(66, 67),]),
        main = "Spacegroup C2",
        type = "l", lty = 1, ylab = "response",
        xlab = expression(paste("2", theta)))
matplot(degelder$thetas, t(degelder$patterns[c(89, 91),]),
        main = "Spacegroup P-1", ylab = "response",
        xlab = expression(paste("2", theta)),
        type = "l", lty = 1)

## compile the function definition of the WCCd dissimilarity function
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

