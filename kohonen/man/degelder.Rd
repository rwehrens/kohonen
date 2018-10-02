\name{degelder}
\alias{degelder}
\docType{data}
\title{Powder pattern data by Rene de Gelder}
\description{
  X-ray powder patterns of 131 crystallographic structures, contributed
  by Rene de Gelder.
}
\usage{data(degelder)}
\format{
  This yields a list with three components: the first component,
  '"patterns"', is a matrix of 131 rows and 441 variables, containing
  the powder patterns; the second component is "thetas", the 2theta values
  at which intensities have been measured. The final component,
  '"properties"', gives information on the crystallographic properties
  of the structures.
}
\source{Rene de Gelder, Institute of Molecules and Materials, Radboud
  University Nijmegen.  
}
\examples{
\dontrun{
data(degelder)
mydata <- list(patterns = degelder$patterns,
               CellVol = log(degelder$properties[,"cell.vol"]))

## custom distance function
require(Rcpp)
sourceCpp(system.file("Distances", "wcc.cpp", package = "kohonen"))
set.seed(7)
powsom <- supersom(data = mydata, grid = somgrid(6, 4, "hexagonal"),
                   dist.fcts = c("WCCd", "sumofsquares"),
                   keep.data = TRUE)
summary(powsom)		   
}
}
\keyword{datasets}
