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
set.seed(1)
geldermap <- wccsom(degelder$patterns, somgrid(6, 6, "hexagonal"))
options(digits = 3)
summary(geldermap, "unit", nr=1, properties = degelder$properties)
}
}
\keyword{datasets}
