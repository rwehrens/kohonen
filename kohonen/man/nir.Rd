\name{nir}
\alias{nir}
\docType{data}
\title{Near-infrared data with temperature effects}
\description{A data object containing near-infrared spectra of ternary
  mixtures of ethanol, water and iso-propanol, measured at five
  different temperatures (30, 40, ..., 70 degrees Centigrade). } 
\references{F. Wulfert , W.Th. Kok, A.K. Smilde: Anal. Chem. 1998, 1761-1767}
\examples{
  data(nir)

  set.seed(3)
  nirnet <- xyf(data = nir$spectra[nir$training,],
                Y = nir$composition[nir$training,],
                xweight=.75,
                grid = somgrid(6, 6, "hexagonal"), rlen=500)
  plot(nirnet, "counts", main="Counts")

  ## Focus on compound 2 (water):
  dev.new(width = 14)
  par(mfrow = c(1,2))
  set.seed(13)
  nirnet <- xyf(data = nir$spectra[nir$training,],
                Y = nir$composition[nir$training, 2],
                grid = somgrid(6, 6, "hexagonal"), rlen=500)
  water.xyf <- predict(nirnet)$prediction
  plot(nirnet, "property", property = water.xyf,
       main="Prediction of water content")
  ## Plot temperatures as circles
  symbols(nirnet$grid$pts[nirnet$unit.classif,] +
          matrix(rnorm(sum(nir$training)*2, sd=.1), ncol=2),
	  circles = (nir$temperature[nir$training] - 20)/250,
	  inches = FALSE, add = TRUE)

  ## Model temperatures
  set.seed(13)
  nirnet2 <- xyf(data = nir$spectra[nir$training,],
                 Y = nir$temperature[nir$training],
                 xweight=.25,
                 grid = somgrid(6, 6, "hexagonal"), rlen=500)
  temp.xyf <- predict(nirnet2)$prediction

  plot(nirnet2, "property", property = temp.xyf,
       palette.name = rainbow,
       main="Prediction of temperatures")
  ## Plot concentrations of water as circles
  symbols(nirnet2$grid$pts[nirnet2$unit.classif,] +
          matrix(rnorm(sum(nir$training)*2, sd=.1), ncol=2),
	  circles = 0.05 + 0.4 * nir$composition[nir$training,2],
	  inches = FALSE, add = TRUE)
}
\keyword{datasets}
