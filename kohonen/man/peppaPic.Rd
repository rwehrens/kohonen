\name{peppaPic}
\alias{peppaPic}
\docType{data}
\title{Synthetic image of a pepper plant with peppers}
\description{A data matrix with four columns representing a 600 by 800
  image of a pepper plant. Each row is a pixel in the image. The first
  column is the class label; the other columns contain the RGB values.
}
\usage{data("peppaPic")}
\source{
http://dx.doi.org/10.4121/uuid:884958f5-b868-46e1-b3d8-a0b5d91b02c0
}
\references{
  This is image 10039 from a set of 10,500 images described in

  Barth R, IJsselmuiden J, Hemming J, and van Henten E (2017). "Data
  Synthesis Methods for Semantic Segmentation in Agriculture. A Capsicum
  annuum Dataset." Submitted.
}
\examples{
data(peppaPic)
head(peppaPic)
## show ground truth per pixel
image(t(matrix(peppaPic[,1], 600, 800))[,600:1], col = rainbow(10))
}
\keyword{datasets}
