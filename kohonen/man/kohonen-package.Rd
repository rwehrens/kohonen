\name{kohonen-package}
\alias{kohonen-package}
\alias{kohonen}
\docType{package}
\title{
\packageTitle{kohonen}
}
\description{
\packageDescription{kohonen}
}
\details{
The kohonen package implements several forms of self-organising maps
(SOMs). Online and batch training algorithms are available; batch
training can also be done in parallel. Multiple data layers may be
presented to the training algorithm, with potentially different distance
measures for each layer. The overall distance is a weighted average of
the layer distances. Layers may be selected through the \code{whatmap}
argument, or by providing a weight of zero. The basic function is
\code{supersom}; \code{som} is simply a wrapper for SOMs using just one
layer (the classical form).

New data may be mapped to a trained SOM using the \code{map.kohonen}
function. Function \code{predict.kohonen} will map data to the SOM, and
will return predictions (i.e., average values for winning units) for
those layers that are not in the new data object.

Several visualisation methods are available in function
\code{plot.kohonen}. 
  
\packageIndices{kohonen}
}
\author{
  \packageAuthor{kohonen}
  
  Maintainer: \packageMaintainer{kohonen}
}
\references{
R. Wehrens and J. Kruisselbrink: Flexible Self-Organising Maps in
kohonen 3.0. Journal of Statistical Software, 87, 7 (2018).}
\keyword{ package }
