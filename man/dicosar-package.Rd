\name{dicosar-package}
\alias{dicosar-package}
\docType{package}
\title{
\packageTitle{dicosar}
}
\description{
\packageDescription{dicosar}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{dicosar}
\packageIndices{dicosar}
Robust and fast two-sample test of equal correlations or correlation matrices
}
\author{
\packageAuthor{dicosar}

Maintainer: \packageMaintainer{dicosar}
}
\references{
He et al. Robust and fast two-sample test of equal correlations with an application to differential co-expression, (2022).
}

\keyword{ package }

\examples{
library(dicosar)
n = 100
x1 = matrix(rnorm(2*n),nrow=n)
x2 = matrix(rnorm(2*n),nrow=n)
dicosar(x1,x2)
}
