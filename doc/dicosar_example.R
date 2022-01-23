## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("lhe17/dicosar")

## ----echo=TRUE----------------------------------------------------------------
library(dicosar)
n = 100
x1 = matrix(rnorm(2*n),nrow=n)
x2 = matrix(rnorm(2*n),nrow=n)

## ----echo=TRUE----------------------------------------------------------------
p = dicosar(x1,x2)
p

## ----echo=TRUE----------------------------------------------------------------
p = dicosar(x1,x2,delta=TRUE)
p

## ----echo=TRUE----------------------------------------------------------------
x1 = matrix(rnorm(4*n),nrow=n)
x2 = matrix(rnorm(4*n),nrow=n)
p = dicosar(x1,x2)
p

