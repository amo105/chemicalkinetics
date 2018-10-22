library(RcppArmadillo)
library(Rcpp)

setwd("E:/Data/Chemometrics/chemopackage")

RcppArmadillo.package.skeleton(name="chemopac", example_code = FALSE,code_files = c("chemopac.R"),force=FALSE)

compileAttributes("chemopac")

file.copy(from="chemo.cpp", to="E:/Data/Chemometrics/chemopackage/chemopac/src/chemopac.cpp")


remove.packages("chemopac")