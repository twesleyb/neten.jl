#!/usr/bin/env Rscript

library(neten)
#library(microbenchmark)

data(butterfly)

netw = neten(butterfly)
#microbenchmark("neten.R"=neten(butterfly),times=10)
