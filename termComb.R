#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Correct arguments must be supplied", call.=FALSE)
}

source("./src/main.R")

categories <- c("BPO", "CCO", "KO", "MFO")

termComb <- function(a, b) {
  for (id in categories) {
    id <<- id
    loadAncestorsAndDescendants(verbose = FALSE)
    calcSimilarityOfTerms(a, b)
  }
}

termComb(args[1], args[2])
