#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop("Correct arguments must be supplied", call.=FALSE)
}

source("./src/main.R")

inputFile <- args[1]
categories <- c("BPO", "CCO", "KO", "MFO")

geneComb <- function(a, b) {
  for (id in categories) {
    id <<- id
    csv <<- inputFile
    
    loadAncestorsAndDescendants(verbose = TRUE)
    loadCsv(verbose = TRUE)
    
    calculateSimilarityOfGenesList(writeResult = TRUE)
  }
}

geneComb()
