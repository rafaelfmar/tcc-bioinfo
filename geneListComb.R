#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Correct arguments must be supplied", call.=FALSE)
}

source("./src/main.R")

inputFolder <- args[1]
outputFile <- args[2]
categories <- c("BPO", "CCO", "KO", "MFO")

geneListComb <- function() {
  for (id in categories) {
    id <<- id
    files <-  list.files(inputFolder, pattern="*.csv", full.names=TRUE)
    
    loadAncestorsAndDescendants(verbose = TRUE)
    
    for (path in files) {
      append <- !(path == files[1])
      csv <<- path
      loadCsv(verbose = TRUE, append = append)
    }
    
    csv <<- outputFile
    calculateSimilarityOfGenesList(clusterResult = TRUE)
  }
}

geneListComb()
