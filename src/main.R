source("./src/util.R")

isA <- 0.4
partOf <- 0.3
cat <- 0.2
sub <- 0.3
path <- 0.4
br <- 0.4

ancestors <- list()
descendants <- list()
genesList <- list()

getSimPath <- function() {
  simPath <- "./data/result_eg_gene_list_in.txt"
  simPath <- paste(simPath, id, "sim", sep = "")
  return(simPath)
}

getPrePath <- function() {
  prePath <- "./data/result_eg_gene_list_in.txt"
  prePath <- paste(prePath, id, "pre", sep = "")
  return(prePath)
}

getClsPath <- function() {
  clsPath <- "./data/result_eg_gene_list_in.txt"
  clsPath <- paste(clsPath, id, "cls", sep = "")
  return(clsPath)
}

loadAncestorsAndDescendants <- function(verbose = TRUE) {
  ancestors <<- list()
  descendants <<- list()

  
  pathAncestors <- paste("./data/", id, "_ancestor.txt", sep = "")
  pathDescendants <- paste("./data/", id, "_descendant.txt", sep = "")
  
  if (verbose) {
    Log(paste('Loading ancestors from ', pathAncestors, sep = ""))
    Log(paste('Loading descendants from ', pathDescendants, sep = ""))
  }

  ancestors <<- dget(pathAncestors)
  descendants <<- dget(pathDescendants)
  
  if (verbose) {
    Log(paste('Ancestors and descendants loaded, number of GO terms: ', length(ancestors) , sep = ""))
  }
}

loadCsv <- function(verbose = TRUE, append = FALSE) {
  if (!append) {
    genesList <<- list()
  }

  if (verbose) {
    Log(paste('Loading genes CSV file from ', csv, sep = ""))
  }
  
  genesCsv <- read.csv(csv, header = F)
  
  for (i in 1:length(genesCsv$V1)) {
    gene <- genesCsv[i,1]
    
    for (j in 2:length(genesCsv[i,])) {
      term <- genesCsv[i,j]
      if (nchar(term) > 0 && !is.null(ancestors[[term]])) {
        genesList[[gene]] <<- c(term, genesList[[gene]])
      }
    }
  }

  if (verbose) {
    Log(paste('CSV loaded, number of genes: ', length(genesList) , sep = ""))
  }
}

getAncestors <- function(term) {
  termAncestors <- c(
    term, 
    ancestors[[term]]$is_a, 
    ancestors[[term]]$part_of, 
    ancestors[[term]]$cat, 
    ancestors[[term]]$sub, 
    ancestors[[term]]$path, 
    ancestors[[term]]$br
  )
  
  for (ancestor in termAncestors) {
    if(ancestor == term){
      next
    }
    ancestorAncestors <- getAncestors(ancestor)
    termAncestors <- c(termAncestors, ancestorAncestors)
  }
  
  termAncestors <- unique(termAncestors)
  
  return(termAncestors)
}

getNumberOfChildren <- function(term) {
  numberOfChildren <- 0
  
  if (!is.null(descendants[[term]])){
    numberOfChildren <- length(descendants[[term]][!is.na(descendants[[term]])])
  }
  
  return(numberOfChildren)
}

getTermsRelation <- function(a, b) {
  relation <- NA
  
  if (b %in% ancestors[[a]]$is_a) {
    relation <- "IS_A"
  } else if (b %in% ancestors[[a]]$part_of) {
    relation <- "PART_OF"
  } else if (b %in% ancestors[[a]]$cat) {
    relation <- "CAT"
  } else if (b %in% ancestors[[a]]$sub) {
    relation <- "SUB"
  } else if (b %in% ancestors[[a]]$path) {
    relation <- "PATH"
  } else if (b %in% ancestors[[a]]$br) {
    relation <- "BR"
  }
  
  return(relation)
}

getLinkWeight <- function(relation) {
  if (is.na(relation) || is.null(relation)) {
    return(0)
  }
  
  linkWeight <- switch (relation,
                        "IS_A" = isA,
                        "PART_OF" = partOf,
                        "CAT" = cat,
                        "SUB" = sub,
                        "PATH" = path,
                        "BR" = br,
                        0
  )
  
  return(linkWeight)
}

getSimilarityMatrix <- function(numLines, names) {
  matrix = c()
  
  matrix <- matrix(numeric(numLines * numLines), nrow = numLines, ncol = numLines, data = NA) 
  rownames(matrix) <- names
  colnames(matrix) <- names
  for(row in 1:nrow(matrix)) {
    for(col in 1:ncol(matrix)) {
      if(row == col) {
        matrix[row, col] <- 1
      }
    }
  }
  
  return(matrix)
}

calcSvalues <- function(term) {
  termAncestors <- getAncestors(term)
  weightMatrix <- matrix(0, nrow = length(termAncestors), ncol = length(termAncestors), dimnames = list(termAncestors, termAncestors))
  Svalues <- list()
  Svalues[[term]] <- 1
  linkWeight <- 0
  numberOfChildren <- 0
  sValue <- 0
  
  for(row in termAncestors) {
    for(col in termAncestors) {
      # Set weight to 0 if col_element == row_element and skip
      if (col == row) {
        weightMatrix[row, col] <- 0
        next
      } 
      
      # Set link weight based on relation; Skip if GO terms aren't linked
      relation <- getTermsRelation(row, col)
      linkWeight <- getLinkWeight(relation)
      if (linkWeight == 0) {
        next
      }
      
      # Calc weight
      numberOfChildren <- getNumberOfChildren(col)
      weight <- round((1 / (0.67 + numberOfChildren) + linkWeight), 3)
      weightMatrix[row, col] <- weight
      
      # Calc S-value
      if (is.null(Svalues[[col]])) {
        Svalues[[col]] <- 0
      }
      
      sValue <- round((Svalues[[row]] * weight), 3)
      
      if(sValue > Svalues[[col]]){
        Svalues[[col]] <- sValue
      }
    }
  }
  
  return(Svalues)
}

calcSimilarityOfTerms <- function(a, b, verbose = TRUE) {
  if (a == b) {
    return(1)
  }

  similarity <- 0
  a.sValues <- calcSvalues(a)
  b.sValues <- calcSvalues(b)
  commonAncestors <- intersect(getAncestors(a), getAncestors(b))
  commonAncestors.sValues <- 0
  
  for (commonAncestor in commonAncestors) {
    commonAncestors.sValues <- commonAncestors.sValues + a.sValues[[commonAncestor]] + b.sValues[[commonAncestor]]
  }
  
  similarity <- commonAncestors.sValues / (sum(unlist(a.sValues)) + sum(unlist(b.sValues)))
  similarity <- round(similarity, 3)
  
  if (verbose) {
    Log(paste( a, b, id, similarity))
  }
  
  return (similarity)
}

calcSimilarityOfGenes <- function(a, b, verbose = TRUE) {
  a.name <- a[1]
  b.name <- b[1]
  a.terms <- a[2:length(a)]
  b.terms <- b[2:length(b)]
  a.sum <- 0
  b.sum <- 0
  similarity <- 0
  
  for(term in a.terms) {
    similarities <- mapply(calcSimilarityOfTerms, a = term, b = b.terms, verbose = FALSE)
    a.sum <- max(unlist(similarities)) + a.sum
  }
  
  for(term in b.terms) {
    similarities <- mapply(calcSimilarityOfTerms, a = term, b = a.terms, verbose = FALSE)
    b.sum <- max(unlist(similarities)) + b.sum
  }
  
  similarity <- (a.sum + b.sum)/(length(a.terms) + length(b.terms))
  
  if (verbose) {
    Log(
      paste(
        "Similarity between", a.name, "and", b.name, id, "=", similarity
      )
    )
  }
  
  similarity <- round(similarity, digits = 3)
  
  return(similarity)
}

clusterGenes <- function() {
  shell(paste("apcluster.exe", getSimPath(), getPrePath(), getClsPath()))
  
  conn <- file(getClsPath(),open="r")
  lines <-readLines(conn)
  cluster <- list()
  dest <- paste("./data/result/", rev(strsplit(csv, "/")[[1]])[1], id, "_result.txt", sep = "")
  
  for (i in 1:length(genesList)){
    gene <- names(genesList)[[i]]
    target <-  names(genesList)[[as.numeric(lines[[i]])]]
    
    if (is.null(cluster[[target]])) {
      cluster[[target]] <- ""
    }
    
    cluster[[target]] <- trimws(paste(gene, cluster[[target]]))
  }
  
  write(id, dest, append = FALSE)
  lapply(cluster, write, dest, append = TRUE)
  close(conn)
}

writeSimResult <- function (sim) {
  dest <- paste("./data/result/", rev(strsplit(csv, "/")[[1]])[1], id, "_result.txt", sep = "")
  fileConn = file(dest)
  writeLines(sim, fileConn)
  close(fileConn)
}

calculateSimilarityOfGenesList <- function(writeResult = FALSE, clusterResult = FALSE) {
  if(length(genesList) == 0) {
    return()
  }
  
  numLines <- length(genesList)
  similarityMatrix <- getSimilarityMatrix(numLines, names(genesList))
  sim <- c()
  
  for (x in 1:(numLines - 1)) {
    for (y in seq(from = (x + 1), to = numLines)) {
      append <- !((x == 1) && (y == 2))
      a <- c(names(genesList)[x], genesList[[x]])
      b <- c(names(genesList)[y], genesList[[y]])
      a.name <- a[1]
      b.name <- b[1]
      
      similarity <- calcSimilarityOfGenes(a, b, FALSE)
      
      if((length(similarity) > 0) && (similarity > 0)) {
        similarityMatrix[a.name, b.name] <- similarity
        similarityMatrix[b.name, a.name] <- similarity
        write(paste(x, y, similarity), getSimPath(), append = append)
      } else {
        write(paste(x, y, "NA"), getSimPath(), append = append)
      }
      
      sim <- c(sim, paste(a.name, b.name, id, similarity))
      
      Log(
        paste(
          a.name, b.name, id, similarity
        )
      )
    }
  }
  
  if (writeResult) {
    writeSimResult(sim)
  }
  
  if (clusterResult) {
    median <- median(similarityMatrix[which(!is.na(similarityMatrix))])
    
    for (x in 1:numLines) {
      append <- (x > 1)
      notNa <- length(na.omit(similarityMatrix[x,]))
      
      if (notNa != 1) {
        write(median, getPrePath(), append = append)
      }
    }
    
    clusterGenes()
  }
}
