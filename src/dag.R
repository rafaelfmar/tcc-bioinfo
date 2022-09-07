alt_ids_BPO <- list()
alt_ids_CCO <- list()
alt_ids_MFO <- list()
ancestors_BPO <- list()
ancestors_CCO <- list()
ancestors_MFO <- list()
descendants_BPO <- list()
descendants_CCO <- list()
descendants_MFO <- list()

addElementToList <- function(namespace, type, id, element) {
  if (type == "ancestor") {
    list <- switch(  
      namespace,  
      "biological_process" = ancestors_BPO,  
      "cellular_component" = ancestors_CCO,  
      "molecular_function" = ancestors_MFO
    ) 
  }
  else if (type == "descendant") {
    list <- switch(  
      namespace,  
      "biological_process" = descendants_BPO,  
      "cellular_component" = descendants_CCO,  
      "molecular_function" = descendants_MFO
    ) 
  }
  
  if(is.null(list[[id]])) {
    list[[id]] <- NA
  }

  list[[id]] <- c(list[[id]][!is.na(list[[id]])], element)

  if (namespace == "biological_process" && type == "ancestor") {
    ancestors_BPO <<- list
  }
  else if (namespace == "biological_process" && type == "descendant") {
    descendants_BPO <<- list
  }
  else if (namespace == "cellular_component" && type == "ancestor") {
    ancestors_CCO <<- list
  }
  else if (namespace == "cellular_component" && type == "descendant") {
    descendants_CCO <<- list
  }
  else if (namespace == "molecular_function" && type == "ancestor") {
    ancestors_MFO <<- list
  }
  else if (namespace == "molecular_function" && type == "descendant") {
    descendants_MFO <<- list
  }    
}

populateAltIds <- function(target, alt_ids) {
  for (id in names(alt_ids)) {
    for (alt_id in alt_ids[[id]]) {
      target[[alt_id]] <- target[[id]]
    }
  }
  
  return(target)
}

writeOutput <- function(list, path) {
  orderedList <- list[order(names(list))]
  dput(orderedList, path)
}

extractData <- function() {
  obo <- "./data/go.obo"
  goAsText <- readChar(obo, file.info(obo)$size)
  goTerms <- strsplit(goAsText, split = "\n\n")[[1]]

  for (i in 2:(length(goTerms) - 1)) {
    goTerm <- goTerms[i]
    lines <- strsplit(goTerm, "\n")[[1]]
    id <- "" 
    namespace <- ""
    ancestor <- ""
    is_a <- c()
    part_of <- c()
    alt_ids <- c()
    is_obsolete <- FALSE
    
    for (j in 2:length(lines)) {
      line <- lines[j]
      
      if (j == 2) {
        id <- strsplit(line, split = " ")[[1]][2]
      } 
      else if (j == 3) {
        next
      } 
      else if (j == 4) {
        namespace <- strsplit(line, split = " ")[[1]][2]
      }
      else if (grepl("is_obsolete: true", line)) {
        is_obsolete <- TRUE
        break
      }
      else if (grepl("alt_id: ", line)) {
        alt_id <- strsplit(line, split = " ")[[1]][2]
        alt_ids <- c(alt_ids, alt_id)
      }
      else if (grepl("is_a: ", line)) {
        ancestor <- strsplit(line, split = " ")[[1]][2]
        is_a <- c(is_a, ancestor)
      }
      else if (grepl("relationship: part_of GO:", line)) {
        ancestor <- strsplit(line, split = " ")[[1]][3]
        part_of <- c(part_of, ancestor)
      }
      
      if (nchar(ancestor) > 0) {
        addElementToList(namespace, "descendant", ancestor, id)
        ancestor <- ""
      }
    }

    if (is_obsolete) {
      next
    }
    else if (namespace == "biological_process") {
      if(is.null(descendants_BPO[[id]])) {
        descendants_BPO[[id]] <<- NA
      }
      alt_ids_BPO[[id]] <<- alt_ids
    }
    else if (namespace == "cellular_component") {
      if(is.null(descendants_CCO[[id]])) {
        descendants_CCO[[id]] <<- NA
      }
      alt_ids_CCO[[id]] <<- alt_ids
    }
    else if (namespace == "molecular_function") {
      if(is.null(descendants_MFO[[id]])) {
        descendants_MFO[[id]] <<- NA
      } 
      alt_ids_MFO[[id]] <<- alt_ids
    }
    
    addElementToList(namespace, "ancestor", id, list(is_a = is_a, part_of = part_of))
  }
  
  ancestors_BPO <<- populateAltIds(ancestors_BPO, alt_ids_BPO)
  ancestors_CCO <<- populateAltIds(ancestors_CCO, alt_ids_CCO)
  ancestors_MFO <<- populateAltIds(ancestors_MFO, alt_ids_MFO)
  descendants_BPO <<- populateAltIds(descendants_BPO, alt_ids_BPO)
  descendants_CCO <<- populateAltIds(descendants_CCO, alt_ids_CCO)
  descendants_MFO <<- populateAltIds(descendants_MFO, alt_ids_MFO)
}

extractData()
writeOutput(ancestors_BPO, "./data/BPO_ancestor.txt")
writeOutput(ancestors_CCO, "./data/CCO_ancestor.txt")
writeOutput(ancestors_MFO, "./data/MFO_ancestor.txt")
writeOutput(descendants_BPO, "./data/BPO_descendant.txt")
writeOutput(descendants_CCO, "./data/CCO_descendant.txt")
writeOutput(descendants_MFO, "./data/MFO_descendant.txt")
