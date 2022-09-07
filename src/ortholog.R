#install.packages("rjson")
rm(list = ls())
library("rjson")

ancestors <- list()
descendants <- list()
ignoredCategories <- c(
  "09190 Not Included in Pathway or Brite"
)
ignoredSubcategories <- c(
  "09112 Not included in regular maps"
)

getCategoryId <- function(category) {
  categoryId <- strsplit(category$name, " ")[[1]][1]
  categoryId <- paste('CAT:', categoryId, sep = '')
  return(categoryId)
}

getSubcategoryId <- function(subcategory) {
  subcategoryId <- strsplit(subcategory$name, " ")[[1]][1]
  subcategoryId <- paste('SUB:', subcategoryId, sep = '')
  return(subcategoryId)
}

getPathwayId <- function(pathway) {
  pathwayId <- toupper(tail(strsplit(pathway$name, split=" ")[[1]],1))
  pathwayId <- gsub('[[]','',pathwayId)
  pathwayId <- gsub('[]]','',pathwayId)
  return(pathwayId)
}

getKoTermId <- function(koTerm) {
  koTermId <- strsplit(koTerm$name, " ")[[1]][1]
  koTermId <- paste('KO:', koTermId, sep = '')
  return(koTermId)
}

addAncestor <- function(koTermId, ancestorToAdd) {
  id <- strsplit(koTermId, ':')[[1]][2]
  ancestor <- strsplit(ancestorToAdd, ':')[[1]][2]
  category <- tolower(strsplit(ancestorToAdd, ':')[[1]][1])
  
  ancestors[[id]][[category]] <<- unique(c(ancestors[[id]][[category]], ancestor))
}

addDescendant <- function(koTermId, descendantToAdd) {
  id <- strsplit(koTermId, ':')[[1]][2]
  descendant <- strsplit(descendantToAdd, ':')[[1]][2]
  
  descendants[[id]] <<- unique(c(descendants[[id]], descendant))
}

orderAncestorsAndDescendantsNames <- function() {
  ancestors <<- ancestors[order(names(ancestors))]
  descendants <<- descendants[order(names(descendants))]
}

extractAncestorsAndDescendants <- function() {
  ko <- fromJSON(file = "./data/ko.json")
  
  for (category in ko$children) {
    if(category$name %in% ignoredCategories) {
      next
    }
    
    categoryId <- getCategoryId(category)
    
    for (subcategory in category$children) {
      if(subcategory$name %in% ignoredSubcategories) {
        next
      }
      
      subcategoryId <- getSubcategoryId(subcategory)
      
      addDescendant(categoryId, subcategoryId)
      addAncestor(subcategoryId, categoryId)
      
      for (pathway in subcategory$children) {
        pathwayId <- getPathwayId(pathway)
        
        addDescendant(subcategoryId, pathwayId)
        addAncestor(pathwayId, subcategoryId)
        
        for (koTerm in pathway$children) {
          koTermId <- getKoTermId(koTerm)

          addDescendant(pathwayId, koTermId)
          addAncestor(koTermId, pathwayId)
        }
      }
    }
  }
}

writeOutput <- function() {
  dput(ancestors, "./data/KO_ancestor.txt")
  dput(descendants, "./data/KO_descendant.txt")
}

extractAncestorsAndDescendants()
orderAncestorsAndDescendantsNames()
writeOutput()
