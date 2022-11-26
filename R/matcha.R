#' Match sample order between count matrix and metadata
#' 
#' Match sample order between count matrix and metadata
#' 
#' Ceeb
#'
#' @param counts data frame, counts matrix
#' @param meta data frame, sample information 
#' 
#' @import dplyr
#' 
#' @return data frames with sample names in right order
#' 
#' @examples
#' matcha()
#' 
#' @export
matcha <- function(counts = counts, meta = meta) {
  
  # In case fix names did something weird    
  
  colnames(counts) <- gsub("^X","",colnames(counts),perl = T)
  
  
  
  if (is.character(counts[[1]])) {
  
  # Check if gene names are first column instead of rownmaes
    cat("Gene names in first column. Removing and assigning to rownames...\n")
    
    rownames(counts) <- counts[,1]
    counts <- counts[-1]
    
  }
  
  # Same number of samples?
  
  if (length(counts) > nrow(meta)) {
    
    cat("More samples in count matrix than metadata. Removing unmatched samples.")
    
    counts <- counts[,colnames(counts) %in% meta$id]
    
  }
  
  if (length(counts) < nrow(meta)) {
    
    cat("More samples in metadata than counts matrix. Removing unmatched samples.")
    
    meta <- meta[meta %in% colnames(counts),]
  
  }
  
  # Check if sample names are in write order
  
  if (any(!meta$id == colnames(counts)) == TRUE) {
    
    cat("Sample are in wrong order. Fixing...\n")
    
    meta <- mutate_all(meta,as.character)
    meta <- meta[order(meta$condition),]
    counts <- counts[,meta$id]
    
  }  

  # Meta as factors
  meta <- mutate_all(meta,as.factor)
  stopifnot(colnames(counts) == meta$id)
  
  cat("Matching complete!\n")
  
  .GlobalEnv$meta <- meta
  .GlobalEnv$counts <- counts
  }
