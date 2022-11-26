#' Converts human genes to mouse orthologs, and vice-versa
#' 
#' Orthologs converter
#' 
#' ceeb
#'
#' @param geneIDs character vector of gene IDs
#' @param from organism from which to convert gene IDs
#' @param to organism to which to convert geneIDs
#' @param identifier type of gene identifier
#' 
#' @import crayon
#' @import biomaRt
#' @import clusterProfiler 
#' 
#' @return Data frame of mapped gene IDs and annotations
#' 
#' @examples
#' getOrthologs()
#'
#' @export
getOrthologs <- function(geneIDs, from = c('mouse','human'), to = c('human','mouse'), identifier = 'entrezgene_id') {
  
  if (identifier == "ensembl") {
    cat(yellow("Converting ensembl to entrez...\n"))
    
    if (from ==  "mouse") {
      org <- "org.Mm.eg.db"   
    } else if (from ==  'human') {
      org <- "org.Hs.eg.db"   
    }
    
    geneIDs <-  bitr(geneIDs, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org)$ENTREZID
    
    }
  
  # Get marts
  
  if (!exists("marts")) {

    
  cat(yellow("Getting Mart...\n"))
  
    martspecies = list(human = 'hsapiens_gene_ensembl', 
                  mouse = 'mmusculus_gene_ensembl')
  
  marts = list()
  
  for (i in names(martspecies)) {
    
    marts[[i]] <- tryCatch(useMart(biomart = 'ENSEMBL_MART_ENSEMBL',martspecies[i])) #get mart
    if (inherits(marts[[i]], "error")) {
      next 
        marts[[i]] <-  tryCatch(useMart(biomart = 'ensembl', martspecies[i]))
      }
    }
  .GlobalEnv$marts <- marts
  }
  # Map orthologs
  
  cat(blue("Mapping orthologs...\n"))

  orthoGenes <- getLDS(attributes = 'entrezgene_id',
                    filters = identifier, 
                    values = geneIDs, 
                    mart = marts[[from]], 
                    attributesL = 'entrezgene_id',
                    martL = marts[[to]], 
                    uniqueRows = TRUE)

  colnames(orthoGenes) <- c("ENTREZID","ENTREZID2")
  write.csv(orthoGenes, paste0(from,"2",to,"_orthologs.csv"))
  .GlobalEnv$orthoGenes <- orthoGenes
  }



















