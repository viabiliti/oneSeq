#' Performs Genewise Trimmed Means of M-values normalisation of RNAseq read counts
#' 
#' Essy geTMM normalisation
#'
#' ceeb
#'
#' @param counts data frame count/mRNA expression matrix with samples in columns and genes in rows
#' @param species character accepts either mouse or human
#' @param expname character name of the experiment
#' @param id character gene id type
#' 
#' @import dplyr
#' @import biomaRt
#' @import clusterProfiler 
#' @import EDASeq
#' @import edgeR
#' 
#' @return data frame of normalised mRNA counts
#' 
#' @examples 
#' geTMMnorm()
#' 
#' @export
geTMMnorm <- function(counts = counts, species = NULL, expname = 'norm', meta=meta, id = "ENSEMBL") {
  
  cat(yellow("Preparing data for normalisation...\n"))

  if (!any(grepl(colnames(meta),pattern = "condition"))){
    stop("No column names 'condition' in meta dataframe!\n")
  }
  
  if (is.null(species) | species ==  "mouse") {
    martspecies <- 'mmusculus_gene_ensembl'
    orgdb <- 'mmu'
    org <- "org.Mm.eg.db"
  } else if (species ==  'human') {
    martspecies  <- c("hsapiens_gene_ensembl")
    orgdb <- "hg38"
    org <- "org.Hs.eg.db"
  }
  
  if (id == 'hgnc') {
  

  ens <- bitr(rownames(counts), fromType = "SYMBOL", toType = "ENSEMBL",
                                           OrgDb = org) # symbol to ensembl

  counts2 <- data.frame(SYMBOL = rownames(counts),counts)
  join <- join <- inner_join(counts2,ens, by = 'SYMBOL')
  join2 <- join[!duplicated(join[,max(length(join))]),] %>% as.data.frame()
  rownames(join2) <- join2$ENSEMBL
  counts <- join2[,-c(1,max(length(join)))]
  
  }
  
  ##  Length
  
  cat(yellow("Fetching gene length and GC content information...\n"))

  ens_char <- as.character(rownames(counts)) #ensembl IDs and character
  gc_length <- getGeneLengthAndGCContent(ens_char,orgdb, mode = "org.db") #get gc and length information
  
  gc_length <- gc_length %>% as.data.frame() # as df
  
  gc <- data.frame(ensembl = rownames(gc_length),gc = gc_length %>% 
                     dplyr::select(gc) %>% 
                     dplyr::transmute(gc = gc*100)) # as a percentage
  
  len <- data.frame(ensembl = rownames(gc_length),length = gc_length$length) #because we need per kbp
  len[is.na(len$length),2] <- median(len[!is.na(len$length),][,2]) # make NAs median length
  
  ## GeTMM
  
  cat(blue("Normalising counts...\n"))
  counts <- mutate_all(counts,as.numeric)
  rpk <- (10e3*counts/len$length) # reads per 1000 bp
  y <- DGEList(counts,group = factor(meta$condition),remove.zeros = TRUE)
  y <- calcNormFactors(y) 
  norm_counts <- cpm(y) %>% as.data.frame()
  cat(blue("Writing normalised count matrix to csv...\n"))
  write.csv(norm_counts,paste0(expname,"_geTMM_counts.csv"))
  .GlobalEnv$geTMMcounts <- norm_counts
  
  }


# require(EDASeq)
# require(dplyr)
# require(edgeR)
# require(crayon)

