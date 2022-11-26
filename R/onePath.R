#' Pipleine for Geneset Enrichement and Pathway analysis
#' 
#' Leverages the wonferful setRank and gsCluster packages for GSE and Pathway analysis
#'
#' ceeb
#'
#' @param sigDEG list of data frames output of lazyDGE function or DE = 'noiseq' of lazySeq function
#' @param mods list of WGCNA modules outputed from oneNet function 
#' @param species character, accepts mouse and human
#' @param orthogenes Data frame output of 'getOrthologs' function
#' @param orthospecies character, takes 'mouse' or 'human'
#' @param setRank logical, perform setRank GSEA?
#' @param CTDbase logical, include Chemical and Toxigenomics Data base in GSEA?
#' @param chemDis logical, if CTDbase, include chemical-gene-disease interactions in GSEA?
#' @param chemIX logical, if CTDbase, include chemical-gene interactions in GSEA?
#' @param msigDB logical, include msigDB in GSEA?
#' @param msigCat list of desired msigDB sub-categories 
#' @param ppiGSEA logical, map GSEA results to protein-protein interaction network?
#' @param clusterGS logical, perform hierarchal clustering on GSEA results?
#' @param expname character, name of experiment
#' @param seed integer, random seed
#' @param directory character, project path
#' @param cores integer, number of cpu cores to use
#' @param cytoscape logical, export networks to cytoscape?
#' 
#' @import EDASeq
#' @import clusterProfiler
#' @import magrittr
#' @import ggplot2
#' @import data.table
#' @import readr
#' @import dplyr
#' @import genefilter
#' @import biomaRt
#' @import cluster
#' @import flashClust
#' @import foreach
#' @import doParallel
#' @import mixOmics
#' @import tweeDEseq
#' @import swamp
#' @import flashClust
#' @import colorspace
#' @import factoextra
#' @import ggsci
#' @import calibrate
#' @import gtools
#' @import stringr
#' @import stringi
#' @import textclean
#' @import msigdbr
#' @import gmt
#' @import htmlwidgets
#' @import enrichplot
#' @import crayon
#' @import igraph
#' @import tidyr
#' @import huxtable
#' @import plyr
#' @import purrr
#' @import RCy3
#' 
#' @examples
#' onePath()
#' 
#' @export
lazyPaths <- function(sigDEGs = NULL,
                       mods = NULL,
                       species = NULL,
                       orthogenes = NULL,
                       orthospecies = NULL,
                       setRank = NULL,
                       CTDbase = NULL,
                       chemDis = NULL,
                       chemIX = NULL,
                       msigDB = NULL,
                       msigCat = NULL,
                       minGenes = 50,
                       useRank = TRUE,
                       geneSetPCut = 0.01,
                       geneSetFDRCut = 0.05,
                       ppiGSEA = NULL,
                       clusterGS = NULL,
                       cytoscape = NULL,
                       expname =  NULL,
                       seed = NULL,
                       directory = NULL,
                     cores = NULL) {
  
    # crayons
    wild <- red $ bold
    yay <- green $ bold
    calc <- cyan $ bold
    prep <- silver $ bold
    fyi <- white $ bold
  
    cat(prep("Preparing data for gene pathway analysis...\n"))
    
    setwd(directory)
    
    if (is.null(species) | species ==  "mouse") {
      spec <- 'Mus musculus'
      org <- "org.Mm.eg.db"
      symb <- "mgi_symbol"
    } else if (species ==  'human') {
      org <- "org.Hs.eg.db"
      spec  <- 'Homo sapiens'
      symb <- "hgnc"
    }
    
    ## Prepare geneinputs
    
    if (!is.null(setRank)) {
      
      if (is.null(mods)) {
       
        if (is.null(sigDEGs)) {
  
          tmp <- list.files(pattern = "^sigDEG.rds", recursive = T)
    
      if (!isEmpty(tmp)) {
    
          sigDEGs <- read_rds(tmp)
    
    } else {
      
    
      stop("Please provide a data frame or list of dataframes output from DGEA!")
      
      }
  
        rm(tmp)
    
        }
    
    # if (!is.null(mods)) {
    #   
    #   sigDEGs <- mods
    # }  
    
    cat(prep("Converting ENSEMBL to ENTREZ...\n"))
    
    #sigDEGs = allGenes
    list <- lapply(sigDEGs, function(x) bitr(x$ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID",
                                            OrgDb = org)) # ensembl to entrez
    
    sigList2 <- sigDEGs
    
    for (i in seq_along(sigList2)) {
      colnames(sigList2[[i]])[1] <- "ENSEMBL"
    }
    
    sigList2 <- lapply(sigList2, function(x)
      lapply(x, function(y) as.character(y))
      %>% as.data.frame())
    
    # Remove those that didn't map
    
    entrez_sigList <- list()
    for (i in seq_along(sigList2)) {
      entrez_sigList[[i]] <- inner_join(list[[i]],sigList2[[i]], by = "ENSEMBL")
    }
    
    names(sigList2) <- names(sigDEGs)
    names(entrez_sigList) <- names(sigList2)
    
    entrez_sigList <-  entrez_sigList[sapply(entrez_sigList, function(x) dim(x)[1]) > 0]#remove empty
    .GlobalEnv$entrez_sigList <- entrez_sigList
    write_rds(entrez_sigList,"entrez_sigDEGs.rds")
    
    # Orth df
    
    if (!is.null(orthogenes)) {
      
    cat(prep("Preparing orthologs data frame...\n"))
    
    tmp <- list.files(pattern = 'orth_entrez_DEG.rds', recursive = T)
    if (!isEmpty(tmp)){
      entrez_DEG <- read_rds(tmp)
      rm(tmp)
    } else {
    
    cat(prep(paste0("Substituting ",species," genes ", "for ", orthospecies, " genes...\n")))
    
      entrez_DEG <- mutate_all(entrez_DEG, as.character)
    
    joins <- list()
    for (i in names(entrez_sigList)) {
      joins[[i]] <- inner_join(entrez_DEG, dplyr::select(entrez_sigList[[i]],c(ENTREZID,log2FC,condition,pval.adjus)), by = "ENTREZID")
      colnames(joins[[i]])[c(1,2)] <- c("ENTREZID","ENTREZID2")
    
      }  
    
    # Annotations
    
    if (orthospecies == "mouse") {
      martspecies <- 'mmusculus_gene_ensembl'
      orgdb <- 'mmu'
      symb = "mgi_symbol"
      spec <- 'Mus musculus'
    } else if (orthospecies ==  "human") {
      martspecies  <- c("hsapiens_gene_ensembl")
      orgdb <- "hg38"
      symb = "hgnc_symbol"
      spec <- 'Homo sapiens'
    }
    
    cat(prep("Collecting gene annotations...\n"))
    cat(magenta("Getting Mart...\n"))
    
    for (i in 1:1) {
      mart = tryCatch(useMart(biomart = 'ENSEMBL_MART_ENSEMBL', martspecies[i])) #get mart
      if (inherits(mart, "error")) {
        next 
        mart = tryCatch(useMart(biomart = 'ensembl', martspecies[i]))
      }}
    
    ids <- list()
    
    for (i in names(joins)) {
    
    ids[[i]] <- getBM(attributes = c('entrezgene_id',
                                     'ensembl_gene_id',                                
                                      symb,
                                      'gene_biotype',
                                      'ensembl_peptide_id',
                                      'wikigene_description'),
                      filters = "entrezgene_id",
                      values = joins[[i]]$ENTREZID2,
                      mart = mart) 
    }
        
    names(ids) <- names(joins)
  
    anno <- list()
    
    for (i in names(ids)) {
      colnames(ids[[i]])[3] <- "ENTREZID2"
      ids[[i]]$ENTREZID2 <- as.character(ids[[i]]$ENTREZID2)
      ids[[i]] <- inner_join(ids[[i]],joins[[i]], by = "ENTREZID2")
      anno[[i]] <- ids[[i]][!duplicated(ids[[i]]$ENTREZID2),] 
      anno[[i]] <- anno[[i]] %>% mutate(select_biotypes = gsub("protein_coding",">protein_coding",anno[[i]]$gene_biotype)) #">" anchor for later filtering
      anno[[i]]$gene_biotype <- forcats::fct_explicit_na(anno[[i]]$gene_biotype, "unannotated_element") #nice funtion to replace NAs with a factor
      anno[[i]]$select_biotypes <- forcats::fct_explicit_na(anno[[i]]$select_biotypes, ">unannotated_element") #nice function to replace NAs with a factor
      colnames(anno[[i]])[3] <- "ENTREZID"
      colnames(anno[[i]]) <- gsub(pattern = c(".x"), replacement = "", x = colnames(anno[[i]]))
    }
    
    .GlobalEnv$entrez_DEG <- anno
    write_rds(entrez_DEG,"orth_entrez_DEG.rds")
    species <- orthospecies
      }
    }
    ## Prepare network annotations
    
    setwd(directory)
   
    cat(prep("Preparing Network Annotations...\n"))
    
    
      if (!is.null(orthogenes)) {
      tmp <- list.files(pattern = 'orth_entrez_DEG.rds', recursive = T)
      entrez_DEG <- read_rds(tmp)
      species <- orthospecies
      
    } 
  
    # Extract information
   
    if (species ==  "mouse") {
      
      geneListDes <- lapply(entrez_DEG, function(x) data.frame(geneID = x$entrezgene_id,
                                                                   symbol = x$mgi_symbol,
                                                                   p = x$pval.adjus,
                                                                   logFC = x$log2FC
      ))
    } else if (species ==  'human') {
      
      
      geneListDes <- lapply(entrez_DEG, function(x) data.frame(geneID = x$entrezgene_id,
                                                                   symbol = x$hgnc_symbol,
                                                                   p = x$pval.adjus,
                                                                   logFC = x$log2FC
      ))
    }
    
    geneListDes  <- lapply(geneListDes, function(x) x[!duplicated(x$geneID),]) #remove duplicates
    allGenes <- unique(unlist(lapply(geneListDes, function(x) x$geneID)))
    
    
    if (!dir.exists('GSEA')) {
      dir.create('GSEA')
    }
    
    setwd(paste0(directory,"/GSEA"))
    
    write_rds(geneListDes,"geneListDes.rds")
    
    .GlobalEnv$geneListDes <- geneListDes
    
    cat(prep("Preparing Gene Input...\n"))
    
    
    # Sort by pval
    
    pvals <- lapply(entrez_DEG, function(x) x$pval.adjus)
    
    for (i in seq_along(pvals)) {
      names(pvals[[i]]) <- as.character(entrez_DEG[[i]]$ENTREZID)
    }
    
    pvals <- lapply(pvals, function(x) sort(x, decreasing = F)) # sort pvals
    geneinput <- lapply(pvals, function(x) names(x)) #Entrez in order of pvals
    geneinput <- geneinput[sapply(geneinput, function(x) length(x)[1]) >= minGenes] #remove <100
    
    write_rds(geneinput, "geneinput.rds")
    .GlobalEnv$geneinput <- geneinput
    geneListDes <- geneListDes[names(geneinput)]
    test  <- lapply(geneListDes, function(x) lapply(x[3:4], function(y) as.numeric(y)))
    geneListDes2 <- list()
    for (i in seq_along(test)) {
      geneListDes2[[i]] <- data.frame(geneListDes[[i]][,c(1:2)],test[[i]],stringsAsFactors = F)
    }
    
    names(geneListDes2) <- names(geneinput)
    write_rds(geneListDes2, "geneListDes2.rds")
    
    # Fold-Change
    
    fc <- lapply(entrez_DEG, function(x) x$log2FC)
    
    for (i in seq_along(fc)) {
      
      names(fc[[i]]) <- as.character(entrez_DEG[[i]]$ENTREZID)
      
      } 
  
    geneinput <- lapply(list, function(x) x$ENTREZ)

#### IF MODS ####
    
    } else {
      
      list <- lapply(mods, function(x) bitr(x$ENSEMBL, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"),
                                               OrgDb = org)) # ensembl to entrez
      
      modGeneList <- lapply(list, function(x) data.frame(geneID = x$ENTREZID,
                                                               symbol = x$SYMBOL
                                                              
      ))
                                                               
      modGeneList  <- lapply(modGeneList, function(x) x[!duplicated(x$geneID),]) #remove duplicates
      allGenes <- unique(unlist(lapply(modGeneList, function(x) x$geneID)))
      geneinput <- lapply(modGeneList, function(x) x$geneID)
      saveRDS(geneinput,"geneinput.rds")
      write_rds(modGeneList, "modGeneList.rds")
      
    }
  
  # Background genes
    
    setwd(directory)
    
    if (!is.null(orthogenes)) {
      
      allGenes <- orthogenes$ENTREZID2
      
    } else if (!exists("allGenes")) {
    
      tmp <- list.files(pattern = "allGenes.rds", recursive = T)
      
      if (!is.null(tmp)) {
        
        allGenes <- read_rds(tmp)      
        
        }
    }
      
    cat(fyi("Loading last checkpoint...\n"))
    
    if (!is.null(orthospecies)) {
      tmp <- list.files(pattern = 'orth_entrez_DEG.rds', recursive = T)
      entrez_DEG <- read_rds(tmp)
      species <- orthospecies
      spec <- orthospecies
    }
    
    ## Already done?  
    
    tmp <- list.files(pattern = "SetRankResult.rds", recursive = T) 
    if (!isEmpty(tmp)) {
      sraList <- read_rds(tmp)
    } 
    # Descriptions if input
    tmp <- list.files(pattern = "setRank_collection.rds", recursive = T) 
    if (!isEmpty(tmp)) {
      collection <- read_rds(tmp)  
      .GlobalEnv$collection <- collection
    } else {
      
    ## Data bases
    
    setwd(paste0(directory,"/GSEA"))
    
    if (msigDB) {
    
    cat(prep("\nPreparing Molecular Signatures Databases...\n"))
    
    .GlobalEnv$msigCat
    names(msigCat) <- msigCat  
    msigDBList <- list()
    
    for (i in names(msigCat)) {
      msigDBList[[i]] <- msigdbr(species = spec,category = msigCat[[i]]) %>% 
        dplyr::select(gs_id,entrez_gene,gs_name, gs_cat, gs_description)
    }
    
    names(msigDBList) <- msigCat
    .GlobalEnv$msigDBList <- msigDBList
    
    colist <- list()
    for (i in names(msigDBList)) {
      colnames(msigDBList[[i]]) <- c("termID","geneID","termName","dbName","description")
      colist[[i]] <- mutate_all(msigDBList[[i]], function(x) as.character(x))
      replace(colist[[i]]$description,which(colist[[i]]$description == ""),"no description")
    }
  }
    
    ## CTDbase
    
    if (CTDbase) {
      
    cat(prep("Preparing Comparative Toxicogenomics Databases...\n" %+% wild("Downloading...\n")))
    
    urlList <- list("http://ctdbase.org/reports/CTD_chemicals.csv.gz", #download dataases from ctdbase.org
                    "http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz",
                    "http://ctdbase.org/reports/CTD_chemicals_diseases.csv.gz",
                    "http://ctdbase.org/reports/CTD_genes_pathways.csv.gz")
    fileNameList <- list("CTD_chemicals.csv.gz", #call them something
                         "CTD_chemicals_ixns.csv.gz",
                         "CTD_chem_disease_ixns.csv.gz",
                         "CTD_chem_pathways_ixns.csv.gz")
    # Read - Really dumb format
    CTD_list <- list()
    for(i in 1:length(urlList)) {
      download.file(url = urlList[[i]],destfile = fileNameList[[i]])
      CTD_list[[i]] <-  read.csv(fileNameList[[i]],
                                 skip = 27,
                                 blank.lines.skip = T,
                                 skipNul = T,
                                 encoding = "UTF-8",
                                 check.names = T,
                                 na.strings = "#",
                                 strip.white = T)
      CTD_list[[i]] <- CTD_list[[i]][-c(1),]
      CTD_list[[i]] <- lapply(CTD_list[[i]], function(x) gsub("MESH:","",x)) %>%
      as.data.frame()
    }
    
    names(CTD_list) <- fileNameList #rename
    .GlobalEnv$CTD_list <- CTD_list
    # Chemical-gene interactions
    
    if (chemIX) {
      chem_ixs <- data.frame(termID = CTD_list[[2]]$ChemicalID,
                             geneID = CTD_list[[2]]$GeneID,
                             termName = CTD_list[[2]]$X..ChemicalName,
                             dbName = rep("CTD",nrow(CTD_list[[2]])),
                             description = CTD_list[[2]]$Interaction)
      chem_ixs <- mutate_all(chem_ixs, function(x) as.character(x))
      chem_ixs <- chem_ixs[,c("termID","geneID","termName","dbName","description")] %>% as.data.frame()
      chem_ixs$termName <- gsub("\\.","_",chem_ixs$termName)
      #Remove
      .GlobalEnv$chem_ixs <- chem_ixs
      #rm(CTD_list,chemids)
      cat(prep("Clearing RAM...\n"))
      gc()
      
      if (msigDB && CTDbase) {
        colist$CTDix <- chem_ixs
        setwd(paste0(directory,"/GSEA"))
      } else {
        setwd(paste0(directory,"/GSEA"))
        colist <- chem_ixs
        #colist <- colist %>% as.data.frame()
        write_rds(colist,"colist_ctd.rds")
        .GlobalEnv$colist <- colist  
      }
    }
    
    if (chemDis) {
      
    ## Chemical disesse interactions
    chemids <- data.frame(ChemicalID = CTD_list[[2]]$ChemicalID,GeneID = CTD_list[[2]]$GeneID)
    chemids <- chemids[!duplicated(chemids$ChemicalID),]
    
    dis <- data.frame(ChemicalID = CTD_list[[3]]$ChemicalID,
                      ChemicalName = CTD_list[[3]]$X..ChemicalName,
                      description = CTD_list[[3]]$DiseaseName,
                      termID = CTD_list[[3]]$DiseaseID)
    
    #chem_dis2 <- merge(chemids,dis)
    dis <- inner_join(dis,chemids, by = "ChemicalID")
    #chem_dis2 <- merge(dis,chemids, by = "ChemicalID")
    dis <- dis[-1]
    colnames(dis)
    colnames(dis) <- c("description","termName","termID","geneID")
    dis <- data.frame(dis,dbName = rep("CTDis",nrow(dis)))
    dis <- mutate_all(dis, function(x) as.character(x))
    dis <- dis[,c("termID","geneID","termName","dbName","description")] %>% as.data.frame()#reorder columns
    dis$termName <- gsub("\\.","_", dis$termName)
    .GlobalEnv$dis <- dis
    
    gc()
    
    if (msigDB && CTDbase) {
      colist$CTDdis <- dis
      setwd(paste0(directory,"/GSEA"))
      colist <- rbindlist(colist)
    } else {
      setwd(paste0(directory,"/GSEA"))
      colist <- dis
      colist <- colist %>% as.data.frame()
      write_rds(colist,"colist_ctd.rds")
      .GlobalEnv$colist <- colist  
        }
      }
    }
    
    # Write
    
    colist <- lapply(colist, function(x) as.data.frame(x))
    colist <- rbindlist(colist)
    write_rds(colist,"colist.rds")
    
    ## Build Collection
    
    cat(prep("Building Set Collection. Will take it's sweet time...\n"))
    
    setwd(directory)
    options(mc.cores = cores) 
    gc()
    
    if (msigDB == FALSE && CTDbase == TRUE) {
    tmp <- list.files(pattern = "colist_ctd.rds", recursive = T)
    colist <- read_rds(tmp)
    } else if (msigDB == TRUE) {
    tmp <- list.files(pattern = "colist.rds", recursive = T)
    colist <- read_rds(tmp)
    }
    
    setwd(paste0(directory,"/GSEA"))
    
    #colist <- colist %>% as.data.frame()
    #colnames(colist) <- gsub(pattern = "H.",replacement = "",x = colnames(colist))
  
    start_time <- Sys.time()
    
    collection <-  buildSetCollection(... = colist,referenceSet = allGenes, maxSetSize = 500)
    
    end_time <- Sys.time()
    
    capture.output((end_time - start_time),file = "time2buildcollection.txt")
    
    cat(blue("Set Collection Complete. Congratulation on your patience\n"))
    
    write_rds(collection,"setRank_collection.rds")
    }

    ## Run setRank Analysis
    
    tmp <- list.files(pattern = "SetRankResult.rds", recursive = T) 

    if (isEmpty(tmp)) {
    
    setwd(paste0(directory,"/GSEA"))
    
    cat(blue("Running SetRank GeneSet Enrichment Analysis...\n"))
    
    options(mc.cores = cores) 
    sraList <- list()  
    gc()
    
    start_time <- Sys.time()
    
    if (!is.null(mods)) {
      
      useRank = FALSE
      
    }
    
    # Loop GSEA
    
    for (i in seq_along(geneinput)) {
      
      sraList[[i]] <- setRankAnalysis(geneIDs = geneinput[[i]],
                                      use.ranks = useRank,
                                      setCollection = collection,
                                      setPCutoff = geneSetPCut,
                                      fdrCutoff = geneSetFDRCut)
      }
    
    end_time <- Sys.time()
    
    capture.output((end_time - start_time),file = "time4setRankanalysis.txt")
    
    
    cat(yay("SetRank GSEA Complete!\n"))
    names(sraList) <- names(geneinput)
    .GlobalEnv$sraList <- sraList
    cat(prep("Writing setRank Output...\n"))
    write_rds(sraList,"SetRankResult.rds")
    } else {
    setwd(directory)
    tmp <- list.files(pattern = "SetRankResult.rds", recursive = T)
    sraList <- read_rds(tmp)
    }
  
    ## Export GeneSet networks
    
    setwd(directory)
  
    tmp <- list.files(pattern = "setrank", recursive = T)
    
    if (isEmpty(tmp)) {
    
    cat(blue("Preparing to Export GeneSet Networks...\n"))  
    
    if (msigDB == FALSE && CTDbase == TRUE && !exists("colist")) {
      tmp <- list.files(pattern = "colist_ctd.rds", recursive = T)
      colist <- read_rds(tmp)
    } else if (msigDB == TRUE && !exists("colist")) {
      tmp <- list.files(pattern = "colist.rds", recursive = T)
      colist <- read_rds(tmp)
    }
    
    # if (!exists("colist")) {
    # tmp <- list.files(pattern = "colist", recursive = T)
    # colist <- read_rds(tmp) %>% as.data.frame()
    # .GlobalEnv$colist <- colist
    # }
    # 
    if (!exists("geneinput")) {
    tmp <- list.files(pattern = "geneinput.rds", recursive = T)
    #tmp <- tmp[-1]
    geneinput <- read_rds(tmp)
    .GlobalEnv$geneinput <- geneinput
    }
    
    setwd(paste0(directory,"/GSEA"))
    collection <- read_rds("setRank_collection.rds")
    
    # Define labels    
    if (!is.null(orthogenes)) {
    if (orthospecies == 'human') {
        org <- "org.Hs.eg.db"
        spec  <- 'Homo sapiens'
        symb <- "hgnc_symbol"
    } else {
        spec <- 'Mus musculus'
        org <- "org.Mm.eg.db"
        symb <- "mgi_symbol"
      }
    } else {
      
    if (species == 'human') {
      org <- "org.Hs.eg.db"
      spec  <- 'Homo sapiens'
      symb <- "hgnc_symbol"
    } else if (species == 'mouse') {
      spec <- 'Mus musculus'
      org <- "org.Mm.eg.db"
      symb <- "mgi_symbol"
      }
    }
    # Create converter
    
    entrez2symbol <-  createIDConverter(org,"ENTREZID","SYMBOL")
    
    if (!dir.exists("Basic_geneSet_networks")) {
        dir.create("Basic_geneSet_networks")
    }
    
    setwd(paste0(directory,"/GSEA/Basic_geneSet_networks"))
    
    cat(calc("Exporting GeneSet Networks...\n"))
    
  ## Remove empty results
  
    cat(fyi("Removing empty results\n"))
      
  sraTest <- sraList[sapply(sraList, function(x) length(V(x)$pSetRank)[1]) > 0] #remove empty element
    
  ## Export 
    
  for (i in seq_along(geneinput)) {
      
      exportMultipleResults(sraTest, 
                            geneinput[i], 
                            collection,
                            IDConverter = entrez2symbol, 
                            outputPath = getwd())
      
    }
    
    cat(yay("Export Complete!\n"))
  
    ## Pathways to csv
    
    cat(prep("Writing pathways...\n"))  
    
    
    
    tbs <- list.files(pattern = "_pathways.txt",recursive = T) # read pathways
    sraResList <- list()
    
    for (i in 1:length(tbs)) {
      sraResList[[i]] <- read.table(tbs[i],sep = "\t",header = T)
    }
    
    names(sraResList) <- tbs
    
    huxlist <- list()
    formats <- list()
    
    for (i in seq_along(sraResList)) {
      
      huxlist[[i]] <- sraResList[[i]] %>% as_hux()
      formats[[i]] <-  huxlist[[i]] %>% theme_article()
      quick_html(formats[[i]],file = paste0(names(sraResList[i]),"_setRank.html"))
      write.csv(formats[[i]],file = paste0(names(sraResList[i]),"_setRank.csv"))
      
        }
      }
    }
  
   ## PPI Maps
  
  if (ppiGSEA) {
      
    setwd(paste0(directory,"/GSEA"))
    
      if (!dir.exists("PPI_geneSet_networks")) {
        dir.create("PPI_geneSet_networks")
        setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
      }
  
    tmp <- list.files(pattern = "gene_net_styles.viz", recursive = T) # checkpoint
    
    if (isEmpty(tmp)) {
      setwd(directory)
      
      if (!exists("geneinput")) {
      tmp <- list.files(pattern = 'geneinput.rds', recursive = T)
      geneinput <- read_rds(tmp)
      }
      
      setwd(directory)
      
      cat(prep("Reading Node and Edge Lists...\n"))
      
      tbs <- list.files(pattern = "_pathways.txt",recursive = T) # read pathways
      
      sraResList <- list()
      for (i in 1:length(tbs)) {
        sraResList[[i]] <- read.table(tbs[i],sep = "\t",header = T)
      }
      names(sraResList) <- tbs
      .GlobalEnv$sraResList <- sraResList
      mems <- list.files(pattern = "_membership.txt",recursive = T) # Read genemapping
      sraMemList <- list()
      for (i in 1:length(mems)) {
        sraMemList[[i]] <- read.table(mems[i],sep = "\t",header = T)
      }
      names(sraMemList) <- mems
      .GlobalEnv$sraMemList <- sraMemList
      
      # Create interaction network
      setwd(paste0(directory))
      tmp <- list.files(pattern = "ppi.rds", recursive = T)
    
    if (!isEmpty(tmp)) {
        ppi <- read_rds(tmp)
        .GlobalEnv$ppi <- ppi
    } else {
        cat(blue("Constructing PPI network...\n"))
    
      # Mouse
    if (species ==  "mouse") {
        
      cat(prep("Downloading...\n"))
        mippie_ppi <- read_tsv(url("http://cbdm-01.zdv.uni-mainz.de/~galanisl/mippie/downloads/mippie_ppi_v1_0.tsv"))[,c(1,2)]
        mippie_pro <- read_tsv(url("http://cbdm-01.zdv.uni-mainz.de/~galanisl/mippie/downloads/mippie_proteins_v1_0.tsv"))
        colnames(mippie_pro)[2] <- c("symbol")
        ppi <- graph_from_data_frame(mippie_ppi,vertices = mippie_pro)
        write_rds(ppi, 'ppi.rds')
      }
        # Human
        
    if (species ==  "human") {
        cat(prep("Downloading...\n"))
        h <- read_tsv(url("http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt"),col_names = F)
        h <- h[order(h$X2),]
        h$symbol <- gsub(pattern = "[_].*",replacement = "",h$X1) #make symbols
        h <- mutate_all(h, function(x) as.character(x)) #as character
        h_ppi <- data.frame(EntrezA = h$X2,EntrezB = h$X4) #make ppi frame
        h_pro <- data.frame(entrez = h_ppi$EntrezA,symbol = h$symbol) #make vertex frame
        h_pro <- h_pro[-which(duplicated(h_pro$symbol)),] # get rid of duplicates
        h_pro <- h_pro[-which(duplicated(h_pro$entrez)),]
        h_ppi <- h_ppi[h_ppi$EntrezA %in% h_pro$entrez,] #choose the ones which appe
        h_ppi <- h_ppi[h_ppi$EntrezB %in% h_pro$entrez,]
        ppi <- graph_from_data_frame(d = h_ppi,vertices = h_pro) # we made it folks
        write_rds(ppi, 'ppi.rds')
        .GlobalEnv$ppi <- ppi}
      cat(yay("PPI Network Complete!\n"))
      
    }
      
    cat(blue("Mapping GSEA results to ppi networks...\n"))
    
    setwd(paste0(directory,"/GSEA"))
    if (!dir.exists("PPI_geneSet_networks")) {
      dir.create("PPI_geneSet_networks")
      setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
    }
    
    tmp <- list.files(pattern = "gene_net_styles.viz")
    
    if (isEmpty(tmp)) {
      
      setwd(paste0(directory,"/GSEA"))
      sraList <- read_rds("SetRankResult.rds")
      collection <- read_rds("setRank_collection.rds")
      #.GlobalEnv$sraList <- sraList
      #names(sraList) <- names(geneListDes2)
      .GlobalEnv$collection <- collection
      
      
      ## Export netowrks
      
      gc()
    
      cat(blue("Exporting Geneset Networks...\n"))
    
    if (is.null(mods)) {
        
      if(!exists("geneListDes2")) {
        
        geneListDes2 <- read_rds("geneListDes2.rds")
      .GlobalEnv$geneListDes2 <- geneListDes2
      
      }
      
      setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
      
      exportGeneNets(geneListDes2,
                     sraList,
                     collection,
                     ppi,
                     geneSetIDs = NULL,
                     fields = c(geneID = "geneID",
                                symbol = "symbol",
                                logFC = "logFC",
                                p = "p"),
                     outDir = getwd())
    
    } else if (!is.null(mods)) {
    
      
      if (!exists("modGeneList")) {
        
        tmp <- list.files(paste0(directory,"/modGeneList.rds"),recursive = T)
        modGeneList <- read_rds(paste0(directory,"/modGeneList.rds"))
        .GlobalEnv$modGeneList <- modGeneList }
      
        g <- lapply(modGeneList, function(x) {
        
        y <- data.frame(x,logFC = rep(0,nrow(x)),p = rep(0,nrow(x)))
        y
      })
      
      setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
        
      exportGeneNets(g,
                     sraList,
                     collection,
                     ppi,
                     geneSetIDs = NULL,
                     fields = c(geneID = "geneID", 
                                symbol = "symbol",
                                logFC = "logFC",
                                p = "p"),
                     outDir = getwd())
      
      } 
      
      cat(yay("Networks exported!\n"))
      }
    }
    
    # Plot networks DGE
    
    if (is.null(mods)) {
    
    setwd(directory)
    
    tmp <- list.files(pattern = 'ppi_gsea_minimal_network_plots_log2FC.pdf', recursive = T)
    
    if (isEmpty(tmp)) {
    
      cat(blue("Plotting DGE mapped GSEA-PPI Networks...\n"))
    
    netlist <- list()
    nets <- list.files(pattern = ".xml",recursive = T)
    nets <- nets[!grepl(pattern = "viz",nets)]
    nets <- nets[!grepl(pattern = "setrank.xml",nets,ignore.case = F,fixed = T)]
    nets <- nets[!grepl(pattern = "Basic_geneSet_networks",nets,ignore.case = F,fixed = F)]
    for (i in 1:length(nets)) {
      netlist[[i]] <- read_graph(format = "graphml",
                                 file = nets[i])
      
    }
    
    #.GlobalEnv$nets <- nets
    #.GlobalEnv$netlist <- netlist
    
    ## Edges and nodes
    tmp <- list.files(pattern = "geneListDes2.rds", recursive = T)
    #tmp <- tmp[-1]
    geneListDes2 <- read_rds(tmp)
    tmp <- list.files(pattern = "geneinput.rds", recursive = T)
    geneinput <- read_rds(tmp)
    .GlobalEnv$geneinput <- geneinput
    
  
    
    if (!is.null(orthogenes)) {
    
      tmp <- list.files(pattern = "orth_entrez_DEG.rds", recursive = T)
      entrez_sigList <- read_rds(tmp)
      entrez_sigList <- entrez_sigList[names(geneinput)]
      .GlobalEnv$entrez_sigList <- entrez_sigList  
      
      for (i in seq_along(entrez_sigList)) {
        colnames(entrez_sigList[[i]])[7] <- "mouseENTREZID"
      }
      
    } else {
    
    tmp <- list.files(pattern = "entrez_sigDEGs.rds", recursive = T)[1]
    entrez_sigList <- read_rds(tmp)
    entrez_sigList <- entrez_sigList[names(geneinput)]
    
    }
    
    tmp1 <- list.files(pattern="nodelist4net.rds", recursive = T)
    tmp2 <- list.files(pattern="edgelist4net.rds", recursive = T)
    tmp3 <- list.files(pattern="split_nodes.rds", recursive = T)
    
    if (!isEmpty(tmp1) && !isEmpty(tmp2) && !isEmpty(tmp3)) {
      nodelist4net <- readRDS(tmp1)
      edgelist4net <- readRDS(tmp2)
      split_nodes <- readRDS(tmp3)
      
      } else {
    
      
    # Read edges and nodes
    edges <- list()
    nodes <- list()
    
    for (i in 1:length(netlist)) {
      .GlobalEnv$edges[[i]] <- igraph::as_data_frame(netlist[[i]], what = "edges")
      .GlobalEnv$nodes[[i]] <- igraph::as_data_frame(netlist[[i]], what = "vertices")
      edges[[i]] <-  igraph::as_data_frame(netlist[[i]], what = "edges")
      nodes[[i]] <-  igraph::as_data_frame(netlist[[i]], what = "vertices")
    }
    
    names(edges) <- nets 
    names(nodes) <- nets
    names(edges) <- gsub("\\.net.*","",names(edges)) #names
    names(edges) <- gsub("GSEA/PPI_geneSet_networks/","",names(edges)) #names
    names(nodes) <- gsub("\\.net.*","",names(nodes)) #names
    names(nodes) <- gsub("GSEA/PPI_geneSet_networks/","",names(nodes)) #names
    names(nodes)
    nodes <- nodes[!grepl(pattern = "setRankresu",fixed = T,ignore.case = F,x = names(nodes))]
    edges <- edges[!grepl(pattern = "setRankresu",fixed = T,ignore.case = F,x = names(edges))]
      
      
    # Split into edges df    
    
    split_edges <- list()
    for (i in 1:length(edges)) {
      split_edges[[i]] <- data.frame(edges[[i]][,c(1:2)],stringsAsFactors = F)
    } 
    for (i in seq_along(edges)){
      split_edges[[i]] <- data.frame(edges[[i]][1],edges[[i]][2],stringsAsFactors = F)
      colnames(split_edges[[i]]) <- c("ENTREZID","TO")
    }                       
    
    # Join to DGEA description df
    
    entrez_sigList <- entrez_sigList[names(geneinput)]
    lengthvec <- seq(from = 1,by = 1,to = length(geneListDes2))
    e_vec <- rep(seq_along(lengthvec),length(split_edges))
    sig_vec <- rep(seq_along(split_edges),length(geneListDes2)) %>% sort()
    
    
    
    edge_joins <- list()
    for (i in 1:length(e_vec)) {
      edge_joins[[i]] <- inner_join(split_edges[[sig_vec[i]]], entrez_sigList[[e_vec[i]]], by = "ENTREZID")
      names(edge_joins)[i] <- paste0(names(edges[sig_vec[i]]),"_", names(entrez_sigList[e_vec[i]]))
    }
    
    edge_joins <- edge_joins[sapply(edge_joins, function(x) dim(x)[1]) > 0] #remove empty elements
    
    #Same with nodez
    lengthvec <- seq(from = 1,by = 1,to = length(geneListDes2)) #1,2,3,4,5
    conds <- names(geneListDes2)
    names(conds) <- conds
    cond_vec <- conds[e_vec] #cond names
    names(cond_vec) <- cond_vec
  
    
    # Remove semi-exclusive, i.e, paths which only appear in one condition
    
    nodes <- nodes[which(sapply(nodes, function(c) ncol(c)) == max(sapply(nodes, function(x) (ncol(x)))))]
    
    e_vec <- rep(seq_along(lengthvec),length(nodes))
    sig_vec <- rep(seq_along(nodes),length(geneListDes2)) %>% sort()
    cond_vec <- conds[e_vec] #cond names
    names(cond_vec) <- cond_vec
    
    split_nodes <- list()
    
    for (i in 1:length(e_vec)) {
      split_nodes[[i]] <- data.frame(nodes[[sig_vec[i]]][,c(1:2)],
                                     nodes[[sig_vec[i]]][,grepl(pattern = cond_vec[i],
                                                                x = colnames(nodes[[1]]),fixed = TRUE)],
                                     sp = rep(species,nrow(nodes[[sig_vec[[i]]]])),stringsAsFactors = F)
      names(split_nodes)[i] <- paste0(names(nodes[sig_vec[i]]), "_", names(entrez_sigList[e_vec[i]]))
      split_nodes[[i]][,paste0(names(cond_vec[i]),".p") < 0.05]
      colnames(split_nodes[[i]])[c(1,2)] <- c("ENTREZID", "Symbol")
      split_nodes[[i]] <- split_nodes[[i]][!duplicated(split_nodes[[i]]$ENTREZID),]
      
    }
    
    split_nodes <- split_nodes[sapply(split_nodes, function(x) dim(x)[1]) > 0] #remove empty element
    
    # Make networks
    
    edge_joins <- edge_joins[names(edge_joins) %in% names(split_nodes)] # intersect
    split_nodes <- split_nodes[names(split_nodes) %in% names(edge_joins)] 
    
    splitzies <- list()
    for (i in seq_along(split_nodes)) {
      splitzies[[i]] <- left_join(split_nodes[[i]],edge_joins[[i]], by = "ENTREZID")
      splitzies[[i]] <- split_nodes[[i]][!duplicated(split_nodes[[i]]$ENTREZID),]
    }
    
    
    write_rds(split_nodes, "split_nodes.rds")
    
    # Remove dups, make sure genes are in both node and edgelist
    
    edgelist4net <- list()
    nodelist4net <- list()
    edge_joins2 <- list()
    
    for (i in seq_along(edge_joins)) {
      edge_joins2[[i]] <- edge_joins[[i]][!duplicated(edge_joins[[i]]$ENTREZID),]
      edge_joins2[[i]] <- edge_joins2[[i]][!duplicated(edge_joins2[[i]]$TO),]
      split_nodes[[i]] <- split_nodes[[i]][!duplicated(split_nodes[[i]]$ENTREZID),]
      edge_joins2[[i]] <- edge_joins2[[i]][edge_joins2[[i]]$ENTREZID %in% split_nodes[[i]]$ENTREZID,] # filter
      edgelist4net[[i]] <- edge_joins2[[i]][edge_joins2[[i]]$TO %in% split_nodes[[i]]$ENTREZID,]
      nodelist4net[[i]] <- inner_join(split_nodes[[i]],edgelist4net[[i]][-2], by = "ENTREZID")   # filter
    }
    
    names(edgelist4net) <- names(edge_joins)
    names(nodelist4net) <- names(edge_joins)
    write_rds(nodelist4net,"nodelist4net.rds")
    write_rds(edgelist4net,"edgelist4net.rds")
    }
  
    ## Make network from data frame
    
    gsea_nets <- list()    
    for (i in seq_along(edgelist4net)) {
      gsea_nets[[i]] <- graph_from_data_frame(d = edgelist4net[[i]],
                                              vertices = split_nodes[[i]],
                                              directed = F) %>% igraph::simplify(remove.loops = T, remove.multiple = F) %>% 
        delete.vertices(which(igraph::degree(.) == 0))
    } 
    names(gsea_nets) <- names(edgelist4net) #names
      
    # Map log2FC to color
    np <- diverge_hcl(n = 500, "Blue-Red 3") # colour ramp
    ll <- length(edgelist4net)
    lp <- ceiling(length(edgelist4net)/2)
    
    ceil <- list()
    for (i in 1:length(gsea_nets)) {
      ceil[[i]] <- max(abs(as.numeric((E(gsea_nets[[i]])$log2FC)))) # ceil
      ceil[[i]][ceil[[i]] == -Inf] <- 0
    }
    
    
    # If for some reason log2FC is not in edges 
    if (any(is.na(unlist(ceil)))) {
    
      ceil <- list()
      for (i in 1:length(gsea_nets)) {
        ceil[[i]] <- max(abs(as.numeric((V(gsea_nets[[i]])$log2FC)))) # ceil
        ceil[[i]][ceil[[i]] == -Inf] <- 0
      }
    } 
      
    names(ceil) <- names(edgelist4net)
    sf <-  sapply(ceil, max)
    sf <- round(max(abs(sf)),digits = 2)
    
    node_cols <- list()
    for (i in 1:length(gsea_nets)) {
      node_cols[[i]] <- (as.numeric(E(gsea_nets[[i]])$log2FC) + sf) / (2*sf) * 500 # mid poin
    }
    # Same as above
    
    if (any(is.na(unlist(node_cols)))) {
    node_cols <- list()
    for (i in 1:length(gsea_nets)) {
        node_cols[[i]] <- (as.numeric(V(gsea_nets[[i]])$log2FC) + sf) / (2*sf) * 500 # mid poin
      }
    }
     
    names(node_cols) <- names(gsea_nets)
    node_cols1 <- node_cols[sapply(node_cols, function(x) length(x) > 0)]
    node_cols1 <- sapply(node_cols1, function(x) round(x,digits = 0))
    gsea_nets1 <- gsea_nets[names(node_cols1)]
    
    
    saveRDS(gsea_nets1,"DGE_GSEA_PPI_Networks.rds")

    # Plots
    setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
    
    pdf("ppi_gsea_minimal_network_plots_log2FC.pdf", width = 11.7, height = 8.3)
    par(mfrow = c(1, 2),cex.main = 0.75)
    
    for (i in seq_along(gsea_nets1)) {
      plot(gsea_nets1[[i]],
           vertex.color = np[node_cols1[[i]]],
           vertex.frame.color = np[node_cols1[[i]]],
           vertex.label = paste0(V(gsea_nets1[[i]])$Symbol,"\n",round(as.numeric(V(gsea_nets1[[i]])$log2FC),digits = 2)),
           vertex.label.color = "black",
           vertex.label.font = 2,
           vertex.size = 20,
           edge.width  = 2,
           edge.color = "gray50",
           edge.curved = 0.1,
           frame = FALSE,
           main = names(gsea_nets1[i]),
           layout = layout_nicely) }
    graphics.off()
    
      }
    
    if (cytoscape && isEmpty(list.files(pattern = "*DGE_ppi_cytoscape.csv"))) {
      setwd(directory)
      
      if (!exists(gsea_nets1)) {
      
      tmp <- list.files(pattern = "DGE_GSEA_PPI_Networks.rds")
        
      if (!isEmpty(tmp)) {
        
        gsea_nets1 <- readRDS("DGE_GSEA_PPI_Networks.rds")  
      } else {
        
        cat(wild("Can't find iGraph networks!\n"))
          
        }
      
      }
      
      setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
      
      if (!dir.exists("Cytoscape")) {
        dir.create("Cytoscape")
        setwd(paste0(directory,"/GSEA/PPI_geneSet_networks/Cytoscape"))
      }
      
      cat(prep("Exporting graphs to cytoscape...\n"))
      # Pause function
      readkey <- function()
      {
        cat("Open cytoscape, then on your keyboard and in the console, Press [Enter] to continue")
        line <- readline()
      }
      readkey()
      cytoscapePing()
      for (i in seq_along(gsea_nets1)) {
        createNetworkFromIgraph(gsea_nets1[[i]])
      }
      write.csv(names(gsea_nets1), "order_of_conditions_DGE_ppi_cytoscape.csv")
      cat(yay("Graphs exported!\n"))
    
     }  
    
    } else if (!is.null(mods)) {
      
      
    setwd(directory) 
    
    tmp <- list.files(pattern = 'PPI_Modules_networks.pdf', recursive = T)
    
    if (isEmpty(tmp)) {
      
      cat(blue("Plotting DGE mapped GSEA-PPI Networks...\n"))
      
      netlist <- list()
      
      nets <- list.files(pattern = ".xml",recursive = T)
      nets <- nets[!grepl(pattern = "viz",nets)]
      nets <- nets[!grepl(pattern = "setrank.xml",nets,ignore.case = F,fixed = T)]
      nets <- nets[!grepl(pattern = "Basic_geneSet_networks",nets,ignore.case = F,fixed = F)]
      for (i in 1:length(nets)) {
        netlist[[i]] <- read_graph(format = "graphml",
                                   file = nets[i])
        
      }
      
      #.GlobalEnv$nets <- nets
      #.GlobalEnv$netlist <- netlist
        
    if (!exists("modGeneList"))  {
    tmp <- list.files(pattern = "modGeneList.rds", recursive = T)
    modGeneList <- read_rds(tmp)
    tmp <- list.files(pattern = "geneinput.rds", recursive = T)
    geneinput <- read_rds(tmp)
    modGeneList <- modGeneList[names(geneinput)]
    }
    #.GlobalEnv$modGeneList <- modGeneList
    
    tmp1 <- list.files(pattern="nodelist4net.rds", recursive = T)
    tmp2 <- list.files(pattern="edgelist4net.rds", recursive = T)
    tmp3 <- list.files(pattern="split_nodes.rds", recursive = T)
    
    if (!isEmpty(tmp1) && !isEmpty(tmp2) && !isEmpty(tmp3)) {
      nodelist4net <- readRDS(tmp1)
      edgelist4net <- readRDS(tmp2)
      split_nodes <- readRDS(tmp3)
      
    } else {
      
      
      # Read edges and nodes
      edges <- list()
      nodes <- list()
      
      for (i in 1:length(netlist)) {
        .GlobalEnv$edges[[i]] <- igraph::as_data_frame(netlist[[i]], what = "edges")
        .GlobalEnv$nodes[[i]] <- igraph::as_data_frame(netlist[[i]], what = "vertices")
        edges[[i]] <-  igraph::as_data_frame(netlist[[i]], what = "edges")
        nodes[[i]] <-  igraph::as_data_frame(netlist[[i]], what = "vertices")
      }
      
      names(edges) <- nets 
      names(nodes) <- nets
      names(edges) <- gsub("\\.net.*","",names(edges)) #names
      names(edges) <- gsub("GSEA/PPI_geneSet_networks/","",names(edges)) #names
      names(nodes) <- gsub("\\.net.*","",names(nodes)) #names
      names(nodes) <- gsub("GSEA/PPI_geneSet_networks/","",names(nodes)) #names
      #names(edges)
      #names(nodes)
      nodes <- nodes[!grepl(pattern = "setRankresu",fixed = T,ignore.case = F,x = names(nodes))]
      edges <- edges[!grepl(pattern = "setRankresu",fixed = T,ignore.case = F,x = names(edges))]
      
      
      # Split into edges df    
      
      split_edges <- list()
      
      for (i in 1:length(edges)) {
        split_edges[[i]] <- data.frame(edges[[i]][,c(1:2)],stringsAsFactors = F)
      } 
      
      for (i in seq_along(edges)){
        split_edges[[i]] <- data.frame(edges[[i]][1],edges[[i]][2],stringsAsFactors = F)
        colnames(split_edges[[i]]) <- c("geneID","TO")
      }                       
      
      # Join to DGEA description df
      
      modGeneList <- modGeneList[names(geneinput)]
      lengthvec <- seq(from = 1,by = 1,to = length(modGeneList))
      e_vec <- rep(seq_along(lengthvec),length(split_edges))
      sig_vec <- rep(seq_along(split_edges),length(modGeneList)) %>% sort()
    
      edge_joins <- list()
      for (i in 1:length(e_vec)) {
        edge_joins[[i]] <- inner_join(split_edges[[sig_vec[i]]], modGeneList[[e_vec[i]]], by = "geneID")
        names(edge_joins)[i] <- paste0(names(edges[sig_vec[i]]),"_", names(modGeneList[e_vec[i]]))
      }
      
      edge_joins <- edge_joins[sapply(edge_joins, function(x) dim(x)[1]) > 0] #remove empty elements
      
      # Same with nodez
      
      lengthvec <- seq(from = 1,by = 1,to = length(modGeneList)) #1,2,3,4,5
      conds <- names(modGeneList)
      names(conds) <- conds
      cond_vec <- conds[e_vec] # cond names
      names(cond_vec) <- cond_vec
      
      # Remove semi-exclusive, i.e, paths which only appear in one condition
      
      nodes <- nodes[which(sapply(nodes, function(c) ncol(c)) == max(sapply(nodes, function(x) (ncol(x)))))]
      
      e_vec <- rep(seq_along(lengthvec),length(nodes))
      sig_vec <- rep(seq_along(nodes),length(modGeneList)) %>% sort()
      cond_vec <- conds[e_vec] #cond names
      names(cond_vec) <- cond_vec
      
      
      
      split_nodes <- list()
      
      for (i in 1:length(e_vec)) {
        split_nodes[[i]] <- data.frame(nodes[[sig_vec[i]]][,c(1:2)],
                                       nodes[[sig_vec[i]]][,grepl(pattern = cond_vec[i],
                                                                  x = colnames(nodes[[1]]),fixed = TRUE)],
                                       sp = rep(species,nrow(nodes[[sig_vec[[i]]]])),stringsAsFactors = F)
        
        names(split_nodes)[i] <- paste0(names(nodes[sig_vec[i]]), "_", names(modGeneList[e_vec[i]]))
        
        split_nodes[[i]][,paste0(names(cond_vec[i]),".p") < 0.05]
        
        colnames(split_nodes[[i]])[c(1,2)] <- c("geneID", "Symbol")
        
        split_nodes[[i]] <- split_nodes[[i]][!duplicated(split_nodes[[i]]$geneID),]
        
      }
      
      split_nodes <- split_nodes[sapply(split_nodes, function(x) dim(x)[1]) > 0] #remove empty element
      
      # Make networks
      
      edge_joins <- edge_joins[names(edge_joins) %in% names(split_nodes)] # intersect
      split_nodes <- split_nodes[names(split_nodes) %in% names(edge_joins)] 
      
      splitzies <- list()
      
      for (i in seq_along(split_nodes)) {
        splitzies[[i]] <- left_join(split_nodes[[i]],edge_joins[[i]], by = "geneID")
        splitzies[[i]] <- split_nodes[[i]][!duplicated(split_nodes[[i]]$geneID),]
      }
      
      
      write_rds(split_nodes, "split_nodes.rds")
      
      # Remove dups, make sure genes are in both node and edgelist
      
      edgelist4net <- list()
      nodelist4net <- list()
      edge_joins2 <- list()
      
      for (i in seq_along(edge_joins)) {
        edge_joins2[[i]] <- edge_joins[[i]][!duplicated(edge_joins[[i]]$geneID),]
        edge_joins2[[i]] <- edge_joins2[[i]][!duplicated(edge_joins2[[i]]$TO),]
        split_nodes[[i]] <- split_nodes[[i]][!duplicated(split_nodes[[i]]$geneID),]
        edge_joins2[[i]] <- edge_joins2[[i]][edge_joins2[[i]]$geneID %in% split_nodes[[i]]$geneID,] # filter
        edgelist4net[[i]] <- edge_joins2[[i]][edge_joins2[[i]]$TO %in% split_nodes[[i]]$geneID,]
        nodelist4net[[i]] <- inner_join(split_nodes[[i]],edgelist4net[[i]][-2], by = "geneID")   # filter
      }
      
      names(edgelist4net) <- names(edge_joins)
      names(nodelist4net) <- names(edge_joins)
      write_rds(nodelist4net,"nodelist4net.rds")
      write_rds(edgelist4net,"edgelist4net.rds")
    }
    
    ## Make network from data frame
  
    gsea_nets <- list()    
    
    for (i in seq_along(edgelist4net)) {
      gsea_nets[[i]] <- graph_from_data_frame(d = edgelist4net[[i]],
                                              vertices = split_nodes[[i]],
                                              directed = F) %>% igraph::simplify(remove.loops = T, 
                                                                                 remove.multiple = F) %>% 
        delete.vertices(which(igraph::degree(.) == 0))
    } 
    
    names(gsea_nets) <- names(edgelist4net) #names
    
    ## Remove empty symbols
    
    symlist <- list()
    idx <- list()
    
    for (i in seq_along(gsea_nets)) {
      
      symlist[[i]] <-  paste0(V(gsea_nets[[i]])$Symbol)
      #idx[[i]] <- !which(length(symlist[[i]]) < 1)
      idx[[i]] <- ifelse(length(symlist[[i]]) < 1, FALSE,TRUE )
    }
    
    u <- unlist(idx)
    gsea_nets <- gsea_nets[u]
    
    
    saveRDS(gsea_nets,"Mods_PPI_GSEA_Networks.rds")
    
    
    ####### Maybe add joining to noiseq_df 
    
    # # Map log2FC to color
    # 
    # np <- diverge_hcl(n = 500, "Blue-Red 3") # colour ramp
    # ll <- length(edgelist4net)
    # lp <- ceiling(length(edgelist4net)/2)
    # 
    # ceil <- list()
    # for (i in 1:length(gsea_nets)) {
    #   ceil[[i]] <- max(abs(as.numeric((E(gsea_nets[[i]])$log2FC)))) # ceil
    #   ceil[[i]][ceil[[i]] == -Inf] <- 0
    # }
    # 
    # 
    # # If for some reason log2FC is not in edges 
    # if (any(is.na(unlist(ceil)))) {
    #   
    #   ceil <- list()
    #   for (i in 1:length(gsea_nets)) {
    #     ceil[[i]] <- max(abs(as.numeric((V(gsea_nets[[i]])$log2FC)))) # ceil
    #     ceil[[i]][ceil[[i]] == -Inf] <- 0
    #   }
    # } 
    # 
    # names(ceil) <- names(edgelist4net)
    # sf <-  sapply(ceil, max)
    # sf <- round(max(abs(sf)),digits = 2)
    # 
    # node_cols <- list()
    # for (i in 1:length(gsea_nets)) {
    #   node_cols[[i]] <- (as.numeric(E(gsea_nets[[i]])$log2FC) + sf) / (2*sf) * 500 # mid poin
    # }
    # # Same as above
    # 
    # if (any(is.na(unlist(node_cols)))) {
    #   node_cols <- list()
    #   for (i in 1:length(gsea_nets)) {
    #     node_cols[[i]] <- (as.numeric(V(gsea_nets[[i]])$log2FC) + sf) / (2*sf) * 500 # mid poin
    #   }
    # }
    # 
    # names(node_cols) <- names(gsea_nets)
    # node_cols1 <- node_cols[sapply(node_cols, function(x) length(x) > 0)]
    # node_cols1 <- sapply(node_cols1, function(x) round(x,digits = 0))
    # gsea_nets1 <- gsea_nets[names(node_cols1)]
    
    
    # Plots
    setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
    
    pdf("PPI_Modules_networks.pdf", width = 11.7, height = 8.3)
    par(mfrow = c(1, 2),cex.main = 0.75)
    
    for (i in seq_along(gsea_nets)) {
      plot(gsea_nets[[i]],
           #vertex.color = np[node_cols1[[i]]],
           #vertex.frame.color = np[node_cols1[[i]]],
           vertex.label = paste0(V(gsea_nets[[i]])$Symbol),
           vertex.label.color = "black",
           vertex.label.font = 2,
           vertex.size = 20,
           edge.width  = 2,
           edge.color = "gray50",
           edge.curved = 0.1,
           frame = FALSE,
           main = names(gsea_nets[i]),
           layout = layout_nicely) 
      }
    
    dev.off()
    
        }  
      }  
    
    if (cytoscape && isEmpty(list.files(pattern = "*mods_ppi_cytoscape.csv"))) {
      
      setwd(directory)
      
      if (!exists(gsea_nets1)) {
        
        tmp <- list.files(pattern = "Mods_PPI_GSEA_Networks.rds")
        
        if (!isEmpty(tmp)) {
          
          gsea_nets1 <- readRDS("Mods_PPI_GSEA_Networks.rds")  
        
          } else {
          
          cat(wild("Can't find iGraph networks!\n"))
          
        }
        
      }
      setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
      
      if (!dir.exists("Cytoscape")) {
        dir.create("Cytoscape")
        setwd(paste0(directory,"/GSEA/PPI_geneSet_networks/Cytoscape"))
      }
      cat(prep("Exporting graphs to cytoscape...\n"))
      # Pause function
      readkey <- function()
      {
        cat("Open cytoscape, then on your keyboard and in the console, Press [Enter] to continue")
        line <- readline()
      }
      readkey()
      cytoscapePing()
      for (i in seq_along(gsea_nets1)) {
        createNetworkFromIgraph(gsea_nets1[[i]])
      }
      write.csv(names(gsea_nets1), "order_of_conditions_mods_ppi_cytoscape.csv")
      cat(yay("Graphs exported!\n"))
      
    }  
    
    
#########################################  
## Gene Set Clustering
    
    if (clusterGS) {
    
    setwd(directory)
    tmp <- list.files(pattern = "gs_clustList.rds", recursive = T)
  
    if (isEmpty(tmp)) {
    
    cat(prep("Reading gene sets for clustering..."))  
    
    if (!is.null(orthogenes)) {
      
      tmp <- list.files(pattern = "orth_entrez_DEG.rds", recursive = T)
      entrez_sigList <- read_rds(tmp)
      entrez_sigList <- entrez_sigList[names(geneinput)]
      
      for (i in seq_along(entrez_sigList)) {
        colnames(entrez_sigList[[i]])[7] <- "mouseENTREZID"
        .GlobalEnv$entrez_sigList <- entrez_sigList
      }
      
      if (orthospecies == 'human') {
        
        symb = 'hgnc_symbol'
        
      } else {
        
        symb = 'mgi_symbol'
      }
    } else {
      
      tmp <- list.files(pattern = "entrez_sigDEGs.rds", recursive = T)
      entrez_sigList <- read_rds(tmp)
      #entrez_sigList <- entrez_sigList[names(geneinput)]
      #.GlobalEnv$entrez_sigList <- entrez_sigList
      
      if (species == "human") {
        
        symb = 'hgnc_symbol'
        
      } else {
        
        symb = 'mgi_symbol'
        
      }
    }
    
    netlist <- list()
    nets <- list.files(pattern = ".xml",recursive = T) 
    nets <- nets[!grepl(pattern = "viz",nets)]
    nets <- nets[!grepl(pattern = "setrank.xml",nets,ignore.case = F,fixed = T)]
    
    for (i in 1:length(nets)) {
      netlist[[i]] <- read_graph(format = "graphml", # read as graphs
                                 file = nets[i])
      
    }
    #.GlobalEnv$nets <- nets
    #.GlobalEnv$netlist <- netlist
    
    ## Prepare gene membership and interaction
    
    setwd(paste0(directory,"/GSEA/Basic_geneSet_networks"))
    paths <- read.csv("pathways.txt", sep = "\t")
    paths <- dplyr::select(paths,-c("score"))
    tmp <- list.files(pattern = "membership")
    memList <- list()
    for (i in seq_along(tmp)) {
      memList[[i]] <- read.csv(tmp[i], sep = "\t")
    }
    names(memList) <- gsub(pattern = "_membership.txt", "", tmp)
    
    ## Melt and split paths
    
    pathMelt <- data.table::melt(paths)
    pathSplit <- split(pathMelt, pathMelt$variable)
    
    for (i in names(pathSplit)) {
      colnames(pathSplit[[i]]) <- c("name","description","condition","score")
    }
    
    memListT <- list()
    
    for (i in seq_along(memList)) {
      memListT[[i]] <- t(memList[[i]]) %>% as.data.frame()
      colnames(memListT[[i]]) <- memListT[[i]][1,]
      memListT[[i]] <- memListT[[i]][-1,]
      memListT[[i]] <- data.frame(name = rownames(memListT[[i]]), memListT[[i]])
    }
    
    
    names(memListT) <- gsub(pattern = "_membership.txt", "", tmp)
    memListT <- memListT[sapply(memListT, function(x) dim(x)[1]) > 0] #remove empt
    pathSplit <- pathSplit[names(memListT)] 
    
    pathJoins <- list()
    
    for (i in seq_along(memListT)) {
      pathJoins[[i]] <- left_join(pathSplit[[i]],memListT[[i]], by = "name")
    }
    
    names(pathJoins) <- names(memListT)
    
    idx <- list()
    pj <- list()
    
    for (i in 1:length(pathJoins)) {
      pj[[i]] <- t(pathJoins[[i]]) %>% as.data.frame()      
      pj[[i]] <- data.frame(symb = rownames(pj[[i]]),pj[[i]])
      #idx[[i]] <- lapply(pj[[i]], function(x) which(x == "X"))
    }
      
    z <- list()
    x <- list()
    for (i in seq_along(pj)) {
      
      z[[i]] <- pj[[i]][c(1,3),] %>% as.data.frame()
      x[[i]] <- lapply(pj[[i]][-c(1,3),], function(y) stri_trans_totitle(y)) %>% as.data.frame()
      x[[i]] <- rbind(z[[i]],x[[i]])
      x[[i]][c(1,2,3),] <- x[[i]][c(1,3,2),]
      rownames(x[[i]])  <- x[[i]][,1]
    }
    pj <- x
    rm(x)
  
  # X marks the gene
    
    idx <- list()
    for (i in 1:length(pj)) {
    
      idx[[i]] <- lapply(pj[[i]], function(x) which(x == "X"))
    
    }
    
    ## Map genes using membership 
    
    setwd(directory) 
    
    tmp <- list.files(pattern = "genesetMap.rds", recursive = T)
    
    if (!isEmpty(tmp)) {
      pjt <- read_rds(tmp)
    } else {
      
      cat(prep('Mapping genes to geneSets...\n' %+% 
                 wild('This might take a while...\n') %+% 
                 fyi('Suggestion: go for a long walk outside..\n')))
      
      pjt <- list()
      for (i in 1:length(pj)) {
        for (j in 2:length(pj[[i]])) {
          pj[[i]][idx[[i]][j][[1]],j] <- pj[[i]][idx[[i]][j][[1]],1]
          pjt[[i]] <- t(pj[[i]]) %>% as.data.frame()   
          pjt[[i]] <- pjt[[i]][-c(1),]
          pjt[[i]][pjt[[i]] == "."] <- NA
          rownames(pjt[[i]]) <- stri_trans_totitle(rownames(pjt[[i]]))
          #pjt[[i]] <- stri_trans_totitle(rownames(pjt[[i]]))
          }}
  
    }
    
    # pjt <- lapply(pjt, function(x) {
    #   x <- lapply(x, function(y) stri_trans_totitle(y)) %>% as.data.frame()
    #   x
    #   
    # })
    # 
    names(pjt) <- names(memListT)
    pjt <- pjt[lengths(pjt) > minGenes] # cutoff of min genes
    write_rds(pjt,"genesetMap.rds")
    
    #pjt <- readRDS("genesetMap.rds")
    
    #.GlobalEnv$pjt <- pjt  
    
    cat(yay("Gene mapping complete!\n"))
    cat(prep("Preparing dataframe for GeneSet Clustering...\n" %+% 
               wild("This too, will take a while!\n")))
    
    
    ##################### start mods #####################
    
    modGeneList <- readRDS("modGeneList.rds") #mods
    
    dlist <- modGeneList[names(pjt)] # subset frames used for analysis

    
    dlist <- lapply(dlist, function(x) {
      
      x$symbol <- stri_trans_totitle(x$symbol)
      x
      
    })
    
    symb = "symbol"
    
    dp <- list()
    
    for (i in seq_along(pjt)){
      dp[[i]] <- t(pjt[[i]]) %>% as.data.frame()   
      colnames(dp[[i]]) <- dp[[i]][1,]
      dp[[i]] <- data.frame(symb = rownames(dp[[i]]),dp[[i]])
      colnames(dp[[i]])[1] <- symb
      dp[[i]] <- left_join(dlist[[i]],dp[[i]], by = symb)
    }
    
    
    pp <- list()
    
    tmp <- list.files(pattern = "_pathways.txt",recursive = T)
    
    
    for (i in 1:length(tmp)) {
      pp[[i]] <- read.csv(tmp[[i]], sep = "\t")    
    }
    
    names(pp) <- gsub("_pathways.txt", "", tmp)
    names(pp) <- gsub(pattern = ".*\\/",replacement = "",names(pp))
    pp <- pp[names(pjt)]
    stopifnot(names(pp) == names(pjt))
    
    
    tdp <- list()
    ups <- list()
    downs <- list()
    
    for (i in seq_along(dp)){
      tdp[[i]] <- t(dp[[i]]) %>% as.data.frame()
      tdp[[i]] <- tdp[[i]][-1,]
      #ups[[i]]  <- apply(tdp[[i]],1, function(x) paste(x[x %in% up[[i]]],sep = ","))
      #ups[[i]][length(ups[[i]]) == 0] <- "banana"
      #tdp[[i]]$Up_regulated <- ups[[i]]
      #tdp[[i]]$Up_regulated <- gsub('[c.\\)\\(\\"]' ,"",tdp[[i]]$Up_regulated)
      #downs[[i]]  <- apply(tdp[[i]],1, function(x) paste(x[x %in% down[[i]]],sep = ","))
      #tdp[[i]]$Down_regulated <- downs[[i]]
      #tdp[[i]]$Down_regulated <- gsub('[c.\\)\\(\\"]' ,"",tdp[[i]]$Down_regulated)
      tdp[[i]] <- data.frame(rownames(tdp[[i]]),tdp[[i]], stringsAsFactors = FALSE)
      colnames(tdp[[i]])[1] <- "name"
      
    }
    
    ppp <- list()

    for (i in seq_along(pp)){

      ppp[[i]] <- inner_join(pp[[i]], tdp[[i]], by = "name")

    }

    names(ppp) <- names(pp)
  
    clust_list <- list()
    
    for (i in seq_along(ppp)) {
      clust_list[[i]] <- data.frame(ID = ppp[[i]]$name,
                                    Term_Description = ppp[[i]]$description,
                                    Up_regulated = ppp[[i]]$Up_regulated,
                                    Down_regulated = ppp[[i]]$Down_regulated,
                                    highest_p = ppp[[i]]$correctedPValue,
                                    lowest_p = ppp[[i]]$adjustedPValue
      )
    }

    ## read nodelists first
    
    
    facs <- str_extract(pattern = "[[:alnum:]\\.]+$", string = names(nodelist4net)) # extract cond names. Use . sep!
    facs2 <- paste0("_",facs)
    dots <- paste0(".",facs)
    names(nodelist4net) <- mgsub(pattern = facs2,replacement = dots,names(nodelist4net))
    names(nodelist4net) <-  gsub(pattern = "\\.(?=C)",replacement = "_", names(nodelist4net),perl = T) # if . sepatatin 'CTD'
    names(nodelist4net)[1]
    
    
    
    nSplit <- split(nodelist4net,factor(facs)) 
    eSplit <-   split(edgelist4net,factor(facs))
    
    ## Drop incomplete results
    
    gsnames <- list()
    
    for (i in seq_along(nSplit)){
      gsnames[[i]] <- data.frame(ENTREZID = unique(str_extract(pattern = ".*[_]", string = names(nSplit[[i]]))))
      gsnames[[i]]$ENTREZID <-   gsub(pattern = "_$", replacement = "", x = gsnames[[i]]$ENTREZID) # fix
    }
    
    
    namebind <- rbindlist(gsnames)
    
    id <- table(namebind$ENTREZID) < length(facs2) # find incomplete results
    id2 <- table(namebind$ENTREZID) > length(facs2)
    
    uni <- names(table(namebind$ENTREZID)[id])
    uni2 <- names(table(namebind$ENTREZID)[id2])
    
    nSplit <- split(nodelist4net,factor(facs)) 
    eSplit <-   split(edgelist4net,factor(facs))
    
    if (length(uni) > 0 || length(uni2) > 0) {
      
      newnSplit <- list()
      neweSplit <- list()
      exc <- list()
      for (i in seq_along(nSplit)){
        exc <-  data.frame(semi_exclusive = names(nSplit[[i]][grepl(paste(uni,collapse="|"),x = names(nSplit[[i]]))]))
        newnSplit[[i]] <- nSplit[[i]][!grepl(paste(uni,collapse="|"),x = names(nSplit[[i]]))]
        neweSplit[[i]] <- eSplit[[i]][!grepl(paste(uni,collapse="|"),x = names(eSplit[[i]]))]
        #newnSplit[[i]] <- newnSplit[[i]][!grepl(paste(uni2,collapse="|"),x = names(newnSplit[[i]]))]
        #neweSplit[[i]] <- neweSplit[[i]][!grepl(paste(uni2,collapse="|"),x = names(neweSplit[[i]]))]
        
        
      }
      
      names(newnSplit) <- names(nSplit)
      names(neweSplit) <- names(eSplit)
      nSplit <- newnSplit 
      eSplit <- neweSplit 
      write.csv(exc, file = "exclusive_pathways.csv", row.names = FALSE)
      huxdf <- exc %>% as_hux()
      format <- huxdf %>% theme_article()
      quick_html(format, "exclusive_pathways.html")
      rm(neweSplit)
      rm(newnSplit)
    } 
    
    eSplit <- lapply(eSplit, function(x) x[sapply(x, function(y) nrow(y) > 0)])
    nSplit <- lapply(nSplit, function(x) x[sapply(x, function(y) nrow(y) > 0)])
    
    # Remove condition label from pathways
    
    # gsnames  <- data.frame(ENTREZID = unique(str_extract(pattern = ".*[_]", string = names(nSplit[[1]])))) # get gsnames
    # gsnames$ENTREZID <- gsub(pattern = "_$", replacement = "", x = gsnames$ENTREZID) # fix
    # 
    # 
    # # Rename
    # 
    # nSplit <- lapply(nSplit, function(x){
    #   names(x) <- gsnames$ENTREZID 
    #   x
    # })  
    # 
    
    
    ## Remove unique columns
    
    #setwd(directory)
    #tmp <- list.files(pattern = "geneListDes2.rds", recursive = T)
    #geneListDes2 <- read_rds(tmp)
    
    
    comps <- names(modGeneList)
    nSplit <- lapply(nSplit, function(x) lapply(x, function(y){
      z <- !grepl(paste(comps,collapse = "|"),colnames(y),perl = T)
      g <- y[,z] %>% as.data.frame()
      g
    }))
    
    
    # eSplit2 <- lapply(eSplit, function(x){
    #   names(x) <- gsnames$ENTREZID 
    #   x # return x
    #   
    # })  
    # 
    
    # Bind into data frame
    
    nodeBound <- lapply(nSplit, function(x) rbindlist(x))
    nodeBound <-  lapply(nodeBound, function(x) as.data.frame(x))
    finalNodes <- lapply(nodeBound, function(x) plyr::rbind.fill(x,gsnames)) # append with GS nodes
    edgeBound <- lapply(eSplit, function(x) rbindlist(x)) # bind edge
    
    for (i in seq_along(edgeBound)) {
      
      colnames(edgeBound[[i]])[3] <- "Symbol"
      
    }
    # Make Edges data frame
    grid <- expand.grid(1:length(nodeBound),gsnames$ENTREZID)
    names <- as.character(grid$Var2) %>% as.data.frame()
    grid$Var1 <- grid$Var1 %>% as.character() %>% as.numeric()
    grid$Var2 <- grid$Var2 %>% as.numeric()
    
    GSnodeIX <- list() 
    
    for (i in 1:nrow(grid)) {
      
      GSnodeIX[[i]] <- data.frame(geneID = rep(names[[i,1]],nrow(nSplit[[grid[[1,1]]]][[grid[[1,2]]]])),
                                  TO = nSplit[[grid[[1,1]]]][[grid[[1,2]]]][1],
                                  check.names = F,
                                  stringsAsFactors = F)
      colnames(GSnodeIX[[i]])[2] <- "TO"
    }
    
    names(GSnodeIX) <- names$.
    
    
    conds <- unique(str_extract(pattern = "[[:alnum:]\\.]+$", string = names(nodeBound)))
    GSnodeSplit <- split(GSnodeIX, factor(conds))
    GSnodeBound <- lapply(GSnodeSplit, function(x) rbindlist(x)) # bind lists
    
    
    for (i in seq_along(GSnodeSplit)) {
      GSedgesList <- lapply(edgeBound, function(x) plyr::rbind.fill(x,GSnodeBound[[i]])) # add GS nodes
    }
    
    # 
    # # Gather genesets, pvals and member/rep status
    # 
    # statdf <- list()
    # for (i in names(pp)) {
    #   statdf[[i]] <- pp[[i]][,-c(3,4)] %>% as.data.frame() # with number label this time
    # }
    # 
    # # Split into status
    # statSplit <- lapply(statdf, function(x) split(x,x$Status)) # Split into rep/mem per condition
    # 
    # # Make GS Edges
    # statJoin <- list()
    # for (i in seq_along(statSplit)) {
    #   
    #   statJoin[[i]] <- left_join(statSplit[[i]][[1]],statSplit[[i]][[2]], by = "Cluster") %>% 
    #     dplyr::select(c(Term_Description.x,Term_Description.y)) # select from and to columns
    #   colnames(statJoin[[i]]) <- c("ENTREZID","TO") # rename to match edge df
    #   
    # }
    
    
    #names(statJoin) <- names(statdf) 
    #names(GSedgesList) <- names(statdf)
    #names(finalNodes) <- names(statdf)
    #GSedgesList <- GSedgesList[names(statdf)]
    finalEdges <- GSedgesList
    #finalEdges <- map2(GSedgesList,statJoin, plyr::rbind.fill) # append GS edges to gene edges
    
    
    # Somtimes a condition/keyword appears in different comparisons. Remove excess columns.
    
    namesvec <- paste0("\\<",names(finalNodes),"\\>") # word boundaries
    
    finalNodes <- lapply(finalNodes, function(x){
      x <- x[,!grepl(pattern = paste(namesvec,collapse = "|"), x = paste0(colnames(x)),fixed = FALSE)] 
      x
    })
    
    
    finalNodes <- finalNodes[names(finalEdges)]
    
    # Ensure that all nodes are in edge list  
    
    if (is.null(species) | species ==  "mouse") {
      org <- "org.Mm.eg.db"
    } else if (species == 'human') {
      org <- "org.Hs.eg.db"
    }
    
    edgeNodes <- list()
    tmp <- list()
    
    for (i in seq_along(finalEdges)){
      
      tmp[[i]] <- rbind(data.frame(geneID=unique(finalEdges[[i]]$geneID)), 
                        data.frame(geneID=unique(finalEdges[[i]]$TO)))
      
      edgeNodes[[i]] <- left_join(tmp[[i]],finalNodes[[i]], by = "geneID")
      
      edgeNodes[[i]] <- subset(edgeNodes[[i]],!duplicated(edgeNodes[[i]]$geneID))
      
      colnames(edgeNodes[[i]]) <- gsub(pattern = "\\.y", "",colnames(edgeNodes[[i]])) # delete join artifact
      
      edgeNodes[[i]]$sp <- rep(species,nrow(edgeNodes[[i]])) # fill NAs in species columns
      
      edgeNodes[[i]][is.na(edgeNodes[[i]]$Symbol),2] <- edgeNodes[[i]][is.na(edgeNodes[[i]]$Symbol),1]
      
      #edgeNodes[[i]]$Symbol <- bitr(edgeNodes[[i]]$geneID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org)$SYMBOL
    }
    
    finalNodes <- edgeNodes
    names(finalNodes) <- names(finalEdges)
    #names(finalEdges) <- names(statJoin)
    #.GlobalEnv$finalEdges <- finalEdges
    #.GlobalEnv$finalNodes <- finalNodes
    
    cat(prep("Making Graphs...\n"))
    
    repNets <- list()    
    for (i in seq_along(finalEdges)) {
      repNets[[i]] <- graph_from_data_frame(d = finalEdges[[i]],
                                            vertices = finalNodes[[i]],
                                            directed = F) %>% 
        igraph::simplify(remove.loops = T, remove.multiple = T) %>% 
        delete.vertices(which(igraph::degree(.) == 0))
    }
    
    # If the network is sparse slash doesnt exist
    
    if (any(lengths(repNets) == 0)) {    
      
      repNets <- list()    
      for (i in seq_along(finalEdges)) {
        repNets[[i]] <- graph_from_data_frame(d = finalEdges[[i]],
                                              vertices = finalNodes[[i]],
                                              directed = F) 
      }
    }
    names(repNets) <- names(finalEdges)
    
    setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
    
    pdf("wgcna_mods_ppi_integrated_networks.pdf", width = 11.7, height = 8.3)
    par(mfrow = c(1, 2),cex.main = 0.75)
    
    for (i in seq_along(repNets)) {
      plot(repNets[[i]],
           #vertex.color = np[node_cols1[[i]]],
           #vertex.frame.color = np[node_cols1[[i]]],
           vertex.label = paste0(V(repNets[[i]])$geneID),
           vertex.label.color = "black",
           vertex.label.font = 2,
           vertex.size = 20,
           edge.width  = 2,
           edge.color = "gray50",
           edge.curved = 0.1,
           frame = FALSE,
           main = names(repNets[i]),
           layout = layout_nicely) 
    }
    
    dev.off()
    
    #### testing mods cont..
    
    cat(prep("Plotting Adjacency Matrix Heatmap...\n"))
    netm <- list()
    for (i in seq_along(repNets)) {
      pdf(paste0(names(repNets[i]),"_repnets_adj_heatmap.pdf"))
      netm[[i]] <- get.adjacency(repNets[[i]], sparse=F) %>% 
        as.data.frame()
      netm[[i]] <- netm[[i]][!rowSums(netm[[i]]) == 0]
      colnames(netm[[i]]) <- V(repNets[[i]])$Symbol
      rownames(netm[[i]]) <- V(repNets[[i]])$Symbol
      palf <- diverge_hcl(50,palette = "Berlin")
      heatmap(as.matrix(netm[[i]]), Rowv = NA, Colv = NA, col = palf,
              scale="none", margins=c(10,10)) 
      graphics.off()
    }
    cat(yay("Adjacency heatmaps exported!\n"))
    ## Write tables
    cat(prep("Writing eigen-centrality and betweeness tables...\n"))
    ec <- list()
    be <- list()
    ecbe <- list()
    orderlist <- list()
    eiglist <- list()
    betlist <- list()
    huxlistE <- list()
    huxlistB <- list()
    formatsE <- list()
    formatsB <- list()
    
    for (i in seq_along(repNets)) {
      ec[[i]] <- data.frame(eigen_centrality(repNets[[i]], directed = FALSE, scale = TRUE,
                                             weights = NULL, options = arpack_defaults))
      be[[i]] <- data.frame(ec[[i]],betweenness(repNets[[i]], v = V(repNets[[i]]), directed = TRUE, weights = NULL,                                                nobigint = TRUE, normalized = FALSE))
      ecbe[[i]] <- be[[i]][,c(1,23)]
      colnames(ecbe[[i]])[c(1,2)] <- c("Eigen-Centrality","Betweeness")
      ecbe[[i]]$Gene <- V(repNets[[i]])$Symbol
      ecbe[[i]] <- ecbe[[i]][,c(3,1,2)]
      eiglist[[i]] <- ecbe[[i]][order(ecbe[[i]]$"Eigen-Centrality",ecbe[[i]]$Betweeness, decreasing = TRUE),]
      betlist[[i]] <- ecbe[[i]][order(ecbe[[i]]$Betweeness,ecbe[[i]]$"Eigen-Centrality", decreasing = TRUE),]
      huxlistE[[i]] <- eiglist[[i]] %>% as_hux()
      huxlistB[[i]] <- betlist[[i]] %>% as_hux()
      formatsE[[i]] <-  huxlistE[[i]] %>% theme_article()
      formatsB[[i]] <- huxlistB[[i]] %>% theme_article()
      quick_html(formatsE[[i]],file = paste0(names(repNets[i]),"_repnets_eigen_centrality_sorted.html"))
      write.csv(formatsE[[i]],file = paste0(names(repNets[i]),"_repnets_eigen_centrality_sorted.csv"))
      quick_html(formatsB[[i]],file = paste0(names(repNets[i]),"_repnets_betweeness_sorted.html"))
      write.csv(formatsB[[i]],file = paste0(names(repNets[i]),"_repnets_betweeness_sorted.csv"))
    }
    
    # Cytoscape export
    
    if (cytoscape) {
      cat(prep("Exporting graphs to cytoscape...\n"))
      # Pause function
      readkey <- function()
      {
        cat("Open cytoscape, then on your keyboard and in the console, Press [Enter] to continue")
        line <- readline()
      }
      readkey()
      cytoscapePing()
      for (i in seq_along(repNets)) {
        createNetworkFromIgraph(repNets[[i]],network.name = names(repNets[1]))
      }
      write.csv(names(finalEdges), "order_of_integrated_mods_conditions_cytoscape.csv")
      cat(yay("Graphs exported!\n"))
      cat(cyan("Analysis complete!\n"))
        }  
      }
    }  

   ##################### end mods ########################## 
  
    ## If for some reason the hgnc symcols are all caps...
    
    dlist <- lapply(dlist, function(x) {
     
      x$hgnc_symbol <- stri_trans_totitle(x$hgnc_symbol)
      x
    })
    
    
    dp <- list()
    
    for (i in seq_along(pjt)){
      dp[[i]] <- t(pjt[[i]]) %>% as.data.frame()   
      colnames(dp[[i]]) <- dp[[i]][1,]
      dp[[i]] <- data.frame(symb = rownames(dp[[i]]),dp[[i]])
      colnames(dp[[i]])[1] <- symb
      dp[[i]] <- left_join(dlist[[i]],dp[[i]], by = symb)
    }
    
    # Map fold changes
    
    up <- list()      
    down <- list()
    for (i in seq_along(dp)) {
      up[[i]] <-  dp[[i]][dp[[i]]$log2FC > 0,][,symb]
      #up[[i]] <-  up[[i]][complete.cases(up[[i]])]
      #up[[i]][is.na(up[[i]][i])] <-  c('')
      #up[[i]] <- paste0(up[[i]][!is.na(up[[i]])], collapse = ",")
      down[[i]] <-  dp[[i]][dp[[i]]$log2FC <= 0,][,symb]
      #down[[i]] <-  down[[i]][complete.cases(down[[i]])]
      #down[[i]] <- paste0(down[[i]][!is.na(down[[i]])], collapse = ",")
    
    }
    
    ## Wot...
    
    up <- lapply(up, function(x) {
      x[length(x) == 0] <- 'Laptm'
      x
    })
    
    down <- lapply(down, function(x) {
      x[length(x) == 0] <- 'Laptm'
      x
    })
    
    tmp <- list.files(pattern = "geneinput", recursive = T)[1]
    geneinput <- read_rds(tmp)
    geneinput <- geneinput[names(pjt)]
    names(dp) <- names(pjt)
  
    tdp <- list()
    ups <- list()
    downs <- list()
    
    for (i in seq_along(dp)){
      tdp[[i]] <- t(dp[[i]]) %>% as.data.frame()
      tdp[[i]] <- tdp[[i]][-1,]
      ups[[i]]  <- apply(tdp[[i]],1, function(x) paste(x[x %in% up[[i]]],sep = ","))
      ups[[i]][length(ups[[i]]) == 0] <- "banana"
      tdp[[i]]$Up_regulated <- ups[[i]]
      tdp[[i]]$Up_regulated <- gsub('[c.\\)\\(\\"]' ,"",tdp[[i]]$Up_regulated)
      downs[[i]]  <- apply(tdp[[i]],1, function(x) paste(x[x %in% down[[i]]],sep = ","))
      tdp[[i]]$Down_regulated <- downs[[i]]
      tdp[[i]]$Down_regulated <- gsub('[c.\\)\\(\\"]' ,"",tdp[[i]]$Down_regulated)
      tdp[[i]] <- data.frame(rownames(tdp[[i]]),tdp[[i]], stringsAsFactors = FALSE)
      colnames(tdp[[i]])[1] <- "name"
      
    }
    
    names(tdp) <- names(dlist)

    
    # Join paths  
    
    pp <- list()
    tmp <- list.files(pattern = "_pathways.txt",recursive = T)
    for (i in 1:length(tmp)) {
      pp[[i]] <- read.csv(tmp[[i]], sep = "\t")    
    }
    names(pp) <- gsub("_pathways.txt", "", tmp)
    names(pp) <- gsub(pattern = ".*\\/",replacement = "",names(pp))
    pp <- pp[names(pjt)]
    stopifnot(names(pp) == names(pjt))
    
    
    ppp <- list()
    
    for (i in seq_along(pp)){
      
      ppp[[i]] <- inner_join(pp[[i]], tdp[[i]], by = "name")
      
    }
    
    names(ppp) <- names(pp)
  
    
    ## Make df for clustering
    
    clust_list <- list()  
    
    for (i in seq_along(ppp)) {
      clust_list[[i]] <- data.frame(ID = ppp[[i]]$name,
                                    Term_Description = ppp[[i]]$description,
                                    Up_regulated = ppp[[i]]$Up_regulated,
                                    Down_regulated = ppp[[i]]$Down_regulated,
                                    highest_p = ppp[[i]]$correctedPValue,
                                    lowest_p = ppp[[i]]$adjustedPValue
                                    )
                                    #Fold_Enrichment = -log10(ppp[[i]]$adjustedPValue)
    
    
    # Loop this to swap hi/lo p
    
    
    idx <- list()
    for (i in seq_along(clust_list)) {
      idx[[i]] <- which(clust_list[[i]]$lowest_p - clust_list[[i]]$highest_p > 0)
      clust_list[[i]][,c(5,6)][idx[[i]],] <- clust_list[[1]][,c(6,5)][idx[[i]],]
      #clust_list[[i]][clust_list[[i]]== "banana"] <- NA 
      }
    }
    
    
    ## Cluster    
    
    setwd(paste0(directory,"/GSEA"))
    if(!dir.exists("Representative_GSEA")){
      dir.create("Representative_GSEA")
    }
    setwd("Representative_GSEA")
    
    names(clust_list) <- gsub("_","\\.",names(ppp))
    
    cat(calc("Now clustering enriched terms...\n"))
    
    gs_clustList <- list()
    
    pdf("gsea_cluster_dendogams.pdf")
    
    for (i in seq_along(clust_list)) {
      gs_clustList[[i]] <- cluster_enriched_terms(clust_list[[i]],
                                                  method = "hierarchical",
                                                  plot_clusters_graph = F,
                                                  use_description = TRUE)
    }
    graphics.off()
    cat(yay("Clustering complete!\n"))
    
    names(gs_clustList) <- names(ppp)
    write_rds(gs_clustList, "gs_clustList.rds")
    .GlobalEnv$gs_clustList <- gs_clustList
    
    # Huxtable
    
    cat(prep("Writing tables...\n"))
    gsl <- list()
    for (i in names(gs_clustList)) {
      gsl[[i]] <- gs_clustList[[i]][,-c(3,4,7)] %>% as_hux()
    }
    formats = list() 
    for (i in names(gsl)) {
      formats[[i]] = gsl[[i]] %>% set_background_color(everywhere,row = gsl[[i]]$Status == "Representative",value = "lightblue") %>% 
        theme_article()
    }
    for (i in names(formats)) {
      quick_html(formats[[i]],file = paste0(expname,"_",i,"_gsea_clusters_tables.html"))
    }
    
    for (i in names(formats)) {
      write.csv(formats[[i]],file = paste0(expname,"_",i,"_gsea_clusters_tables.csv"))
    }
    
    ### Representative only
    
    gsl2 <- list()
    for (i in names(gs_clustList)) {
      gsl2[[i]] <-  subset(gs_clustList[[i]],gs_clustList[[i]]$Status == "Representative" )
      gsl2[[i]] <- gsl2[[i]][,-c(3,4,7)] %>% as_hux()
    }
    formats = list() 
    for (i in names(gsl2)) {
      formats[[i]] = gsl2[[i]] %>% theme_article()
    }
    for (i in names(formats)) {
      quick_html(formats[[i]],file = paste0(expname,"_",i,"_gsea_representative_only_GS.html"))
    }
    .GlobalEnv$rep_gs <- gsl2
    cat(yay("Tables exported!\n"))
    
    } else {
      
      cat(fyi('Gene set clustering result already exists!\n'))
    }
    
    
    
    ## PPi-GSclutser Plots
    
    cat(prep("Preparing to plot universal PPI-GeneSet Cluster network...\n"))
    
    setwd(directory)
    
    tmp <- list.files(pattern = "nodelist4net.rds", recursive = T)
    tmp2 <- list.files(pattern = "edgelist4net.rds", recursive = T)
    tmp3 <- list.files(pattern = "gs_clustList.rds", recursive = T)
    
    if (!isEmpty(tmp)) {
      nodelist4net <- read_rds(tmp)
      edgelist4net <- read_rds(tmp2)
      gs_clustList <- read_rds(tmp3)
      .GlobalEnv$edgelist4net <- edgelist4net
      .GlobalEnv$nodelist4net <- nodelist4net
      
    } else {
      
      stop("Please run ppiGSEA before gseaCluster!\n")
      
    }
  
    # Split into conditions
  
    
    if (mods) {
      
    facs <- str_extract(pattern = "[[:alnum:]\\.]+$", string = names(nodelist4net)) # extract cond names. Use . sep!
    facs2 <- paste0("_",facs)
    dots <- paste0(".",facs)
    names(nodelist4net) <- mgsub(pattern = facs2,replacement = dots,names(nodelist4net))
    names(nodelist4net) <-  gsub(pattern = "\\.(?=C)",replacement = "_", names(nodelist4net),perl = T) # if . sepatatin 'CTD'
    names(nodelist4net)[1]

    } else {
  
    facs <- str_extract(pattern = "[[:alnum:]\\.]+$", string = names(nodelist4net)) # extract cond names. Use . sep!
    names(nodelist4net) <-  gsub(pattern = "\\.(?=C)",replacement = "_", names(nodelist4net),perl = T) # if . sepatatin 'CTD'
    
    }
    
    nSplit <- split(nodelist4net,factor(facs)) 
    eSplit <-   split(edgelist4net,factor(facs))
    
    ## Drop incomplete results
    
    gsnames <- list()
    
    for (i in seq_along(nSplit)){
      gsnames[[i]] <- data.frame(ENTREZID = unique(str_extract(pattern = ".*[_]", string = names(nSplit[[i]]))))
      gsnames[[i]]$ENTREZID <-   gsub(pattern = "_$", replacement = "", x = gsnames[[i]]$ENTREZID) # fix
    }
    
    
    namebind <- rbindlist(gsnames)
  
    id <- table(namebind$ENTREZID) < length(facs2) # find incomplete results
    id2 <- table(namebind$ENTREZID) > length(facs2)
    
    uni <- names(table(namebind$ENTREZID)[id])
    uni2 <- names(table(namebind$ENTREZID)[id2])
    
    nSplit <- split(nodelist4net,factor(facs)) 
    eSplit <-   split(edgelist4net,factor(facs))
    
    if (length(uni) > 0 || length(uni2) > 0) {
      
    newnSplit <- list()
    neweSplit <- list()
    exc <- list()
    for (i in seq_along(nSplit)){
      exc <-  data.frame(semi_exclusive = names(nSplit[[i]][grepl(paste(uni,collapse="|"),x = names(nSplit[[i]]))]))
      newnSplit[[i]] <- nSplit[[i]][!grepl(paste(uni,collapse="|"),x = names(nSplit[[i]]))]
      neweSplit[[i]] <- eSplit[[i]][!grepl(paste(uni,collapse="|"),x = names(eSplit[[i]]))]
      #newnSplit[[i]] <- newnSplit[[i]][!grepl(paste(uni2,collapse="|"),x = names(newnSplit[[i]]))]
      #neweSplit[[i]] <- neweSplit[[i]][!grepl(paste(uni2,collapse="|"),x = names(neweSplit[[i]]))]
      
      
      }
    
    names(newnSplit) <- names(nSplit)
    names(neweSplit) <- names(eSplit)
    nSplit <- newnSplit 
    eSplit <- neweSplit 
    write.csv(exc, file = "semi_exclusive_pathways.csv", row.names = FALSE)
    huxdf <- exc %>% as_hux()
    format <- huxdf %>% theme_article()
    quick_html(format, "semi_exclusive_pathways.html")
    rm(neweSplit)
    rm(newnSplit)
  } 

    
    eSplit <- lapply(eSplit, function(x) x[sapply(x, function(y) nrow(y) > 0)])
    nSplit <- lapply(nSplit, function(x) x[sapply(x, function(y) nrow(y) > 0)])
    
    # Remove condition label from pathways
    
    gsnames  <- data.frame(ENTREZID = unique(str_extract(pattern = ".*[_]", string = names(nSplit[[1]])))) # get gsnames
    gsnames$ENTREZID <- gsub(pattern = "_$", replacement = "", x = gsnames$ENTREZID) # fix
  
    
    # Rename
    
    nSplit <- lapply(nSplit, function(x){
      names(x) <- gsnames$ENTREZID 
      x
    })  
    
    ## Remove unique columns
  
    setwd(directory)
    tmp <- list.files(pattern = "geneListDes2.rds", recursive = T)
    geneListDes2 <- read_rds(tmp)
    comps <- names(geneListDes2)
    
    nSplit <- lapply(nSplit, function(x) lapply(x, function(y){
      z <- !grepl(paste(comps,collapse = "|"),colnames(y),perl = T)
      g <- y[,z] %>% as.data.frame()
      g
      }))
  
    
    eSplit <- lapply(eSplit, function(x){
      names(x) <- gsnames$ENTREZID 
      x # return x
      
    })  
  
  
    # Bind into data frame
  
    nodeBound <- lapply(nSplit, function(x) rbindlist(x))
    nodeBound <-  lapply(nodeBound, function(x) as.data.frame(x))
    finalNodes <- lapply(nodeBound, function(x) plyr::rbind.fill(x,gsnames)) # append with GS nodes
    edgeBound <- lapply(eSplit, function(x) rbindlist(x)) # bind edge
    
    # Make Edges data frame
    grid <- expand.grid(1:length(nodeBound),gsnames$ENTREZID)
    names <- as.character(grid$Var2) %>% as.data.frame()
    grid$Var1 <- grid$Var1 %>% as.character() %>% as.numeric()
    grid$Var2 <- grid$Var2 %>% as.numeric()
    
    GSnodeIX <- list() 
    
    for (i in 1:nrow(grid)) {
      
      GSnodeIX[[i]] <- data.frame(ENTREZID = rep(names[[i,1]],nrow(nSplit[[grid[[1,1]]]][[grid[[1,2]]]])),
                                  TO = nSplit[[grid[[1,1]]]][[grid[[1,2]]]][1],
                                  check.names = F,
                                  stringsAsFactors = F)
      colnames(GSnodeIX[[i]])[2] <- "TO"
    }
  
    names(GSnodeIX) <- names$.
    
    conds <- unique(str_extract(pattern = "[[:alnum:]\\.]+$", string = names(nodeBound)))
    GSnodeSplit <- split(GSnodeIX, factor(conds))
    GSnodeBound <- lapply(GSnodeSplit, function(x) rbindlist(x)) # bind lists
    
    for (i in seq_along(GSnodeSplit)) {
      GSedgesList <- lapply(edgeBound, function(x) plyr::rbind.fill(x,GSnodeBound[[i]])) # add GS nodes
    }
    
    # Gather genesets, pvals and member/rep status
    statdf <- list()
    for (i in names(gs_clustList)) {
      statdf[[i]] <- gs_clustList[[i]][,-c(3,4)] %>% as.data.frame() # with number label this time
    }
    
    # Split into status
    statSplit <- lapply(statdf, function(x) split(x,x$Status)) # Split into rep/mem per condition
    
    # Make GS Edges
    statJoin <- list()
    for (i in seq_along(statSplit)) {
      
      statJoin[[i]] <- left_join(statSplit[[i]][[1]],statSplit[[i]][[2]], by = "Cluster") %>% 
        dplyr::select(c(Term_Description.x,Term_Description.y)) # select from and to columns
      colnames(statJoin[[i]]) <- c("ENTREZID","TO") # rename to match edge df
      
    }
  
    names(statJoin) <- names(statdf) 
    names(GSedgesList) <- names(statdf)
    names(finalNodes) <- names(statdf)
    GSedgesList <- GSedgesList[names(statdf)]
    finalEdges <- map2(GSedgesList,statJoin, plyr::rbind.fill) # append GS edges to gene edges
    
    
    # Somtimes a condition/keyword appears in different comparisons. Remove excess columns.
    
    namesvec <- paste0("\\<",names(finalNodes),"\\>") # word boundaries
    
    finalNodes <- lapply(finalNodes, function(x){
      x <- x[,!grepl(pattern = paste(namesvec,collapse = "|"), x = paste0(colnames(x)),fixed = FALSE)] 
      x
    })
    
    finalNodes <- finalNodes[names(finalEdges)]
    
    # Ensure that all nodes are in edge list  
    
    if (is.null(species) | species ==  "mouse") {
      org <- "org.Mm.eg.db"
    } else if (species == 'human') {
      org <- "org.Hs.eg.db"
    }
    edgeNodes <- list()
    tmp <- list()
    
    for (i in seq_along(finalEdges)){
      tmp[[i]] <- rbind(data.frame(ENTREZID=unique(finalEdges[[i]]$ENTREZID)), 
                        data.frame(ENTREZID=unique(finalEdges[[i]]$TO)))
      
      edgeNodes[[i]] <- left_join(tmp[[i]],finalNodes[[i]], by = "ENTREZID")
      
      edgeNodes[[i]] <- subset(edgeNodes[[i]],!duplicated(edgeNodes[[i]]$ENTREZID))
      
      colnames(edgeNodes[[i]]) <- gsub(pattern = "\\.y", "",colnames(edgeNodes[[i]])) # delete join artifact
      
      edgeNodes[[i]]$sp <- rep(species,nrow(edgeNodes[[i]])) # fill NAs in species columns
      
      edgeNodes[[i]][is.na(edgeNodes[[i]]$Symbol),2] <- edgeNodes[[i]][is.na(edgeNodes[[i]]$Symbol),1]
      
      #edgeNodes[[i]]$Symbol <- bitr(edgeNodes[[i]]$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org)$SYMBOL
    }
    
  
    finalNodes <- edgeNodes
    names(finalNodes) <- names(statJoin)
    names(finalEdges) <- names(statJoin)
    .GlobalEnv$finalEdges <- finalEdges
    .GlobalEnv$finalNodes <- finalNodes
    
    cat(prep("Making Graphs...\n"))
    
    
    repNets <- list()    
    for (i in seq_along(finalEdges)) {
      repNets[[i]] <- graph_from_data_frame(d = finalEdges[[i]],
                                            vertices = finalNodes[[i]],
                                            directed = F) %>% 
        igraph::simplify(remove.loops = T, remove.multiple = T) %>% 
        delete.vertices(which(igraph::degree(.) == 0))
    }
  
    # If the network is sparse slash doesnt exist
    
    if (any(lengths(repNets) == 0)) {    
    
      repNets <- list()    
    for (i in seq_along(finalEdges)) {
      repNets[[i]] <- graph_from_data_frame(d = finalEdges[[i]],
                                            vertices = finalNodes[[i]],
                                            directed = F) 
      }
    }
    names(repNets) <- names(finalEdges)
    
    
    .GlobalEnv$repNets <- repNets
    
    # Map log2FC to color
    n <- diverge_hcl(n = 500, "Blue-Red 3") # colour ramp
    ll <- length(finalEdges)
    lp <- ceiling(length(finalEdges)/2)
    
    fclist <- list()
    
    for (i in seq_along(repNets)){
      fclist[[i]] <- any(abs(as.numeric((E(repNets[[i]]))$log2FC)) <= 1 )
      fc <- unlist(fclist)
      #rm(fclist)
    } 
    
    if (any(fc) == TRUE) {
      
      ceil <- list()
      
      for (i in 1:length(repNets)) {
        ceil[[i]] <- abs(as.numeric((E(repNets[[i]])$log2FC))) # ceil
        ceil[[i]][ceil[[i]] == -Inf] <- 0
        }
    } else {
      
      ceil <- list()
      
      for (i in 1:length(repNets)) {
        ceil[[i]] <- max(abs(as.numeric((E(repNets[[i]])$log2FC)))) # ceil
        ceil[[i]][ceil[[i]] == -Inf] <- 0
      }
    }
      
    # If for some reason log2FC is not in edges 
    if (any(unlist(ceil) == 0 | any(is.na(ceil)))) {
    
      fclist <- list()
      
      for (i in seq_along(repNets)){
        fclist[[i]] <- any(abs(as.numeric((V(repNets[[i]]))$log2FC)) <= 1 )
        fc <- unlist(fclist)
        #rm(fclist)
      } 
      
    if (any(fc) == TRUE) {
        
      ceil <- list()
      for (i in 1:length(repNets)) {
        ceil[[i]] <- abs(as.numeric((V(repNets[[i]])$log2FC))) # ceil
        ceil[[i]][ceil[[i]] == -Inf] <- 0
        }
      } else {
    
      ceil <- list()
      for (i in 1:length(repNets)) {
        ceil[[i]] <- max(abs(as.numeric((V(repNets[[i]])$log2FC)))) # ceil
        ceil[[i]][ceil[[i]] == -Inf] <- 0
          }
        }
      } 
    
    names(ceil) <- names(finalEdges)
    
    sf <-  sapply(ceil, function(x) {
      x <- x[complete.cases(x)] 
      x <- max(x)}
      )
    sf <- round(max(abs(sf)),digits = 2)
    
    
    nodeColors <- list()
    for (i in 1:length(repNets)) {
      nodeColors[[i]] <- (as.numeric(E(repNets[[i]])$log2FC) + sf) / (2*sf) * 500 # mid poin
    }
    # Same as above
    
    if (any(is.na(unlist(nodeColors)))) {
      nodeColors <- list()
      for (i in 1:length(repNets)) {
        nodeColors[[i]] <- (as.numeric(V(repNets[[i]])$log2FC) + sf) / (2*sf) * 500 # mid poin
      }
    }
    
    
    # Vertex label size
    for (i in seq_along(repNets)) {
      V(repNets[[i]])$label.cex <- 0.5
      l <- layout_in_circle(repNets[[i]],order = igraph::degree(repNets[[i]]))
      #l <- layout.auto(repNets[[i]])
    }
    
    # Plot em
    
    cat(prep("Plotting Representative GSEA Networks...\n"))
    
    setwd(paste0(directory,"/GSEA"))
    if(!dir.exists("Representative_GSEA")){
      dir.create("Representative_GSEA")
      setwd("Representative_GSEA")
    }
    
    #par(mfrow = c(1, 2),cex.main = 0.75)
    
    
    par(cex.main = 0.75)
  
    for (i in 1:length(repNets)) {
    pdf(paste0(names(repNets[i]),"_representative_gsea_network_plots.pdf"), width = 11.7, height = 8.3)
    plot(repNets[[i]],
         vertex.color = n[nodeColors[[i]]],
         vertex.frame.color = n[nodeColors[[i]]],
         vertex.label = V(repNets[[i]])$Symbol,
         vertex.label.cex = 0.5,
         vertex.frame.width = 3,
         vertex.label.color = "black",
         vertex.label.font = 1,
         vertex.size = 8,
         vertex.frame.color = V(repNets[[i]])$sp,
         edge.width  = 2,
         edge.color = "grey",
         edge.curved = 0.2,
         layout = l,
         frame = FALSE,
         main = names(repNets[i]))
    
    # legend("bottomleft", 
    #        legend = levels(as.factor(V(repNets[[i]])$sp)), 
    #        col = pals[[i]] , 
    #        bty = "n", 
    #        pch = 20, 
    #        pt.cex = 3, 
    #        cex = 1.5, 
    #        text.col = pals[[i]], 
    #        horiz = FALSE, 
    #        inset = c(0.1, 0.1))
    
    graphics.off()
  }
  cat(yay("Plots Exported!\n"))
  
  # Heatmap of adjacency matrix
  
  cat(prep("Plotting Adjacency Matrix Heatmap...\n"))
  netm <- list()
  for (i in seq_along(repNets)) {
    pdf(paste0(names(repNets[i]),"_repnets_adj_heatmap.pdf"))
    netm[[i]] <- get.adjacency(repNets[[i]], sparse=F) %>% 
      as.data.frame()
    netm[[i]] <- netm[[i]][!rowSums(netm[[i]]) == 0]
    colnames(netm[[i]]) <- V(repNets[[i]])$Symbol
    rownames(netm[[i]]) <- V(repNets[[i]])$Symbol
    palf <- diverge_hcl(50,palette = "Berlin")
    heatmap(as.matrix(netm[[i]]), Rowv = NA, Colv = NA, col = palf,
            scale="none", margins=c(10,10)) 
    graphics.off()
  }
  cat(yay("Adjacency heatmaps exported!\n"))
  ## Write tables
  cat(prep("Writing eigen-centrality and betweeness tables...\n"))
  ec <- list()
  be <- list()
  ecbe <- list()
  orderlist <- list()
  eiglist <- list()
  betlist <- list()
  huxlistE <- list()
  huxlistB <- list()
  formatsE <- list()
  formatsB <- list()
  
  for (i in seq_along(repNets)) {
    ec[[i]] <- data.frame(eigen_centrality(repNets[[i]], directed = FALSE, scale = TRUE,
                                           weights = NULL, options = arpack_defaults))
    be[[i]] <- data.frame(ec[[i]],betweenness(repNets[[i]], v = V(repNets[[i]]), directed = TRUE, weights = NULL,                                                nobigint = TRUE, normalized = FALSE))
    ecbe[[i]] <- be[[i]][,c(1,23)]
    colnames(ecbe[[i]])[c(1,2)] <- c("Eigen-Centrality","Betweeness")
    ecbe[[i]]$Gene <- V(repNets[[i]])$Symbol
    ecbe[[i]] <- ecbe[[i]][,c(3,1,2)]
    eiglist[[i]] <- ecbe[[i]][order(ecbe[[i]]$"Eigen-Centrality",ecbe[[i]]$Betweeness, decreasing = TRUE),]
    betlist[[i]] <- ecbe[[i]][order(ecbe[[i]]$Betweeness,ecbe[[i]]$"Eigen-Centrality", decreasing = TRUE),]
    huxlistE[[i]] <- eiglist[[i]] %>% as_hux()
    huxlistB[[i]] <- betlist[[i]] %>% as_hux()
    formatsE[[i]] <-  huxlistE[[i]] %>% theme_article()
    formatsB[[i]] <- huxlistB[[i]] %>% theme_article()
    #quick_html(formatsE[[i]],file = paste0(names(repNets[i]),"_repnets_eigen_centrality_sorted.html"))
    #write.csv(formatsE[[i]],file = paste0(names(repNets[i]),"_repnets_eigen_centrality_sorted.csv"))
    quick_html(formatsB[[i]],file = paste0(names(repNets[i]),"_repnets_betweeness_sorted.html"))
    write.csv(formatsB[[i]],file = paste0(names(repNets[i]),"_repnets_betweeness_sorted.csv"))
  }
  
  # Cytoscape export
  
  if (cytoscape) {
    cat(prep("Exporting graphs to cytoscape...\n"))
    # Pause function
    readkey <- function()
    {
      cat("Open cytoscape, then on your keyboard and in the console, Press [Enter] to continue")
      line <- readline()
    }
    readkey()
    cytoscapePing()
    for (i in seq_along(repNets)) {
      createNetworkFromIgraph(repNets[[i]])
    }
    write.csv(names(statJoin), "order_of_conditions_cytoscape.csv")
    cat(yay("Graphs exported!\n"))
    cat(cyan("Analysis complete!\n"))
      }  
    }
      


