#' Differential gene expressionn analysis pipeline
#' 
#' DGEA leveraging the wonderful noiSeq package
#' 
#' ceeb
#'
#' @param  counts data frame Raw bulk mRNA seq counts genes as rows, columns as samples. 
#' @param  meta data frame Metadata table with samples as rows and variables as columns. Compulsory columns include 'id' with sample names corresponding to those in counts, 'condition' with relevant tests and controls, 'priority' with integers starting from 0 to 1 for indicating preferred direction of pairwise comparison and 'class' indicating whether the sample is a control or disease/test. 
#' @param  species character Takes 'mouse' or 'human'. 
#' @param  expname character Name of expreriment. E.g., "Exp1_DGEA"
#' @param  cores numeric Specify number of CPU cores.
#' @param  seed numeric Seed from which randomness grows.
#' @param  directory character Name of working folder.
#' @param  CPMfilter logical If TRUE, raw counts will be filtered based of counts per milion (CPM) cutoff
#' @param  CPMcutoff numeric Cutoff for CPM filtering e.g., 2
#' @param  normalise logical, character If TRUE, raw counts will be geTMM normalised. If "CQN", raw counts will be quantile normalised
#' @param  orthologs logical If TRUE, maps genes to orthologs of given species and continues analysis with these orthologous genes.  
#' @param  orthogenes character list List of character vector of orthologous ENRTEZID genes (mouse or human)
#' @param  orthospecies character Takes 'mouse' or 'human'
#' @param  batch logical If TRUE, data will undergo batch correction
#' @param  dropOut logical If TRUE, outliers will be dropped based on PCA distance.
#' @param  svaCor logical If TRUE, performs surrogatew variable analysis and correction.
#' @param  svaFormula character Formula for SVA e.g. id ~ condition + gender
#' @param  svaCovar character Covariates of interest whose influence should not be removed e.g., condition
#' @param  conCovar character Potentially confounding covariates whose influence might want to be removed e.g., gender
#' @param  includeZeros logical If TRUE. zero values will be nudged by 1x10^-16 to facilitate PCA. Otherwise, rows with zero values will be removed.
#' @param  writeOut logical If TRUE, tables generated during analysis will be written to csv.
#' @param  QCplots logical IF TRUE,  various quality control plots will be generated.
#' @param  DE logical If TRUE, performs noisSeq differential gene expression analysis
#' @param  GSEA logical If TRUE, performs setRank network based pathway analysis of differentially expressed genes.
#' @param  msigDB logical nclude msigDB in GSEA?
#' @param  msigCat list of character vectors list of desired msigDB sub-categories e.g., c("H", "C1")
#' @param  CTDbase logical Include Chemical and Toxigenomics Data base in GSEA?
#' @param  chemDis logical if CTDbase, include chemical-gene-disease interactions in GSEA?
#' @param  chemIX logical if CTDbase, include chemical-gene interactions in GSEA?
#' @param  ppiGSEA logical map GSEA results to HIPPIE protein-protein interaction network?
#' @param  clusterGS logical perform hierarchal clustering on GSEA results?
#' @param  cytoscape logical export networks to cytoscape?
#' 
#' @import EDASeq
#' @import DESeq2
#' @import edgeR
#' @import rrcovHD
#' @import clusterProfiler
#' @import limma
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
#' @import CoExpNets
#' @import swamp
#' @import outliers
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
#' oneDiff()
#' 
#' @export
oneDiff <- function(counts = counts,
                    meta = meta,
                    species = NULL,
                    expname =  NULL,
                    normalise = NULL,
                    CPMfilter = NULL,
                    CPMcutoff = NULL,
                    orthologs = FALSE,
                    dropOut = NULL,
                    svaFormula = NULL,
                    svaCor = NULL,
                    svaCovar = NULL,
                    conCovar = NULL,
                    batch = FALSE,
                    includeZeros = TRUE,
                    writeOut = NULL,
                    QCplots = NULL,
                    DE = NULL,
                    GSEA = NULL,
                    orthogenes = NULL,
                    orthospecies = NULL,
                    clusterGS = NULL,
                    CTDbase = NULL,
                    chemDis = NULL,
                    chemIX = NULL,
                    msigCat = NULL,
                    msigDB = NULL,
                    cytoscape = NULL,
                    seed = NULL,
                    directory = NULL,
                    cores = NULL,
                    ppiGSEA = NULL) {

  options(stringsAsFactors = FALSE)

  # crayons
  wild <- red $ bold
  yay <- green $ bold
  calc <- cyan $ bold
  prep <- silver $ bold
  fyi <- white $ bold
  
    
  ## Check if everything is in order

  if (DE == "noiseq") {
    
    if (any(!colnames(counts) %in% meta$id)) { # see cat 
      stop("Not all samples in meta are in counts matrix!")}
    meta <- meta[order(meta$condition),]
    meta  <- mutate_all(meta, function(x) as.character(x)) #need to convert to character
    counts <- counts[,meta$id] #order by sample name
    if (any(!colnames(counts) == meta$id)) { # see cat
      stop("Sample order is not the same between counts matrix and meta!")}
    #rownames(meta) <- paste0(meta$condition,"_",colnames(counts)) # Create new varabole: condition_sampleID
    colnames(counts) <- rownames(meta) # rename samples in count matrix
    .GlobalEnv$meta <- meta # save new meta to global
    .GlobalEnv$counts <- counts
    csums <- colSums(counts)
  }
  
  ##CQN Normalisation
  
  if (normalise ==  "CQN") {
    cat(green("CQN Normalisation underway...\n"))
    counts_CQN <- tweeDEseq::normalizeCounts(counts, method = "CQN")
    if (writeOut) {
      cat("Writing CQN Normalised Counts")}
    if (!dir.exists("CountMatrices")) { # Create folder
      dir.create("CountMatrices")}
    write.csv(counts_CQN,file.path(countspath,paste0(expname,"CQN_counts.csv")),row.names = FALSE)
    
    if (QCplots) {  
      ## Density Plot
      dnames <-  dimnames(counts)
      colors <-  sequential_hcl(ncol(counts), "Viridis")
      if (!dir.exists("QC_plots")) { # Create folder
        dir.create("QC_plots")}
      qcpath <- paste0(getwd(),"/QC_plots") # save path
      pdf(file = paste0(qcpath,"/Density_plots_CQN_normalisation.pdf"),onefile = T)
      cat("Plotting Densities...\n")
      plot(density(x = as.numeric(unlist(counts))),main = "Un-Normalised",
           xlim = c(0,20000),xlab = "Exp")
      mask <-  sample(1:length(dnames[[2]]))
      ret <-  lapply(mask[2:length(mask)],function(x) {
        i = which(mask = x)
        lines(density(as.numeric(counts[,x])),col = colors[i])
      })
      plot(density(x = as.numeric(unlist(counts_CQN))),main = "CQN Normalised",
           xlim = c(0,20000),xlab = "Exp")
      mask <-  sample(1:length(dnames[[2]]))
      ret <-  lapply(mask[2:length(mask)],function(x) {
        i = which(mask = x)
        lines(density(as.numeric(counts[,x])),col = colors[i])
      })
      graphics.off()
      counts <- counts_CQN # reassign to counts
    }
  }
  
  ## Gene-wise Trimmed Means of M-values normalisation
  
  if (normalise == "geTMM" | normalise == TRUE) {
    
    cat(green("GeTMM Normalisation Underway...\n"))
    
    if (is.null(species) | species ==  "mouse") {
      martspecies <- 'mmusculus_gene_ensembl'
      orgdb <- 'mmu'
    } else if (species ==  'human') {
      martspecies  <- c("hsapiens_gene_ensembl")
      orgdb <- "hg38"
    }
    
    ens_char <- as.character(rownames(counts)) #ensembl IDs and character
    gc_length <- getGeneLengthAndGCContent(ens_char,orgdb, mode = "org.db") #get gc and length information
    gc_length <- gc_length %>% as.data.frame() # as df
    gc <- data.frame(ensembl = rownames(gc_length),gc = gc_length %>% 
                       dplyr::select(gc) %>% 
                       dplyr::transmute(gc = gc*100)) # as a percentage
    len <- data.frame(ensembl = rownames(gc_length),length = gc_length$length) #because we need per kbp
    len[is.na(len$length),2] <- median(len[!is.na(len$length),][,2]) # make NAs median length
    
    ## GeTMM
    rpk <- (10e3*counts/len$length) # reads per 1000 bp
    y <- DGEList(counts=rpk,group = factor(meta$condition),remove.zeros = TRUE)
    y <- calcNormFactors(y) 
    norm_counts <- cpm(y) %>% as.data.frame()
    counts <- norm_counts
    write.csv(norm_counts,paste0(expname,"_geTMM_normalized_countMatrix.csv"))
    
  }
  
  ## CPM Filter
  
  if (CPMfilter) {
    meta  <- mutate_all(meta, function(x) as.factor(x)) # convert to factor
    counts_filt <-  filtered.data(counts,
                                  factor = meta$condition,
                                  norm = TRUE,
                                  method = 1,
                                  cv.cutoff = 95,
                                  cpm = CPMcutoff, 
                                  p.adj = "fdr")
    .GlobalEnv$norm_filt_counts <- counts_filt
    counts <- counts_filt
    
    write.csv(counts_filt,paste0(expname,"_GeTMM_",CPMcutoff,"CPMfilt_Counts.csv"))
    
    ## Background genes for gsea    
    
    if (!dir.exists("GSEA")) {
      dir.create("GSEA")
      setwd(paste0(directory,"/GSEA"))
    }
    if (is.null(species) | species ==  "mouse") {
      org <- "org.Mm.eg.db"   
    } else if (species ==  'human') {
      org <- "org.Hs.eg.db"   
    }
    
    allGenes <- bitr(rownames(counts), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org)$ENTREZID 
    .GlobalEnv$allGenes <- allGenes
    write_rds(allGenes, "allGenes.rds")
    setwd(directory)
  }
  
  if (orthologs) {
    
    if (is.null(species) | species ==  "mouse") {
      org <- "org.Mm.eg.db"   
      from = "mouse"
      to = "human"
    } else if (species ==  'human') {
      org <- "org.Hs.eg.db"   
      from = "human"
      to = "mouse"
    }
    
    geneIDs <-  bitr(rownames(counts), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org)$ENTREZID
    
    # Get marts
    cat(yellow("Getting Marts...\n"))
    
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
    
    cat(blue("Mapping orthologs...\n"))
    identifier <- 'entrezgene_id'
    
    tmp <- getLDS(attributes = identifier,
                  filters = identifier, 
                  values = geneIDs , 
                  mart = marts[[from]], 
                  attributesL = identifier,
                  martL = marts[[to]], 
                  uniqueRows = TRUE)
    
    orthoGenes <- data.frame(ENTREZID = as.character(unique(tmp[,2])))
    write.csv(orthoGenes, paste0(from,"2",to,"ortholog_entrezid.csv"))
    .GlobalEnv$orthoGenes <- orthoGenes
  }
  
  
  
  
  ## Outliers
  
  if (dropOut == "grubbs") {
    
    cat(green("Detecting Outliers using Grubb's method...\n"))
    
    meta <- mutate_all(meta, function(x) as.factor(x))
    colours <- list()
    finalcolours <- list()
    meta2 <- meta %>% dplyr::select(-c(id,priority,class)) %>% colnames()
    names(meta2) <- meta2
    
    for(i in names(meta2)) {
      
      pdf(paste0(expname,"_",i,"_MDS.pdf"))  
      colours[[i]] <-  sequential_hcl(length(unique(meta[,i])),palette = "Viridis")
      finalcolours[[i]] <-  colours[[i]][as.numeric(meta[,i])]
      plotMDS(counts,col = finalcolours[[i]],
              main = paste0("MDS using ",i),pch = 16)
      legend("topright",fill = colours[[i]],
             legend = levels(meta[,i]))
    }
    graphics.off() 
    
    
    ## Grubbs test and Dendogram 
    
    outlist <- NULL
    
    sampdist <-  as.matrix(dist(CoExpNets::trasposeDataFrame(counts,F)))
    result <-  grubbs.test(apply(sampdist,1,sum))
    
    if (result$p.value < 0.05) {
      cat(paste0("Offenders:",names(result$p.value < 0.05)),"\n")
      capture.output(names(result$p.value < 0.05),file = "theOutsiders.txt")
      pdf("outliers_dendograms.pdf")
      flashres <- flashClust(dist(CoExpNets::trasposeDataFrame(counts,F)), method = "average") #%>% as.dendrogram()
      fviz_dend(flashres, cex = 0.8, lwd = 0.8, k = length(levels(meta$condition)) + 1, 
                rect = TRUE, 
                k_colors = "jco", 
                rect_border = "jco", 
                rect_fill = TRUE,
                type = "rectangle",
                repel = TRUE)
    }
    graphics.off()
    meta <- mutate_all(meta, function(x) as.factor(x))
    ## Outlier Removal
    
    if (any(length(result$p.value < 0) <=  0)) {
      
      cat(green("No Outliers to boot!\n")) }
    
    if (any(length(result$p.value < 0) > 0)) {
      
      cat(red("We don't take kindly to outliers 'round here. Boot 'em!\n"))
      
      drop <- c(names(result$p.value < 0.05)) # Drop samples sig in Grubbs test
      counts_osr <-  counts[,!(names(counts) %in% drop)] # from counts matrix  
      meta_osr <- meta[!(rownames(meta) %in% drop),] # from metadata
      
    }
    
    if (any(!colnames(counts_osr) ==  rownames(meta_osr))) { 
      
      stop("Sample names in metadata don't match those in expression matrix, yo!")}
    
    cat(yellow("Writing Outlier Sample Removed Metadata...\n"))
    
    write.csv(counts_osr,paste0(expname,"_normalised",CPMcutoff,"CPM_","counts_OSR.csv"))
    write.csv(meta_osr,paste0(expname,"_normalised",CPMcutoff,"CPM_","meta_OSR.csv"))
    
    .GlobalEnv$counts_osr <- counts_osr
    .GlobalEnv$meta_osr <- meta_osr  #Normalize samples
    
    counts <- counts_osr
    meta <- meta_osr
  }
  
  
  if (dropOut == "PCAdist") {
    
    
    cat(yellow("Detecting outliers using PCAdist method...\n"))
    
    pdf("PCAdist_outlier_plots.pdf")
    
    obj <- OutlierPCDist(x = t(counts),grouping = as.factor(meta$condition))
    
    classes <- list()
    idxs <- list()
    f <- levels(factor(meta$condition))
    for (i in 1:length(levels(factor(meta$condition)))){
      classes[[i]] <- getFlag(obj,i) %>% as.data.frame()
      idxs[[i]] <- which(classes[[i]][1] == 0) %>% as.data.frame()
      plot(obj, class = i)
    }
    
    graphics.off()
    
    idxs <- rbindlist(idxs,use.names = T)
    outs <- idxs[!duplicated(idxs$.),] 
    outs <- outs$.
    counts <- counts[,-c(outs)] %>% as.data.frame()
    meta <- meta[-c(outs),] %>% as.data.frame()
    .GlobalEnv$meta_osr <- meta
    .GlobalEnv$norm_osr_counts <- counts
    write.csv(meta_osr,paste0(expname,'_meta_osr.csv'))
  }
  
  ## Batch Correction
  
  if (batch) {
    
    meta2 <- meta %>% select(-c(id,priority,class)) %>% colnames()
    names(meta2) <- meta2
    cat(green("Batch correction underway...\n"))
    meta <- mutate_all(meta, function(x) as.factor(x)) # all as factor
    rownames(meta) <- colnames(counts) #rename rownames as counts
    counts <- sva::ComBat(dat = counts,batch = meta$batch,
                          mod = model.matrix(~1,data = meta[,names(meta2)]),
                          par.prior = T)
    covariate <- "batch"
    mask <- sample(1:ncol(counts),ncol(counts))
    colors = sequential_hcl(length(unique(as.numeric(meta[,covariate]))),"Viridis")
    finalcolors = colors[as.numeric(meta[,covariate])][mask]
    
    if (!dir.exists("QC_plots")) { # Create folder
      dir.create("QC_plots")}
    qcpath <- paste0(getwd(),"/QC_plots") # save path
    pdf(paste0(qcpath,"/Batch_corrected_mds.pdf"))
    plotMDS(counts[,mask],col = finalcolors,
            main = paste0("Batch Corrected MDS"))
    legend("topleft",fill = colors,
           legend = levels(meta[,covariate]))
    graphics.off()
    
    # Plot covariate contribution
    
    if (!dir.exists("QC_plots")) { # Create folder
      dir.create("QC_plots")}
    qcpath <- paste0(getwd(),"/QC_plots") # save path
    #for(i in names(princeList)) {
    pdf(paste0(qcpath,"/","prince_batch_correction_heatmap.pdf"))
    pcres <-  prince(as.matrix(counts),meta[,colnames(meta) !=  "condition"],top = ncol(counts))
    CoExpNets::princePlot(prince = pcres,main = "Batch Uncorrected")
    pcres <-  prince(as.matrix(counts),meta[,colnames(meta) !=  "condition"],top = ncol(counts))
    CoExpNets::princePlot(prince = pcres,main = "Batch Corrected")
    graphics.off()
    
    .GlobalEnv$combat_counts <- counts
    counts <- counts
  }
  
  ## SVA correction
  
  if (svaCor) {
    
    meta <- mutate_all(meta, function(x) as.factor(x)) # all as factor
    newmeta <- meta[,c(svaCovar)] %>% as.data.frame() # only covara of interest
    colnames(newmeta) <- svaCovar
    newmeta <- mutate_all(newmeta, function(x) as.factor(x)) # all as factor
    newmeta <- lapply(newmeta, function(x) droplevels(x)) %>% as.data.frame() # drop any hijacking levels
    .GlobalEnv$newmeta <- newmeta
    
    mm = model.matrix(svaFormula, data = newmeta) # formula
    nullm = model.matrix(~ 1,data = newmeta) # null formula
    
    counts <- counts %>% as.data.frame() 
    
    counts1 <- counts[!apply(counts <= 0,1,any,na.rm = TRUE),] # make sure there's nothing <= 0 or NAs
    counts0 <- counts[apply(counts <= 0,1,any,na.rm = TRUE),] # make sure there's nothing <= 0 or NAs
    
    ## Choose number of surrogate variables
    
    cat(calc("\nChoosing number of SVs...\n"))
    
    svas <- sva::svaseq(dat = as.matrix(counts1),
                        mod = mm,
                        mod0 = nullm,
                        B = 40)
    
    
    .GlobalEnv$svas <- svas
    
    # Test for latent variables which correlate to condition of interest and
    
    numeric.meta <- mutate_all(newmeta, function(x) as.numeric(x)) #as numeric
    
    linp <-  matrix(ncol = svas$n.sv,
                    nrow = ncol(numeric.meta)) #NA matrix
    
    rownames(linp) <-  colnames(numeric.meta) #rename rows
    colnames(linp) <-  paste0("SV",1:svas$n.sv) #rename columns
    
    linp[] <-  0 #0 matrix
    
    for(cov in 1:ncol(numeric.meta)) {
      for(sva in 1:svas$n.sv) {
        if (svas$n.sv == 1)
          axis = svas$sv
        else
          axis = svas$sv[,sva]
        linp[cov,sva] = cor.test(as.numeric(numeric.meta[,cov]),axis)$p.value #correlation pvals
      }}
    
    # Plot
    
    smallest = -10
    linp10 <- log10(linp) 
    linp10 <- replace(linp10, linp10 <= smallest, smallest)  # recode smallest to -10
    tonote <- signif(linp, 1) # 1 sig figs
    ## Corrplots - Which Covaririates are these SVs correlated to??
    
    pdf("SVA_Covariate_Correlation.pdf")
    if (nrow(linp10) == 1) {
      corrplot::corrplot(as.matrix(linp10),is.corr = F,
                         col =  sequential_hcl(10, "Viridis"),
                         method = "square")
      
    } else { 
      gplots::heatmap.2(linp10, Colv = F, Rowv = F, dendrogram = "none",
                        trace = "none", symbreaks = F, symkey = F, breaks = seq(-20, 0, length.out = 100),
                        key = T, colsep = NULL, rowsep = NULL, sepcolor = "black",
                        sepwidth = c(0.05, 0.05), main = "Correlation - SVs-Covariates",
                        labCol = colnames(linp10),
                        labRow = colnames(newmeta),
                        xlab = "Surrogate variables")
    }
    graphics.off()
    
    # Whip out corrected counts matrix
    
    meta.rs <-  newmeta[,match(colnames(numeric.meta),colnames(newmeta))] %>% as.data.frame()
    
    cat(prep("Discarding SVs correlated to variables of interest...\n"))
    
    # Drop significant correlations to variables of interest
    
    if (is.null(conCovar)) {
      conCovar <- c()
    }
    
    linp_t <- t(linp) %>% as.data.frame()
    idx_drop <- sapply(dplyr::select(linp_t,-conCovar), function(x) which(x < 0.05))
    idx_keep <- sapply(dplyr::select(linp_t,conCovar), function(x) which(x < 0.05))
    
    if (length(idx_keep) >= 1){
      idx_keep2 <- unlist(idx_keep)[-c(unlist(idx_drop))]
    } else {
      idx_keep2 <- c(1:nrow(linp_t))[-unlist(idx_drop)]
    }
    
    if (length(idx_keep2) <= 1 ) {
      
      cat(yellow("Not enough significant surrogate variables!\n"))
      cat(green("Continuing...\n"))
      
    } else {
      
      cat(blue("Generating SV corrected count matrix...\n"))
      
      svs <- svas$sv %>% as.data.frame()
      svs <- svs[,c(idx_keep2)] %>% as.data.frame()
      Y <- as.matrix(t(counts1))
      W <- cbind(mm,svs) %>% as.matrix()
      alpha <- (solve(t(W) %*% W) %*% t(W)) %*% Y
      SVAcor_counts <- t(Y - W[,-c(1:ncol(mm))] %*% alpha[-c(1:ncol(mm)),]) %>% as.data.frame()
      SVAcor_counts <- lapply(SVAcor_counts, function(x) round(x,digits = 3)) %>% as.data.frame()
      SVAcor_counts[SVAcor_counts == 0] <- 1e-16 # get rid of 0s
      rownames(SVAcor_counts) <-  rownames(counts1)
      
      cat(yay("SVA corrected matrix complete!\n"))
      
      # Pals and such
      
      levelpals <-  lapply(newmeta, function(x) sequential_hcl(length(levels(x)),palette = "viridis"))
      leg <- lapply(newmeta, function(x) colnames(x))
      names(leg) <- colnames(newmeta)
      groups <- lapply(newmeta, function(x) as.factor(x))
      linp_t <- t(linp) %>% as.data.frame()
      
      if (length(idx_keep2) >= 1 ) {
        
        cat(yellow("Plotting PCAs...\n"))
      
        # PCA
        
        sva_pca <- prcomp(t(SVAcor_counts),scale = TRUE, center = T)
        res_pca <- prcomp(t(counts1), scale = TRUE, center = T) 
        
        #Plots
        
        pdf("SVA_PCA.pdf")
        
        for (i in seq_along(leg)) {
          
          print(fviz_pca_ind(res_pca, 
                             axes = c(1, 2),
                             geom = c("point"),
                             label = "all", 
                             palette = c(levelpals[[i]]),
                             invisible = "none", 
                             labelsize = 4,
                             pointsize = 4, 
                             habillage = "none",
                             addEllipses = F, 
                             ellipse.level = 0.95, 
                             col.ind = groups[[i]], 
                             alpha.ind = 1,
                             repel = T, 
                             legend.title = names(leg[i])))
          
          print(fviz_pca_ind(sva_pca, 
                             axes = c(1, 2),
                             geom = c("point"),
                             label = "all", 
                             palette = levelpals[[i]],
                             invisible = "none", 
                             labelsize = 4,
                             pointsize = 4, 
                             habillage = "none",
                             addEllipses = F, 
                             ellipse.level = 0.95, 
                             col.ind = groups[[i]], 
                             alpha.ind = 1,
                             repel = T, 
                             legend.title = names(leg[i])))
          
        }
        graphics.off()
        
        if (includeZeros) {
        
        colnames(SVAcor_counts) <- gsub('X','', colnames(SVAcor_counts))
        SVAcor_counts <- rbind(SVAcor_counts,counts0)
        SVAcor_counts[SVAcor_counts == 0] <- 1e-16 # planck nudge
        
        }
        cat(yay("PCA plots exported!...\n"))
        counts <- SVAcor_counts 
        .GlobalEnv$SVAcor_counts <- SVAcor_counts
        write_rds(counts,"norm_filt_osr_svaCor_counts.rds")
        cat(yay("SVAcor complete!\n"))
        }
      }
    }
  
    if (DE == "noiseq") {
    
    setwd(directory)
    
    if (!dir.exists("NoiSeq")) { # Create folder
      
      dir.create("Noiseq")
      
      setwd(paste0(directory, "/Noiseq"))
      
    }
    
    # Annotations        
    
    if (is.null(species) | species ==  "mouse") {
      martspecies <- 'mmusculus_gene_ensembl'
      orgdb <- 'mmu'
      symb = "mgi_symbol"
      .GlobalEnv$symb
    } else if (species ==  'human') {
      martspecies  <- c("hsapiens_gene_ensembl")
      orgdb <- "hg38"
      symb = "hgnc_symbol"
      .GlobalEnv$symb <- symb
    }
    
    cat(yellow("Collecting  gene annotations...\n"))
    
    genes <- rownames(counts) #gene names 
    
    cat(magenta("Getting Mart...\n"))
    
    for (i in 1:1) {
      mart = tryCatch(useMart(biomart = 'ENSEMBL_MART_ENSEMBL',martspecies[i])) #get mart
      if (inherits(mart, "error")) {
        next 
        mart = tryCatch(useMart(biomart = 'ensembl',martspecies[i]))
      }}
    
    ids <- getBM(attributes = c('ensembl_gene_id',                                
                                'entrezgene_id',
                                symb,
                                'chromosome_name',
                                'strand',
                                'start_position',
                                'end_position',
                                'gene_biotype',
                                'ensembl_peptide_id',
                                'wikigene_description',
                                'transcription_start_site'),
                 filters = "ensembl_gene_id", values = genes,
                 mart = mart) 
    
    ## Gene Annotations
    
    anno <- left_join(data.frame(ensembl_gene_id = rownames(counts)),ids, by = "ensembl_gene_id") #left join
    anno <- anno[!duplicated(anno$ensembl_gene_id),] #remove weird generted duplicate
    anno <- anno %>% mutate(select_biotypes = gsub("protein_coding",">protein_coding",anno$gene_biotype)) #">" anchor for later filtering
    anno$gene_biotype <- forcats::fct_explicit_na(anno$gene_biotype, "unannotated_element") #nice funtion to replace NAs with a factor
    anno$select_biotypes <- forcats::fct_explicit_na(anno$select_biotypes, ">unannotated_element") #nice function to replace NAs with a factor
    .GlobalEnv$annotations <- anno
    ## GC, length and chromosomes
    ens_char <- as.character(rownames(counts)) #ensembl IDs and character
    cat(magenta("Getting length and GC content\n"))
    gc_length <- getGeneLengthAndGCContent(ens_char,orgdb, mode = "org.db") #get gc and length information
    gc_length <- gc_length %>% as.data.frame() # as df
    gc <- data.frame(ensembl = rownames(gc_length),gc = gc_length %>% 
                       dplyr::select(gc) %>% 
                       dplyr::transmute(gc = gc*100)) # as a percentage
    len <- data.frame(ensembl = rownames(gc_length),length = gc_length$length)
    bio <- anno[,c(1,7,11)] # biotype df
    chrom <- anno[,c(3,5,6)] #chromosome info
    rownames(chrom) <- bio$ensembl_gene_id #rownames as gene names
    ##Order
    gc <- gc[order(bio$ensembl_gene_id),] %>% as.data.frame()
    len <- len[order(bio$ensembl_gene_id),]
    counts <- counts[order(match(rownames(counts),bio$ensembl_gene_id)), , drop = FALSE] %>% as.data.frame()
    rownames(counts) =  bio$ensembl_gene_id 
    
    annolist <- list(gc = gc,length = len,chrom = chrom,biotypes = bio,gene_annotations = anno)
    
    for(i in names(annolist))
      
      write.csv(i,paste0(expname,"_",i,".csv"))
    
    .GlobalEnv$annoList <- annolist
    
    setwd(directory)
    
    
    ## Make eSet Object
    
    tmp <- list.files(pattern = "_noiseq_results.rds", recursive = TRUE)
    if(!isEmpty(tmp)){
      res <- read_rds(tmp)
    } else { 
      
      
      setwd(paste0(directory,"/Noiseq"))
      
      conds <- data.frame(conds = meta$condition)
      conds$conds <- conds$conds %>% as.factor
      
      noise <- NOISeq::readData(data = counts,
                                length = len, 
                                gc = gc, 
                                biotype = bio[,c(1:2)],
                                chromosome = chrom, 
                                factors = conds)
      
      write_rds(noise,paste0(expname,"_Noiseq_eSet.rds"))
      
      ## Make comparison for looping
      
      meta$pcond <- paste0(meta$priority,meta$class,meta$condition)
      meta <- mutate_all(meta,function(x) as.factor(x))
      
      l <- length(levels(meta$condition)) #number of levels in condition
      levs <- levels(meta$pcond)
      
      combos <- combinations(n = l, r = 2, v = 1:l, repeats.allowed = FALSE) %>% 
        as.data.frame() # encode unique combinationns
      combos <- mutate_all(combos, function(x) as.factor(x)) # convert to factors
      n <- data.frame(V1 = factor(c(1:l)),V2 = c(levs), stringsAsFactors = FALSE) # Make df for mapping levels to combo matrix
      facs <- left_join(combos,n, by = "V1") #leftjoin
      colnames(facs)[c(1,2)] <- c("V","V1") #adjust rownames
      factors <- left_join(facs,n, by  = "V1") #leftjoin again for remaining levels
      factors <- factors[,c(3,4)] #remove uneeded columns
      colnames(factors) <- c("d1","d2") #change colnames
      factors <- factors[,c(1,2)]
      
      ## Order to desired pairwise comparison
      idx1 <- which(str_detect(string = factors[,2], pattern = "disease"))   
      factors[c(idx1),c(1,2)] <- factors[c(idx1),c(2,1)]
      idx2 <- which(str_detect(string = factors[,2], pattern = as.character(max(as.numeric(as.character(meta$priority))))))  
      factors[c(idx2),c(1,2)] <- factors[c(idx2),c(2,1)]
      
      factors <- lapply(factors, function(x) stri_replace_all_regex(x, "^[0-9]", replacement = "", vectorize_all = FALSE)) %>% 
        as.data.frame()
      factors <- lapply(factors, function(x) mgsub_fixed(x,pattern = c("disease","control"),
                                                         replacement = c(""), 
      )) %>% 
        as.data.frame()
      
      write.csv(factors,"DGEA_pairwise_comparisons.csv")
      
      # Pause to check comparisons
      
      readkeyz <- function()
        
      {
        cat("If you wish to remove or modify any comparison, do so in the csv called 'DGEA_pairwise_comparisons' which is found
          in the 'Noiseq' folder and save it. Then, or otherwise, on your keyboard and in the console, Press [Enter] to continue")
        line <- readline()
      }
      readkeyz()
      
      factors <-  read.csv("DGEA_pairwise_comparisons.csv", row.names = 1)
      facnames <- paste0(factors$d1,"vs",factors$d2) #make names for list
      .GlobalEnv$factors <- factors
      .GlobalEnv$facnames <- facnames
      
      st <- list()
      for (i in 1:nrow(factors)) {
        
        st[[i]] <- c(factors[i,1],factors[i,2])
      }
      names(
        st) <- facnames #rename elements
      
      .GlobalEnv$st <- st
      
      ## DGEA
      
      cat(green("\nNoiseq underway...\n"))
      gc()
      rm(mm)
      res <-  list() #result list
      for (i in 1:length(st)) {
        
        set.seed(seed)
        res[[i]] <- noiseqbio(input = noise,
                              factor = "conds",
                              norm = "n",
                              conditions = st[[i]], #loop across each comparison
                              plot = F,
                              random.seed = seed,
                              a0per = 0.9, #!a0per
                              r = 200,
                              filter = 0) 
        
      }
      
      gc()
      
      # Write result
      
      names(res) <- names(st) #rename list elements
      .GlobalEnv$noiseq_result <- res
      res <- noiseq_result 
      write_rds(res,paste0(expname,"_noiseq_results.rds"))
    }
    ## Extract all genes
    
    
    cat(blue("Extracting results...\n"))
    setwd(paste0(directory,"/Noiseq"))
    
    factors <-  read.csv("DGEA_pairwise_comparisons.csv", row.names = 1)
    facnames <- paste0(factors$d1,"vs",factors$d2) #make names for list
    
    if (!dir.exists("Tables")) {
      dir.create("Tables") }
    setwd(paste0(getwd(),"/Tables"))
    
    
    all_list <- list() #de list
    for (i in seq_along(res))
      
      all_list[[i]] = degenes(res[[i]], q = 0, M = NULL) #1-q = pval.adjust
    
    names(all_list) <- facnames #rename
    
    ## Make dataframe for biotypes selectio
    anno = annotations
    bio_sel <- data.frame(ensembl_gene_id = anno$ensembl_gene_id,select_biotypes = anno$select_biotypes)
    
    joins <- list()
    for (i in seq_along(all_list)) {
      joins[[i]] <- data.frame(ensembl_gene_id = rownames(all_list[[i]]),all_list[[i]],stringsAsFactors = FALSE)
      joins[[i]] <- left_join(joins[[i]],anno, by = "ensembl_gene_id")
      joins[[i]] <- data.frame(joins[[i]], condition = rep(facnames[i],nrow(joins[[i]])),stringsAsFactors = FALSE)
      joins[[i]] <- joins[[i]] %>% 
        dplyr::mutate(hiprob = prob-10e-17) %>% #cant be having -log(0)
        dplyr::mutate(pval.adjus = 1-hiprob )
    }
    
    names(joins) <- facnames
    for (i in names(joins)) {
      write.csv(joins[[i]],paste0(expname,i,"_AllGenes_NoiSeq.csv"))}
    ## If there's duplicates for some reason
    joins_deduped <- list()
    for (i in 1:length(joins)) {
      joins_deduped[[i]] <- joins[[i]][!duplicated(joins[[i]]$ensembl_gene_id),]
    }
    names(joins_deduped) <- facnames #rename list elements
    ## Extract and write only significant genes
    
    sigDEG <- list()
    for (i in seq_along(joins_deduped)) {
      joins_deduped[[i]] <- joins_deduped[[i]][order(joins_deduped[[i]]$hiprob,abs(joins_deduped[[i]]$log2FC),decreasing = TRUE),] #order
      sigDEG[[i]] <- joins_deduped[[i]][joins_deduped[[i]]$pval.adjus <=  5e-2,] #sig cutoff fixed at 0.05
      #sigDEG[[i]] <- sigDEG[[i]][complete.cases(sigDEG[[i]]),] #remove nas
    }
    names(sigDEG) <- facnames #rename list elements
    .GlobalEnv$sigDEG <- sigDEG
    write_rds(sigDEG,"sigDEGs.rds")
    ## Write
    
    for (i in names(sigDEG)) {
      write.csv(sigDEG[[i]], paste0(expname,"_",i,"_sigGenes_",nrow(sigDEG[[i]]),".csv"))
    }
    
    ## Extract non-procoding
    
    nonpro <- lapply(sigDEG, function(x) x[!grepl(x$select_biotypes,pattern = ">"),])
    
    for (i in names(nonpro)) {
      write.csv(nonpro[[i]], paste0(expname,"_",i,"_sigGenes_nonProteinCoding_",nrow(nonpro[[i]]),".csv"))
    }
    
    sList <- list()  
    for (i in names(nonpro)) {
      sList[[i]] <- nonpro[[i]] %>% dplyr::select(ensembl_gene_id,symb,gene_biotype,wikigene_description,log2FC) %>%
        as_hux()
    }
    sList <- lapply(sList, function(x) x[order(x$log2FC,decreasing = T),])
    names(sList) <- names(nonpro)
    formats = list() 
    for (i in names(sList)) {
      formats[[i]] = sList[[i]] %>% theme_article()
    }
    for (i in names(formats)) {
      quick_html(formats[[i]],file = paste0(i,"_nonProcoding.html"))
    }
    
    ## In the darkness bind them
    
    sig_bound <- rbindlist(sigDEG, fill = T)
    write.csv(sig_bound, paste0(expname,"_MasterSummary_NoiSeq_DEGs.csv"))
    nonpro_bound <- rbindlist(nonpro, fill = T)
    write.csv(nonpro_bound, paste0(expname,"_MasterSummary_NoiSeq_NonProCoding_DEGs.csv"))
    setwd(directory)
    
    ## Volcano Plots
    
    cat(yellow("Generating Volcanos...\n"))
    
    vol <- list()
    
    for (i in 1:length(joins_deduped)) {
      vol[[i]] <- data.frame(logFC = joins_deduped[[i]]$log2FC,Adjusted.Pvalue = joins_deduped[[i]]$pval.adjus,
                             name =joins_deduped[[i]]$ensembl_gene_id) #df
      vol[[i]] <- vol[[i]][order(abs(vol[[i]]$logFC),decreasing = TRUE),] #order largest to smallest
      vol[[i]] <- vol[[i]][vol[[i]]$logFC != 0,] #remove 0 foldchange genes
      vol[[i]] <- drop_empty_row(vol[[i]]) #remove any empty rows
      vol[[i]] <- inner_join(vol[[i]],data.frame(name = anno$ensembl_gene_id,desc = anno$wikigene_description),
                             by = "name") 
      
    }
    names(vol) <- facnames #rename list elements
    
    # Plot em
    
    setwd(paste0(directory,"/Noiseq"))
    
    if (!dir.exists("Volcanos")) {
      dir.create("Volcanos") }
    setwd(paste0(getwd(),"/Volcanos"))
    
    pdf(paste0(expname,"volcanos_withnames.pdf")) # coloured by logfc and named based on lowest p-val logfc > 6 or logfc > 10
    for (i in seq_along(vol)) {
      with(vol[[i]], plot(logFC, -log10(Adjusted.Pvalue), pch = 20, main = paste(mgsub(x = names(vol[i]),c("vs"), c(" vs "))))) #rename main title
      abline(h = -log10(5e-2),col = "red",lty = 2) #significane boundaries
      with(subset(vol[[i]], Adjusted.Pvalue >=5e-2 & abs(logFC)>=0), points(logFC, -log10(Adjusted.Pvalue), 
                                                                            pch = 20, col = "#252525"))
      with(subset(vol[[i]], logFC > 0 & Adjusted.Pvalue <=  5e-2), points(logFC, -log10(Adjusted.Pvalue), 
                                                                          pch = 20, col = "#c62828"))
      with(subset(vol[[i]], logFC < 0 & Adjusted.Pvalue <= 5e-2), points(logFC, -log10(Adjusted.Pvalue), 
                                                                         pch = 20, col = "#1565c0"))
      with(subset(vol[[i]], Adjusted.Pvalue <= 10e-16 & abs(logFC) >= 6 | abs(logFC) >= 8), textxy(logFC, -log10(Adjusted.Pvalue), labs = desc,
                                                                                                   cex = .5))
    }
    graphics.off()
    
    pdf("GSE126848_noiseqDE_volcanos_nonames.pdf") # Same as above but no gene names
    for (i in seq_along(vol)) {
      with(vol[[i]], plot(logFC, -log10(Adjusted.Pvalue), pch = 20, main = paste(mgsub(x = names(vol[i]),c("vs"), c(" vs ")))))
      abline(h = -log10(5e-2),col = "red",lty = 2)
      with(subset(vol[[i]], Adjusted.Pvalue >=5e-2 & abs(logFC)>=0), points(logFC, -log10(Adjusted.Pvalue), 
                                                                            pch = 20, col = "#252525"))
      with(subset(vol[[i]], logFC > 0 & Adjusted.Pvalue <=  5e-2), points(logFC, -log10(Adjusted.Pvalue), 
                                                                          pch = 20, col = "#c62828"))
      with(subset(vol[[i]], logFC < 0 & Adjusted.Pvalue <= 5e-2), points(logFC, -log10(Adjusted.Pvalue), 
                                                                         pch = 20, col = "#1565c0"))
    }
    graphics.off()
    setwd(directory) #reset wd
  }
  
  ## GSEA
  
  if (GSEA) {
    
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
    
    tmp <- list.files(pattern = "^sigDEGs.rds", recursive = T)
    sigDEGs <- read_rds(tmp)
    
    cat(yellow("Converting ENSEMBL to ENTREZ...\n"))
    
    cat(prep("Converting ENSEMBL to ENTREZ...\n"))
    
    list <- lapply(sigDEGs, function(x) bitr(x$ensembl_gene_id, fromType = "ENSEMBL", toType = "ENTREZID",
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
    } else {
      entrez_DEG <- entrez_sigList
    }
    
    
    # Extract information
    
    if (species ==  "mouse") {
      
      geneListDes <- lapply(entrez_DEG, function(x) data.frame(geneID = x$entrezgene_id,
                                                               symbol = x$mgi_symbol,
                                                               p = x$pval.adjus,
                                                               logFC = x$log2FC
      ))
    } else if (species ==  'human') {
      geneListDes <- lapply(entrez_DEG, function(x) data.frame(geneID = x$ENTREZID,
                                                               symbol = x$hgnc_symbol,
                                                               p = x$pval.adjus,
                                                               logFC = x$log2FC
      ))
    }
    
    geneListDes  <- lapply(geneListDes, function(x) x[!duplicated(x$geneID),]) #remove duplicates
    
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
    geneinput <- geneinput[sapply(geneinput, function(x) length(x)[1]) > 100] #remove <100
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
    
    # Background genes
    
    setwd(directory)
    
    if (!is.null(orthogenes)) {
      
      allGenes <- orthogenes$ENTREZID2
      
    } else {
      
      tmp <- list.files(pattern = "allGenes.rds", recursive = T)
      allGenes <- read_rds(tmp[1])      
      
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
      
      for (i in seq_along(geneinput)) {
        sraList[[i]] <- setRankAnalysis(geneIDs = geneinput[[i]],
                                        use.ranks = T,
                                        setCollection = collection,
                                        setPCutoff = 0.01,
                                        fdrCutoff = 0.05)}
      
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
      
      if (msigDB == FALSE && CTDbase == TRUE) {
        tmp <- list.files(pattern = "colist_ctd.rds", recursive = T)
        colist <- read_rds(tmp)
      } else if (msigDB == TRUE) {
        tmp <- list.files(pattern = "colist.rds", recursive = T)
        colist <- read_rds(tmp)
      }
      
      tmp <- list.files(pattern = "colist", recursive = T)
      colist <- read_rds(tmp) %>% as.data.frame()
      .GlobalEnv$colist <- colist
      
      tmp <- list.files(pattern = "geneinput.rds", recursive = T)
      #tmp <- tmp[-1]
      geneinput <- read_rds(tmp)
      .GlobalEnv$geneinput <- geneinput
      setwd(paste0(directory,"/GSEA"))
      collection <- read_rds("setRank_collection.rds")
      
      
      if (!is.null(orthospecies)) {
        if (orthospecies == 'human') { 
          org <- "org.Hs.eg.db"
          spec  <- 'Homo sapiens'
          symb <- "hgnc_symbol"
        } else {
          spec <- 'Mus musculus'
          org <- "org.Mm.eg.db"
          symb <- "mgi_symbol"
        }
      }
      
      entrez2symbol <-  createIDConverter(org,"ENTREZID","SYMBOL")
      
      if (!dir.exists("Basic_geneSet_networks")) {
        dir.create("Basic_geneSet_networks")
      }
      
      setwd(paste0(directory,"/GSEA/Basic_geneSet_networks"))
      
      cat(calc("Exporting GeneSet Networks...\n"))
      
      for (i in seq_along(geneinput)) {
        
        exportMultipleResults(sraList, 
                              geneinput[i], 
                              collection,
                              IDConverter = entrez2symbol, 
                              outputPath = getwd())
        
      } 
      
      cat(yay("Export Complete!\n"))
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
        tmp <- list.files(pattern = "gene_net_styles.viz", recursive = T)
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
            ppi <- graph_from_data_frame(mippie_ppi,vertices = mippie_pro)}
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
            .GlobalEnv$ppi <- ppi}
          cat(yay("PPI Network Complete!\n"))
        }
        
        cat(blue("Mapping GSEA results to ppi networks...\n"))
        
        setwd(paste0(directory,"/GSEA"))
        if (!dir.exists("PPI_geneSet_networks")) {
          dir.create("PPI_geneSet_networks")
          setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
        }
        
        if (isEmpty(list.files(pattern = "gene_net_styles.viz"))){
          
          setwd(paste0(directory,"/GSEA"))
          geneListDes2 <- read_rds("geneListDes2.rds")
          sraList <- read_rds("SetRankResult.rds")
          collection <- read_rds("setRank_collection.rds")
          .GlobalEnv$sraList <- sraList
          names(sraList) <- names(geneListDes2)
          .GlobalEnv$geneListDes2 <- geneListDes2
          .GlobalEnv$collection <- collection
          
          setwd(paste0(directory,"/GSEA"))
          if(!dir.exists("PPI_geneSet_networks"))
            dir.create("PPI_geneSet_networks")
          setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
          gc()
          cat(blue("Exporting Geneset Networks..."))
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
          
          cat(yay("Networks exported!\n"))
        }
        # Plot networks
        setwd(directory)
        
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
        
        .GlobalEnv$nets <- nets
        .GlobalEnv$netlist <- netlist
        
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
          
          tmp <- list.files(pattern = "entrez_sigDEGs.rds", recursive = T)
          entrez_sigList <- read_rds(tmp)
          entrez_sigList <- entrez_sigList[names(geneinput)]
          .GlobalEnv$entrez_sigList <- entrez_sigList
          
        }
        
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
        
        
        #edge_list <- unlist(split_edges, recursive = FALSE) #one flag, one flag was all that was needed. I will burn my white flag. 
        ## Same with nodez
        
        lengthvec <- seq(from = 1,by = 1,to = length(geneListDes2)) #1,2,3,4,5
        conds <- names(geneListDes2)
        names(conds) <- conds
        cond_vec <- conds[e_vec] #cond names
        names(cond_vec) <- cond_vec
        
        
        # Remove semi-exclusive  
        
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
        
        ## Make networks
        
        edge_joins <- edge_joins[names(edge_joins) %in% names(split_nodes)] # intersect
        split_nodes <- split_nodes[names(split_nodes) %in% names(edge_joins)] 
        
        
        for (i in seq_along(split_nodes)) {
          split_nodes[[i]] <- inner_join(split_nodes[[i]],edge_joins[[i]], by = "ENTREZID") 
        }
        
        
        ## Remove dups, make sure genes are in both node and edgelist
        
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
        
        ## Make network from data frame
        
        gsea_nets <- list()    
        for (i in seq_along(edgelist4net)) {
          gsea_nets[[i]] <- graph_from_data_frame(d = edgelist4net[[i]],
                                                  vertices = split_nodes[[i]],
                                                  directed = F) %>% igraph::simplify(remove.loops = T, remove.multiple = F) %>% 
            delete.vertices(which(igraph::degree(.) == 0))
        } 
        names(gsea_nets) <- names(edgelist4net) #names
        
        # Plots
        
        np <- diverge_hcl(n = 500, "Blue-Red 3") # colour ramp
        ll <- length(entrez_sigList)
        lp <- ceiling(length(entrez_sigList)/2)
        
        ceil <- list()
        for (i in 1:length(gsea_nets)) {
          ceil[[i]] <- max(abs(as.numeric((E(gsea_nets[[i]])$log2FC)))) # ceil
          ceil[[i]][ceil[[i]] == -Inf] <- 0
        }
        
        names(ceil) <- names(edgelist4net)
        sf <-  sapply(ceil, max)
        sf <- round(max(abs(sf)),digits = 2)
        node_cols <- list()
        
        for (i in 1:length(gsea_nets)) {
          node_cols[[i]] <- (as.numeric(E(gsea_nets[[i]])$log2FC) + sf) / (2*sf) * 500 # mid poin
        }
        
        names(node_cols) <- names(gsea_nets)
        
        node_cols1 <- node_cols[sapply(node_cols, function(x) length(x) > 0)]
        node_cols1 <- sapply(node_cols1, function(x) round(x,digits = 0))
        gsea_nets1 <- gsea_nets[names(node_cols1)]
        
        # Plots
        
        setwd(paste0(directory,"/GSEA/PPI_geneSet_networks"))
        
        pdf("gsea_network_plots_log2FC.pdf", width = 11.7, height = 8.3)
        
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
               layout = layout.fruchterman.reingold) }
        graphics.off()
      }  
    }
    
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
        entrez_sigList <- entrez_sigList[names(geneinput)]
        .GlobalEnv$entrez_sigList <- entrez_sigList
        
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
      1:length(nets)
      for (i in 1:length(nets)) {
        netlist[[i]] <- read_graph(format = "graphml", # read as graphs
                                   file = nets[i])
        
      }
      .GlobalEnv$nets <- nets
      .GlobalEnv$netlist <- netlist
      
      ## Prepare gene membership and interactions
      
      
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
          }}
        
        write_rds(pjt,"genesetMap.rds")
        .GlobalEnv$pjt <- pjt
      }
      
      cat(yay("Gene mapping complete!\n"))
      cat(prep("Preparing dataframe for GeneSet Clustering...\n" %+% 
                 wild("This too, will take a while!\n")))
           
      
      names(pjt) <- names(geneinput)
      
      dlist <- entrez_sigList[names(pjt)] # subset frames used for analysis
      
      ## Join to add gene set names
      
      dp <- list()
      for (i in seq_along(pjt)){
        dp[[i]] <- t(pjt[[i]]) %>% as.data.frame()   
        colnames(dp[[i]]) <- dp[[i]][1,]
        dp[[i]] <- data.frame(symb = rownames(dp[[i]]),dp[[i]])
        colnames(dp[[i]])[1] <- symb
        dp[[i]] <- left_join(dlist[[i]],dp[[i]], by = symb)
      }
      
      
      
      ## Map fold changes 
      
      up <- list()      
      down <- list()
      
      for (i in seq_along(dp)){
        up[[i]] <-  dp[[i]][dp[[i]]$log2FC > 0,][,symb]
        #up[[i]] <-  up[[i]][complete.cases(up[[i]])]
        #up[[i]] <- paste0(up[[i]][!is.na(up[[i]])], collapse = ",")
        down[[i]] <-  dp[[i]][dp[[i]]$log2FC < 0,][,symb]
        #down[[i]] <-  down[[i]][complete.cases(down[[i]])]
        #down[[i]] <- paste0(down[[i]][!is.na(down[[i]])], collapse = ",")
      }
      
      tdp <- list()
      ups <- list()
      downs <- list()
      
      for (i in seq_along(dp)){
        tdp[[i]] <- t(dp[[i]]) %>% as.data.frame()
        tdp[[i]] <- tdp[[i]][-1,]
        ups[[i]]  <- apply(tdp[[i]][-c(1:4)],1, function(x) paste(x[x %in% up[[i]]],sep = ","))
        downs[[i]]  <- apply(tdp[[i]][-c(1:4)],1, function(x) paste(x[x %in% down[[i]]],sep = ","))
        tdp[[i]]$Up_regulated <- ups[[i]]
        tdp[[i]]$Down_regulated <- downs[[i]]
        tdp[[i]]$Up_regulated <- gsub('[c.\\)\\(\\"]' ,"",tdp[[i]]$Up_regulated)
        tdp[[i]]$Down_regulated <- gsub('[c.\\)\\(\\"]' ,"",tdp[[i]]$Down_regulated)
        tdp[[i]] <- data.frame(rownames(tdp[[i]]),tdp[[i]])
        colnames(tdp[[i]])[1] <- "name"
      }
      names(tdp) <- names(geneinput)
      
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
      for (i in seq_along(pjt)){
        ppp[[i]] <- inner_join(pp[[i]], tdp[[i]], by = "name")
        
      }
      
      
      
      ## Make df for clustering
      
      clust_list <- list()  
      for (i in seq_along(ppp)) {
        clust_list[[i]] <- data.frame(ID = ppp[[i]]$name,
                                      Term_Description = ppp[[i]]$description,
                                      Up_regulated = ppp[[i]]$Up_regulated,
                                      Down_regulated = ppp[[i]]$Down_regulated,
                                      highest_p = ppp[[i]]$correctedPValue,
                                      lowest_p = ppp[[i]]$adjustedPValue
                                      #Fold_Enrichment = -log10(ppp[[i]]$adjustedPValue)
        )}
      ##Loop this to swap hi/lo p
      idx <- list()
      for (i in seq_along(clust_list)) {
        idx[[i]] <- which(clust_list[[i]]$lowest_p - clust_list[[i]]$highest_p > 0)
        clust_list[[i]][,c(5,6)][idx[[i]],] <- clust_list[[1]][,c(6,5)][idx[[i]],]
      }
    
      ## Cluster    
      
      setwd(paste0(directory,"/GSEA"))
      if(!dir.exists("Representative_GSEA")){
        dir.create("Representative_GSEA")
      }
      setwd("Representative_GSEA")
      
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
      
      names(gs_clustList) <- names(geneinput)
      write_rds(gs_clustList, "gs_clustList.rds")
      .GlobalEnv$gs_clustList <- gs_clustList
      
      #Huxtable
      
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
        quick_html(formats[[i]],file = paste0(i,"_gsea_clusters_tables.html"))
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
        quick_html(formats[[i]],file = paste0(i,"_gsea_representative_only_GS.html"))
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
      
      facs <- str_extract(pattern = "[[:alpha:]\\.]+$", string = names(nodelist4net)) # extract cond names. Use . sep!
      nSplit <- split(nodelist4net,factor(facs)) 
      eSplit <-   split(edgelist4net,factor(facs))
      
      ## Drop incomplete results
      
      gsnames <- list()
      
      for (i in seq_along(nSplit)){
        gsnames[[i]] <- data.frame(ENTREZID = unique(str_extract(pattern = ".*[_]", string = names(nSplit[[i]]))))
        gsnames[[i]]$ENTREZID <-   gsub(pattern = "_$", replacement = "", x = gsnames[[i]]$ENTREZID) # fix
      }
      
      
      namebind <- rbindlist(gsnames)
      
      id <- table(namebind$ENTREZID) < length(levels(factor(facs))) # find incomplete results
      uni <- names(table(namebind$ENTREZID)[id])
      
      if (length(uni) > 0) {
        
        newnSplit <- list()
        neweSplit <- list()
        exc <- list()
        for (i in seq_along(nSplit)){
          exc <-  data.frame(semi_exclusive = names(nSplit[[i]][grepl(paste(uni,collapse="|"),x = names(nSplit[[i]]))]))
          newnSplit[[i]] <- nSplit[[i]][!grepl(paste(uni,collapse="|"),x = names(nSplit[[i]]))]
          neweSplit[[i]] <- eSplit[[i]][!grepl(paste(uni,collapse="|"),x = names(eSplit[[i]]))]
        }
        names(newnSplit) <- names(nSplit)
        names(neweSplit) <- names(eSplit)
        nSplit2 <- newnSplit 
        eSplit2 <- neweSplit 
        write.csv(exc, file = "semi_exclusive_pathways.csv", row.names = FALSE)
        huxdf <- exc %>% as_hux()
        format <- huxdf %>% theme_article()
        quick_html(format, "semi_exclusive_pathways.html")
        rm(neweSplit)
        rm(newnSplit)
      }
      
      # Remove condition label from pathways
      
      gsnames  <- data.frame(ENTREZID = unique(str_extract(pattern = ".*[_]", string = names(nSplit2[[1]])))) # get gsnames
      
      gsnames$ENTREZID <- gsub(pattern = "_$", replacement = "", x = gsnames$ENTREZID) # fix
      
      nSplit3 <- lapply(nSplit2, function(x){
        names(x) <- gsnames$ENTREZID 
        x
        })  
      
      eSplit3 <- lapply(eSplit2, function(x){
        names(x) <- gsnames$ENTREZID 
        x # return x
    
          })  
      
      
      # Bind into data frame
      
      nodeBound <- lapply(nSplit3, function(x) rbindlist(x))
      nodeBound <-  lapply(nodeBound, function(x) as.data.frame(x))
      finalNodes <- lapply(nodeBound, function(x) plyr::rbind.fill(x,gsnames)) # append with GS nodes
      edgeBound <- lapply(eSplit3, function(x) rbindlist(x)) # bind edge
      nSplit <- nSplit3
      eSplit <- eSplit3
      rm(eSplit3,eSplit2,nSplit3,nSplit2)
      
      # Make Edges data frame
      
      grid <- expand.grid(1:length(nSplit),gsnames$ENTREZID)
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
      
      names(GSnodeIX) <- names(nSplit[[1]])
      conds <- unique(str_extract(pattern = "[[:alpha:]\\.]+$", string = names(nodelist4net)))
      
      
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
      GSedgesList <- GSedgesList[names(statdf)]
      finalEdges <- map2(GSedgesList,statJoin, plyr::rbind.fill) # append GS edges to gene edges
      
      finalNodes <- finalNodes[names(statdf)]
      
      # Somtimes a condition/keyword appears in different comparisons. Remove excess columns.
      
      namesvec <- paste0("\\<",names(finalNodes),"\\>") # word boundaries
      
      finalNodes <- lapply(finalNodes, function(x){
        x <- x[,!grepl(pattern = paste(namesvec,collapse = "|"), x = paste0(colnames(x)),fixed = FALSE)] 
        x
      })
      
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
      names(finalNodes) <- conds
      names(finalEdges) <- conds
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
      names(repNets) <- names(finalEdges)
      .GlobalEnv$repNets <- repNets
      # Plot
      
      n <- diverge_hcl(n = 500, "Blue-Red 3") # colour ramp for fold change
      ll <- length(repNets)
      lp <- ceiling(length(repNets)/2)
      ceil <- list()
      for (i in 1:length(repNets)) {
        ceil[[i]] <- max(abs(as.numeric((V(repNets[[i]])$log2FC))),na.rm = T) # ceil
        ceil[[i]][ceil[[i]] == -Inf] <- 0
      }
      
      names(ceil) <- names(finalEdges)
      ceil <- ceil[sapply(ceil, function(x) x > 0)]
      sf <-  sapply(ceil, max)
      sf <- round(max(abs(sf)),digits = 2)
      nodeColors <- list()
      for (i in 1:length(repNets)) {
        nodeColors[[i]] <- (as.numeric(V(repNets[[i]])$log2FC) + sf) / (2*sf) * 500 # mid poin
      }
      pals <- list()
      spPals <- list()
      for (i in seq_along(repNets)) {
        pals[[i]] <- sequential_hcl(length(levels(as.factor(V(repNets[[i]])$sp))),'YlGnBu')
        spPals[[i]] <- pals[as.numeric(as.factor(V(repNets[[i]])$sp))]
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
        quick_html(formatsE[[i]],file = paste0(names(repNets[i]),"repnets_eigen_centrality_sorted.html"))
        write.csv(formatsE[[i]],file = paste0(names(repNets[i]),"repnets_eigen_centrality_sorted.csv"))
        quick_html(formatsB[[i]],file = paste0(names(repNets[i]),"repnets_betweeness_sorted.html"))
        write.csv(formatsB[[i]],file = paste0(names(repNets[i]),"repnets_betweeness_sorted.csv"))
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
        cat(yay("Graphs exported!\n"))
        cat(cyan("Analysis complete!\n"))
          }  
        }
      }
    }
          
