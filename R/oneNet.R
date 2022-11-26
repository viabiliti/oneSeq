#' Automated WGCNA analysis pipeline for Bulk RNAseq data
#' 
#' Run and forget WGCNA analysis 
#' 
#' ceeb
#' 
#' @param counts data frame  Raw bulk mRNA seq counts genes as rows, columns as samples. 
#' @param meta data frame Metadata table with samples as rows and variables as columns. Compulsory columns include 'id' with sample names corresponding to those in counts, 'condition' with relevant tests and controls, 'priority' with integers starting from 0 to 1 for indicating preferred direction of pairwise comparison and 'class' indicating whether the sample is a control or disease/test. 
#' @param species character Takes 'mouse' or 'human'. 
#' @param expname character Name of expreriment. E.g., "Exp1_DGEA"
#' @param cores numeric Specify number of CPU cores.
#' @param seed numeric Seed from which randomness grows.
#' @param directory character Name of working folder.
#' @param CPMfilter logical If TRUE, raw counts will be filtered based of counts per milion (CPM) cutoff
#' @param CPMcutoff numeric Cutoff for CPM filtering e.g., 2
#' @param normalise logical, character If TRUE, raw counts will be geTMM normalised. If "CQN", raw counts will be quantile normalised
#' @param orthologs logical If TRUE, maps genes to orthologs of given species and continues analysis with these orthologous genes.  
#' @param orthogenes character list List of character vector of orthologous ENRTEZID genes (mouse or human)
#' @param orthospecies character Takes 'mouse' or 'human'
#' @param batch logical If TRUE, data will undergo batch correction
#' @param dropOut logical If TRUE, outliers will be dropped based on PCA distance.
#' @param svaCor logical If TRUE, performs surrogatew variable analysis and correction.
#' @param svaFormula character Formula for SVA e.g. id ~ condition + gender
#' @param svaCovar character Covariates of interest whose influence should not be removed e.g., condition
#' @param conCovar character Potentially confounding covariates whose influence might want to be removed e.g., gender
#' @param includeZeros logical If TRUE. zero values will be nudged by 1x10^-16 to facilitate PCA. Otherwise, rows with zero values will be removed.
#' @param RsqCut numeric Desired Rsquared cutoff for scale-free fit e.g, 0.9
#' @param normFormu character Formula used for variance stabilising normalisation e.g., conditon ~ id
#' @param choosePower logical If TRUE, will attempt to choose integer exponent which results in scale-free gene expression covariance network. If initial correlation network is not scale free, will iterativelty remove genes based on mean expression cutoff. If "bicor" algorithm is chosen initially and the function fails to produce a scale free network with R^2 > 0.75, 
#' @param corFunc character The chosen correlation algorithm e.g., "bicor" for bi-correlation or "cor" for pearson's correlation. See WGCNA vignette for more details. 
#' @param meanExprStep numeric If ChoosePower is TRUE and the intial correlation network is not scale free, the percentage of mean expression cutoff at which genes are shaved off at each iteration.
#' @param netType character The type of correlation network e.g., "signed", "unsigned".
#' @param spQN logical If TRUE, performs spatial quantile normalisation on scale-free correlation network. Attempt to correlation bias of higher expressed genes. 
#' @param spQNplots logical If TRUE, plots all quality control plots from spQN
#' @param QCplots logical If TRUE, plots all quality control plots associated with analysis
#' @param wgcna logical If TRUE, performs weighted gene co-expression neteork analysis. See WGCNA vigette for more details.
#' @param cutDepth numeric Depth at which to cut heirarchal clustering tree default: 3. See WGCNA vigette for more details.
#' @param kMEmerge logical If TRUE, merges clusters based on gene module membership. See WGCNA vigette for more details.
#' @param kMEclustThres numeric The percentage genes of modules to be included for merging e.g., 50. See WGCNA vigette for more details.
#' @param kMEclustMrgPcnt numeric The percentage of genes above the threshold percentage in a given module e.g., 75. See WGCNA vigette for more details.
#' 
#' @import WGCNA
#' @import EDASeq
#' @import DESeq2
#' @import edgeR
#' @import spqn
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
#' @return List of WGCNA modules
#' 
#' @examples
#' oneNet()
#' 
#' @export
oneNet <- function(counts = counts,
                    meta = meta,
                    species = NULL,
                    directory = getwd(),
                    normalise = NULL,
                    expname =  "oneNet_exp",
                    seed = 42069,
                    cores = NULL,
                    includeZeros = FALSE,
                    CPMfilter = NULL,
                    CPMcutoff = NULL,
                    dropOut = NULL,
                    batch = NULL,
                    svaCor = NULL,
                    svaCovar = NULL,
                    conCovar = NULL,
                    svaFormula = NULL,
                    RsqCut = NULL,
                    normFormula = NULL,
                    choosePower = NULL,
                    corFunc = 'bicor',
                    meanExprStep = 0.01,
                    netType = "signed",
                    #levelOrder = NULL,
                    spQN = NULL,
                    spQNplots = FALSE,
                    QCplots = NULL,
                    wgcna = NULL,
                    cutDepth = 3,
                    kMEmerge = NULL,
                    kMEclustThresh = 50, 
                    kMEclustMrgPcnt = 50,
                    MergeCutHeight = 0.2) {

  # crayonst
  wild <- red $ bold
  yay <- green $ bold
  calc <- cyan $ bold
  prep <- silver $ bold
  fyi <- white $ bold

  # WGCNA options
  
  options(stringsAsFactors = FALSE)
  enableWGCNAThreads()
  ALLOW_WGCNA_THREADS = cores 
  
  # Match Samples
    
  #colnames(counts) <- gsub("^X","",colnames(counts),perl = T)
  
  if (is.character(counts[[1]])) {
    
    # Check if gene names are first column instead of rownmaes
    cat("Gene names in first column. Removing and assigning to rownames...\n")
    rownames(counts) <- counts[,1]
    counts <- counts[-1]
    
  }
  
  # Same number of samples?
  
  if (length(counts) < nrow(meta)) {
    
    cat("More samples in metadata than counts matrix. Removing unmatched samples.")
    
    readkey <- function()
    {
      cat("If you're cool with the above, then on your keyboard and in the console, Press [Enter] to continue")
      line <- readline()
    }
    readkey()
    meta <- meta[meta$id %in% colnames(counts),]
    
   }
    
  if (length(counts) > nrow(meta)) {
      
    cat("More samples in count matrix than metadata. Removing unmatched samples.")
    
    readkey <- function()
    {
      cat("If you're cool with the above, then on your keyboard and in the console, Press [Enter] to continue")
      line <- readline()
    }
    readkey()
    
      counts <- counts[,colnames(counts) %in% meta$id]
      
    }
    # Check if sample names are in write order
    
  if (any(!meta$id == colnames(counts))) {
      
      cat("Sample are in wrong order. Fixing...\n")
      
      meta <- mutate_all(meta,as.character)
      meta <- meta[order(meta$condition),]
      counts <- counts[,meta$id]
      }      
    # Meta as factors
    meta <- mutate_all(meta,as.factor)
    stopifnot(colnames(counts) == meta$id)
    
    cat("Sample matching complete!\n")
  
    # all genes unfiltered
    
    setwd(directory)
    if (is.null(species) | species ==  "mouse") {
      org <- "org.Mm.eg.db"   
    } else if (species ==  'human') {
      org <- "org.Hs.eg.db"   
    }
    allGenes <- bitr(rownames(counts), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org)$ENTREZID 
    write_rds(allGenes, "allGenes.rds")

    ## Filter 
    
    if (CPMfilter) {
      
      cat(yellow("Filtering data...\n"))
      
      setwd(directory)
      
      cpm <- 10^6*t(t(counts)/colSums(counts)) %>% as.data.frame()
      filtidx <- which(rowSums(cpm)/ncol(cpm) <= CPMcutoff)
      counts_filt <- counts[-c(filtidx),] %>% as.data.frame()
      write.csv(counts_filt,paste0(expname,"_unormalized",CPMcutoff,"CPM_Counts.csv"))
      
      # All genes
      
      if (is.null(species) | species ==  "mouse") {
        org <- "org.Mm.eg.db"   
      } else if (species ==  'human') {
        org <- "org.Hs.eg.db"   
      }
      allGenes <- bitr(rownames(counts_filt), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org)$ENTREZID 
      write_rds(allGenes, "allGenes.rds")
      setwd(directory)
      }
      
    if (normalise) {
    
    cat(yellow("Performing Variance Stabilsing Transformation\n"))
    meta <-  mutate_all(meta, function(x) as.factor(x)) #all as factor
    counts <- mutate_all(counts, function(x) as.integer(x)) #all as integer
    dso <- DESeqDataSetFromMatrix(counts,meta,design = normFormula) # DEseq object with predefined formula
    dso <- estimateSizeFactors(dso) #size factors
    
    # Dispersion parameters
    
    cat(blue("Optimising Dispersion Parameters...\n"))
    ## Parametric
    par <- estimateDispersions(dso, fitType = "parametric")
    parresid <-   mcols(par)$dispGeneEst - mcols(par)$dispFit
    parabsresid <- abs(parresid) %>% as.data.frame()
    parabsresid <- parabsresid[complete.cases(parabsresid),]
    parmedabs <- median(parabsresid)
    ## Local 
    loc <- estimateDispersions(dso, fitType = "local")
    locresid <-   mcols(loc)$dispGeneEst - mcols(loc)$dispFit
    locabsresid <- abs(locresid) %>% as.data.frame()
    locabsresid <- locabsresid[complete.cases(locabsresid),]
    locmedabs <- median(locabsresid)
    
    ## Choose best fit
    if (any(locmedabs < parresid)) {
      cat(cyan("Using Local Fit...\n"))
      dds <- loc
    }else{
      cat(cyan("Using Parametric Fit...\n"))
      dds <- par
    }
    
    
    vst_dds <- DESeq2::varianceStabilizingTransformation(dds,blind = FALSE) # VST
    vst_counts  <- assay(vst_dds) %>% as.data.frame()
    cat(yellow("Writing VST Normalised Counts\n"))
    
    if (CPMfilter) {
      
      counts <- vst_counts[-c(filtidx),] %>% as.data.frame()
      write.csv(counts_filt,paste0(expname,"_VST_",CPMcutoff,"CPM_Counts.csv"))
      counts <- vst_counts
      .GlobalEnv$vst_counts_filt <- vst_counts                
      
    } else {
    
    counts <- vst_counts
    .GlobalEnv$vst_counts <- vst_counts          
    write.csv(vst_counts,paste0(expname,"_VST_counts.csv"))
    cat(green("Normalisation complete!\n"))
    
    }
    
    if (QCplots) {
    
    ## Density Plots
    dnames <-  dimnames(counts)
    colors <-  sequential_hcl(ncol(counts), "Viridis")
    if (!dir.exists("QC_plots")) { # Create folder
      dir.create("QC_plots")}
    
    qcpath <- paste0(directory,"/QC_plots") # save path
    setwd(qcpath)
  
    pdf('Density_plots_VST.pdf')
    par(mfrow = c(2,1))
  
    plot(density(x = as.numeric(unlist(counts))),main = "Un-normalised",
         xlim = c(0,20000),xlab = "Exp")
    mask <-  sample(1:length(dnames[[2]]))
    ret <-  lapply(mask[2:length(mask)], function(x) {
      i = which(mask == x)
      lines(density(as.numeric(counts[,x])),col = colors[i])
    })
    
    plot(density(x = as.numeric(unlist(vst_counts))),main = "VST Normalised",
         xlim = c(0,20000),xlab = "Exp")
    mask <-  sample(1:length(dnames[[2]]))
    ret <-  lapply(mask[2:length(mask)],function(x) {
      i = which(mask == x)
      lines(density(as.numeric(counts[,x])),col = colors[i])
    })
    graphics.off()
  }
  
    cat(green("Normalisation Successful!\n"))
   
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
  
  ## PCAdist Outlier Removal
    
  if (dropOut == "PCAdist" || dropOut == TRUE) {
    
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
    write.csv(norm_osr_counts,paste0(expname,'_VST_CPM_OSR_counts'))
    
    cat(green("Outliers removed!\n"))
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
    
    cat(prep("\nDiscarding SVs correlated to variables of interest...\n"))
    
    # Drop significant correlations to variables of interest
    
    if (is.null(conCovar)) {
      conCovar <- c()
    }
    
    linp_t <- t(linp) %>% as.data.frame()
    idx_drop <- sapply(dplyr::select(linp_t,-conCovar), function(x) which(x < 0.05))
    idx_keep <- sapply(dplyr::select(linp_t,conCovar), function(x) which(x < 0.05))
    
    if (length(idx_drop) > 1) {
    
    if (length(idx_keep) >= 1){
      idx_keep2 <- unlist(idx_keep)[-c(unlist(idx_drop))]
    } else if ((length(idx_keep) < 1)) {
      idx_keep2 <- c(1:nrow(linp_t))[-c(unlist(idx_drop))]
      } else {
      idx_keep2 <- c()
        }
      }
    
    if (length(idx_keep2) <= 1 ) {
      
      cat(wild("Not enough significant surrogate variables!\n"))
      cat(green("Continuing...\n"))
      
    } else {
      
    cat(prep("Generating SV corrected count matrix...\n"))
    
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
    
     if (length(idx_keep2) > 1 ) {
    
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


    ## Choose power  
  
    if (choosePower) {
      
    cat(yellow("Picking soft-thresholding power...\n"))
      
    # Setup 
    options(stringsAsFactors = FALSE)
    enableWGCNAThreads()
    ALLOW_WGCNA_THREADS = cores 
    dropsList <- list()
    res <- list()
    beta_list <- list()
    splits <- list()
    Rcut <- RsqCut
    beta_plts <- list()
    expMeanCut = 1
    gc()   
    
    # Choose beta
    tcounts <-  t(counts) %>% as.data.frame() # transpose
    .GlobalEnv$tcounts <- tcounts
    
    beta <- pickSoftThreshold(tcounts,
                              dataIsExpr = TRUE, 
                              corFnc = corFunc, 
                              networkType = netType,
                              RsquaredCut = Rcut,
                              moreNetworkConcepts = F,
                              verbose = T) 

    # Loop
  
    beta_list[[1]] <- beta
    names(beta_list) <- "picks"
    
    # Plot df
    
    pwrvsR <- list()
    
    for (i in seq_along(beta_list)) {
      
      pwrvsR[[i]] <- data.frame(beta_list[[i]]$fitIndices) %>%
        dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq) 
      }
    
    .GlobalEnv$pwrvsR <- pwrvsR
    .GlobalEnv$beta_list <- beta_list
    
    if (max(pwrvsR[[1]]$model_fit) >=  Rcut) {
      
    cat(paste(green("\nLucky! Congratulations on your Scale-Free Network @R^2 =",max(pwrvsR[[1]]$model_fit),"!\n")))
  
    wgcna_in <- tcounts
    pwr <- beta_list$pick$powerEstimate
    write_rds(pwr,"pwr.rds")
    write_rds(wgcna_in,"wgcna_IN_matrix.rds")
    #.GlobalEnv$beta_list <- beta_list
    #.GlobalEnv$wgcna_in <- wgcna_in
    
    # Plots
    
    for (i in seq_along(pwrvsR)) {
    
    ## B vs R^2
      beta_plts[[i]] <-  ggplot(pwrvsR[[i]], aes(x = Power, y = model_fit, label = Power)) +
      # Plot the points
      geom_point() +
      # We'll put the Power labels slightly above the data points
      geom_text(nudge_y = 0.1) +
      # We will plot what WGCNA recommends as an R^2 cutoff
      geom_hline(yintercept = Rcut, col = "purple") +
      # Just in case our values are low, we want to make sure we can still see the 0.80 level
      ylim(c(min(pwrvsR[[i]]$model_fit), 1.05)) +
      # We can add more sensible labels for our axis
      xlab("Soft Threshold (power)") +
      ylab("Scale Free Topology Model Fit, signed R^2") +
      ggtitle("Scale independence") +
      theme_classic()
    }
    ## Print
    pdf("beta_vs_fit_picks.pdf")
    print(beta_plts[[1]])
    graphics.off()
    gc()   
    
    }
    
    if (!exists('wgcna_in') && max(pwrvsR[[1]]$model_fit) <  Rcut ) {  
    
    cat(yellow("Shaving by mean expression...\n"))  
    expMeanCut = 1
    
    while (max(pwrvsR[[1]]$model_fit) <  Rcut && expMeanCut >= 0.75) { 

    expMeanCut <- expMeanCut - meanExprStep # drop 5% from bottom until R^2 < 0.8
    means <-  apply(tcounts, 2, mean) %>% as.data.frame()
    dropidx <- which(means$. >= quantile(means$., expMeanCut, na.rm = T))
    drops <- tcounts[,-c(dropidx)] %>% as.data.frame()
    dropsList <- list()
    dropsList[[1]] <- drops %>% as.data.frame()
    names(dropsList) <- 'picks'
    
    beta_list <- list()
    
    for (i in names(dropsList)) {
    
    cat(magenta(paste("Trying r^2 =", Rcut ,"\n")))
      
    cat(cyan(paste("Trying Mean Expression Quantile Cutoff = ", expMeanCut ,"\n")))
    
    beta_list[[i]] <- pickSoftThreshold(dropsList[[i]], # run on both picks and rejects
                                        dataIsExpr = TRUE, 
                                        corFnc = corFunc, # bicor looks more legit
                                        networkType = netType, # this is the sign of the correlation
                                        RsquaredCut =  Rcut, # scale-free model cut-off
                                        moreNetworkConcepts = F,
                                        verbose = T)    
    
    
    pwrvsR <- list()
    
    for (i in seq_along(beta_list)) {
      
      pwrvsR[[i]] <- data.frame(beta_list[[i]]$fitIndices) %>%
        dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq) }
    
    .GlobalEnv$pwrvsR <- pwrvsR
    .GlobalEnv$expMeanCut <- expMeanCut
    .GlobalEnv$Rcut <-  Rcut 
    .GlobalEnv$dropsList <- dropsList
    
    if (expMeanCut == 0.75 && max(pwrvsR[[1]]$model_fit) < Rcut ) {
      
      expMeanCut <- 1
      Rcut  <- Rcut - 0.05
    
      }
    
    if (expMeanCut == 0.75 && Rcut  == 0.85) {
      
      expMeanCut <- 1
      Rcut  <- Rcut - 0.05
          
          }
        }
      }
    }
      
      
    if (max(pwrvsR[[1]]$model_fit) >= Rcut && expMeanCut >= 0.75) {
      
    cat(paste(green("\nCongratulations on your Scale-Free Network @R^2 = ", Rcut, "!\n")))
    
    wgcna_in <- dropsList$picks %>% as.data.frame()
    
    .GlobalEnv$beta_list <- beta_list
    .GlobalEnv$wgcna_in <- wgcna_in
    .GlobalEnv$dropsList <- dropsList
    pwr <- beta_list$picks$powerEstimate
    write_rds(pwr,"pwr.rds")
    write_rds(wgcna_in,"wgcna_IN_matrix.rds")
    
    for (i in seq_along(beta_list)) {
  
    ## B vs R^2
      beta_plts[[i]] <-  ggplot(pwrvsR[[i]], aes(x = Power, y = model_fit, label = Power)) +
      # Plot the points
      geom_point() +
      # We'll put the Power labels slightly above the data points
      geom_text(nudge_y = 0.1) +
      # We will plot what WGCNA recommends as an R^2 cutoff
      geom_hline(yintercept = Rcut, col = "purple") +
      # Just in case our values are low, we want to make sure we can still see the 0.80 level
      ylim(c(min(pwrvsR[[i]]$model_fit), 1.05)) +
      # We can add more sensible labels for our axis
      xlab("Soft Threshold (power)") +
      ylab("Scale Free Topology Model Fit, signed R^2") +
      ggtitle("Scale independence") +
        theme_classic()
      }
      ## Print
    
      pdf("beta_vs_fit_picks.pdf")
      print(beta_plts[[1]])
      graphics.off()
      
    } 
    
    ## Try Pearson's  
    
    if (!exists('wgcna_in') && corFunc == 'bicor' && max(pwrvsR[[1]]$model_fit) < Rcut) { 
  
      cat(yellow('Trying Pearsons correlation...\n'))
        
      corFunc = 'cor'
      Rcut = 0.9
      expMeanCut = 1
      
      if (max(pwrvsR[[1]]$model_fit) <  Rcut ) {  
    
      while (max(pwrvsR[[1]]$model_fit) <  Rcut  && expMeanCut >= 0.75) { 
    
      expMeanCut <- expMeanCut - meanExprStep # drop 5% from bottom until R^2 < 0.8
      
      means <-  apply(tcounts, 2, mean) %>% as.data.frame()
      dropidx <- which(means$. >= quantile(means$., expMeanCut, na.rm = T))
      drops <- counts[,-c(dropidx)] %>% as.data.frame()
      dropsList <- list()
      dropsList[[1]] <- drops
      names(dropsList) <- 'picks'
      #splits <- split(counts,counts$mean >= quantile(counts$mean, expMeanCut, na.rm = T)) # split by cutoff
      
      
      for (i in names(dropsList)) {
        
        cat(magenta(paste("Trying r^2 =", Rcut ,"\n")))
        
        cat(cyan(paste("Trying Mean Expression Quantile Cutoff = ", expMeanCut ,"\n")))
        
        beta_list[[i]] <- pickSoftThreshold(dropsList[[i]], # run on both picks and rejects
                                            dataIsExpr = TRUE, 
                                            corFnc = corFunc, # bicor looks more legit
                                            networkType = netType, # this is the sign of the correlation
                                            RsquaredCut =  Rcut, # scale-free model cut-off
                                            moreNetworkConcepts = F,
                                            verbose = T)    
        .GlobalEnv$expMeanCut <- expMeanCut
        .GlobalEnv$Rcut <-  Rcut 
        
        pwrvsR <- list()
        for (i in seq_along(beta_list)) {
          
          pwrvsR[[i]] <- data.frame(beta_list[[i]]$fitIndices) %>%
            dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq) }
        
        .GlobalEnv$pwrvsR <- pwrvsR
        
          
      if (expMeanCut == 0.75 && max(pwrvsR[[1]]$model_fit) < Rcut ) {
        
        expMeanCut <- 1
        Rcut  <- Rcut - 0.05
        
      }
      
      if (expMeanCut == 0.75 && Rcut  == 0.85) {
        
        expMeanCut <- 1
        Rcut  <- Rcut - 0.05  ## add rest        
        
                }
              }
            }
          }
        }
      
      if (!exists('wgcna_in') &&  max(pwrvsR[[1]]$model_fit) >= Rcut && expMeanCut >= 0.75 && corFunc == "cor") {
        
        cat(paste(green("\nCongratulations on your Scale-Free Network @R^2 = ", Rcut, "!\n")))
        
        wgcna_in <- dropsList$picks
        .GlobalEnv$beta_list <- beta_list
        .GlobalEnv$dropsList <- dropsList
        pwr <- beta_list$picks$powerEstimate
        write_rds(pwr,"pwr.rds")
        
        for (i in seq_along(beta_list)) {
          
        # B vs R^2
        beta_plts[[i]] <-  ggplot(pwrvsR[[i]], aes(x = Power, y = model_fit, label = Power)) +
        # Plot the points
        geom_point() +
        # We'll put the Power labels slightly above the data points
        geom_text(nudge_y = 0.1) +
        # We will plot what WGCNA recommends as an R^2 cutoff
        geom_hline(yintercept = Rcut, col = "purple") +
        # Just in case our values are low, we want to make sure we can still see the 0.80 level
        ylim(c(min(pwrvsR[[i]]$model_fit), 1.05)) +
        # We can add more sensible labels for our axis
        xlab("Soft Threshold (power)") +
        ylab("Scale Free Topology Model Fit, signed R^2") +
        ggtitle("Scale independence") +
        theme_classic()
        }
        ## Print
        
        pdf("beta_vs_fit_picks.pdf")
        print(beta_plts[[1]])
        graphics.off()
        
        pdf("beta_vs_fit_rejects.pdf")
        print(beta_plts[[2]])
        
        graphics.off()
          
        } 
      }

      ## Spatial Quantile Normalisation
    
      if (spQN) {
      
      setwd(directory)
        
      cat(prep("Preparing for Quantile Normalisation underway...\n"))
      
      if (!exists('wgcna_in') && !isEmpty(list.files(pattern = "wgcna_IN_matrix.rds",recursive = T))) {
        
        wgcna_in <- read_rds("wgcna_IN_matrix.rds")
        
      } 
        
      # Correlation matrix
      
      cat(blue("\nCalculating correlations...\n"))  
      
      # Order high to low
      
      wgcna_in <- wgcna_in %>% as.data.frame()
      avs <- rowMeans(t(wgcna_in)) %>% as.numeric()
      names(avs) <- colnames(wgcna_in)
      avs <- avs[order(avs,decreasing = T)]
      tmp <- avs %>% as.data.frame()
      wgcna_in <- wgcna_in[,rownames(tmp)] %>% as.data.frame()
      stopifnot(colnames(wgcna_in) == rownames(tmp))
      rm(tmp)
      
      if (corFunc == "bicor") {
      
        bi <- bicor(wgcna_in,
                    verbose = T) # absolute
    
      .GlobalEnv$bi <- bi
      
      } else {
        
        bi <- WGCNA::cor(wgcna_in, verbose = T)
        
      }
      cat(blue("\nCalculating IQR\n"))
        
      # IQR
      
      IQR_list <- get_IQR_condition_exp(bi,avs)
      .GlobalEnv$IQR_list <- IQR_list
      
      # Optional Plots
      if (spQNplots) {
        
      cat(green("\nQC Plots: this will take a while...\n"))
      
      if (!dir.exists("QC_plots")) { # Create folder
        dir.create("QC_plots")}
      qcpath <- paste0(getwd(),"/QC_plots") # save path 
      
      pdf(paste0(qcpath,"/IQR_vs_minaverage_spQN_unnormalised.pdf"))
      #pdf("IQR_vs_minaverage_spQN_unnormalised.pdf")  
      plot(rep(IQR_list$grp_mean, times = 1:10), # compare submatrices
           IQR_unlist,
           xlab = "min(average(pseudocounts))", ylab = "IQR", cex.lab = 1.5, cex.axis = 1.2, col = qualitative_hcl(1, "Dark2"))
      dev.off()
      pdf(paste0(qcpath,"/IQR_vs_minaverage_spQN_unnormalised.pdf"))
      #pdf("IQR_vs_minaverage_spQN_unnormalised.pdf")    
      par(mfrow = c(3,3)) 
      for(j in c(1:8,10)) {
        qqplot_condition_exp(bi, avs, j, j) } #comapare distribution to rest
      graphics.off()
      IQR_unlist <- unlist(lapply(1:10, function(ii) IQR_list$IQR_cor_mat[ii, ii:10])) # arbitrary 10-split
      pdf(paste0(qcpath,"/IQR_vs_minaverage_spQN_normalised.pdf"))
      plot(rep(IQR_spqn_list$grp_mean, times = 1:10),IQR_unlist,xlab = "min(average(VST Counts))", 
           ylab = "IQR", cex.lab = 1.5, cex.axis = 1.2, col = qualitative_hcl(1,"Dark2"))
      dev.off()
      }
      # Mandatory Plots
      
      cat(cyan("Plotting spQN QC plots...\n"))
      
      pdf("signal_Cor_plots.pdf")  
      
      print(plot_signal_condition_exp(bi, avs, signal = 0))
      print(plot_signal_condition_exp(bi, avs, signal = 0.001))
      
      graphics.off()
      
      #IQR
      
      pdf("Unormalised_Cor_Bin_plot.pdf")  
      
      print(plot_IQR_condition_exp(IQR_list)) #plot it
      
      graphics.off()
      
      ## Run Normalisation
      
      cat(yellow("\nspQN Underway...\n"))
      
      sz <- round(ncol(bi)/60*2) # Fixed nGrp at 60, product of it and grp sizw must be at least ncol/nrow matrix. Made it 2x
      spqn_mat <- normalize_correlation(cor_mat = as.matrix(bi), ave_exp = avs, ngrp = 60, size_grp = sz, ref_grp = 58)
      cat(green("\nspQN Complete!...\n"))
      write_rds(spqn_mat, "spqn_matrix.rds")
      
      # IQR 
      IQR_spqn_list <- get_IQR_condition_exp(spqn_mat, avs)
      .GlobalEnv$IQR_spqn_list <- IQR_spqn_list
      IQR_unlist <- unlist(lapply(1:10, function(ii) IQR_spqn_list$IQR_cor_mat[ii, ii:10]))
      
      # Optiional Plots
      
      if (spQNplots) {
        
        cat(cyan("Plotting QC plots of normalised data\n"))
        pdf(paste0(qcpath,"/IQR_vs_minaverage_spQN_normalised.pdf"))
        plot(rep(IQR_spqn_list$grp_mean, times = 1:10),IQR_unlist,xlab = "min(average(log2RPKM))", 
             ylab = "IQR", cex.lab = 1.5, cex.axis = 1.2, col = qualitative_hcl(1,"Dark2"))
        graphics.off()
        pdf(paste0(qcpath,"/qqplots_spQN_normalised.pdf"))
        par(mfrow = c(3,3))
        for(j in c(1:8,10)) {
          qqplot_condition_exp(spqn_mat, avs, j, j)
        }
        graphics.off()
      }
    
      # Normed Plots    
      
      pdf("spQN_Normalised_bckrnd_sig.pdf")
      
      print(plot_signal_condition_exp(spqn_mat, avs, signal = 0))
      print(plot_signal_condition_exp(spqn_mat, avs, signal = 0.001))
      
      graphics.off()
      
      pdf("spQN_normalised_IQR.pdf")
      
      print(plot_IQR_condition_exp(IQR_spqn_list))
      
      graphics.off()
      
      # write
      
      write_rds(spqn_mat,"spqn_matrix.rds")
      
      if(!exists("pwr")){
        
        pwr <-   read_rds("pwr.rds")
      
        } else {
        
          stop("Please run choosePower for input to adjacency function!")
      
          }
      
      cat(blue("Making adjacency matrix\n"))
      
      # Adjacency Matrices
      
      adj <- adjacency.fromSimilarity(as.matrix(spqn_mat),
                                      type = netType,
                                      power = pwr)
      
      cat(green("Adjacency matrix complete!\n"))
      write_rds(adj, "adjacency_matrix.rds")
      .GlobalEnv$adj <- adj
      
      }
    
      ## Adjacency 
    
      if (is.null(spQN || spQN == FALSE)) {
      
      if (!isEmpty(list.files(pattern = "adjacency_matrix.rds"))) { 
          adj <- read_rds("adjacency_matrix.rds") 
      
        } 
      
      if (!exists("wgcna_in") && !isEmpty(list.files(pattern = "wgcna_in.rds"))) {
        stop("Input matrix doesn't exist!\n Run 'choosePower'.\n")
      } else {
        
        wgcna_in <- read_rds("wgcna_in.rds")
      }    
      
      gc()   
      
      cat(blue("Making adjacency matrix\n"))
      
      adj <- adjacency(as.matrix(wgcna_in),
                       type = netType,
                       power = pwr,
                       corFnc = corFunc)
      
  
      cat(green("Adjacency matrix complete!\n"))
      .GlobalEnv$adj <- adj
      
      write_rds(adj, "adjacency_matrix.rds")
      gc()   
      
      }
    
      ## TOM
    
      gc()   
      cat(yellow("\nCalculating Topological Overlap Measure...\n"))
      
      if (!exists("adj")) {
        if(!isEmpty(list.files(pattern = "adjacency_matrix.rds"))){
          adj <- read_rds("adjacency_matrix.rds")
        } else {
          stop("Adjacency matrix not found!")
        }
      }
      
      if (wgcna) {
      
      if(!exists("wgcna_IN_matrix.rds")){
        wgcna_in <- read_rds("wgcna_IN_matrix.rds")
      }
      tmp <- list.files(pattern = "TOM.rds", recursive = T)
      if(!isEmpty(tmp)){
        TOM <- read_rds(tmp)
      } else {
        
        TOM <-  TOMsimilarity(as.matrix(adj),
                              TOMType = "signed",
                              TOMDenom = "mean",
                              suppressNegativeTOM = F,
                              suppressTOMForZeroAdjacencies = F,
                              verbose = T) #topo overlap
        gc()   
        write_rds(TOM, "TOM_spqn.rds")
      }
      # Distance matrix
      
      dissTOM <- 1 - TOM
      
      cat(green("Toplogical overlap measure calculation complete!\n"))
      
      # HClust distances - the ones which are different the same will aggregate
      
      cat(yellow("Clustering TOM distance matrix...\n"))
      
      clustTree = hclust(as.dist(dissTOM), method = "average") #standard hclust
      
      cat(green("Clustering Complete!\n"))
      
      # Cut dat tree
      cat(yellow("Cutting branches using dynamic tree cut...\n"))
      
      dynamicMods = cutreeDynamic(dendro = clustTree, 
                                  distM = dissTOM,
                                  deepSplit = cutDepth, 
                                  pamRespectsDendro = FALSE, 
                                  minClusterSize = 20,
                                  method = "hybrid",
                                  verbose = T)
      
      # Assign colours
      
      dynamicColors = labels2colors(dynamicMods)
      
      cat(green("Clustering Complete!\n") %+% blue("Number of Modules:",(length(unique(dynamicColors))),"\n"))
      
      ## Plots
      
      cat(yellow("Plotting dendograms...\n"))
      
      pdf("dissTOM_hclust.pdf")
      
      p <- plot(clustTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
      print(p)
      graphics.off()
      
      pdf("dissTOM_dynamic_treecut.pdf")
      
      plotDendroAndColors(clustTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                          main = "Gene dendrogram and module colors")
      graphics.off()
      
      cat(green("Dendograms exported!\n"))
      gc()   
      
      ## Eigengenes
      
      cat(yellow("Clustering Eigengenes...\n"))
      
      MEList <-  moduleEigengenes(wgcna_in, colors = dynamicColors)
      MEs <-  MEList$eigengenes
      
      # Calculate dissimilarity and Cluster
      
      MEDiss = 1 - cor(MEs,method = "pearson")
      METree = hclust(as.dist(MEDiss), method = "average");
      
      ## Plot
    
      pdf("Eigengene_clustering.pdf")
      
      p <- plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
      print(p)
      graphics.off()
      
      cat(green("Clustering complete!"))
      
      ## Merge close    
      
      cat(prep("Merging close modules..."))
      
      merge = mergeCloseModules(wgcna_in, 
                                dynamicColors, 
                                cutHeight = MergeCutHeight, 
                                verbose = 4,
                                MEs = MEs)    
      
      mergeColors = merge$colors #colours
      mergedMEs = merge$newMEs #MEs
      newMEs <- moduleEigengenes(wgcna_in,colors = mergeColors)
      mergedMEs <- newMEs$eigengenes 
      colorOrder = c("grey", standardColors(100))
      moduleLabels = match(mergeColors, colorOrder)-1
      .GlobalEnv$mergeColors <- mergeColors
    
    
    ## kME merge
    
      if (kMEmerge) {
        
        
        #kMEclustThresh = 20
        #kMEclustMrgPcnt = 90
        
        kmerge <- moduleMergeUsingKME(
          wgcna_in, mergeColors, ME = mergedMEs,
          threshPercent = kMEclustThresh, 
          mergePercent = kMEclustMrgPcnt,
          reassignModules = TRUE,
          convertGrey = TRUE,
          omitColors = "grey",
          reassignScale = 1,
          threshNumber = NULL)
    
    # New cols and eigens  
        
        mergeCols <- kmerge$moduleColors
        newMEs <- moduleEigengenes(wgcna_in,colors = mergeCols)
        clustME <- newMEs$eigengenes
        colorOrder = c("grey", standardColors(100))
        moduleLabels = match(mergeColors, colorOrder)-1
        cat(green("Number of Modules:",length(unique(mergeCols)),"\n"))
        
        .GlobalEnv$mergeCols <- mergeCols
        .GlobalEnv$newMEs <- newMEs
        .GlobalEnv$clustME <- clustME
        .GlobalEnv$colorOrder <- colorOrder
        .GlobalEnv$moduleLabels <- moduleLabels
        
        cat(green("Merging complete!\n"))
      }
      
      ## Trait Cors  
      
      cat(blue("Correlating traits and making summary tables...\n"))
      
      nSamples = nrow(clustME)
      
      tmp <- list.files(pattern = "meta_osr.csv")
      if (!isEmpty(tmp)) {
        meta <- read.csv(tmp, row.names = 1)
      } 
      
      # Make sure all samples in order
      
      dat <- as.data.frame(wgcna_in)
      meta <- meta[rownames(dat) %in% meta$id,]
      meta <- mutate_all(meta,as.character)
      tmeta <- t(meta) %>% as.data.frame()
      tdat <- t(dat) %>% as.data.frame()
      colnames(tmeta) <-  meta$id
      tmeta <- tmeta[,colnames(tdat)]
      meta <- t(tmeta) %>% as.data.frame()
      stopifnot(rownames(dat) == meta$id)
      
      
      bintraits <- binarizeCategoricalVariable(meta$condition)
      
      modTraitCor = cor(clustME, bintraits, use = "p")
      
      modTraitP = corPvalueFisher(modTraitCor, nSamples)
      
      #plot    
      
      textMatrix = paste(signif(modTraitCor, 2), "\n(",
                         signif(modTraitP, 1), ")", sep = "")
      dim(textMatrix) = dim(modTraitCor)
      par(mar = c(10, 5, 3, 3))
    
      #.GlobalEnv$bintraits <- bintraits
      #.GlobalEnv$modTraitCor <- modTraitCor
      #.GlobalEnv$modTraitP <- modTraitP
      #.GlobalEnv$textMatrix <- textMatrix
      
      pdf("module_trait_corrplot.pdf")
      
      par(mar = c(5.1, 9.1, 4.1, 2.1))
      
      labeledHeatmap(Matrix = modTraitCor, 
                     xLabels = names(as.data.frame(bintraits)),
                     yLabels = names(clustME), 
                     ySymbols = names(clustME),
                     colorLabels = mergeCols,
                     colors = blueWhiteRed(100),
                     textMatrix=textMatrix,
                     setStdMargins = FALSE, 
                     cex.text = 0.5,
                     zlim = c(-1,1),
                     #textAdj = 2,
                     main = paste("Module-trait relationships"))
      
      graphics.off()
      
      
      ## Extract Modules and info
      
      tophubs <- chooseTopHubInEachModule(
        wgcna_in,
        mergeCols,
        omitColors = "grey",
        power = pwr,
        type = netType)
      write.csv(tophubs, "tophub_per_mods.csv")
      
      # Write kME summary 
      
      # Split into modules
      
      color_counts <- t(wgcna_in) %>% as.data.frame()
      color_counts$modules <- mergeCols # Add color labels
      mod_list <- split(data.frame(ENSEMBL=rownames(color_counts)), color_counts$modules) # Extract gene
      write_rds(mod_list, "modlist.rds")
      
      geneModuleMembership <- as.data.frame(cor(wgcna_in, clustME, use = "p")) # kME
      MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) # pval of Km
      colnames(geneModuleMembership) <- gsub("ME","kME",colnames(geneModuleMembership)) # change colnames
      colnames(MMPvalue) <- gsub("kME","pval.kME",colnames(MMPvalue)) # change colname
      kmedf <- cbind(color_counts$modules,geneModuleMembership,MMPvalue)
      write.csv(kmedf,"kME_summary.csv")
      
      cat(green("Export Complete!\n"))
     
      # convert to ENTREZ

      cat(prep("Converting ENSEMBL to ENTREZ...\n"))

      allgenes <- readRDS("allgenes.rds")
      list <- lapply(allgenes, function(x) bitr(x$ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID",
                                               OrgDb = org)) # ensembl to entre

      for (i in seq_along(mod_list)) {
        colnames(mod_list[[i]])[1] <- "ENSEMBL"
      }

      mod_list <- lapply(mod_list, function(x)
        lapply(x, function(y) as.character(y))
        %>% as.data.frame())

      # Remove those that didn't map

      join <- list()
      
      for (i in seq_along(mod_list)) {
        join[[i]] <- inner_join(list[[i]],mod_list[[i]], by = "ENSEMBL")
      }
      join <-  join[sapply(join, function(x) dim(x)[1]) > 0] #remove empty
      #.GlobalEnv$join <- join
      write_rds(join,"entrez_mods.rds")
      }
    }
      






