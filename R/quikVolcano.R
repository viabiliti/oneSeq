#' Generate volcano plots from differentially expressed genes table
#' 
#' Simple function for plotting volcanoes in base R
#' 
#' ceeb
#'
#' @param DEGlist list of data frames output of lazyDGE function or DE = 'noiseq' of lazySeq function
#' @param id character, name of symbol column
#' @param Log2FC character, name of log fold change column
#' @param pval character, name of pval column
#' @param cutoff p-value cutoff
#' @param expname name of experiment
#' 
#' @examples
#' quikVolcano()
#' 
#' @return pdf with volcano plots
#' 
#' @export
quikVolcano <- function(DEGList = NULL,
                        id = NULL,
                        Log2FC = NULL,
                        pval = NULL,
                        cutoff = 0.05,
                        expname = NULL) {
  vol <- list()
  for (i in seq_along(DEGList)) {
  vol[[i]] <- data.frame(log2FC = DEGList[[i]]$log2FC,adjusted.pval = DEGList[[i]]$pval,
                         name = DEGList[[i]][id],stringsAsFactors = F) #df
  
  vol[[i]]$log2FC <-  as.numeric(vol[[i]]$log2FC)
  vol[[i]]$adjusted.pval <-  as.numeric(vol[[i]]$adjusted.pval)
  vol[[i]] <- vol[[i]][order(abs(vol[[i]]$log2FC),decreasing = TRUE),] #order largest to smallest
  vol[[i]] <- vol[[i]][vol[[i]]$log2FC != 0,] #remove 0 foldchange genes
  vol[[i]] <- drop_empty_row(vol[[i]]) #remove any empty rows
}
names(vol) <- names(DEGList) # rename list elements

# Create directory
if (!dir.exists("Volcanos")) {
  dir.create("Volcanos") }
setwd(paste0(getwd(),"/Volcanos"))

# Plot
pdf(paste0(expname,"_volcanos.pdf"))

for (i in seq_along(vol)) {
  
  with(vol[[i]], plot(log2FC, -log10(adjusted.pval), pch = 20, main = paste(mgsub(x = names(vol[i]),c("vs"), c(" vs "))))) 
  abline(h = -log10(cutoff),col = "red",lty = 2) #significane boundaries
  with(subset(vol[[i]], adjusted.pval >=cutoff & abs(log2FC)>=0), points(log2FC, -log10(adjusted.pval), 
                                                                          pch = 20, col = "#252525"))
  with(subset(vol[[i]], log2FC > 0 & adjusted.pval <=  cutoff), points(log2FC, -log10(adjusted.pval), 
                                                                        pch = 20, col = "#c62828"))
  with(subset(vol[[i]], log2FC < 0 & adjusted.pval <= cutoff), points(log2FC, -log10(adjusted.pval), 
                                                                       pch = 20, col = "#1565c0"))
  with(subset(vol[[i]], adjusted.pval <= 10e-12 & abs(log2FC) >= 6 | abs(log2FC) >= 8), textxy(log2FC, -log10(adjusted.pval), labs = desc,
                                                                                               cex = .5))
}
graphics.off()
}