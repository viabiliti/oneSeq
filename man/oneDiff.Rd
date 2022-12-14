% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneDiff.R
\name{oneDiff}
\alias{oneDiff}
\title{Differential gene expressionn analysis pipeline}
\usage{
oneDiff(
  counts = counts,
  meta = meta,
  species = NULL,
  expname = NULL,
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
  ppiGSEA = NULL
)
}
\arguments{
\item{counts}{data frame Raw bulk mRNA seq counts genes as rows, columns as samples.}

\item{meta}{data frame Metadata table with samples as rows and variables as columns. Compulsory columns include 'id' with sample names corresponding to those in counts, 'condition' with relevant tests and controls, 'priority' with integers starting from 0 to 1 for indicating preferred direction of pairwise comparison and 'class' indicating whether the sample is a control or disease/test.}

\item{species}{character Takes 'mouse' or 'human'.}

\item{expname}{character Name of expreriment. E.g., "Exp1_DGEA"}

\item{normalise}{logical, character If TRUE, raw counts will be geTMM normalised. If "CQN", raw counts will be quantile normalised}

\item{CPMfilter}{logical If TRUE, raw counts will be filtered based of counts per milion (CPM) cutoff}

\item{CPMcutoff}{numeric Cutoff for CPM filtering e.g., 2}

\item{orthologs}{logical If TRUE, maps genes to orthologs of given species and continues analysis with these orthologous genes.}

\item{dropOut}{logical If TRUE, outliers will be dropped based on PCA distance.}

\item{svaFormula}{character Formula for SVA e.g. id ~ condition + gender}

\item{svaCor}{logical If TRUE, performs surrogatew variable analysis and correction.}

\item{svaCovar}{character Covariates of interest whose influence should not be removed e.g., condition}

\item{conCovar}{character Potentially confounding covariates whose influence might want to be removed e.g., gender}

\item{batch}{logical If TRUE, data will undergo batch correction}

\item{includeZeros}{logical If TRUE. zero values will be nudged by 1x10^-16 to facilitate PCA. Otherwise, rows with zero values will be removed.}

\item{writeOut}{logical If TRUE, tables generated during analysis will be written to csv.}

\item{QCplots}{logical IF TRUE,  various quality control plots will be generated.}

\item{DE}{logical If TRUE, performs noisSeq differential gene expression analysis}

\item{GSEA}{logical If TRUE, performs setRank network based pathway analysis of differentially expressed genes.}

\item{orthogenes}{character list List of character vector of orthologous ENRTEZID genes (mouse or human)}

\item{orthospecies}{character Takes 'mouse' or 'human'}

\item{clusterGS}{logical perform hierarchal clustering on GSEA results?}

\item{CTDbase}{logical Include Chemical and Toxigenomics Data base in GSEA?}

\item{chemDis}{logical if CTDbase, include chemical-gene-disease interactions in GSEA?}

\item{chemIX}{logical if CTDbase, include chemical-gene interactions in GSEA?}

\item{msigCat}{list of character vectors list of desired msigDB sub-categories e.g., c("H", "C1")}

\item{msigDB}{logical nclude msigDB in GSEA?}

\item{cytoscape}{logical export networks to cytoscape?}

\item{seed}{numeric Seed from which randomness grows.}

\item{directory}{character Name of working folder.}

\item{cores}{numeric Specify number of CPU cores.}

\item{ppiGSEA}{logical map GSEA results to HIPPIE protein-protein interaction network?}
}
\description{
DGEA leveraging the wonderful noiSeq package
}
\details{
ceeb
}
\examples{
oneDiff()

}
