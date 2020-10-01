###########################
### Tuesday, 10th Aug 2020
###   Sam Old
###
### Data Preparation Script
###   for OG4953 Shiny App
###
###########################

## platform       x86_64-apple-darwin17.0
## arch           x86_64
## os             darwin17.0
## system         x86_64, darwin17.0
## major          4
## minor          0.2
## svn rev        78730
## language       R
## version.string R version 4.0.2 (2020-06-22)
## nickname       Taking Off Again


# Load packages and data --------------------------------------------------

setwd("/Users/siold/Proj/RNA-Shiny/Shiny_folder/ProjName/")
library(biomaRt)    # v2.44.1
library(DESeq2)     # v1.28.1
library(magrittr)   # v1.5
library(tidyverse)  # v1.3.0

raw_counts <- read.csv("https://raw.githubusercontent.com/fronchese/Mayer_et_al_2020/main/data/RNASeq_raw_counts_gene_STAR_FR_2019-11-01.csv",
                       row.names = 1, header = T)
metadata <- read.csv("https://raw.githubusercontent.com/fronchese/Mayer_et_al_2020/main/data/metadata.csv",
                     row.names = 1, header = T)


# Reformat metadata and remove sample that failed QC ----------------------

removed <- "14"
metadata <- metadata %>%
  dplyr::filter(!samfileID %in% removed)
metadata$Sample.name <- gsub("-", ".", metadata$Sample.name)
metadata$Marker <- as.factor(gsub("-", "_", metadata$Marker))

counts <- raw_counts %>%
  arrange(rownames(.)) %>%
  dplyr::select(metadata$Sample.name)


# Perform DESeq2 ----------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = metadata,
                              ~ Marker) %>%
  DESeq()

B6m24p  <- c("B61.m.24p",  "B62.m.24p",  "B63.m.24p")
B6m24n  <- c("B61.m.24n",                "B63.m.24n" )
B6s103p <- c("B61.s.103p", "B62.s.103p", "B63.s.103p")
B6s103n <- c("B61.s.103n", "B62.s.103n", "B63.s.103n")
B6e11bp <- c("B61.e.11bp", "B62.e.11bp", "B63.e.11bp")
B6eTN   <- c("B61.e.TN",   "B62.e.TN",   "B63.e.TN")
KOm24p  <- c("KO1.m.24p",  "KO2.m.24p",  "KO3.m.24p")
KOm24n  <- c("KO1.m.24n",  "KO2.m.24n",  "KO3.m.24n" )
KOs103p <- c("KO1.s.103p", "KO2.s.103p", "KO3.s.103p")
KOs103n <- c("KO1.s.103n", "KO2.s.103n", "KO3.s.103n")
KOe11bp <- c("KO1.e.11bp", "KO2.e.11bp", "KO3.e.11bp")
KOeTN   <- c("KO1.e.TN",   "KO2.e.TN",   "KO3.e.TN")

comparisons <- c(
  "KOm24pvsB6m24p",
  "KOm24nvsB6m24n",
  "KOs103pvsB6s103p",
  "KOs103nvsB6s103n",
  "KOe11bpvsB6e11bp",
  "KOeTNvsB6eTN",
  
  "B6eTNvsB6e11bp",
  "KOeTNvsKOe11bp",
  "B6m24nvsB6m24p",
  "KOm24nvsKOm24p",
  "B6s103nvsB6s103p",
  "KOs103nvsKOs103p",
  
  "KOe11bpvsB6eTN",
  "KOeTNvsB6e11bp",
  
  "B6eTNvsB6m24n",
  "B6eTNvsB6m24p",
  "B6eTNvsB6s103n",
  "B6eTNvsB6s103p",
  "B6e11bpvsB6m24n",
  "B6e11bpvsB6m24p",
  "B6e11bpvsB6s103n",
  "B6e11bpvsB6s103p",
  
  "B6m24pvsB6s103p",
  "B6m24pvsB6s103n",
  "B6m24nvsB6s103p",
  "B6m24nvsB6s103n"
)
length(unique(comparisons)) == length(comparisons) # should print TRUE

calcComp <- function(comparisons = NULL, comp_dds = NULL) {
  
  # Objectives are to:
  #   1) perform DESeq2 contrasts
  #   2) return adjusted p-values
  #   3) return log2 fold change
  #   4) return normalised counts
  
  require(DESeq2)
  
  if (!length(unique(comparisons)) == length(comparisons)) {
    return("Comparisons are not unique.")
  }
  
  all_pval <- all_l2fc <- tibble(gene_name = rownames(counts))
  
  for (i in 1:length(comparisons)) {
    ss <- unlist(strsplit(comparisons[i], "vs"))
    base = ss[1]
    comp = ss[2]
    
    cat("Comparison", base, "vs", comp, "[", i, "/", length(comparisons), "]\n")
    
    out <- comp_dds %>%
      DESeq2::results(contrast = c("Marker",
                                   as.character(unique(metadata[get(comp), "Marker"])),
                                   as.character(unique(metadata[get(base), "Marker"])))) %>%
      as.data.frame() %>%
      tibble::add_column(gene_name = rownames(.)) %>%
      tibble::as_tibble() %>%
      dplyr::select(gene_name, log2.FC = log2FoldChange, adj.P.val = padj) %>%
      tidyr::replace_na(list(log2.FC = 0, adj.P.val = 1)) %>%
      dplyr::mutate(adj.P.val = replace(adj.P.val, adj.P.val < 1e-16, 1e-16))
    
    all_pval <- all_pval %>%
      dplyr::left_join(out %>%
                         dplyr::select(gene_name, !!comparisons[i] := adj.P.val), by = "gene_name")
    all_l2fc <- all_l2fc %>%
      dplyr::left_join(out %>%
                         dplyr::select(gene_name, !!comparisons[i] := log2.FC), by = "gene_name")
  }
  return(list(all_pval, all_l2fc))
}

out_list <- calcComp(comparisons, comp_dds = dds)
all_pval <- out_list[[1]]
all_l2fc <- out_list[[2]]

all_normalised <- dds %>%
  counts(normalized = TRUE) %>%
  tibble::as_tibble() %>%
  tibble::add_column(gene_name = rownames(dds)) %>%
  dplyr::select(gene_name, metadata$Sample.name) %>%
  dplyr::arrange(gene_name)


# Create VSTpk data table -------------------------------------------------

mart <- useEnsembl(biomart = "ensembl",
                   dataset = "mmusculus_gene_ensembl",
                   version = 91)
transcripts <- getBM(attributes = c("external_gene_name", "ensembl_gene_id_version", "transcript_length", "description", "gene_biotype"),
                     filters = "external_gene_name",
                     values = rownames(counts),
                     bmHeader = FALSE,
                     mart = mart)

# Check correct version of ensembl database was used
table(rownames(counts) %in% transcripts$external_gene_name)

# Keep only the longest transcript for VSTpk gene normalisation
longest_transcripts <- transcripts %>%
  dplyr::group_by(external_gene_name) %>%
  dplyr::slice(which.max(transcript_length)) %>%
  dplyr::arrange(external_gene_name)

# Prepare VST Matrix
dds.vst <- varianceStabilizingTransformation(dds, fitType = "local")
dds.vst.mat <- assays(dds.vst)[[1]];

# log rules to calculate VSTpk normalisation. log1p is used to estimate log2 while accounting for zero counts
vst.logklengths <- log1p(longest_transcripts[,"transcript_length"]/1000)/log1p(2);
rownames(vst.logklengths) <- longest_transcripts$external_gene_name

# Create and fill final VSTpk matrix
dds.vstpk.mat <- matrix( data = NA, nrow = nrow(dds.vst.mat), ncol = ncol(dds.vst.mat), dimnames = dimnames(dds.vst.mat))
for (i in rownames(dds.vstpk.mat)) {
  dds.vstpk.mat[i,] <- dds.vst.mat[i,] - vst.logklengths[i,]
}

all_vstpk <- as_tibble(dds.vstpk.mat) %>%
  add_column(gene_name = rownames(dds.vstpk.mat)) %>%
  dplyr::select(gene_name, all_of(colnames(dds.vstpk.mat))) %>%
  arrange(gene_name)


# Write tidyverse compatible files to disk --------------------------------

setwd("/Users/siold/Proj/FR_RNA-Seq/shiny_august2020/FR19/data/")
write_csv(all_normalised, "allnormalised.csv")
write_csv(all_vstpk, "vstpk_data.csv")
write_csv(all_pval, "all_pval.csv")
write_csv(all_l2fc, "all_l2fc.csv")
write_csv(metadata, "metadata.csv")
write_csv(longest_transcripts, "gene_info.csv")
