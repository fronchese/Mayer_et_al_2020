#!/usr/bin/env Rscript

library(tidyverse);
library(Seurat);

matrix.fileName <-
    "M23_Combined_H2GYLDRXY_1_210203_FD09251586_Other_CGAGGCTG_R_210203_DAVGAL_INDEXLIBNOVASEQ_M001_Expression_Data.st";
tags.fileName <-
    "M23_H2GYLDRXY_1_210203_FD09251586_Other_CGAGGCTG_R_210203_DAVGAL_INDEXLIBNOVASEQ_M001_Sample_Tag_Calls.csv";

## cells as columns and features as rows

read_table2(matrix.fileName, comment="#") %>%
    pivot_wider(id_cols=Gene, names_from=Cell_Index, values_from=RSEC_Adjusted_Molecules) %>%
    column_to_rownames("Gene") %>%
    as.matrix() -> out.counts.mat

read_csv(tags.fileName, comment="#") -> cellTags;

## sanity check: make sure tag order is the same as the cell order
all(cellTags$Cell_Index == colnames(out.counts.mat));

colnames(out.counts.mat) <- paste0(cellTags$Sample_Name, "_", cellTags$Cell_Index);

pdf("out.pdf");
plot(density(log(na.omit(out.counts.mat))));
invisible(dev.off());

out.counts.mat[is.na(out.counts.mat)] <- 0;

## Okay, that's all set up; now onto Seurat
## [https://satijalab.org/seurat/articles/pbmc3k_tutorial.html]
## [https://satijalab.org/seurat/articles/sctransform_vignette.html]

c57.s6KO.all <- CreateSeuratObject(counts = out.counts.mat, project = "c57.s6KO",
                                   min.cells = 3, min.features = 200);

c57.s6KO.all[["orig.ident"]] <- Idents(c57.s6KO.all);

## Mouse mt genes start with lower case 'mt'
## [additional pattern to exclude Ribosomal and tRNA counts
c57.s6KO.all[["percent.mt"]] <- PercentageFeatureSet(c57.s6KO.all, pattern = "^mt-[^TR]");

VlnPlot(c57.s6KO.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
ggsave(sprintf("M23_FeatureCount_all_7B_%s.png", Sys.Date()), width=11, height=8);

plot1 <- FeatureScatter(c57.s6KO.all, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    geom_hline(yintercept=15);
plot2 <- FeatureScatter(c57.s6KO.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(yintercept=c(200,3500));
plot1 + plot2;
ggsave(sprintf("M23_FeatureScatter_all_7B_%s.png", Sys.Date()), width=11, height=8);

c57.s6KO.filt <-
    subset(c57.s6KO.all,
           subset = nFeature_RNA > 200 &
               nFeature_RNA < 3500 & percent.mt < 15);

c57.s6KO.filt <- SCTransform(c57.s6KO.filt, vars.to.regress="percent.mt");
## [iteration limit warnings]

c57.s6KO.filt <- RunPCA(c57.s6KO.filt, verbose = FALSE, npcs=100);

DimPlot(c57.s6KO.filt, reduction = "pca");
ggsave(sprintf("M23_PCA_filt_SCT_7B_%s.png", Sys.Date()), width=8, height=8);

png(sprintf("M23_DimHeatmap_filt_SCT_7B_%s.png", Sys.Date()), width=1280*2, height=720*2);
## Note: DimHeatmap is not a ggplot thing
DimHeatmap(c57.s6KO.filt, dims = 1:15, cells = 500, balanced = TRUE);
invisible(dev.off());

ElbowPlot(c57.s6KO.filt, ndims=100);
ggsave(sprintf("M23_Elbowplot_filt_SCT_7b_%s.png", Sys.Date()), width=8, height=8);

## Choose 50 dimensions for reductions
c57.s6KO.filt <- RunUMAP(c57.s6KO.filt, dims = 1:50, n.epochs=1000);

c57.s6KO.filt <- FindNeighbors(c57.s6KO.filt, dims = 1:50);
c57.s6KO.filt <- FindClusters(c57.s6KO.filt, resolution=0.8);
c57.s6KO.filt[["cluster"]] <- Idents(c57.s6KO.filt);

## Dims flipped to make the rabbit easier to see
DimPlot(c57.s6KO.filt, label=FALSE, group.by="orig.ident", dims=c(2,1),
        pt.size=1.5) +
    NoLegend();
ggsave(sprintf("M23_UMAP_origValues_filt_SCT_7B_noLabel_%s.png", Sys.Date()), width=11, height=8);

DimPlot(c57.s6KO.filt, label = TRUE, group.by="orig.ident", dims=c(2,1),
        pt.size=1.5);
ggsave(sprintf("M23_UMAP_origValues_filt_SCT_7B_%s.png", Sys.Date()), width=11, height=8);
DimPlot(c57.s6KO.filt, label = TRUE, group.by="cluster", dims=c(2,1),
        pt.size=1.5);
ggsave(sprintf("M23_UMAP_byCluster_filt_SCT_7B_%s.png", Sys.Date()), width=11, height=8);

## Teams meeting, 2021-Feb-16

## Almost finished Bangert clustering
## Priority - STAT6KO data

## Elbow plot - use 60 PCs
## Jack straw
## Resolution 0.8

## DC2s - express marker called SIRPalpha (main marker)
## MHC2Hi, CD11c intermediate
## XCr1
## ItGax - used for sorting
## CD326 should see some intermediate
## CD326-/Sirp-alpha
## Sub-population of cells
## Johannes to send sort panel

## List of genes that should be assoc with CD11bHi / Low
## List of genes that should be expressed in STAT6
## e.g. CCL17 (CD11BHi; IFNA5, SNED1 - CD11bLow; IFNA5 should disappear in STAT6KO)
## Do DEG on WT vs STAT6, compare with Bulk RNA Seq
## Bulk RNASeq - 300-400 genes, a bit more associated with CD11bHi (higher for some IFN response genes)
## CD11bLo should be low proportion (15%); ~15% CD326
## Not sure that CD11bLo is a subset, or a flavour of DC2s
## [feature plots with these genes]

## UMAP plot with cluster, split WT/STAT6

## UMAP clustering
