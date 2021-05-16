#!/usr/bin/env Rscript

library(tidyverse);
library(Seurat);
library(readxl);
library(patchwork);
library(scales);
library(viridis);
library(cowplot);
library(ggrepel);
library(ggtree);
library(RColorBrewer);
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

out.counts.mat[is.na(out.counts.mat)] <- 0;

## Okay, that's all set up; now onto Seurat
## [https://satijalab.org/seurat/articles/pbmc3k_tutorial.html]
## [https://satijalab.org/seurat/articles/sctransform_vignette.html]

c57.S6KO.7B.all <- CreateSeuratObject(counts = out.counts.mat, project = "c57.s6KO",
                                      min.cells = 3, min.features = 200);

c57.S6KO.7B.all[["orig.ident"]] <- Idents(c57.S6KO.7B.all);

## Mouse mt genes start with lower case 'mt'
## Note: 7B data seems lower for mt counts, even when including tRNA and rRNA
c57.S6KO.7B.all[["percent.mt"]] <- PercentageFeatureSet(c57.S6KO.7B.all, pattern = "^mt-");

VlnPlot(c57.S6KO.7B.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
ggsave(sprintf("M23_FeatureCount_all_7B_%s.png", Sys.Date()), width=11, height=8);

plot1 <- FeatureScatter(c57.S6KO.7B.all, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    geom_hline(yintercept=15) +
plot2 <- FeatureScatter(c57.S6KO.7B.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(yintercept=c(200,4000));
plot1 + plot2;
ggsave(sprintf("M23_FeatureScatter_all_7B_%s.png", Sys.Date()), width=11, height=8);

c57.S6KO.7B.filt <-
    subset(c57.S6KO.7B.all,
           subset = nFeature_RNA > 200 &
               nFeature_RNA < 4000 & percent.mt < 15);

c57.S6KO.7B.filt <- SCTransform(c57.S6KO.7B.filt, vars.to.regress="percent.mt",
                                return.only.var.genes=FALSE);
## [iteration limit warnings]

c57.S6KO.7B.filt <- RunPCA(c57.S6KO.7B.filt, verbose = FALSE, npcs=100);

DimPlot(c57.S6KO.7B.filt, reduction = "pca");
ggsave(sprintf("M23_PCA_filt_SCT_7B_%s.png", Sys.Date()), width=8, height=8);

png(sprintf("M23_DimHeatmap_filt_SCT_7B_%s.png", Sys.Date()), width=1280*2, height=720*2);
## Note: DimHeatmap is not a ggplot thing
DimHeatmap(c57.S6KO.7B.filt, dims = 1:15, cells = 500, balanced = TRUE);
invisible(dev.off());

ElbowPlot(c57.S6KO.7B.filt, ndims=100);
ggsave(sprintf("M23_Elbowplot_filt_SCT_7b_%s.png", Sys.Date()), width=8, height=8);

## Choose 60 dimensions for reductions
for(it in c(1:100,11:100 * 10)){
    c57.S6KO.7B.filt <- RunUMAP(c57.S6KO.7B.filt, dims = 1:60, n.epochs=it);
    ## Dims flipped to make the shape easier to see
    DimPlot(c57.S6KO.7B.filt, label=FALSE, group.by="orig.ident", dims=c(2,1),
            pt.size=1.5) +
        NoLegend();
    ggsave(sprintf("anim/M23_UMAP_origValues_filt_SCT_7B_noLabel_%s_frame%04d.png", Sys.Date(), it), width=11, height=8);
}

## Choose 770 epochs (because it looks like a dog's face)
c57.S6KO.7B.filt <- RunUMAP(c57.S6KO.7B.filt, dims = 1:60, n.epochs=770);

DimPlot(c57.S6KO.7B.filt, label=FALSE, group.by="orig.ident", dims=c(2,1),
        pt.size=1.5) +
    NoLegend();
ggsave(sprintf("M23_UMAP770_origValues_filt_SCT_7B_noLabel_%s.png", Sys.Date()), width=11, height=8);

c57.S6KO.7B.filt <- FindNeighbors(c57.S6KO.7B.filt, dims = 1:60);
c57.S6KO.7B.filt <- FindClusters(c57.S6KO.7B.filt, resolution=0.8);
c57.S6KO.7B.filt[["cluster"]] <- Idents(c57.S6KO.7B.filt);

DimPlot(c57.S6KO.7B.filt, label = TRUE, group.by="orig.ident", dims=c(2,1),
        pt.size=1.5);
ggsave(sprintf("M23_UMAP770_origValues_filt_SCT_7B_%s.png", Sys.Date()), width=11, height=8);
DimPlot(c57.S6KO.7B.filt, label = TRUE, group.by="cluster", dims=c(2,1),
        pt.size=1.5, cols=DiscretePalette(8));
ggsave(sprintf("M23_UMAP770_byCluster_filt_SCT_7B_%s.png", Sys.Date()), width=11, height=8);

DimPlot(c57.S6KO.7B.filt, label = TRUE, group.by="cluster",
        pt.size=1.5, cols=DiscretePalette(8), split.by="orig.ident");
ggsave(sprintf("M23_UMAP770_byCluster_splitByCell_filt_SCT_7B_%s.png", Sys.Date()), width=18, height=6);


### Alevin data processing
## UMAP clustering
library(tximeta);
library(SingleCellExperiment);
library(org.Mm.eg.db);
library(tidyverse);
##library(alevinQC);

##alevinQCReport(basedir="/mnt/ufds/jmayer/salmon_1.4_5_27_8_JM_2021-02-12_L1/alevin",
##               outputFile="alevinReport.html",
##               outputDir="~/bioinf/MIMR-2019-Dec-09-GBIS/FR/SingleCell/MIMR_2021-02-13/alevin_L1");

alevin.base.L1 <- "/mnt/ufds/jmayer/salmon_1.4_5_27_8_JM_2021-02-12_L1";
alevin.base.L2 <- "/mnt/ufds/jmayer/salmon_1.4_5_27_8_JM_2021-02-12_L2";

files.L1 <- file.path(c(alevin.base.L1), "alevin", "quants_mat.gz");
files.L2 <- file.path(c(alevin.base.L2), "alevin", "quants_mat.gz");

se.L1 <- tximeta(files.L1, type="alevin", alevinArgs=list(filterBarcodes=TRUE));
se.L2 <- tximeta(files.L2, type="alevin", alevinArgs=list(filterBarcodes=TRUE));

sce.L1 <- as(se.L1, "SingleCellExperiment");
sce.L2 <- as(se.L2, "SingleCellExperiment");

## Check number of common cells
length(intersect(colnames(sce.L1), colnames(sce.L2)));
length(union(colnames(sce.L1), colnames(sce.L2)));

sce.L1 <- addIds(sce.L1, "SYMBOL");
sce.L2 <- addIds(sce.L2, "SYMBOL");

tags.fileName <-
    "/mnt/ufds/jmayer/sampleTagCount_H2GYLDRXY_12.txt.gz";

tagCounts <- read_table2(tags.fileName, col_names=FALSE) %>%
    pivot_wider(id_cols=X2, names_from=X3, values_from=X1) %>%
    column_to_rownames("X2") %>%
    as.matrix()
## replace NA with 0
tagCounts[is.na(tagCounts)] <- 0;

## determine tag noise level (unused tags)
tagNoiseLevel <- max(rowSums(tagCounts[,-c(1,2)]));
## zero out any counts below 0.5X the noise level
tagCounts[tagCounts < (tagNoiseLevel * 0.5)] <- 0;
## filter out low-count cell tags
tagCounts <- tagCounts[rowSums(tagCounts[, c(1,2)]) > tagNoiseLevel, c(1,2)];

## Check remaining proportions
table(round(tagCounts[,2] / rowSums(tagCounts), 2));

sum((round(tagCounts[,2] / rowSums(tagCounts), 2) > 0) &
    (round(tagCounts[,2] / rowSums(tagCounts), 2) < 1));

tagCounts[(round(tagCounts[,2] / rowSums(tagCounts), 2) > 0) &
          (round(tagCounts[,2] / rowSums(tagCounts), 2) < 1),];

mixedCounts <-
    tagCounts[(round(tagCounts[,2] / rowSums(tagCounts), 2) > 0) &
              (round(tagCounts[,2] / rowSums(tagCounts), 2) <= 0.5),];

mixedCounts[order(mixedCounts[,1]),];

highReadCells <-
    read_table("/mnt/ufds/jmayer/cellReads_high.txt", col_names=FALSE)$X2

allReadCells <-
    read_table("/mnt/ufds/jmayer/cellLabels_allMatching_counts.txt.gz", col_names=FALSE)$X2

allReadData <-
    read_table("/mnt/ufds/jmayer/cellLabels_allMatching_counts.txt.gz", col_names=FALSE)


## Okay... I'm not comfortable with excluding any more of these tags,
## or setting their tag counts to zero

tagProp2 <- (round(tagCounts[,2] / rowSums(tagCounts), 2));

## Create tag labels
tagLabels <- paste0(ifelse(tagProp2 == 0, "c57",
                    ifelse(tagProp2 == 1, "Stat6KO",
                           paste0("multiplet", "_", round(tagCounts[,1] / rowSums(tagCounts), 1)))),
                    "_", rownames(tagCounts));
names(tagLabels) <- rownames(tagCounts);

length(setdiff(highReadCells, names(tagLabels)));
length(setdiff(names(tagLabels), highReadCells));

length(setdiff(names(tagLabels), allReadCells));


lowReadData <- allReadData[allReadData$X2 %in% setdiff(names(tagLabels), highReadCells),];
table(round(lowReadData$X1, -4));

## Check number of common cells
length(intersect(colnames(sce.L1), colnames(sce.L2)));

## Check number of tag labels with common cells
length(intersect(names(tagLabels), intersect(colnames(sce.L1), colnames(sce.L2))));
noTags <- setdiff(union(colnames(sce.L1), colnames(sce.L2)), names(tagLabels));

## fill in untagged cells
tagLabels[noTags] <- paste0("noTag_", noTags);

colnames(sce.L1) <- tagLabels[colnames(sce.L1)];
colnames(sce.L2) <- tagLabels[colnames(sce.L2)];

## check to make sure row names are consistent
all(rownames(sce.L1) == rownames(sce.L2));

## combine counts from L2 and L2
intNames <- intersect(colnames(sce.L1), colnames(sce.L2));
excNames.L1 <- setdiff(colnames(sce.L1), intNames);
excNames.L2 <- setdiff(colnames(sce.L2), intNames);
counts.combined <- cbind(counts(sce.L1[,intNames]) + counts(sce.L2[,intNames]),
                         counts(sce.L1[,excNames.L1]),
                         counts(sce.L2[,excNames.L2]));

## replace gene names
rownames(counts.combined) <-
    paste(mcols(sce.L1)$SYMBOL, mcols(sce.L1)$gene_id, sep="-");

## remove genes with no matching symbol (or multiple symbols)
counts.combined <- counts.combined[!grepl("^NA-", rownames(counts.combined)),];

## tximeta annotation doesn't add 'mt' to front of mt genes
mtNames <- grep("(^LOC|[0-9][0-9]\\-|^LTO|AAdac)", invert=TRUE, value=TRUE,
                grep("^[A-Z][A-Z]", rownames(counts.combined), value=TRUE));
rownames(counts.combined)[rownames(counts.combined) %in% mtNames] <-
    paste0("mt-", rownames(counts.combined)[rownames(counts.combined) %in% mtNames]);

## Okay, that's all set up; now onto Seurat
## [https://satijalab.org/seurat/articles/pbmc3k_tutorial.html]
## [https://satijalab.org/seurat/articles/sctransform_vignette.html]

c57.S6KO.alv.all <- CreateSeuratObject(counts = counts.combined, project = "c57.s6KO",
                                   min.cells = 3, min.features = 200);

c57.S6KO.alv.all[["orig.ident"]] <- Idents(c57.S6KO.alv.all);

## Mouse mt genes start with lower case 'mt'
c57.S6KO.alv.all[["percent.mt"]] <- PercentageFeatureSet(c57.S6KO.alv.all, pattern = "^mt-");

VlnPlot(c57.S6KO.alv.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
ggsave(sprintf("FeatureCount_alevin_all_%s.png", Sys.Date()), width=11, height=8);

plot1 <- FeatureScatter(c57.S6KO.alv.all, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    geom_hline(yintercept=10) +
    geom_vline(xintercept=c(5000, 25000));
plot2 <- FeatureScatter(c57.S6KO.alv.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(yintercept=c(200,4000)) +
    geom_vline(xintercept=c(5000, 25000));
plot3 <- FeatureScatter(c57.S6KO.alv.all, feature1 = "nFeature_RNA", feature2 = "percent.mt") +
    geom_hline(yintercept=10) +
    geom_vline(xintercept=c(200,4000));
plot1 | (plot2 / plot3);
ggsave(sprintf("FeatureScatter_alevin_all_%s.png", Sys.Date()), width=11, height=8);

plot1 <- FeatureScatter(c57.S6KO.alv.all, feature1 = "nCount_RNA", feature2 = "percent.mt");
plot2 <- FeatureScatter(c57.S6KO.alv.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA");
plot3 <- FeatureScatter(c57.S6KO.alv.all, feature1 = "nFeature_RNA", feature2 = "percent.mt");
plot1 | (plot2 / plot3);
ggsave(sprintf("FeatureScatter_noLines_alevin_all_%s.png", Sys.Date()), width=11, height=8);


table(c57.S6KO.alv.filt[["orig.ident"]]);

## Genes that we want to exclude as T/B-cell (or likely T/B-cell)
exclusion.set <- c("Zap70", "Trbc2", "Ms4a4b", "Trac", "Igkc", "Ighm", "Cd79a");
geneIndex <- match(exclusion.set, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.all)));
exclusion.set <- na.omit(rownames(c57.S6KO.alv.all)[geneIndex]);

cell.totalCounts <- colSums(GetAssayData(c57.S6KO.alv.all, "counts", "RNA"));

str(GetAssayData(c57.S6KO.alv.all, "counts", "RNA"));
range(GetAssayData(c57.S6KO.alv.all, "counts", "RNA"));

exclusion.mat <- GetAssayData(c57.S6KO.alv.all, "counts", "RNA")[exclusion.set,];
table(round(colSums(exclusion.mat) / (cell.totalCounts / 10000), 4) > 3);
table(colSums(exclusion.mat) > 3);

## Sventja's threshold: >3 per 10,000 molecules mapped to T/B cells
include.cells <- which(round(colSums(exclusion.mat) / (cell.totalCounts / 10000), 4) <= 0);

c57.S6KO.alv.filt <-
    subset(c57.S6KO.alv.all,
           subset = nFeature_RNA > 200 &
               nFeature_RNA < 4000 &
               percent.mt < 10 &
               nCount_RNA > 5000 &
               nCount_RNA < 25000 &
               orig.ident %in% c("c57", "Stat6KO"),
           cells = include.cells);

## Find variable gene set via cell subsampling
set.geneNames <- NULL;
set.replicate <- NULL;
for(rep in 1:20){
    cat(rep, "\n");
    c57.S6KO.alv.filt.subsamp <- SCTransform(c57.S6KO.alv.filt, vars.to.regress="percent.mt",
                                             ncells=2000, seed=NULL, variable.features.n=4000);
    geneList <- rownames(GetAssayData(c57.S6KO.alv.filt.subsamp, "scale.data", "SCT"));
    set.geneNames <- c(set.geneNames, geneList);
    set.replicate <- c(set.replicate, rep(rep, length(geneList)));
}

tSet <- table(set.geneNames);
idealSet <- setdiff(names(tSet)[which(tSet == max(tSet))], "Stat6-ENSMUSG00000002147.18");
c57.S6KO.alv.filt.bs <- SCTransform(c57.S6KO.alv.filt, vars.to.regress="percent.mt",
                                    residual.features=idealSet,
                                    return.only.var.genes = FALSE);

c57.S6KO.alv.filt.bs <- RunPCA(c57.S6KO.alv.filt.bs, verbose = TRUE, npcs=100);
c57.S6KO.alv.filt.bs <- RunUMAP(c57.S6KO.alv.filt.bs, dims = 1:60, n.epochs=320);
c57.S6KO.alv.filt.bs <- FindNeighbors(c57.S6KO.alv.filt.bs, dims = 1:60);
c57.S6KO.alv.filt.bs <- FindClusters(c57.S6KO.alv.filt.bs, resolution=0.8);
c57.S6KO.alv.filt.bs[["cluster"]] <- Idents(c57.S6KO.alv.filt.bs);

DimPlot(c57.S6KO.alv.filt.bs, label=FALSE, group.by="cluster", pt.size=1.5, split.by="orig.ident") +
    ggtitle(sprintf("Bootstrapped SCTransform without Stat6, cluster resolution 0.8"));
ggsave(sprintf("M23_UMAP320_filt_bootstrapped_SCT_alv_noStat6_%s.png", Sys.Date()), width=12, height=8);

idealSet <- names(tSet)[which(tSet == max(tSet))];
c57.S6KO.alv.filt.bs2 <- SCTransform(c57.S6KO.alv.filt, vars.to.regress="percent.mt",
                                    residual.features=idealSet,
                                    return.only.var.genes = FALSE);

c57.S6KO.alv.filt.bs2 <- RunPCA(c57.S6KO.alv.filt.bs2, verbose = TRUE, npcs=100);
c57.S6KO.alv.filt.bs2 <- RunUMAP(c57.S6KO.alv.filt.bs2, dims = 1:60, n.epochs=320);
c57.S6KO.alv.filt.bs2 <- FindNeighbors(c57.S6KO.alv.filt.bs2, dims = 1:60);
c57.S6KO.alv.filt.bs2 <- FindClusters(c57.S6KO.alv.filt.bs2, resolution=0.8);
c57.S6KO.alv.filt.bs2[["cluster"]] <- Idents(c57.S6KO.alv.filt.bs2);
c57.S6KO.alv.filt.bs[["cluster.bs2"]] <- Idents(c57.S6KO.alv.filt.bs2);
c57.S6KO.alv.filt.bs2[["cluster.bs"]] <- Idents(c57.S6KO.alv.filt.bs);
DimPlot(c57.S6KO.alv.filt.bs2, label=FALSE, group.by="cluster", pt.size=1.5, split.by="orig.ident") +
    ggtitle(sprintf("Bootstrapped SCTransform with Stat6, cluster resolution 0.8"));
ggsave(sprintf("M23_UMAP320_filt_bootstrapped_SCT_alv_withStat6_%s.png", Sys.Date()), width=12, height=8);

DimPlot(c57.S6KO.alv.filt.bs, label=FALSE, group.by="cluster.bs2", pt.size=1.5, split.by="orig.ident") +
    ggtitle(sprintf("Bootstrapped SCTransform with Stat6, cluster resolution 0.8, remapped to without Stat6 UMAP"));
ggsave(sprintf("M23_UMAP320_filt_bootstrapped_SCT_alv_withStat6_remappedWithout_%s.png", Sys.Date()), width=12, height=8);

DimPlot(c57.S6KO.alv.filt.bs2, label=FALSE, group.by="cluster.bs", pt.size=1.5, split.by="orig.ident") +
    ggtitle(sprintf("Bootstrapped SCTransform without Stat6, cluster resolution 0.8, remapped to with Stat6 UMAP"));
ggsave(sprintf("M23_UMAP320_filt_bootstrapped_SCT_alv_withoutStat6_remappedWith_%s.png", Sys.Date()), width=12, height=8);

c57.S6KO.alv.filt.nobs <- SCTransform(c57.S6KO.alv.filt, vars.to.regress="percent.mt",
                                      return.only.var.genes = FALSE);

c57.S6KO.alv.filt.nobs <- RunPCA(c57.S6KO.alv.filt.nobs, verbose = TRUE, npcs=100);
c57.S6KO.alv.filt.nobs <- RunUMAP(c57.S6KO.alv.filt.nobs, dims = 1:60, n.epochs=320);
c57.S6KO.alv.filt.nobs <- FindNeighbors(c57.S6KO.alv.filt.nobs, dims = 1:60);
c57.S6KO.alv.filt.nobs <- FindClusters(c57.S6KO.alv.filt.nobs, resolution=0.8);
c57.S6KO.alv.filt.nobs[["cluster"]] <- Idents(c57.S6KO.alv.filt.nobs);
DimPlot(c57.S6KO.alv.filt.nobs, label=FALSE, group.by="cluster", pt.size=1.5, split.by="orig.ident") +
    ggtitle(sprintf("Non-bootstrapped SCTransform with Stat6, cluster resolution 0.8"));
ggsave(sprintf("M23_UMAP320_filt_nonBootstrapped_SCT_alv_withStat6_%s.png", Sys.Date()), width=12, height=8);

c57.S6KO.alv.filt.noStat6 <- SCTransform(c57.S6KO.alv.filt, vars.to.regress="percent.mt",
                                       variable.features.n = 3001,
                                       return.only.var.genes = FALSE);

c57.S6KO.alv.filt.noStat6 <- RunPCA(c57.S6KO.alv.filt.noStat6, verbose = TRUE, npcs=100,
                                  features=setdiff(VariableFeatures(c57.S6KO.alv.filt.noStat6),
                                                   "Stat6-ENSMUSG00000002147.18"));
c57.S6KO.alv.filt.noStat6 <- RunUMAP(c57.S6KO.alv.filt.noStat6, dims = 1:60, n.epochs=320);
c57.S6KO.alv.filt.noStat6 <- FindNeighbors(c57.S6KO.alv.filt.noStat6, dims = 1:60);




for(cRes in seq(1.5, 1.5, by=0.1)){
    c57.S6KO.alv.filt.noStat6 <- FindClusters(c57.S6KO.alv.filt.noStat6, resolution=cRes);
    c57.S6KO.alv.filt.noStat6[["orig.ident"]] <- sub("c57", "WT", as.character(unlist(c57.S6KO.alv.filt.noStat6[["orig.ident"]])));
    c57.S6KO.alv.filt.noStat6[["cluster"]] <- Idents(c57.S6KO.alv.filt.noStat6);
    clen <- length(unique(unlist(c57.S6KO.alv.filt.noStat6[["cluster"]])));
    if(clen <= 8){
        clusterCols.noStat6 <- brewer.pal(clen, "Dark2");
    } else {
        clusterCols.noStat6 <- viridis(clen);
    }
    DimPlot(c57.S6KO.alv.filt.noStat6, label=FALSE, group.by="cluster", pt.size=1.5, split.by="orig.ident",
            cols=clusterCols.noStat6) +
        ggtitle(sprintf("SCTransform without Stat6, cluster resolution %0.1f", cRes));
    ggsave(sprintf("M23_UMAP320_filt_SCT_alv_withoutStat6_%0.1fres_%s.png", cRes, Sys.Date()), width=12, height=8);
    ## DotPlot
    genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);
    geneIndex <- match(c(genes.OL13, "Zap70", "Stat6"), sub("-ENSM.*$", "", rownames(c57.S6KO.alv.filt.noStat6)));
    listSub <- na.omit(rownames(c57.S6KO.alv.filt.noStat6)[geneIndex]);
    maxRes <- NULL;
    cairo_pdf(sprintf("DotPlot_OL13_byCluster_alv_withoutStat6_%0.1fres_%s.pdf", cRes, Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.filt.noStat6, features=listSub,
                   group.by="cluster",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
            scale_colour_viridis() +
            ylab("Cluster") +
            ggtitle(sprintf("Gene Set - OL13 (res: %0.1f)", cRes)) +
            scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
    ## Cluster counts
    idents.rt <- as.character(unlist(c57.S6KO.alv.filt.noStat6[["cluster"]]));
    res.cluster <- table(list(idents.rt, as.character(unlist(c57.S6KO.alv.filt.noStat6[["orig.ident"]]))));
    res.cluster <- res.cluster[order(as.numeric(rownames(res.cluster))),];
    write.csv(res.cluster, sprintf("cluster_counts_withoutStat6_%0.1fres_%s.csv", cRes, Sys.Date()));
}


## Choose resolution=1.5 as the target
c57.S6KO.alv.filt.noStat6 <- FindClusters(c57.S6KO.alv.filt.noStat6, resolution=1.5);
c57.S6KO.alv.filt.noStat6[["orig.ident"]] <- sub("c57", "WT", as.character(unlist(c57.S6KO.alv.filt.noStat6[["orig.ident"]])));
c57.S6KO.alv.filt.noStat6[["cluster"]] <- Idents(c57.S6KO.alv.filt.noStat6);

## Confirm that Stat6 is still in the results
grep("Stat6", rownames(c57.S6KO.alv.filt.noStat6), value=TRUE);

loading.mat <-
    apply(c57.S6KO.alv.filt.noStat6[["pca"]]@feature.loadings,2,function(x){x <- x[order(-abs(x))]; names(x)})[,1:60];

loading.vals <-
    apply(c57.S6KO.alv.filt.noStat6[["pca"]]@feature.loadings,2,function(x){x <- x[order(-abs(x))]; x})[,1:60];


name.tbl <- as_tibble(loading.mat[1:30,]) %>%
    mutate(PCrank=1:30) %>%
    pivot_longer(names_to="PC", values_to="gene", cols=-PCrank)

value.tbl <- as_tibble(loading.vals[1:30,]) %>%
    mutate(PCrank=1:30) %>%
    pivot_longer(names_to="PC", values_to="loading", cols=-PCrank)

## Stat6 is no longer in the PCA
##stat6.pos <- which(loading.mat == "Stat6-ENSMUSG00000002147.18", arr.ind=TRUE);
##stat6.val <- loading.vals[stat6.pos];

merged.tbl <- name.tbl %>%
    inner_join(value.tbl) %>%
    ## bind_rows(tibble(PCrank=stat6.pos[,1],
    ##                  PC=colnames(loading.mat)[stat6.pos[,2]],
    ##                  gene="Stat6-ENSMUSG00000002147.18",
    ##                  loading=stat6.val)) %>%
    arrange(PC, PCrank)

write_csv(merged.tbl[,c("PC","PCrank", "gene", "loading")], sprintf("PC_loadings_mostExtreme30_%s.csv", Sys.Date()));


max.loadings <- apply(c57.S6KO.alv.filt.noStat6[["pca"]]@feature.loadings,2,function(x){nameX <- names(x)[which(x == max(x))]; max(x)});
median.loadings <- apply(c57.S6KO.alv.filt.noStat6[["pca"]]@feature.loadings,2,function(x){nameX <- names(x)[which(x == median(x))]; median(x)});
##Stat6.loadings <- c57.S6KO.alv.filt[["pca"]]@feature.loadings["Stat6-ENSMUSG00000002147.18",];

head(cbind(max.loadings, median.loadings))

DimPlot(c57.S6KO.alv.filt.noStat6, reduction = "pca");
ggsave(sprintf("PCA_filt_byCluster_SCT_alv_withoutStat6_%s.png", Sys.Date()), width=8, height=8);

DimPlot(c57.S6KO.alv.filt.noStat6, reduction = "pca", group.by="orig.ident");
ggsave(sprintf("PCA_filt_byCellType_SCT_alv_withoutStat6_%s.png", Sys.Date()), width=8, height=8);

png(sprintf("DimHeatmap_filt_SCT_alv_withoutStat6_%s.png", Sys.Date()), width=1280*2, height=720*2);
## Note: DimHeatmap is not a ggplot thing
DimHeatmap(c57.S6KO.alv.filt.noStat6, dims = 1:15, cells = 500, balanced = TRUE);
invisible(dev.off());

ElbowPlot(c57.S6KO.alv.filt.noStat6, ndims=100);
ggsave(sprintf("Elbowplot_filt_SCT_alv_withoutStat6_%s.png", Sys.Date()), width=8, height=8);

## ## Choose 60 dimensions for reductions
## for(it in c(1:100,11:100 * 10)){
##     c57.S6KO.alv.filt <- RunUMAP(c57.S6KO.alv.filt, dims = 1:60, n.epochs=it);
##     ## Dims flipped to make the shape easier to see
##     DimPlot(c57.S6KO.alv.filt, label=FALSE, group.by="orig.ident", dims=c(2,1),
##             pt.size=1.5) +
##         NoLegend();
##     ggsave(sprintf("anim_alv/M23_UMAP_origValues_filt_SCT_7B_noLabel_%s_frame%04d.png", Sys.Date(), it), width=11, height=8);
## }

## ## Choose 320, 60 dimensions
## c57.S6KO.alv.filt <- RunUMAP(c57.S6KO.alv.filt, dims = 1:60, n.epochs=320);
## c57.S6KO.alv.filt.noStat6 <- RunUMAP(c57.S6KO.alv.filt.noStat6, dims = 1:60, n.epochs=320);

## c57.S6KO.alv.filt <- FindNeighbors(c57.S6KO.alv.filt, dims = 1:60);
## c57.S6KO.alv.filt <- FindClusters(c57.S6KO.alv.filt, resolution=0.8);
## c57.S6KO.alv.filt[["cluster"]] <- Idents(c57.S6KO.alv.filt);

## c57.S6KO.alv.filt.noStat6 <- FindNeighbors(c57.S6KO.alv.filt.noStat6, dims = 1:60);
## c57.S6KO.alv.filt.noStat6 <- FindClusters(c57.S6KO.alv.filt.noStat6, resolution=0.8);
## c57.S6KO.alv.filt.noStat6[["cluster"]] <- Idents(c57.S6KO.alv.filt.noStat6);
## c57.S6KO.alv.filt[["cluster.noStat6"]] <- Idents(c57.S6KO.alv.filt.noStat6);

## ## Sanity check to make sure that rows are ordered the same [yes]
## all(rownames(c57.S6KO.alv.filt.noStat6[["cluster"]]) == rownames(c57.S6KO.alv.filt[["cluster"]]));

## combined.mat <- cbind(withStat6=c57.S6KO.alv.filt[["cluster"]], noStat6=c57.S6KO.alv.filt.noStat6[["cluster"]]);
## table(withStat6=combined.mat[,1], noStat6=combined.mat[,2]);

## ##                                           *   *   * *   *
## clusterCols.orig <- brewer.pal(8, "Dark2")[c(1,5,3,6,7,8,4,2)];
## ##clusterCols.orig <- brewer.pal(8, "Dark2");
## names(clusterCols.orig) <- unique(sort(as.character(unlist(c57.S6KO.alv.noMx[["cluster"]]))));

## ## Comparison DimPlots to test exclusion of Stat6
## DimPlot(c57.S6KO.alv.filt, label=FALSE, group.by="cluster", pt.size=1.5, cols=clusterCols.orig) +
## ggtitle("SCTransform with Stat6");
## ggsave(sprintf("M23_UMAP320_origValues_filt_SCT_alv_%s.png", Sys.Date()), width=11, height=8);
## DimPlot(c57.S6KO.alv.filt, label=FALSE, group.by="cluster", pt.size=1.5, cols=clusterCols.orig, split.by="orig.ident") +
## ggtitle("SCTransform with Stat6");
## ggsave(sprintf("M23_UMAP320_origValues_filt_split_SCT_alv_%s.png", Sys.Date()), width=11, height=8);

## clusterCols.noStat6 <- clusterCols.orig[c(3,1,2,4,6,5,7,8)];
## names(clusterCols.noStat6) <- c(unique(sort(as.character(unlist(c57.S6KO.alv.filt.noStat6[["cluster"]])))), "Other");

## DimPlot(c57.S6KO.alv.filt, label=FALSE, group.by="cluster.noStat6", pt.size=1.5, cols=clusterCols.noStat6) +
## ggtitle("SCTransform without Stat6 (mapped to 'with Stat6' UMAP)");
## ggsave(sprintf("M23_UMAP320_origValues_filt_Stat6Remapped_SCT_alv_%s.png", Sys.Date()), width=11, height=8);
## DimPlot(c57.S6KO.alv.filt, label=FALSE, group.by="cluster.noStat6", pt.size=1.5, cols=clusterCols.noStat6, split.by="orig.ident") +
## ggtitle("SCTransform without Stat6 (mapped to 'with Stat6' UMAP)");
## ggsave(sprintf("M23_UMAP320_origValues_filt_Stat6Remapped_split_SCT_alv_%s.png", Sys.Date()), width=11, height=8);

## DimPlot(c57.S6KO.alv.filt.noStat6, label=FALSE, group.by="cluster", pt.size=1.5, cols=clusterCols.noStat6) +
## ggtitle("SCTransform without Stat6");
## ggsave(sprintf("M23_UMAP320_origValues_filt_noStat6_SCT_alv_%s.png", Sys.Date()), width=11, height=8);

## DimPlot(c57.S6KO.alv.filt, label = TRUE, group.by="orig.ident", cols=DiscretePalette(4), pt.size=1.5);
## ggsave(sprintf("M23_UMAP320_origValues_alevin_filt_SCT_%s.png", Sys.Date()), width=11, height=8);
## DimPlot(c57.S6KO.alv.filt, label=FALSE, group.by="cluster", cols=DiscretePalette(8), pt.size=1.5);
## ggsave(sprintf("M23_UMAP320_byCluster_alevin_filt_SCT_%s.png", Sys.Date()), width=11, height=8);
## DimPlot(c57.S6KO.alv.filt, label = TRUE, group.by="cluster", cols=DiscretePalette(8), pt.size=1.5);
## ggsave(sprintf("M23_UMAP320_byCluster_alevin_filt_SCT_%s.png", Sys.Date()), width=11, height=8);

## DimPlot(c57.S6KO.alv.filt, label = TRUE, group.by="cluster", cols=DiscretePalette(8),
##         pt.size=1.5, split.by="orig.ident") +
##         geom_hline(yintercept=3.8);
## ggsave(sprintf("M23_UMAP320_byCluster_splitByCell_alevin_filt_SCT_%s.png", Sys.Date()), width=18, height=6);

## c57.S6KO.alv.filt[["shape"]] <- factor(16);
## { ## Alevin / Ighm + Cd79a
##     listSub <- c("Ighmbp2", "Cd79a", "Ms4a1", "Trbc2", "Ms4a4b", "Tcra", "Tcrb", "Tcrd", "Cd274", "Zap70");
##     geneIndex <- match(listSub, sub("-.*$", "", rownames(c57.S6KO.alv.filt)));
##     listSub <- na.omit(rownames(c57.S6KO.alv.filt)[geneIndex]);
##     selectPoss <- seq(1, length(listSub), by=4);
##     maxRes <- NULL;
##     cairo_pdf(sprintf("FeaturePlot_OL4_alv_filt_%s.pdf", Sys.Date()), onefile=TRUE,
##               width=8, height=11, pointsize=8);
##     for(listStart in selectPoss){
##         cat(listStart, "\n");
##         geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
##         res <- list();
##         for(gene in geneSubset){
##             res <- c(res, FeaturePlot(c57.S6KO.alv.filt, col=c("lightgrey", "#e31836"),
##                                       split.by="orig.ident", features=gene,
##                                       order=TRUE, pt.size=1, shape.by="shape",
##                                       combine=FALSE));
##         }
##         res <- lapply(res, function(x){
##             x$labels$title <- sub("-.*$", "", x$labels$title);
##             x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
##         });
##         resPre <- res;
##         if(length(res) < 8){
##             plotsToDo <- (8 - length(res)) / 2;
##             plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
##             plotsToAdd <- lapply(plotsToAdd, function(x){
##                 for(n in names(x$labels)){
##                     x$labels[[n]] <- "";
##                 }
##                 x + xlab("") + ylab("") + guides(col=FALSE) +
##                     scale_shape_manual(values=NA) +
##                     ggtitle("") +
##                     theme(panel.border=element_blank(),
##                           axis.line=element_blank(),
##                           axis.text=element_blank(),
##                           axis.text.y.right=element_blank(),
##                           axis.title.y.right=element_blank(),
##                           axis.ticks=element_blank());
##             });
##             for(plotSeq in seq_len(plotsToDo)){
##                 res <- c(res, plotsToAdd);
##             }
##         }
##         suppressWarnings(print(wrap_plots(res, ncol=2)));
##     }
##     invisible(dev.off());
## }

## c57.S6KO.alv.filt[["shape"]] <- factor(16);
## { ## Alevin / Cluster 7 markers
##     listSub <- c("Cd74", "Ms4a1", "Ms4a4b", "Zap70");
##     geneIndex <- match(listSub, sub("-.*$", "", rownames(c57.S6KO.alv.filt)));
##     listSub <- na.omit(rownames(c57.S6KO.alv.filt)[geneIndex]);
##     selectPoss <- seq(1, length(listSub), by=4);
##     maxRes <- NULL;
##     cairo_pdf(sprintf("FeaturePlot_Cd74-Ms4a1-Ms4a4b-Zap70_alv_filt_%s.pdf", Sys.Date()), onefile=TRUE,
##               width=8, height=11, pointsize=8);
##     for(listStart in selectPoss){
##         cat(listStart, "\n");
##         geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
##         res <- list();
##         for(gene in geneSubset){
##             res <- c(res, FeaturePlot(c57.S6KO.alv.filt, col=c("lightgrey", "#e31836"),
##                                       split.by="orig.ident", features=gene,
##                                       order=TRUE, pt.size=1, shape.by="shape",
##                                       combine=FALSE));
##         }
##         res <- lapply(res, function(x){
##             x$labels$title <- sub("-.*$", "", x$labels$title);
##             x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
##         });
##         resPre <- res;
##         if(length(res) < 8){
##             plotsToDo <- (8 - length(res)) / 2;
##             plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
##             plotsToAdd <- lapply(plotsToAdd, function(x){
##                 for(n in names(x$labels)){
##                     x$labels[[n]] <- "";
##                 }
##                 x + xlab("") + ylab("") + guides(col=FALSE) +
##                     scale_shape_manual(values=NA) +
##                     ggtitle("") +
##                     theme(panel.border=element_blank(),
##                           axis.line=element_blank(),
##                           axis.text=element_blank(),
##                           axis.text.y.right=element_blank(),
##                           axis.title.y.right=element_blank(),
##                           axis.ticks=element_blank());
##             });
##             for(plotSeq in seq_len(plotsToDo)){
##                 res <- c(res, plotsToAdd);
##             }
##         }
##         suppressWarnings(print(wrap_plots(res, ncol=2)));
##     }
##     invisible(dev.off());
## }

## ## look at extent of cluster 7
## umap.mx <- c57.S6KO.alv.filt[["umap"]]@cell.embeddings;
## range(umap.mx[c57.S6KO.alv.filt[["cluster"]] == 7,"UMAP_2"]);

## 2021-Mar-03 choose UMAP_2 < 3.8 as a cutoff

## now we have two UMAPped / clustered datasets
## c57.S6KO.alv.filt
## c57.S6KO.7B.filt

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

## ## UMAP plot with cluster, split WT/STAT6

## clustersToCheck <- as.character(unique(unlist(c57.S6KO.alv.filt[["cluster"]])));

## {
##     marker.list <- NULL;
##     for(id1 in clustersToCheck){
##         cat(id1,"\n");
##         test.markers <- FindMarkers(c57.S6KO.alv.filt, group.by="cluster", test.use="DESeq2",
##                                     ident.1 = id1, ident.2 = setdiff(clustersToCheck, id1));
##         test.markers <- test.markers[order(test.markers$avg_log2FC),];
##         marker.list <- rbind(marker.list, cbind(cluster=id1, head(test.markers, 10)));
##         marker.list <- rbind(marker.list, cbind(cluster=id1, tail(test.markers, 10)));
##     }
##     marker.list.adjFiltered <- unique(subset(marker.list, p_val_adj < 0.1));
##     write.csv(marker.list.adjFiltered,
##               file=sprintf("ClusterDifferentiatingMarkers_all_alevin_%s.csv",
##                            format(Sys.Date())))
##     DoHeatmap(c57.S6KO.alv.filt, features=rownames(marker.list.adjFiltered), slot="data", group.by="cluster") +
##         scale_fill_viridis();
##     ggsave(sprintf("Heatmap_byCluster_clusterDifferentiating_M23_alv_filt_SCT_%s.png", Sys.Date()), width=11, height=8);
##     ## same for SevenBridges
##     clustersToCheck <- as.character(unique(unlist(c57.S6KO.7B.filt[["cluster"]])));
##     marker.list.7B <- NULL;
##     for(id1 in clustersToCheck){
##         cat(id1,"\n");
##         test.markers <- FindMarkers(c57.S6KO.7B.filt, group.by="cluster", test.use="DESeq2",
##                                     ident.1 = id1, ident.2 = setdiff(clustersToCheck, id1));
##         test.markers <- test.markers[order(test.markers$avg_log2FC),];
##         marker.list.7B <- rbind(marker.list.7B, cbind(cluster=id1, head(test.markers, 10)));
##         marker.list.7B <- rbind(marker.list.7B, cbind(cluster=id1, tail(test.markers, 10)));
##     }
##     marker.list.adjFiltered.7B <- unique(subset(marker.list.7B, p_val_adj < 0.1));
##     write.csv(marker.list.adjFiltered.7B,
##               file=sprintf("ClusterDifferentiatingMarkers_all_7B_%s.csv",
##                            format(Sys.Date())));
##     DefaultAssay(c57.S6KO.7B.filt) <- "RNA";
##     DoHeatmap(c57.S6KO.7B.filt, features=rownames(marker.list.adjFiltered.7B), slot="data", group.by="cluster") +
##         scale_fill_viridis();
##     ggsave(sprintf("Heatmap_byCluster_clusterDifferentiating_M23_7B_filt_SCT_%s.png", Sys.Date()), width=11, height=8);
## }

## DefaultAssay(c57.S6KO.7B.filt) <- "RNA";
## DoHeatmap(c57.S6KO.7B.filt, features="Stat6", group.by="orig.ident", slot="data");

## ggsave(sprintf("Heatmap_Stat6_M23_7B_filt_SCT_%s.png", Sys.Date()), width=11, height=8);

##geneLists <- read_excel("List of genes associated with WT and STAT6 cells.xlsx");

## DefaultAssay(c57.S6KO.7B.filt) <- "SCT";
## c57.S6KO.7B.filt[["shape"]] <- factor(16);
## for(listName in unique(geneLists$Category)){
##     cat(listName, "\n");
##     listSub <- geneLists$Gene[geneLists$Category == listName];
##     listName <- gsub(" ", "-", listName);
##     DoHeatmap(c57.S6KO.7B.filt, features=listSub, slot="data", group.by="orig.ident") +
##         scale_fill_viridis();
##     ggsave(sprintf("Heatmap_OL3-%s_M23_7B_filt_SCT_%s.png", listName, Sys.Date()), width=11, height=8);
##     DoHeatmap(c57.S6KO.7B.filt, features=listSub, slot="data", group.by="cluster") +
##         scale_fill_viridis();
##     ggsave(sprintf("Heatmap_byCluster_OL3-%s_M23_7B_filt_SCT_%s.png", listName, Sys.Date()), width=11, height=8);
##     for(listStart in seq(1, length(listSub), by=9)){
##         cat(listStart, "\n");
##         geneSubset <- na.omit(listSub[listStart:(listStart+8)]);
##         res <- FeaturePlot(c57.S6KO.7B.filt, slot="data", col=c("lightgrey", "#e31836"),
##                            features=geneSubset, order=TRUE, pt.size=1, shape.by="shape", combine=FALSE);
##         res <- lapply(res, function(x){
##             x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
##         });
##         print(wrap_plots(res, ncol=3));
##         ggsave(sprintf("FeaturePlot_OL3-%s_M23_7B_filt_SCT_%02d_%s.png", listName, listStart, Sys.Date()), width=11, height=8);
##     }
## }

## DefaultAssay(c57.S6KO.alv.filt) <- "SCT";
## c57.S6KO.alv.filt[["shape"]] <- factor(16);
## for(listName in unique(geneLists$Category)){
##     cat(listName, "\n");
##     listSub <- geneLists$Gene[geneLists$Category == listName];
##     if(!all(listSub %in% sub("-.*$", "", rownames(c57.S6KO.alv.filt)))){
##         cat("Not all genes found; these were not: ",
##             listSub[!(listSub %in% sub("-.*$", "", rownames(c57.S6KO.alv.filt)))],
##             "\n");
##     }
##     geneIndex <- match(listSub, sub("-.*$", "", rownames(c57.S6KO.alv.filt)));
##     listSub <- na.omit(rownames(c57.S6KO.alv.filt)[geneIndex]);
##     listName <- gsub(" ", "-", listName);
##     DoHeatmap(c57.S6KO.alv.filt, features=listSub, slot="data", group.by="orig.ident") +
##         scale_fill_viridis();
##     ggsave(sprintf("Heatmap_OL3-%s_M23_alv_filt_SCT_%s.png", listName, Sys.Date()), width=11, height=8);
##     DoHeatmap(c57.S6KO.alv.filt, features=listSub, slot="data", group.by="cluster") +
##         scale_fill_viridis();
##     ggsave(sprintf("Heatmap_byCluster_OL3-%s_M23_alv_filt_SCT_%s.png", listName, Sys.Date()), width=11, height=8);
##     for(listStart in seq(1, length(listSub), by=9)){
##         cat(listStart, "\n");
##         geneSubset <- na.omit(listSub[listStart:(listStart+8)]);
##         res <- FeaturePlot(c57.S6KO.alv.filt, slot="data", col=c("lightgrey", "#e31836"),
##                            features=geneSubset, order=TRUE, pt.size=1, shape.by="shape", combine=FALSE);
##         res <- lapply(res, function(x){
##             x$labels$title <- sub("-.*$", "", x$labels$title);
##             x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
##         });
##         print(wrap_plots(res, ncol=3));
##         ggsave(sprintf("FeaturePlot_OL3-%s_M23_alv_filt_SCT_%02d_%s.png", listName, listStart, Sys.Date()), width=11, height=8);
##     }
## }

## 2020-Feb-19
## Discussion with FR
## - find most distinguishing clusters
## - do DE of c57 vs Stat6KO in cells within those clusters

## ## Look at c57 vs Stat6KO proportions in each cluster
## tc.alv <- table(unlist(c57.S6KO.alv.filt[["orig.ident"]]), unlist(c57.S6KO.alv.filt[["cluster"]]));
## tc.alv;
## round(t(t(tc.alv) / colSums(tc.alv)), 2);
## abs(round(t(t(tc.alv) / colSums(tc.alv)), 2) - .5);

## ## Clusters 6, 2, 3 are the most distinguishing
## c57.S6KO.alv.distSub <-
##     subset(c57.S6KO.alv.filt, subset=cluster %in% c("6", "2", "3"));
## ncol(c57.S6KO.alv.distSub); ## 926 cells
## table(c57.S6KO.alv.distSub[["orig.ident"]]); ## 372 / 508

## ## find distinguishing markers
## dist.markers.alv <- FindMarkers(c57.S6KO.alv.distSub, group.by="orig.ident", test.use="DESeq2",
##                                 ident.1 = "c57", ident.2 = "Stat6KO");

## ## write to file
## dist.markers.alv <- dist.markers.alv[order(-dist.markers.alv$avg_log2FC),]
## dist.markers.alv <- dist.markers.alv[!is.na(dist.markers.alv$p_val_adj),]
## write.csv(dist.markers.alv[dist.markers.alv$p_val_adj < 0.1,],
##           file=sprintf("CellDistinguishingMarkers_cluster-2-3-6_alv_%s.csv",
##           format(Sys.Date())));

## ## look at differentiating cluster vs all non-differentiating clusters
## {
##     c57.S6KO.alv.filt[["pop"]] <-
##         ifelse(unlist(c57.S6KO.alv.filt[["cluster"]]) %in% c("6", "2", "3"),
##                as.character(unlist(c57.S6KO.alv.filt[["cluster"]])), "other");
##     marker.list <- NULL;
##     for(cli in c("6", "2", "3")){
##         clCells <- which(unlist(c57.S6KO.alv.filt[["pop"]]) == cli);
##         cat(cli,"\n");
##         test.markers <- FindMarkers(c57.S6KO.alv.filt, test.use="DESeq2", group.by="pop",
##                                     ident.1 = cli, ident.2 = "other");
##         test.markers <- test.markers[order(-test.markers$avg_log2FC),];
##         test.markers <- test.markers[!is.na(test.markers$p_val_adj),];
##         test.markers <- test.markers[test.markers$p_val_adj < 0.1,];
##         markersToAdd <- unique(rbind(head(test.markers, 20), tail(test.markers, 20)));
##         marker.list <- rbind(marker.list,
##                              cbind(cluster=cli, marker=rownames(markersToAdd),
##                                    markersToAdd));
##     }
##     write_csv(marker.list,
##               file=sprintf("ClusterDistinguishingMarkers_cluster-2-3-6_alv_%s.csv",
##                            format(Sys.Date())));
## }


## tc.7B <- table(unlist(c57.S6KO.7B.filt[["orig.ident"]]), unlist(c57.S6KO.7B.filt[["cluster"]]));
## tc.7B;
## round(t(t(tc.7B) / colSums(tc.7B)), 2);
## abs(round(t(t(tc.7B) / colSums(tc.7B)), 2) - .5);

## ## Clusters 1, 5 are the most distinguishing
## c57.S6KO.7B.distSub <-
##     subset(c57.S6KO.7B.filt, subset=cluster %in% c("1", "5"));
## ncol(c57.S6KO.7B.distSub); ## 871 cells
## table(c57.S6KO.7B.distSub[["orig.ident"]]); ## 349 / 514

## ## find distinguishing markers
## dist.markers.7B <- FindMarkers(c57.S6KO.7B.distSub, group.by="orig.ident", test.use="DESeq2",
##                                ident.1 = "C57", ident.2 = "STAT6KO");

## ## write to file
## dist.markers.7B <- dist.markers.7B[order(-dist.markers.7B$avg_log2FC),]
## dist.markers.7B <- dist.markers.7B[!is.na(dist.markers.7B$p_val_adj),]
## write.csv(dist.markers.7B[dist.markers.7B$p_val_adj < 0.1,],
##           file=sprintf("CellDistinguishingMarkers_cluster-1-5_7B_%s.csv",
##           format(Sys.Date())));

## ## look at differentiating cluster vs all non-differentiating clusters
## {
##     c57.S6KO.7B.filt[["pop"]] <-
##         ifelse(unlist(c57.S6KO.7B.filt[["cluster"]]) %in% c("1", "5"),
##                as.character(unlist(c57.S6KO.7B.filt[["cluster"]])), "other");
##     marker.list.7B <- NULL;
##     for(cli in c("1", "5")){
##         clCells <- which(unlist(c57.S6KO.7B.filt[["pop"]]) == cli);
##         cat(cli,"\n");
##         test.markers <- FindMarkers(c57.S6KO.7B.filt, test.use="DESeq2", group.by="pop",
##                                     ident.1 = cli, ident.2 = "other");
##         test.markers <- test.markers[order(-test.markers$avg_log2FC),];
##         test.markers <- test.markers[!is.na(test.markers$p_val_adj),];
##         test.markers <- test.markers[test.markers$p_val_adj < 0.1,];
##         markersToAdd <- unique(rbind(head(test.markers, 20), tail(test.markers, 20)));
##         marker.list.7B <- rbind(marker.list.7B,
##                              cbind(cluster=cli, marker=rownames(markersToAdd),
##                                    markersToAdd));
##     }
##     write_csv(marker.list.7B,
##               file=sprintf("ClusterDistinguishingMarkers_cluster-1-5_7B_%s.csv",
##                            format(Sys.Date())));
## }

## No need to remove additional cells / doublets, because we've already done that
c57.S6KO.alv.noMx <- c57.S6KO.alv.filt.noStat6;

## Refactor clusters
c57.S6KO.alv.noMx[["orig.ident"]] <- factor(unlist(c57.S6KO.alv.noMx[["orig.ident"]]));
cLevs <- unique(as.character(unlist(c57.S6KO.alv.noMx[["cluster"]])));
cLevs <- cLevs[order(as.numeric(cLevs))];
c57.S6KO.alv.noMx[["cluster"]] <- factor(unlist(c57.S6KO.alv.noMx[["cluster"]]), levels=cLevs);
c57.S6KO.alv.noMx[["clusterREV"]] <- factor(unlist(c57.S6KO.alv.noMx[["cluster"]]), levels=rev(cLevs));
## c57.S6KO.7B.noMx <-
##     subset(c57.S6KO.7B.filt, subset=orig.ident %in% c("STAT6KO", "C57"));
## ## Refactor clusters
## c57.S6KO.7B.noMx[["orig.ident"]] <- factor(unlist(c57.S6KO.7B.noMx[["orig.ident"]]));

## Redo cluster counting
tc.alv <- table(unlist(c57.S6KO.alv.noMx[["orig.ident"]]),
                unlist(c57.S6KO.alv.noMx[["cluster"]]));
write_csv(tc.alv %>% as.data.frame(),
          sprintf("ClusterCounts_filt_SCT_alv_withoutStat6_%s.csv", Sys.Date()));
write_csv(round(t(t(tc.alv) / colSums(tc.alv)), 2) %>% as.data.frame(),
          sprintf("ClusterProps_filt_SCT_alv_withoutStat6_%s.csv", Sys.Date()));
abs(round(t(t(tc.alv) / colSums(tc.alv)), 2) - .5);

## ## Redo cluster counting
## tc.7B <- table(unlist(c57.S6KO.7B.noMx[["orig.ident"]]),
##                 unlist(c57.S6KO.7B.noMx[["cluster"]]));
## tc.7B;
## round(t(t(tc.7B) / colSums(tc.7B)), 2);
## abs(round(t(t(tc.7B) / colSums(tc.7B)), 2) - .5);

## ## Carry out cluster-based WT vs STAT6 DE testing (Seven Bridges)
## {
##     marker.list.7B <- NULL;
##     Idents(c57.S6KO.7B.noMx) <- c57.S6KO.7B.noMx[["cluster"]];
##     for(cli in unique(unlist(c57.S6KO.7B.noMx[["cluster"]]))){
##         cat(cli,"\n");
##         test.markers <- FindMarkers(c57.S6KO.7B.noMx, test.use="DESeq2",
##                                     subset.ident=cli,
##                                     group.by="orig.ident",
##                                     ident.1="C57", ident.2="STAT6KO",);
##         test.markers <- test.markers[order(-test.markers$avg_log2FC),];
##         test.markers <- test.markers[!is.na(test.markers$p_val_adj),];
##         test.markers <- test.markers[test.markers$p_val_adj < 0.1,];
##         markersToAdd <- test.markers; #unique(rbind(head(test.markers, 20), tail(test.markers, 20)));
##         if(nrow(markersToAdd) > 0){
##             marker.list.7B <- rbind(marker.list.7B,
##                                     cbind(cluster=cli, marker=rownames(markersToAdd),
##                                           markersToAdd));
##         }
##     }
##     write_csv(marker.list.7B,
##               file=sprintf("ClusterDE_WT-STAT6KO_noMx_7B_%s.csv",
##                            format(Sys.Date())));
## }

## Carry out cluster-based WT vs STAT6 DE testing (Alevin)
{
    marker.list.alv <- NULL;
    marker.list.all.alv <- NULL;
    Idents(c57.S6KO.alv.noMx) <- c57.S6KO.alv.noMx[["cluster"]];
    for(cli in unique(unlist(c57.S6KO.alv.noMx[["cluster"]]))){
        cat(cli,"\n");
        test.markers <- FindMarkers(c57.S6KO.alv.noMx, test.use="DESeq2",
                                    subset.ident=cli,
                                    group.by="orig.ident",
                                    ident.1="Stat6KO", ident.2="WT");
        test.markers <- test.markers[order(-test.markers$avg_log2FC),];
        marker.list.all.alv <- rbind(marker.list.all.alv, cbind(gene=rownames(test.markers), cluster=cli, test.markers));
        test.markers <- test.markers[!is.na(test.markers$p_val_adj),];
        test.markers <- test.markers[test.markers$p_val_adj < 0.1,];
        markersToAdd <- test.markers; #unique(rbind(head(test.markers, 20), tail(test.markers, 20)));
        if(nrow(markersToAdd) > 0){
            marker.list.alv <- rbind(marker.list.alv,
                                     cbind(cluster=cli, marker=rownames(markersToAdd),
                                           markersToAdd));
        }
    }
    write_csv(marker.list.alv,
              file=sprintf("ClusterDE_Stat6KO-WT_filt_SCT_alv_withoutStat6_%s.csv",
                           format(Sys.Date())));
    write_csv(marker.list.all.alv,
              file=sprintf("ClusterDE_all_Stat6KO-WT_filt_SCT_alv_withoutStat6_%s.csv",
                           format(Sys.Date())));
}

## Carry out cluster-based WT vs STAT6 DE testing (Alevin) with OL labels
{
    marker.list.alv <- NULL;
    marker.list.all.alv <- NULL;
    Idents(c57.S6KO.alv.noMx) <- c57.S6KO.alv.noMx[["OLcluster"]];
    for(cli in unique(unlist(c57.S6KO.alv.noMx[["OLcluster"]]))){
        cat(cli,"\n");
        test.markers <- FindMarkers(c57.S6KO.alv.noMx, test.use="DESeq2",
                                    subset.ident=cli,
                                    group.by="orig.ident",
                                    ident.1="Stat6KO", ident.2="WT");
        test.markers <- test.markers[order(-test.markers$avg_log2FC),];
        marker.list.all.alv <- rbind(marker.list.all.alv, cbind(gene=rownames(test.markers), cluster=cli, test.markers));
        test.markers <- test.markers[!is.na(test.markers$p_val_adj),];
        test.markers <- test.markers[test.markers$p_val_adj < 0.1,];
        markersToAdd <- test.markers; #unique(rbind(head(test.markers, 20), tail(test.markers, 20)));
        if(nrow(markersToAdd) > 0){
            marker.list.alv <- rbind(marker.list.alv,
                                     cbind(cluster=cli, marker=rownames(markersToAdd),
                                           markersToAdd));
        }
    }
    write_csv(marker.list.alv,
              file=sprintf("ClusterDE_OLLabels_STAT6KO-WT_noMx_alv_%s.csv",
                           format(Sys.Date())));
    write_csv(marker.list.all.alv,
              file=sprintf("ClusterDE_all_OLLabels_STAT6KO-WT_noMx_alv_%s.csv",
                           format(Sys.Date())));
}


## Carry out whole population WT vs STAT6 DE testing (Alevin)
{
    marker.list.alv <- NULL;
    marker.list.all.alv <- NULL;
    test.markers <- FindMarkers(c57.S6KO.alv.noMx, test.use="DESeq2",
                                group.by="orig.ident",
                                ident.1="Stat6KO", ident.2="WT");
    test.markers <- test.markers[order(-test.markers$avg_log2FC),];
    marker.list.all.alv <- rbind(marker.list.all.alv, cbind(gene=rownames(test.markers), cluster=cli, test.markers));
    test.markers <- test.markers[!is.na(test.markers$p_val_adj),];
    test.markers <- test.markers[test.markers$p_val_adj < 0.1,];
    markersToAdd <- test.markers; #unique(rbind(head(test.markers, 20), tail(test.markers, 20)));
    if(nrow(markersToAdd) > 0){
        marker.list.alv <- rbind(marker.list.alv,
                                 cbind(cluster=cli, marker=rownames(markersToAdd),
                                       markersToAdd));
    }
    write_csv(marker.list.alv,
              file=sprintf("WholeGroupDE_Stat6KO-WT_filt_SCT_alv_withoutStat6_%s.csv",
                           format(Sys.Date())));
    write_csv(marker.list.all.alv,
              file=sprintf("WholeGroupDE_all_Stat6KO-WT_filt_SCT_alv_withoutStat6_%s.csv",
                           format(Sys.Date())));
}

## Write data to file
##saveRDS(c57.S6KO.alv.noMx, file="FR_2.3k_320epochs_noStat6_Alevin.rds");

c57.S6KO.alv.noMx <- readRDS(file="FR_2.3k_320epochs_noStat6_Alevin.rds");

## Feature plots for Ccl22, split by cell type

## ccl22.gene <- rownames(c57.S6KO.alv.noMx)[which(sub("-.*$", "", rownames(c57.S6KO.alv.noMx)) %in% c("Ccl22", "Sned1"))];
## print(FeaturePlot(c57.S6KO.alv.noMx, features=ccl22.gene, split.by="orig.ident"));

## Gene list plots, Alevin
geneLists <- read_excel("List of genes associated with WT and STAT6 cells.xlsx");

c57.S6KO.alv.noMx[["shape"]] <- factor(16);

for(listName in unique(geneLists$Category)){
    cat(listName, "\n");
    listSub <- geneLists$Gene[geneLists$Category == listName];
    geneIndex <- match(listSub, sub("-.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    listName <- gsub(" ", "-", listName);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    cairo_pdf(sprintf("FeaturePlot_%s_alv_withoutStat6_%s.pdf", listName, Sys.Date()), onefile=TRUE,
              width=8, height=11, pointsize=8);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
                                      split.by="orig.ident", features=gene,
                                      order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        resPre <- res;
        if(length(res) < 8){
            plotsToDo <- (8 - length(res)) / 2;
            plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
            plotsToAdd <- lapply(plotsToAdd, function(x){
                for(n in names(x$labels)){
                    x$labels[[n]] <- "";
                }
                x + xlab("") + ylab("") + guides(col=FALSE) +
                    scale_shape_manual(values=NA) +
                    ggtitle("") +
                    theme(panel.border=element_blank(),
                          axis.line=element_blank(),
                          axis.text=element_blank(),
                          axis.text.y.right=element_blank(),
                          axis.title.y.right=element_blank(),
                          axis.ticks=element_blank());
            });
            for(plotSeq in seq_len(plotsToDo)){
                res <- c(res, plotsToAdd);
           }
        }
        suppressWarnings(print(wrap_plots(res, ncol=2)));
    }
    invisible(dev.off());
}

## DefaultAssay(c57.S6KO.7B.noMx) <- "SCT";
## for(listName in unique(geneLists$Category)){
##     cat(listName, "\n");
##     listSub <- geneLists$Gene[geneLists$Category == listName];
##     geneIndex <- match(listSub, sub("-.*$", "", rownames(c57.S6KO.7B.noMx)));
##     listSub <- na.omit(rownames(c57.S6KO.7B.noMx)[geneIndex]);
##     listName <- gsub(" ", "-", listName);
##     selectPoss <- seq(1, length(listSub), by=4);
##     maxRes <- NULL;
##     cairo_pdf(sprintf("FeaturePlot_OL3_%s_7B_NoMx_%s.pdf", listName, Sys.Date()), onefile=TRUE,
##               width=8, height=11, pointsize=8);
##     for(listStart in selectPoss){
##         cat(listStart, "\n");
##         geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
##         res <- list();
##         for(gene in geneSubset){
##             res <- c(res, FeaturePlot(c57.S6KO.7B.noMx, col=c("lightgrey", "#e31836"),
##                                       split.by="orig.ident", features=gene,
##                                       order=TRUE, pt.size=1, shape.by="shape",
##                                       combine=FALSE));
##         }
##         res <- lapply(res, function(x){
##             x$labels$title <- sub("-.*$", "", x$labels$title);
##             x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
##         });
##         resPre <- res;
##         if(length(res) < 8){
##             plotsToDo <- (8 - length(res)) / 2;
##             plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
##             plotsToAdd <- lapply(plotsToAdd, function(x){
##                 for(n in names(x$labels)){
##                     x$labels[[n]] <- "";
##                 }
##                 x + xlab("") + ylab("") + guides(col=FALSE) +
##                     scale_shape_manual(values=NA) +
##                     ggtitle("") +
##                     theme(panel.border=element_blank(),
##                           axis.line=element_blank(),
##                           axis.text=element_blank(),
##                           axis.text.y.right=element_blank(),
##                           axis.title.y.right=element_blank(),
##                           axis.ticks=element_blank());
##             });
##             for(plotSeq in seq_len(plotsToDo)){
##                 res <- c(res, plotsToAdd);
##            }
##         }
##         suppressWarnings(print(wrap_plots(res, ncol=2)));
##     }
##     invisible(dev.off());
## }

## pdf("out.pdf");
## DotPlot(c57.S6KO.alv.noMx,
##         split.by="orig.ident", features=na.omit(listSub));
## invisible(dev.off());

## Redo cluster-differentiating markers without outlier clusters
{
    clustersToCheck <- as.character(unique(unlist(c57.S6KO.alv.noMx[["cluster"]])));
    marker.list <- NULL;
    marker.list.all <- NULL;
    for(id1 in clustersToCheck){
        cat(id1,"\n");
        test.markers <- FindMarkers(c57.S6KO.alv.noMx, group.by="cluster", test.use="DESeq2",
                                    ident.1 = id1, ident.2 = setdiff(clustersToCheck, id1));
        test.markers <- test.markers[order(test.markers$avg_log2FC),];
        marker.list <- rbind(marker.list, cbind(gene=head(rownames(test.markers), 10), cluster=id1, head(test.markers, 10)));
        marker.list <- rbind(marker.list, cbind(gene=tail(rownames(test.markers), 10), cluster=id1, tail(test.markers, 10)));
        marker.list.all <- rbind(marker.list.all, cbind(gene=rownames(test.markers), cluster=id1, test.markers));
    }
    marker.list.adjFiltered <- unique(subset(marker.list, p_val_adj < 0.1));
    marker.list.adjFiltered <- marker.list.adjFiltered[rownames(marker.list.adjFiltered) %in% rownames(c57.S6KO.alv.noMx),];
    write.csv(marker.list.adjFiltered,
              file=sprintf("ClusterDifferentiatingMarkers_filt_SCT_alv_withoutStat6_%s.csv",
                           format(Sys.Date())))
    write.csv(marker.list.all,
              file=sprintf("ClusterDifferentiatingMarkers_all_filt_SCT_alv_withoutStat6_%s.csv",
                           format(Sys.Date())))
    DoHeatmap(c57.S6KO.alv.noMx, features=rownames(marker.list.adjFiltered), slot="data", group.by="cluster") +
        scale_fill_viridis();
    ggsave(sprintf("Heatmap_byCluster_clusterDifferentiating_filt_SCT_alv_withoutStat6_%s.png", Sys.Date()), width=11, height=16);
}


## ## Redo cluster-differentiating markers without outlier clusters
## {
##     clustersToCheck <- as.character(unique(unlist(c57.S6KO.alv.noMx[["cluster"]])));
##     for(groupName in c("WT", "Stat6KO")){
##         marker.list <- NULL;
##         marker.list.all <- NULL;
##         c57.S6KO.alv.sub <- subset(c57.S6KO.alv.noMx, subset=orig.ident == groupName);
##         for(id1 in clustersToCheck){
##             cat(id1,"\n");
##             test.markers <- FindMarkers(c57.S6KO.alv.sub, group.by="cluster", test.use="DESeq2",
##                                         ident.1 = id1, ident.2 = setdiff(clustersToCheck, id1));
##             test.markers <- test.markers[order(test.markers$avg_log2FC),];
##             marker.list <- rbind(marker.list, cbind(gene=head(rownames(test.markers), 10), cluster=id1, head(test.markers, 10)));
##             marker.list <- rbind(marker.list, cbind(gene=tail(rownames(test.markers), 10), cluster=id1, tail(test.markers, 10)));
##             marker.list.all <- rbind(marker.list.all, cbind(gene=rownames(test.markers), cluster=id1, test.markers));
##         }
##         marker.list.adjFiltered <- unique(subset(marker.list, p_val_adj < 0.1));
##         marker.list.adjFiltered <- marker.list.adjFiltered[rownames(marker.list.adjFiltered) %in% rownames(c57.S6KO.alv.sub),];
##         write.csv(marker.list.adjFiltered,
##                   file=sprintf("ClusterDifferentiatingMarkers_%s_noMx_no7_alevin_%s.csv",
##                                groupName, format(Sys.Date())))
##         write.csv(marker.list.all,
##                   file=sprintf("ClusterDifferentiatingMarkers_%s_all_noMx_no7_alevin_%s.csv",
##                                groupName, format(Sys.Date())))
##         DoHeatmap(c57.S6KO.alv.sub, features=rownames(marker.list.adjFiltered), slot="data", group.by="cluster") +
##             scale_fill_viridis();
##         ggsave(sprintf("Heatmap_byCluster_clusterDifferentiating_%s-only_M23_alv_noMx_no7_SCT_%s.png", groupName, Sys.Date()), width=11, height=8);
##     }
## }




## FeaturePlot for DE genes and cluster-differentiating genes

{ ## Alevin - DE Genes
    listSub <-
        read_csv("WholeGroupDE_Stat6KO-WT_filt_SCT_alv_withoutStat6_2021-03-25.csv") %>%
        pull(marker) %>% sort() %>% unique();
    selectPoss <- seq(1,length(listSub), by=4);
    maxRes <- NULL;
    cairo_pdf(sprintf("FeaturePlot_WholeGroupDEGs_STAT6KO-WT_alv_withoutStat6_%s.pdf",
                      Sys.Date()), onefile=TRUE, width=8, height=11, pointsize=8);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
                                      split.by="orig.ident", features=gene,
                                      order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        resPre <- res;
        if(length(res) < 8){
            plotsToDo <- (8 - length(res)) / 2;
            plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
            plotsToAdd <- lapply(plotsToAdd, function(x){
                for(n in names(x$labels)){
                    x$labels[[n]] <- "";
                }
                x + xlab("") + ylab("") +
                    guides(col=FALSE) + scale_shape_manual(values=NA) + ggtitle("") +
                    theme(panel.border=element_blank(), axis.line=element_blank(),
                          axis.text=element_blank(), axis.text.y.right=element_blank(),
                          axis.title.y.right=element_blank(), axis.ticks=element_blank());
            });
            for(plotSeq in seq_len(plotsToDo)){
                res <- c(res, plotsToAdd);
            }
        }
        suppressWarnings(print(wrap_plots(res, ncol=2)));
    }
    invisible(dev.off());
}

{ ## Alevin - Cluster-differentiating genes
    listSub <-
        read_csv("ClusterDifferentiatingMarkers_filt_SCT_alv_withoutStat6_2021-03-25.csv") %>%
        arrange(as.numeric(cluster), abs(avg_log2FC)) %>%
        pull(gene) %>% sort() %>% unique();
    selectPoss <- seq(1,length(listSub), by=4);
    maxRes <- NULL;
    cairo_pdf(sprintf("FeaturePlot_ClusterDMs_SCT_alv_withoutStat6_%s.pdf",
                      Sys.Date()), onefile=TRUE, width=8, height=11, pointsize=8);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
                                      split.by="orig.ident", features=gene,
                                      order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        resPre <- res;
        if(length(res) < 8){
            plotsToDo <- (8 - length(res)) / 2;
            plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
            plotsToAdd <- lapply(plotsToAdd, function(x){
                for(n in names(x$labels)){
                    x$labels[[n]] <- "";
                }
                x + xlab("") + ylab("") +
                    guides(col=FALSE) + scale_shape_manual(values=NA) + ggtitle("") +
                    theme(panel.border=element_blank(), axis.line=element_blank(),
                          axis.text=element_blank(), axis.text.y.right=element_blank(),
                          axis.title.y.right=element_blank(), axis.ticks=element_blank());
            });
            for(plotSeq in seq_len(plotsToDo)){
                res <- c(res, plotsToAdd);
            }
        }
        suppressWarnings(print(wrap_plots(res, ncol=2)));
    }
    invisible(dev.off());
}


## Volcano and MA plots for Stat6 vs WT for each cluster
marker.list.withGeneNames <- read_csv("ClusterDE_all_Stat6KO-WT_filt_SCT_alv_withoutStat6_2021-03-25.csv");

count.mat <- GetAssayData(c57.S6KO.alv.noMx, "counts", "RNA");


geneLogSums <- log2(rowSums(GetAssayData(c57.S6KO.alv.noMx, "counts", "RNA")));

marker.list.withGeneNames$logSum <-geneLogSums[marker.list.withGeneNames$gene];
marker.list.withGeneNames$shortName <- sub("-ENSM.*$", "", marker.list.withGeneNames$gene);

##marker.list.MAdata$p_val_adj <- (marker.list.MAdata$p_val_adj < 0.05);

## ## plot to include p < 0.1 OR log2FC > 1
## for(cl in unique(marker.list.withGeneNames$cluster)){
##     cat(cl,"\n");
##     marker.list.withGeneNames %>%
##         filter(cluster == cl) %>%
##         ggplot() +
##         aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
##         geom_point() +
##         geom_text_repel(data=marker.list.withGeneNames %>% filter(cluster == cl, (p_val_adj < 0.1) | abs(avg_log2FC) > 1))
##     ggsave(sprintf("MAplot_M23_ClusterDE_STAT6KO-WT_%s_SCT_alv_withoutStat6_%s.png", cl, Sys.Date()), width=11, height=8);
##     marker.list.withGeneNames %>%
##         filter(cluster == cl) %>%
##         ggplot() +
##         aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
##         scale_y_log10() +
##         geom_point() +
##         geom_text_repel(data=marker.list.withGeneNames %>% filter(cluster == cl, (p_val_adj < 0.1) | abs(avg_log2FC) > 1))
##     ggsave(sprintf("VolcanoPlot_M23_ClusterDE_STAT6KO-WT_%s_SCT_alv_withoutStat6_%s.png", cl, Sys.Date()), width=11, height=8);
## }

## plot to include p < 0.05 AND abs(log2FC) > 0.59 [log2(1.50247)]
for(cl in unique(marker.list.withGeneNames$cluster)){
    cat(cl,"\n");
    marker.list.withGeneNames %>%
        filter(cluster == cl) %>%
        ggplot() +
        aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
        geom_point() +
        geom_text_repel(data=marker.list.withGeneNames %>% filter(cluster == cl, (p_val_adj < 0.05) & abs(avg_log2FC) > 0.59))
    ggsave(sprintf("MAplot_M23_ClusterDE_0.05_0.59_STAT6KO-WT_%s_SCT_alv_withoutStat6_%s.png", cl, Sys.Date()), width=11, height=8);
    marker.list.withGeneNames %>%
        filter(cluster == cl) %>%
        ggplot() +
        aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
        scale_y_log10() +
        geom_point() +
        geom_text_repel(data=marker.list.withGeneNames %>% filter(cluster == cl, (p_val_adj < 0.05) & abs(avg_log2FC) > 0.59))
    ggsave(sprintf("VolcanoPlot_M23_ClusterDE_0.05_0.59_STAT6KO-WT_%s_SCT_alv_withoutStat6_%s.png", cl, Sys.Date()), width=11, height=8);
}

## { ## SevenBridges - Cluster DE genes
##     listSub <- read_csv("ClusterDE_WT-STAT6KO_noMx_7B_2021-02-22.csv") %>%
##         pull(marker) %>% sort() %>% unique();
##     selectPoss <- seq(1, length(listSub), by=4);
##     maxRes <- NULL;
##     cairo_pdf(sprintf("FeaturePlot_ClusterDEGs_7B_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
##               width=8, height=11, pointsize=8);
##     for(listStart in selectPoss){
##         cat(listStart, "\n");
##         geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
##         res <- list();
##         for(gene in geneSubset){
##             res <- c(res, FeaturePlot(c57.S6KO.7B.noMx, col=c("lightgrey", "#e31836"),
##                                       split.by="orig.ident", features=gene,
##                                       order=TRUE, pt.size=1, shape.by="shape",
##                                       combine=FALSE));
##         }
##         res <- lapply(res, function(x){
##             x$labels$title <- sub("-.*$", "", x$labels$title);
##             x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
##         });
##         resPre <- res;
##         if(length(res) < 8){
##             plotsToDo <- (8 - length(res)) / 2;
##             plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
##             plotsToAdd <- lapply(plotsToAdd, function(x){
##                 for(n in names(x$labels)){
##                     x$labels[[n]] <- "";
##                 }
##                 x + xlab("") + ylab("") + guides(col=FALSE) +
##                     scale_shape_manual(values=NA) +
##                     ggtitle("") +
##                     theme(panel.border=element_blank(),
##                           axis.line=element_blank(),
##                           axis.text=element_blank(),
##                           axis.text.y.right=element_blank(),
##                           axis.title.y.right=element_blank(),
##                           axis.ticks=element_blank());
##             });
##             for(plotSeq in seq_len(plotsToDo)){
##                 res <- c(res, plotsToAdd);
##             }
##         }
##         suppressWarnings(print(wrap_plots(res, ncol=2)));
##     }
##     invisible(dev.off());
## }

## { ## SevenBridges - ClusterDifferentiating markers
##     listSub <-
##         read_csv("ClusterDifferentiatingMarkers_all_7B_2021-02-19.csv") %>%
##         pull(X1) %>% sort() %>% unique();
##     geneIndex <- match(listSub, rownames(c57.S6KO.7B.noMx));
##     listSub <- na.omit(rownames(c57.S6KO.7B.noMx)[geneIndex]);
##     selectPoss <- seq(1,length(listSub), by=4);
##     maxRes <- NULL;
##     cairo_pdf(sprintf("FeaturePlot_ClusterDifferentiatingMarkers_7B_NoMx_%s.pdf",
##                       Sys.Date()), onefile=TRUE, width=8, height=11, pointsize=8);
##     for(listStart in selectPoss){
##         cat(listStart, "\n");
##         geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
##         res <- list();
##         for(gene in geneSubset){
##             res <- c(res, FeaturePlot(c57.S6KO.7B.noMx, col=c("lightgrey", "#e31836"),
##                                       split.by="orig.ident", features=gene,
##                                       order=TRUE, pt.size=1, shape.by="shape",
##                                       combine=FALSE));
##         }
##         res <- lapply(res, function(x){
##             x$labels$title <- sub("-.*$", "", x$labels$title);
##             x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
##         });
##         resPre <- res;
##         if(length(res) < 8){
##             plotsToDo <- (8 - length(res)) / 2;
##             plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
##             plotsToAdd <- lapply(plotsToAdd, function(x){
##                 for(n in names(x$labels)){
##                     x$labels[[n]] <- "";
##                 }
##                 x + xlab("") + ylab("") +
##                     guides(col=FALSE) + scale_shape_manual(values=NA) + ggtitle("") +
##                     theme(panel.border=element_blank(), axis.line=element_blank(),
##                           axis.text=element_blank(), axis.text.y.right=element_blank(),
##                           axis.title.y.right=element_blank(), axis.ticks=element_blank());
##             });
##             for(plotSeq in seq_len(plotsToDo)){
##                 res <- c(res, plotsToAdd);
##             }
##         }
##         suppressWarnings(print(wrap_plots(res, ncol=2)));
##     }
##     invisible(dev.off());
## }

## grep("Tc", rownames(c57.S6KO.alv.noMx), value=TRUE);


## { ## Alevin / Ighm + Cd79a
##     listSub <- c("Ighmbp2", "Cd79a", "Ms4a1", "Trbc2", "Ms4a4b", "Tcra", "Tcrb", "Tcrd", "Cd274", "Zap70");
##     geneIndex <- match(listSub, sub("-.*$", "", rownames(c57.S6KO.alv.noMx)));
##     listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
##     selectPoss <- seq(1, length(listSub), by=4);
##     maxRes <- NULL;
##     cairo_pdf(sprintf("FeaturePlot_OL4_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
##               width=8, height=11, pointsize=8);
##     for(listStart in selectPoss){
##         cat(listStart, "\n");
##         geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
##         res <- list();
##         for(gene in geneSubset){
##             res <- c(res, FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
##                                       split.by="orig.ident", features=gene,
##                                       order=TRUE, pt.size=1, shape.by="shape",
##                                       combine=FALSE));
##         }
##         res <- lapply(res, function(x){
##             x$labels$title <- sub("-.*$", "", x$labels$title);
##             x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
##         });
##         resPre <- res;
##         if(length(res) < 8){
##             plotsToDo <- (8 - length(res)) / 2;
##             plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
##             plotsToAdd <- lapply(plotsToAdd, function(x){
##                 for(n in names(x$labels)){
##                     x$labels[[n]] <- "";
##                 }
##                 x + xlab("") + ylab("") + guides(col=FALSE) +
##                     scale_shape_manual(values=NA) +
##                     ggtitle("") +
##                     theme(panel.border=element_blank(),
##                           axis.line=element_blank(),
##                           axis.text=element_blank(),
##                           axis.text.y.right=element_blank(),
##                           axis.title.y.right=element_blank(),
##                           axis.ticks=element_blank());
##             });
##             for(plotSeq in seq_len(plotsToDo)){
##                 res <- c(res, plotsToAdd);
##             }
##         }
##         suppressWarnings(print(wrap_plots(res, ncol=2)));
##     }
##     invisible(dev.off());
## }

## ## redo UMAP plots with filtered dataset (excluding outliers)
## DimPlot(c57.S6KO.alv.noMx, label = TRUE, group.by="cluster", cols=DiscretePalette(8),
##         pt.size=1.5);
## ggsave(sprintf("M23_UMAP320_byCluster_alevin_noMx_no7_SCT_%s.png", Sys.Date()), width=11, height=8);
## ggsave(sprintf("M23_UMAP320_byCluster_alevin_noMx_no7_SCT_%s.pdf", Sys.Date()), width=11, height=8);

## DimPlot(c57.S6KO.alv.noMx, label = TRUE, group.by="cluster", cols=DiscretePalette(8),
##         split.by="orig.ident", pt.size=1.5);
## ggsave(sprintf("M23_UMAP320_byCluster_split_alevin_noMx_no7_SCT_%s.png", Sys.Date()), width=16, height=8);
## ggsave(sprintf("M23_UMAP320_byCluster_split_alevin_noMx_no7_SCT_%s.pdf", Sys.Date()), width=16, height=8);


## marker.list.all <- read_csv()


## ## Cluster-differentiating Sned1 [at any level] vs Sned1-absent
## {
##     marker.list <- NULL;
##     marker.list.all <- NULL;
##     c57.S6KO.alv.WT <- subset(c57.S6KO.alv.noMx, subset=orig.ident == "WT");
##     c57.S6KO.alv.WT[["Sned1Expr"]] <- ifelse(GetAssayData(c57.S6KO.alv.WT, "counts", "RNA")["Sned1-ENSMUSG00000047793.13",] == 0, "0", "1");
##     geneLogSumsWT <- log2(rowSums(GetAssayData(c57.S6KO.alv.WT, "counts", "RNA")));
##     test.markers <- FindMarkers(c57.S6KO.alv.WT, test.use="DESeq2",
##                                 group.by="Sned1Expr", ident.1 = 0);
##     test.markers$logSum <-geneLogSumsWT[rownames(test.markers)];
##     test.markers$shortName <- sub("-ENSM.*$", "", rownames(test.markers));
##     test.markers <- test.markers[order(test.markers$avg_log2FC),];
##     marker.list <- rbind(marker.list, cbind(gene=head(rownames(test.markers), 10), head(test.markers, 10)));
##     marker.list <- rbind(marker.list, cbind(gene=tail(rownames(test.markers), 10), tail(test.markers, 10)));
##     marker.list.all <- rbind(marker.list.all, cbind(gene=rownames(test.markers), test.markers));
##     marker.list.adjFiltered <- unique(subset(marker.list, p_val_adj < 0.1));
##     marker.list.adjFiltered <- marker.list.adjFiltered[rownames(marker.list.adjFiltered) %in% rownames(c57.S6KO.alv.WT),];
##     write.csv(marker.list.adjFiltered,
##               file=sprintf("Sned1Markers_WT_noMx_no7_alevin_%s.csv",
##                            format(Sys.Date())))
##     write.csv(marker.list.all,
##               file=sprintf("Sned1Markers_WT_all_noMx_no7_alevin_%s.csv",
##                            format(Sys.Date())))
##     DoHeatmap(c57.S6KO.alv.WT, features=rownames(marker.list.adjFiltered), slot="data", group.by="cluster") +
##         scale_fill_viridis();
##     ggsave(sprintf("Heatmap_byCluster_Sned1Markers_WT_alv_noMx_no7_SCT_%s.png", Sys.Date()), width=11, height=8);
##     DoHeatmap(c57.S6KO.alv.WT, features=rownames(marker.list.adjFiltered), slot="data", group.by="Sned1Expr") +
##         scale_fill_viridis();
##     ggsave(sprintf("Heatmap_bySned1_Sned1Markers_WT_alv_noMx_no7_SCT_%s.png", Sys.Date()), width=11, height=8);
##     marker.list.all %>%
##         ggplot() +
##         aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
##         geom_point() +
##         geom_text_repel(data=marker.list.all %>% filter((p_val_adj < 0.05) | abs(avg_log2FC) > log2(1.5)))
##     ggsave(sprintf("MAplot_WT_Sned1Markers_alv_noMx_no7_SCT_%s.png", Sys.Date()), width=11, height=8);
## }


## ##marker.list.withGeneNames$gene <-
## ##    rownames(c57.S6KO.alv.noMx)[match(sub("-.*$", "", marker.list.withGeneNames$gene), sub("-.*$", "", rownames(c57.S6KO.alv.noMx)))];
## marker.list.withGeneNames <- read_csv("ClusterDifferentiatingMarkers_all_noMx_no7_alevin_2021-03-03.csv");

## geneLogSums <- log2(rowSums(GetAssayData(c57.S6KO.alv.noMx, "counts", "RNA")));

## marker.list.withGeneNames$logSum <-geneLogSums[marker.list.withGeneNames$gene];
## marker.list.withGeneNames$shortName <- sub("-ENSM.*$", "", marker.list.withGeneNames$gene);

## ##marker.list.MAdata$p_val_adj <- (marker.list.MAdata$p_val_adj < 0.05);

## for(cl in unique(marker.list.withGeneNames$cluster)){
##     cat(cl,"\n");
##     marker.list.withGeneNames %>%
##         filter(cluster == cl) %>%
##         ggplot() +
##         aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
##         geom_point() +
##         geom_text_repel(data=marker.list.withGeneNames %>% filter(cluster == cl, (p_val_adj < 0.1) | abs(avg_log2FC) > 1))
##     ggsave(sprintf("MAplot_M23_alv_noMx_clusterDifferentiating_%s_no7_SCT_%s.png", cl, Sys.Date()), width=11, height=8);
##     marker.list.withGeneNames %>%
##         filter(cluster == cl) %>%
##         ggplot() +
##         aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
##         scale_y_log10() +
##         geom_point() +
##         geom_text_repel(data=marker.list.withGeneNames %>% filter(cluster == cl, (p_val_adj < 0.1) | abs(avg_log2FC) > 1))
##     ggsave(sprintf("VolcanoPlot_M23_alv_noMx_clusterDifferentiating_%s_no7_SCT_%s.png", cl, Sys.Date()), width=11, height=8);
## }

## cluster.fac <- as.character(unlist(c57.S6KO.alv.noMx[["cluster"]]));
## c57.S6KO.alv.noMx[["clusterREV"]] <- factor(cluster.fac, levels=sort(unique(cluster.fac), decreasing=TRUE));

## OLcluster.fac <- as.character(unlist(c57.S6KO.alv.noMx[["OL.name"]]));
## c57.S6KO.alv.noMx[["OLclusterREV"]] <- factor(OLcluster.fac, levels=sort(unique(OLcluster.fac), decreasing=TRUE));
## c57.S6KO.alv.noMx[["OLcluster"]] <- factor(OLcluster.fac, levels=sort(unique(OLcluster.fac), decreasing=FALSE));

## OLclusterDE.fac <- as.character(paste(unlist(c57.S6KO.alv.noMx[["OL.name"]]), unlist(c57.S6KO.alv.noMx[["orig.ident"]]), sep="-"));
## c57.S6KO.alv.noMx[["OLclusterDEREV"]] <- factor(OLclusterDE.fac, levels=sort(unique(OLclusterDE.fac), decreasing=TRUE));
## c57.S6KO.alv.noMx[["OLclusterDE"]] <- factor(OLclusterDE.fac, levels=sort(unique(OLclusterDE.fac), decreasing=FALSE));

## ## Check marker clusters
## geneList.OL5 <- c("Aldh1a2", "Ccl9", "Ccl17", "Itgax", "Fgl2", "Crispld2", "Il9r", "Sned1");
## geneList.OL5 <- rownames(c57.S6KO.alv.noMx)[match(geneList.OL5, sub("-ENSM.*", "", rownames(c57.S6KO.alv.noMx)))];

## pdf(sprintf("Dotplot_OL5_OLnamed_noMx_alv_%s.pdf", Sys.Date()));
## DefaultAssay(c57.S6KO.alv.noMx) <- "SCT";
## res <- DotPlot(c57.S6KO.alv.noMx, features=geneList.OL5,
##                group.by="OLclusterREV",
##                dot.min=0.0001, scale.by="size", scale=TRUE,
##                col.min=0, col.max=1) +
##     scale_colour_viridis() +
##     ylab("Cluster") +
##     scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
##     theme(axis.text.x=element_text(angle=45, hjust=1))
## ##res$labels$title <- sub("-.*$", "", res$labels$title);
## print(res);
## invisible(dev.off());


## ## carry out specific cluster DE
## {
##     res.all <- NULL;
##     for(cl1 in as.character(0:6)){
##         for(cl2 in as.character(0:6)){
##             if(cl2 <= cl1){
##                 next;
##             }
##             cat(sprintf("%s vs %s\n", cl1, cl2));
##             res <- NULL;
##             test.markers <- FindMarkers(c57.S6KO.alv.noMx, group.by="cluster", test.use="DESeq2",
##                                         ident.1 = cl2, ident.2 = cl1);
##             test.markers <- cbind(test=sprintf("%s-%s", cl2, cl1),
##                                             gene=rownames(test.markers),
##                                             logSum=geneLogSums[rownames(test.markers)], test.markers) %>%
##                 mutate(shortName=sub("-ENSM.*", "", gene));
##             test.markers.all <- test.markers;
##             res.all <- rbind(res.all, test.markers);
##             test.markers <- subset(test.markers, p_val_adj < 0.1);
##             test.markers <- test.markers[order(-test.markers$avg_log2FC),];
##             write_csv(test.markers, file=sprintf("ClusterDE_%svs%s_noMx_alv_%s.csv", cl2, cl1, Sys.Date()));
##             test.markers.all %>%
##                 ggplot() +
##                 aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
##                 geom_point() +
##                 geom_text_repel(data=test.markers.all %>%
##                                     filter((p_val_adj < 0.05) | abs(avg_log2FC) > log2(1.5)));
##             ggsave(sprintf("MAplot_%svs%s_noMx_alv_%s.png", cl2, cl1, Sys.Date()), width=11, height=8);
##             test.markers.all %>%
##                 ggplot() +
##                 aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
##                 scale_y_log10() +
##                 geom_point() +
##                 geom_text_repel(data=test.markers.all %>%
##                                     filter((p_val_adj < 0.05) | abs(avg_log2FC) > log2(1.5)));
##             ggsave(sprintf("VolcanoPlot_%svs%s_noMx_alv_%s.png", cl2, cl1, Sys.Date()), width=11, height=8);
##         }
##     }
##     write_csv(res.all, file=sprintf("ClusterDE_all_XvsX_noMx_alv_%s.csv", Sys.Date()));
## }

## ## mini-bulk test

## ## carry out mini-bulk RNA sequencing
## bulkCat <- as.character(unlist(c57.S6KO.alv.noMx[["cluster"]]));
## bulkCat[bulkCat %in% c("3", "6")] <- "mb36";
## bulkCat[bulkCat %in% c("0", "1", "2", "4", "5")] <- "mb01245";
## c57.S6KO.alv.noMx[["bulkCat"]] <- bulkCat;

## {
##     res <- NULL;
##     res.all <- NULL;
##     cat("mb36; mb01245\n");
##     for(bc in c("mb36", "mb01245")){
##         cat(bc,"\n");
##         c57.S6KO.alv.noMx.distSub <- subset(c57.S6KO.alv.noMx, subset=(bulkCat == bc));
##         test.markers <- FindMarkers(c57.S6KO.alv.noMx.distSub, group.by="orig.ident", test.use="DESeq2",
##                                 ident.1 = "Stat6KO", ident.2 = "c57");
##         res.all <- rbind(res.all, cbind(subCluster=bc, test="Stat6KO-c57",
##                                 gene=rownames(test.markers), logSum=geneLogSums[rownames(test.markers)], test.markers));
##         test.markers <- subset(test.markers, p_val_adj < 0.1);
##         test.markers <- test.markers[order(-test.markers$avg_log2FC),];
##         res <- rbind(res, cbind(subCluster=bc, test="Stat6KO-c57",
##                                 gene=rownames(test.markers), logSum=geneLogSums[rownames(test.markers)], test.markers));
##     }
##     cat("Stat6KO; c57\n");
##     for(bc in c("Stat6KO", "c57")){
##         cat(bc, "\n");
##         c57.S6KO.alv.noMx.distSub <- subset(c57.S6KO.alv.noMx, (orig.ident == bc));
##         test.markers <- FindMarkers(c57.S6KO.alv.noMx.distSub, group.by="bulkCat", test.use="DESeq2",
##                                     ident.1 = "mb01245", ident.2 = "mb36");
##         res.all <- rbind(res.all, cbind(subCluster=bc, test="mb01245-mb36",
##                                 gene=rownames(test.markers), logSum=geneLogSums[rownames(test.markers)], test.markers));
##         test.markers <- subset(test.markers, p_val_adj < 0.1);
##         test.markers <- test.markers[order(-test.markers$avg_log2FC),];
##         res <- rbind(res, cbind(subCluster=bc, test="mb01245-mb36",
##                                 gene=rownames(test.markers), logSum=geneLogSums[rownames(test.markers)], test.markers));
##     }
##     write_csv(res, file=sprintf("ClusterDE_miniBulk_noMx_alv_%s.csv", Sys.Date()));
##     write_csv(res.all, file=sprintf("ClusterDE_all_miniBulk_noMx_alv_%s.csv", Sys.Date()));
## }


## ## MA / Volcano plots for miniBulk experiment

## res.tbl <- read_csv("ClusterDE_all_miniBulk_noMx_alv_2021-03-03.csv") %>%
##     mutate(shortName = sub("-ENSM.*", "", gene));
## for(sub in unique(res.tbl$subCluster)){
##     for(testName in unique(res.tbl$test)){
##         res.filt <- res.tbl %>% filter(subCluster == sub, test == testName);
##         if(nrow(res.filt) == 0){
##             next;
##         }
##         cat(sub, testName,"\n");
##         res.filt %>%
##             ggplot() +
##             aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
##             geom_point() +
##             geom_text_repel(data=res.filt %>%
##                                 filter((p_val_adj < 0.1) | abs(avg_log2FC) > 1))
##         ggsave(sprintf("MAplot_%s_%s_alv_noMx_no7_SCT_%s.png", sub, testName,
##                        Sys.Date()), width=11, height=8);
##         res.filt %>%
##             ggplot() +
##             aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
##             scale_y_log10() +
##             geom_point() +
##             geom_text_repel(data=res.filt %>%
##                                 filter((p_val_adj < 0.1) | abs(avg_log2FC) > 1))
##         ggsave(sprintf("VolcanoPlot_%s_%s_alv_noMx_no7_SCT_%s.png", sub, testName,
##                        Sys.Date()), width=11, height=8);
##     }
## }


## res.tbl <- read_csv("ClusterDE_all_miniBulk_noMx_alv_2021-03-03.csv") %>%
##     mutate(shortName = sub("-ENSM.*", "", gene));
## for(sub in unique(res.tbl$subCluster)){
##     for(testName in unique(res.tbl$test)){
##         res.filt <- res.tbl %>% filter(subCluster == sub, test == testName);
##         if(nrow(res.filt) == 0){
##             next;
##         }
##         cat(sub, testName,"\n");
##         res.filt %>%
##             ggplot() +
##             aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
##             geom_point()
##         ggsave(sprintf("nolabel/MAplot_noLabel_%s_%s_alv_noMx_no7_SCT_%s.png", sub, testName,
##                        Sys.Date()), width=11, height=8);
##         res.filt %>%
##             ggplot() +
##             aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
##             scale_y_log10() +
##             geom_point() +
##         ggsave(sprintf("nolabel/VolcanoPlot_noLabel_%s_%s_alv_noMx_no7_SCT_%s.png", sub, testName,
##                        Sys.Date()), width=11, height=8);
##     }
## }

## for(test in c("1vs4", "1vs6", "3vs6", "2vs5")){
##     cat(test, "\n");
##     if(!file.exists(sprintf("ClusterDE_all_%s_noMx_alv_2021-03-11.csv", test))){
##         next;
##     }
##     res.df <- read_csv(sprintf("ClusterDE_all_%s_noMx_alv_2021-03-11.csv", test)) %>%
##         mutate(shortName = sub("-ENSM.*", "", gene));
##     res.df %>%
##         ggplot() +
##         aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
##         geom_point() +
##         geom_text_repel(data=res.df %>%
##                             filter((p_val_adj < 0.1) | abs(avg_log2FC) > 1))
##     ggsave(sprintf("MAplot_%s_noMx_alv_%s.png", test, Sys.Date()), width=11, height=8);
##     res.df %>%
##         ggplot() +
##         aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
##         scale_y_log10() +
##         geom_point() +
##         geom_text_repel(data=res.df %>%
##                             filter((p_val_adj < 0.1) | abs(avg_log2FC) > 1))
##     ggsave(sprintf("VolcanoPlot_%s_noMx_alv_%s.png", test, Sys.Date()), width=11, height=8);
## }

## count.matrix <- GetAssayData(c57.S6KO.alv.noMx, "counts", "RNA");
## grep("Malat1", rownames(count.matrix), value=TRUE);
## grep("Ccng2", rownames(count.matrix), value=TRUE);


## grep("Ccl5", rownames(count.matrix), value=TRUE);
## grep("Ccl22", rownames(count.matrix), value=TRUE);

## grep("Il13ra1", rownames(count.matrix), value=TRUE);
## grep("Il4ra", rownames(count.matrix), value=TRUE);


## png(sprintf("Il13ra1_vs_Il4ra_%s.png", Sys.Date()), width=1200, height=1200);
## count.matrix <- GetAssayData(c57.S6KO.alv.noMx, "counts", "RNA");
## plot(jitter(count.matrix["Il13ra1-ENSMUSG00000017057.9",]),
##      jitter(count.matrix["Il4ra-ENSMUSG00000030748.9",]));
## invisible(dev.off());

## png(sprintf("Ccl5_vs_Ccl22_%s.png", Sys.Date()), width=1200, height=1200);
## count.matrix <- GetAssayData(c57.S6KO.alv.noMx, "counts", "RNA");
## plot(count.matrix["Ccl5-ENSMUSG00000035042.2",],
##      count.matrix["Ccl22-ENSMUSG00000031779.3",]);
## invisible(dev.off());

## cor.test(count.matrix["Ccl5-ENSMUSG00000035042.2",],
##      count.matrix["Ccl22-ENSMUSG00000031779.3",]);


## png(sprintf("Malat1_vs_percent.mt_orig_%s.png", Sys.Date()), width=1200, height=1200);
## count.matrix <- GetAssayData(c57.S6KO.alv.all, "counts", "RNA");
## plot(count.matrix["Malat1-ENSMUSG00000092341.3",],
##      unlist(c57.S6KO.alv.all[["percent.mt"]]))
## invisible(dev.off());

## png("Malat1_vs_percent.mt_2021-03-03.png", width=1200, height=1200);
## plot(count.matrix["Malat1-ENSMUSG00000092341.3",],
##      unlist(c57.S6KO.alv.noMx[["percent.mt"]]))
## invisible(dev.off());

## png("Ccng2_vs_percent.mt_2021-03-03.png", width=1200, height=1200);
## plot(count.matrix["Ccng2-ENSMUSG00000029385.14",],
##      unlist(c57.S6KO.alv.noMx[["percent.mt"]]))
## invisible(dev.off());

## png("Malat1_vs_transcript_counts_2021-03-03.png", width=1200, height=1200);
## plot(count.matrix["Malat1-ENSMUSG00000092341.3",],
##      colSums(count.matrix));
## invisible(dev.off());

## png(sprintf("Malat1_vs_transcript_counts_orig_%s.png", Sys.Date()), width=1200, height=1200);
## count.matrix <- GetAssayData(c57.S6KO.alv.all, "counts", "RNA");
## plot(count.matrix["Malat1-ENSMUSG00000092341.3",],
##      colSums(count.matrix));
## invisible(dev.off());


## png("transcript_counts_density_2021-03-03.png", width=1200, height=1200);
## plot(density(colSums(count.matrix)));
## invisible(dev.off());

## png(sprintf("transcript_counts_density_orig_%s.png", Sys.Date()), width=1200, height=1200);
## count.matrix <- GetAssayData(c57.S6KO.alv.all, "counts", "RNA");
## plot(density(colSums(count.matrix)));
## invisible(dev.off());

## t(table(colSums(count.matrix) < 5000, unlist(c57.S6KO.alv.noMx[["cluster"]])));


## genes.OL6 <- c("Tlr1", "Tlr2", "Tlr11", "Unc93b1", "Arc");
## geneIndex <- match(genes.OL6, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
## genes.OL6 <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);

## cells.c57 <- c57.S6KO.alv.noMx[["orig.ident"]] == "c57";
## OL6.mat.c57 <- as.matrix(GetAssayData(c57.S6KO.alv.noMx, "counts", "RNA"))[genes.OL6,cells.c57];
## OL6.mat.Stat6KO <- as.matrix(GetAssayData(c57.S6KO.alv.noMx, "counts", "RNA"))[genes.OL6,!cells.c57];

## ## Counts
## t(apply(OL6.mat.c57,1,function(x){table(x > 0)}));
## t(apply(OL6.mat.Stat6KO,1,function(x){table(x > 0)}));

## ## Proportions
## t(apply(OL6.mat.c57,1,function(x){round(table(x > 0) / length(x), 2)}));
## t(apply(OL6.mat.Stat6KO,1,function(x){round(table(x > 0) / length(x), 2)}));

## ## DotPlot
## pdf(sprintf("Dotplot_OL6_noMx_alv_%s.pdf", Sys.Date()));
## DefaultAssay(c57.S6KO.alv.noMx) <- "SCT";
## res <- DotPlot(c57.S6KO.alv.noMx, features=genes.OL6,
##                group.by="clusterREV", split.by="orig.ident",
##                dot.min=0.0001, scale.by="size", scale=TRUE,
##                col.min=0, col.max=1) +
##     scale_colour_viridis(discrete=TRUE) +
##     ylab("Cluster") +
##     scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
##     guides(col=FALSE) +
##     theme(axis.text.x=element_text(angle=45, hjust=1))
## print(res);
## invisible(dev.off());

## { ## Alevin / Ighm + Cd79a
##     genes.OL7 <- c("Fscn1", "Tlr2", "Fgl2");
##     geneIndex <- match(genes.OL7, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
##     listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
##     selectPoss <- seq(1, length(listSub), by=4);
##     maxRes <- NULL;
##     cairo_pdf(sprintf("FeaturePlot_OL7_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
##               width=8, height=11, pointsize=8);
##     for(listStart in selectPoss){
##         cat(listStart, "\n");
##         geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
##         res <- list();
##         for(gene in geneSubset){
##             res <- c(res, FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
##                                       split.by="orig.ident", features=gene,
##                                       order=TRUE, pt.size=1, shape.by="shape",
##                                       combine=FALSE));
##         }
##         res <- lapply(res, function(x){
##             x$labels$title <- sub("-.*$", "", x$labels$title);
##             x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
##         });
##         resPre <- res;
##         if(length(res) < 8){
##             plotsToDo <- (8 - length(res)) / 2;
##             plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
##             plotsToAdd <- lapply(plotsToAdd, function(x){
##                 for(n in names(x$labels)){
##                     x$labels[[n]] <- "";
##                 }
##                 x + xlab("") + ylab("") + guides(col=FALSE) +
##                     scale_shape_manual(values=NA) +
##                     ggtitle("") +
##                     theme(panel.border=element_blank(),
##                           axis.line=element_blank(),
##                           axis.text=element_blank(),
##                           axis.text.y.right=element_blank(),
##                           axis.title.y.right=element_blank(),
##                           axis.ticks=element_blank());
##             });
##             for(plotSeq in seq_len(plotsToDo)){
##                 res <- c(res, plotsToAdd);
##             }
##         }
##         suppressWarnings(print(wrap_plots(res, ncol=2)));
##     }
##     invisible(dev.off());
## }

## ## Update annotation; c57 -> WT
## c57.S6KO.alv.noMx[["orig.ident"]] <-
##     ifelse(unlist(c57.S6KO.alv.noMx[["orig.ident"]]) == "c57", "WT",
##            as.character(unlist(c57.S6KO.alv.noMx[["orig.ident"]])));

## ## Order with WT first
## c57.S6KO.alv.noMx[["orig.ident"]]

## ## Check mapping
## table(c57.S6KO.alv.noMx[["orig.ident"]]);

## ## Cluster Tree for cluster relationships [2021-Mar-04]
## ## [https://stackoverflow.com/questions/62882606/]
## c57.S6KO.alv.noMx <- BuildClusterTree(c57.S6KO.alv.noMx,assay="SCT")

## myPhyTree <- Tool(object=c57.S6KO.alv.noMx, slot = "BuildClusterTree")
## ggtree(myPhyTree) +
## geom_tiplab()+
## theme_tree()
## ggsave(sprintf("tree_alv_NoMx_%s.pdf", Sys.Date()));


## ## Cluster Tree for cluster relationships, split by cell type [2021-Mar-11]
## ## [https://stackoverflow.com/questions/62882606/]
## for(groupName in c("WT", "Stat6KO")){
##     cat(groupName, "\n");
##     c57.S6KO.alv.sub <- subset(c57.S6KO.alv.noMx, subset=orig.ident == groupName);
##     c57.S6KO.alv.sub <- BuildClusterTree(c57.S6KO.alv.sub, assay="SCT");
##     myPhyTree <- Tool(object=c57.S6KO.alv.sub, slot = "BuildClusterTree")
##     ggtree(myPhyTree) +
##         geom_tiplab()+
##         theme_tree() +
##         ggtitle(groupName)
##     ggsave(sprintf("tree_%s_alv_NoMx_%s.pdf", groupName, Sys.Date()), width=8, height=8);
## }

## Update with Cluster Names
names.tbl <- read_xlsx("Cluster Names.xlsx")
OL.names <- factor(unlist(names.tbl[,2]), levels=sort(unlist(names.tbl[,2]), decreasing=TRUE));
names(OL.names) <- as.character(unlist(names.tbl[,1]));

c57.S6KO.alv.noMx[["OL.name"]] <-
    factor(OL.names[unlist(c57.S6KO.alv.noMx[["cluster"]])], levels=OL.names);

## redo UMAP plots with OL Names
DimPlot(c57.S6KO.alv.filt, label = TRUE, group.by="cluster", cols=brewer.pal(8, "Dark2"),
        pt.size=1.5);
ggsave(sprintf("M23_UMAP320_byOLCluster_alevin_filt_SCT_%s.png", Sys.Date()), width=11, height=8);
ggsave(sprintf("M23_UMAP320_byOLCluster_alevin_filt_SCT_%s.pdf", Sys.Date()), width=11, height=8);

## redo UMAP plots with OL Names
DimPlot(c57.S6KO.alv.noMx, label = TRUE, group.by="OL.name", cols=brewer.pal(8, "Dark2"),
        pt.size=1.5);
ggsave(sprintf("M23_UMAP320_byOLCluster_alevin_noMx_no7_SCT_%s.png", Sys.Date()), width=11, height=8);
ggsave(sprintf("M23_UMAP320_byOLCluster_alevin_noMx_no7_SCT_%s.pdf", Sys.Date()), width=11, height=8);

DimPlot(c57.S6KO.alv.noMx, label = TRUE, group.by="OL.name", cols=brewer.pal(8, "Dark2"),
        split.by="orig.ident", pt.size=1.5);
ggsave(sprintf("M23_UMAP320_byOLCluster_split_alevin_noMx_no7_SCT_%s.png", Sys.Date()), width=16, height=8);
ggsave(sprintf("M23_UMAP320_byOLCluster_split_alevin_noMx_no7_SCT_%s.pdf", Sys.Date()), width=16, height=8);

pdf("brewer_pals.pdf");
par(mar=c(3,4,2,2));
display.brewer.all(colorblindFriendly=FALSE);
invisible(dev.off());

{ ## OL-selected Genes
    genes.OL8 <- read_xlsx("Genes for David.xlsx") %>% pull(`Gene Name`);
    geneIndex <- match(genes.OL8, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    cairo_pdf(sprintf("FeaturePlot_OL8_OLnamed_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=8, height=11, pointsize=8);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
                                      split.by="orig.ident", features=gene,
                                      order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        resPre <- res;
        if(length(res) < 8){
            plotsToDo <- (8 - length(res)) / 2;
            plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
            plotsToAdd <- lapply(plotsToAdd, function(x){
                for(n in names(x$labels)){
                    x$labels[[n]] <- "";
                }
                x + xlab("") + ylab("") + guides(col=FALSE) +
                    scale_shape_manual(values=NA) +
                    ggtitle("") +
                    theme(panel.border=element_blank(),
                          axis.line=element_blank(),
                          axis.text=element_blank(),
                          axis.text.y.right=element_blank(),
                          axis.title.y.right=element_blank(),
                          axis.ticks=element_blank());
            });
            for(plotSeq in seq_len(plotsToDo)){
                res <- c(res, plotsToAdd);
            }
        }
        suppressWarnings(print(wrap_plots(res, ncol=2)));
    }
    invisible(dev.off());
}

{
    genes.OL8 <- read_xlsx("Genes for David.xlsx") %>% pull(`Gene Name`);
    geneIndex <- match(genes.OL8, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    cairo_pdf(sprintf("DotPlot_OL8_OLnamed_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=16, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
               group.by="OL.name",
               dot.min=0.0001, scale.by="size", scale=TRUE,
               col.min=0, col.max=1) +
        scale_colour_viridis() +
        ylab("Cluster") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
}

{ ## OL-selected DE Genes
    genes.OL9 <- read_xlsx("STAT6-KO vs WT DEGs all.xlsx", col_names="Gene Name") %>% pull(`Gene Name`);
    geneIndex <- match(genes.OL9, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    cairo_pdf(sprintf("FeaturePlot_OL9_S6-KO-DEGs_OLnamed_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=8, height=11, pointsize=8);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
                                      split.by="orig.ident", features=gene,
                                      order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        resPre <- res;
        if(length(res) < 8){
            plotsToDo <- (8 - length(res)) / 2;
            plotsToAdd <- list(res[[length(res) - 1]], res[[length(res)]]);
            plotsToAdd <- lapply(plotsToAdd, function(x){
                for(n in names(x$labels)){
                    x$labels[[n]] <- "";
                }
                x + xlab("") + ylab("") + guides(col=FALSE) +
                    scale_shape_manual(values=NA) +
                    ggtitle("") +
                    theme(panel.border=element_blank(),
                          axis.line=element_blank(),
                          axis.text=element_blank(),
                          axis.text.y.right=element_blank(),
                          axis.title.y.right=element_blank(),
                          axis.ticks=element_blank());
            });
            for(plotSeq in seq_len(plotsToDo)){
                res <- c(res, plotsToAdd);
            }
        }
        suppressWarnings(print(wrap_plots(res, ncol=2)));
    }
    invisible(dev.off());
}

{
    genes.OL9 <- read_xlsx("STAT6-KO vs WT DEGs all.xlsx", col_names="Gene Name") %>% pull(`Gene Name`);
    geneIndex <- match(genes.OL9, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    cairo_pdf(sprintf("DotPlot_OL9_S6-KO-DEGs_OLnamed_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=16, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
               group.by="OL.name", split.by="orig.ident",
               dot.min=0.0001, scale.by="size", scale=TRUE,
               col.min=0, col.max=1) +
        scale_colour_viridis(discrete=TRUE) +
        guides(col=FALSE) +
        ylab("Cluster") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
}

{
    genes.dotplots.OL10 <- read_xlsx("Dotplots and heatmaps for identification.xlsx");
    colnames(genes.dotplots.OL10) <- c("plot", "Genes");
    DefaultAssay(c57.S6KO.alv.noMx) <- "SCT";
    cairo_pdf(sprintf("DotPlot_OL10_identification_OLnamed_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    for(plotName in unique(genes.dotplots.OL10$plot)){
        cat(plotName, "\n");
        genes.OL10 <- genes.dotplots.OL10 %>%
            filter(plot == plotName) %>%
            pull(Genes);
        geneIndex <- match(genes.OL10, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
        listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
        res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
                       group.by="OLclusterREV",
                       dot.min=0.0001, scale.by="size", scale=TRUE,
                       col.min=0, col.max=1) +
            scale_colour_viridis() +
            ylab("Cluster") +
            ggtitle(plotName) +
            scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            theme(axis.text.x=element_text(angle=45, hjust=1))
        print(res);
    }
    invisible(dev.off());
    cairo_pdf(sprintf("Heatmap_OL10_identification_OLnamed_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    for(plotName in unique(genes.dotplots.OL10$plot)){
        cat(plotName, "\n");
        genes.OL10 <- genes.dotplots.OL10 %>%
            filter(plot == plotName) %>%
            pull(Genes);
        geneIndex <- match(genes.OL10, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
        listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
        DefaultAssay(c57.S6KO.alv.noMx) <- "SCT";
        res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub,
                         group.by="OLcluster") +
            scale_fill_viridis() +
            scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            guides(col=FALSE)
        print(res);
    }
    invisible(dev.off());
}

{
    genes.dotplots.OL10 <- read_xlsx("Dotplots and heatmaps for selected DEGs.xlsx");
    colnames(genes.dotplots.OL10) <- c("plot", "Genes");
    cairo_pdf(sprintf("DotPlot_OL10_%s_DEGs_OLnamed_alv_NoMx_%s.pdf", plotName, Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    for(plotName in unique(genes.dotplots.OL10$plot)){
        genes.OL10 <- genes.dotplots.OL10 %>%
            filter(plot == plotName) %>%
            pull(Genes);
        geneIndex <- match(genes.OL10, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
        listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
        DefaultAssay(c57.S6KO.alv.noMx) <- "RNA";
        res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
                       group.by="OLcluster", split.by="orig.ident",
                       dot.min=0.0001, scale.by="size", scale=TRUE,
                       col.min=0, col.max=1) +
            scale_colour_viridis(discrete=TRUE) +
            guides(col=FALSE) +
            ggtitle(plotName) +
            ylab("Cluster") +
            scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            theme(axis.text.x=element_text(angle=45, hjust=1))
        print(res);
    }
    invisible(dev.off());
    DefaultAssay(c57.S6KO.alv.noMx) <- "SCT";
    cairo_pdf(sprintf("Heatmap_OL10_%s_DEGs_OLnamed_alv_NoMx_%s.pdf", plotName, Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    for(plotName in unique(genes.dotplots.OL10$plot)){
        genes.OL10 <- genes.dotplots.OL10 %>%
            filter(plot == plotName) %>%
            pull(Genes);
        geneIndex <- match(genes.OL10, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
        res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub,
                         group.by=c("OLclusterDE")) +
            scale_fill_viridis() +
            scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            guides(col=FALSE)
        print(res);
    }
    invisible(dev.off());
}




## try to find the longest gene with reasonable expression
head(sort(geneLogSums, decreasing=TRUE), 20)
## Tmsb4x: 768 bp
## Actb: 1920 bp
## Tmem123: 2907 bp
## Eef1a1: 1794 bp
## Ccl5: 532 bp
## B2m: 860 bp
## Actg1: 1952 bp

## Okay, so looks like Tmem123 is the transcript to try out a variant analysis with
## Tmem123-ENSMUSG00000050912.15

## Per-cluster statistics:
nrow(c57.S6KO.alv.noMx);
nrow(c57.S6KO.alv.filt);
nrow(c57.S6KO.alv.all);

## Number of genes with positive counts
## Number of reads
## Number of transcripts
## Number of cells
## rows: Gene; columns: Cluster/Tag + Total data:Log2(total expression) per cluster
cell.reads.tbl <- read_table("/mnt/ufds/jmayer/cellLabels_allMatching_counts.txt.gz", col_names=c("count", "barcode"));
cell.reads.tbl$annotated <- colnames(c57.S6KO.alv.all)[match(cell.reads.tbl$barcode, sub("^.*_", "", colnames(c57.S6KO.alv.all)))];
cell.reads.tbl$cluster <- colnames(c57.S6KO.alv.all)[match(cell.reads.tbl$barcode, sub("^.*_", "", colnames(c57.S6KO.alv.all)))];


head(colnames(c57.S6KO.alv.all));

pdf("FeaturePlot_Irf4_noMx.pdf", width=12, height=6);
geneIndex <- match("Irf4", sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
res <- c(FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
                     split.by="orig.ident", features=listSub,
                     order=TRUE, pt.size=1, shape.by="shape",
                     combine=FALSE));
print(wrap_plots(res, ncol=2));
invisible(dev.off());

{ ## OL-selected DE Genes
    genes.OL11 <- c("Zap70", "Ms4a4b");
    geneIndex <- match(genes.OL11, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    cairo_pdf(sprintf("FeaturePlot_OL11_Zap70-Ms4a4b_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=12, height=6, pointsize=8);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=2)));
    }
    invisible(dev.off());
}

{ ## OL-selected DE Genes
    genes.OL11 <- c("Zap70", "Ms4a4b");
    geneIndex <- match(genes.OL11, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    cairo_pdf(sprintf("FeaturePlot_OL11_Zap70-Ms4a4b_alv_filt_%s.pdf", Sys.Date()), onefile=TRUE,
              width=12, height=6, pointsize=8);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.filt, col=c("lightgrey", "#e31836"),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=2)));
    }
    invisible(dev.off());
}


{ ## OL-selected DE Genes
    genes.OL12 <- c("Il4ra", "Il13ra1");
    geneIndex <- match(genes.OL12, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    cairo_pdf(sprintf("FeaturePlot_OL12_Il4_Il13_S6-KO-DEGs_OLnamed_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8, pointsize=8);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx, col=c("lightgrey", "#e31836"),
                                      split.by="orig.ident", features=gene,
                                      order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        resPre <- res;
        suppressWarnings(print(wrap_plots(res, ncol=2)));
    }
    invisible(dev.off());
}

str(rownames(GetAssayData(c57.S6KO.alv.noMx)));

## Publication figures (2021-Mar-16)

cluster.tbl <- read_xlsx('Cluster annotation for figure.xlsx');
c57.S6KO.alv.noMx[["Annotation"]] <- factor(cluster.tbl$Annotation[match(unlist(c57.S6KO.alv.noMx[["cluster"]]), cluster.tbl$Ident)],
                                            levels=unique(cluster.tbl$Annotation));
c57.S6KO.alv.noMx[["AnnotationREV"]] <- factor(cluster.tbl$Annotation[match(unlist(c57.S6KO.alv.noMx[["cluster"]]), cluster.tbl$Ident)],
                                            levels=rev(unique(cluster.tbl$Annotation)));

table(as.character(unlist(c57.S6KO.alv.noMx[["Annotation"]])), as.character(unlist(c57.S6KO.alv.noMx[["cluster"]])))

library(DescTools);
clusterCols <- brewer.pal(8, "Dark2");
names(clusterCols) <- c(unique(sort(as.character(unlist(c57.S6KO.alv.noMx[["Annotation"]])))), "Other");
clusterCols["C7"] <- clusterCols["C4"];
grey.pos <- which(names(clusterCols) %in% c("C1", "C2", "C3", "C4"));
clusterCols[grey.pos] <- c("grey20", "grey40", "grey65", "grey85");

## Redo UMAP with shaded colours
DimPlot(c57.S6KO.alv.noMx, label = TRUE, group.by="Annotation", cols=clusterCols,
        pt.size=1.5);
ggsave(sprintf("M23_UMAP320_pub_byAnnot_alevin_noMx_no7_SCT_%s.png", Sys.Date()), width=11, height=8);
ggsave(sprintf("M23_UMAP320_pub_byAnnot_alevin_noMx_no7_SCT_%s.pdf", Sys.Date()), width=11, height=8);

DimPlot(c57.S6KO.alv.noMx, label = TRUE, group.by="Annotation", cols=clusterCols,
        split.by="orig.ident", pt.size=1.5);
ggsave(sprintf("M23_UMAP320_pub_byAnnot_split_alevin_noMx_no7_SCT_%s.png", Sys.Date()), width=16, height=8);
ggsave(sprintf("M23_UMAP320_pub_byAnnot_split_alevin_noMx_no7_SCT_%s.pdf", Sys.Date()), width=16, height=8);

## Feature plot and dot plot for gene list
genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);
{ ## OL-selected DE Genes
    geneIndex <- match(genes.OL13, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    subMat <- GetAssayData(c57.S6KO.alv.noMx, slot="scale.data")[listSub,];
    cellOrder <- hclust(dist(t(subMat)))$order;
    cairo_pdf(sprintf("Heatmap_OL13_pub_byAnnot_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub, cells=cellOrder,
                     group.by="Annotation", raster=FALSE, lines.width=20, group.colors=clusterCols) +
        scale_fill_viridis() +
        scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        guides(col=FALSE)
    print(res);
    invisible(dev.off());
    cairo_pdf(sprintf("Dotplot_OL13_pub_byAnnot_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
                   group.by="AnnotationREV", split.by="orig.ident",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
        scale_colour_viridis(discrete=TRUE) +
        guides(col=FALSE) +
        ggtitle("Figure X") +
        ylab("Cluster") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
    cairo_pdf(sprintf("Dotplot_OL13_pub_unsplit_byAnnot_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
                   group.by="AnnotationREV",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
        scale_colour_viridis() +
        ggtitle("Figure X") +
        ylab("Cluster") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
}


## Heatmap for gene list
genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);
{ ## OL-selected DE Genes
    geneIndex <- match(genes.OL13, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    subMat <- GetAssayData(c57.S6KO.alv.noMx, slot="scale.data")[listSub,];
    cellOrder <- hclust(dist(t(subMat)))$order;
    cairo_pdf(sprintf("Heatmap_OL13_pub_unscaled_unsorted_byAnnot_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub, slot="data",
                     group.by="Annotation", raster=FALSE, lines.width=20, group.colors=clusterCols) +
        scale_fill_viridis() +
        scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        guides(col=FALSE)
    print(res);
    invisible(dev.off());
    cairo_pdf(sprintf("Heatmap_OL13_pub_unsorted_byAnnot_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub, slot="scale.data",
                     group.by="Annotation", raster=FALSE, lines.width=20, group.colors=clusterCols) +
        scale_fill_viridis() +
        scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        guides(col=FALSE)
    print(res);
    invisible(dev.off());
}


gene.counts <- GetAssayData(c57.S6KO.alv.noMx, slot="data");
colnames(gene.counts) <- paste0(unlist(c57.S6KO.alv.noMx[["Annotation"]]), "_", sub("^c57", "WT", colnames(gene.counts)));

expression.summaries <- gene.counts %>%
    as_tibble(rownames="gene") %>%
    pivot_longer(names_sep="_", names_to=c("cluster", "label", "cell"), cols=-gene) %>%
    group_by(gene, cluster) %>%
    mutate(Total=round(mean(value), 3)) %>%
    group_by(gene, cluster, label) %>%
    summarise(subMean=round(mean(value), 3), Total=dplyr::first(Total)) %>%
    pivot_wider(id_cols=c("gene", "cluster", "Total"), names_from=label, values_from=subMean) %>%
    pivot_longer(names_to="stat", cols=c("Total", "WT", "Stat6KO"), values_to="expression") %>%
    pivot_wider(id_cols="gene", names_from=c("cluster", "stat"), values_from=expression);

write_csv(expression.summaries, sprintf("Expression_summaries_pub_Cluster_label_byAnnot_alv_NoMx_%s.csv", Sys.Date()));

#### 'Without Stat6' investigation ##

gene.counts.noStat6 <- GetAssayData(c57.S6KO.alv.filt.noStat6, slot="data");
colnames(gene.counts.noStat6) <- paste0(unlist(c57.S6KO.alv.filt.noStat6[["cluster"]]), "_",
                                        sub("^c57", "WT", colnames(gene.counts.noStat6)));

expression.summaries.noStat6 <- gene.counts.noStat6 %>%
    as_tibble(rownames="gene") %>%
    pivot_longer(names_sep="_", names_to=c("cluster", "label", "cell"), cols=-gene) %>%
    group_by(gene, cluster) %>%
    mutate(Total=round(mean(value), 3)) %>%
    group_by(gene, cluster, label) %>%
    filter(label %in% c("WT", "Stat6KO")) %>%
    summarise(subMean=round(mean(value), 3), Total=dplyr::first(Total)) %>%
    pivot_wider(id_cols=c("gene", "cluster", "Total"), names_from=label, values_from=subMean) %>%
    pivot_longer(names_to="stat", cols=c("Total", "WT", "Stat6KO"), values_to="expression") %>%
    pivot_wider(id_cols="gene", names_from=c("cluster", "stat"), values_from=expression);

write_csv(expression.summaries.noStat6, sprintf("Expression_summaries_filtered_noStat6_alv_NoMx_%s.csv", Sys.Date()));

gene.counts.withStat6 <- GetAssayData(c57.S6KO.alv.filt, slot="data");
colnames(gene.counts.withStat6) <- paste0(unlist(c57.S6KO.alv.filt[["cluster"]]), "_",
                                        sub("^c57", "WT", colnames(gene.counts.withStat6)));

expression.summaries.withStat6 <- gene.counts.withStat6 %>%
    as_tibble(rownames="gene") %>%
    pivot_longer(names_sep="_", names_to=c("cluster", "label", "cell"), cols=-gene) %>%
    group_by(gene, cluster) %>%
    mutate(Total=round(mean(value), 3)) %>%
    group_by(gene, cluster, label) %>%
    filter(label %in% c("WT", "Stat6KO")) %>%
    summarise(subMean=round(mean(value), 3), Total=dplyr::first(Total)) %>%
    pivot_wider(id_cols=c("gene", "cluster", "Total"), names_from=label, values_from=subMean) %>%
    pivot_longer(names_to="stat", cols=c("Total", "WT", "Stat6KO"), values_to="expression") %>%
    pivot_wider(id_cols="gene", names_from=c("cluster", "stat"), values_from=expression);

write_csv(expression.summaries.withStat6, sprintf("Expression_summaries_filtered_withStat6_alv_NoMx_%s.csv", Sys.Date()));

c57.S6KO.alv.noMx[["cluster_cell"]] <-
    paste0(unlist(c57.S6KO.alv.noMx[["cluster"]]), "_",
           unlist(c57.S6KO.alv.noMx[["orig.ident"]]));
ccVals <- unique(unlist(c57.S6KO.alv.noMx[["cluster_cell"]]));
ccVals <- ccVals[order(-as.numeric(sub("_.*$", "", ccVals)),
                       xtfrm(sub("^.*_", "", ccVals)))];
c57.S6KO.alv.noMx[["cluster_cell"]] <-
    factor(unlist(c57.S6KO.alv.noMx[["cluster_cell"]]), levels=ccVals);

## Feature plot and dot plot for gene list
genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);
{ ## OL-selected DE Genes
    geneIndex <- match(c(genes.OL13, "Stat6"), sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    subMat <- GetAssayData(c57.S6KO.alv.noMx, slot="scale.data")[listSub,];
    cellOrder <- hclust(dist(t(subMat)))$order;
    cairo_pdf(sprintf("Heatmap_OL13_byCluster_SCT_alv_withoutStat6_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub, cells=cellOrder,
                     group.by="cluster", raster=FALSE, lines.width=20, group.colors=clusterCols) +
        scale_fill_viridis() +
        scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        guides(col=FALSE)
    print(res);
    invisible(dev.off());
    cairo_pdf(sprintf("Dotplot_OL13_split_byCluster_SCT_alv_withoutStat6_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
                   group.by="cluster_cell",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
        scale_colour_viridis() +
        ggtitle("Figure X") +
        ylab("Cluster") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
    cairo_pdf(sprintf("Dotplot_OL13_unsplit_byCluster_SCT_alv_withoutStat6_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
                   group.by="clusterREV",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
        scale_colour_viridis() +
        ggtitle("Figure X") +
        ylab("Cluster") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
}


## Heatmap for gene list
{ ## OL-selected DE Genes
    genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);
    geneIndex <- match(c(genes.OL13, "Stat6"), sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    subMat <- GetAssayData(c57.S6KO.alv.noMx, slot="scale.data")[listSub,];
    ## order within clusters
    clNames <- as.character(unique(unlist(c57.S6KO.alv.noMx[["cluster"]])));
    cellOrder <- NULL;
    for(cl in clNames[order(as.numeric(clNames))]){
        cat(cl,"\n");
        newOrder <- hclust(dist(t(subMat[, unlist(c57.S6KO.alv.noMx[["cluster"]]) == cl])))$order;
        cellOrder <- c(cellOrder, newOrder + length(cellOrder));
    }
    cairo_pdf(sprintf("Heatmap_OL13_noStat6_unscaled_unsorted_byCluster_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub, slot="data",
                     group.by="cluster", raster=FALSE, lines.width=20, group.colors=clusterCols) +
        scale_fill_viridis() +
        scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        guides(col=FALSE)
    print(res);
    invisible(dev.off());
    cairo_pdf(sprintf("Heatmap_OL13_noStat6_unsorted_byCluster_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub, slot="scale.data",
                     group.by="cluster", raster=FALSE, lines.width=20, group.colors=clusterCols) +
        scale_fill_viridis() +
        scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        guides(col=FALSE)
    print(res);
    invisible(dev.off());
}

## Cluster Tree for cluster relationships [2021-Mar-04]
## [https://stackoverflow.com/questions/62882606/]
c57.S6KO.alv.noMx <- BuildClusterTree(c57.S6KO.alv.noMx, assay="SCT")
myPhyTree <- Tool(object=c57.S6KO.alv.noMx, slot = "BuildClusterTree")
ggtree(myPhyTree) +
geom_tiplab()+
theme_tree()
ggsave(sprintf("tree_SCT_alv_withoutStat6_%s.pdf", Sys.Date()));


## Recluster with increasing resolution, until sub-clusters are found
## that fully separate identified clusters from 'withStat6' and
## 'withoutStat6'

c57.S6KO.alv.filt <- FindNeighbors(c57.S6KO.alv.filt, dims = 1:60);
c57.S6KO.alv.filt.clusterTest <- c57.S6KO.alv.filt;

res <- table(list(res0.8=Idents(FindClusters(c57.S6KO.alv.filt.clusterTest, resolution=0.8)),
                  Stat6=Idents(c57.S6KO.alv.filt)))

res <- table(list(res0.8=Idents(FindClusters(c57.S6KO.alv.filt.clusterTest, resolution=0.8)),
                  noStat6=Idents(c57.S6KO.alv.noMx)))

max(apply(res,1,function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[1] / cmp.vals[2]}));

for(rTest in seq(0.4, 1.5, by=0.1)){
    #cat(rTest, "\n");
    idents.rt <- Idents(FindClusters(c57.S6KO.alv.filt.clusterTest, resolution=rTest, n.start=10, verbose=FALSE));
    res.b <- table(list(Idents(c57.S6KO.alv.noMx), Idents(c57.S6KO.alv.filt)));
    res.1 <- table(list(idents.rt, Idents(c57.S6KO.alv.filt)));
    res.2 <- table(list(idents.rt, Idents(c57.S6KO.alv.noMx)));
    res.cluster <- table(list(idents.rt, as.character(unlist(c57.S6KO.alv.filt[["orig.ident"]]))));
    write.csv(res.cluster, sprintf("cluster_counts_res%0.1f_%s.csv", rTest, Sys.Date()));
    maxProp.1.1 <- max(apply(res.1, 1, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[1] / cmp.vals[2]}));
    maxProp.2.1 <- max(apply(res.2, 1, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[1] / cmp.vals[2]}));
    maxProp.1.2 <- max(apply(res.1, 2, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[1] / cmp.vals[2]}));
    maxProp.2.2 <- max(apply(res.2, 2, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[1] / cmp.vals[2]}));
    maxDiff.1.1 <- max(apply(res.1, 1, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[2] - cmp.vals[1]}));
    maxDiff.2.1 <- max(apply(res.2, 1, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[2] - cmp.vals[1]}));
    maxDiff.1.2 <- max(apply(res.1, 2, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[2] - cmp.vals[1]}));
    maxDiff.2.2 <- max(apply(res.2, 2, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[2] - cmp.vals[1]}));
    cat(sprintf("%s: %0.2f (%d) with Stat6; %0.2f (%d) without Stat6\n",
                rTest, maxProp.1.1, maxDiff.1.1, maxProp.2.1, maxDiff.2.1));
    cat(sprintf("%s: %0.2f (%d) with Stat6; %0.2f (%d) without Stat6\n",
                rTest, maxProp.1.2, maxDiff.1.2, maxProp.2.2, maxDiff.2.2));
    c57.S6KO.alv.filt[["cluster.new"]] <- idents.rt;
    c57.S6KO.alv.filt.clusterTest[["cluster.new"]] <- idents.rt;
    DimPlot(c57.S6KO.alv.filt, label=FALSE, group.by="cluster.new", pt.size=1.5, split.by="orig.ident") +
        ggtitle(sprintf("SCTransform with Stat6, cluster resolution %0.1f", rTest));
    ggsave(sprintf("M23_UMAP320_origValues_filt_Stat6Remapped_res%0.1f_SCT_alv_%s.png", rTest, Sys.Date()), width=11, height=8);
    ## DotPlot
    genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);
    geneIndex <- match(c(genes.OL13, "Zap70"), sub("-ENSM.*$", "", rownames(c57.S6KO.alv.filt)));
    listSub <- na.omit(rownames(c57.S6KO.alv.filt)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    subMat <- GetAssayData(c57.S6KO.alv.filt, slot="scale.data")[listSub,];
    cellOrder <- hclust(dist(t(subMat)))$order;
    cairo_pdf(sprintf("DotPlot_OL13_withStat6_unscaled_unsorted_byCluster_res%0.1f_alv_NoMx_%s.pdf", rTest, Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.filt, features=listSub,
                   group.by="cluster.new",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
            scale_colour_viridis() +
            ylab("Cluster") +
            ggtitle(sprintf("%s (res: %0.1f)", plotName, rTest)) +
            scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
    Idents(c57.S6KO.alv.filt.clusterTest) <- idents.rt;
    c57.S6KO.alv.filt.clusterTest <- BuildClusterTree(c57.S6KO.alv.filt.clusterTest, assay="SCT")
    myPhyTree <- Tool(object=c57.S6KO.alv.filt.clusterTest, slot = "BuildClusterTree")
    ggtree(myPhyTree) +
        geom_tiplab()+
        theme_tree() +
        ggtitle(sprintf("Dendrogram, with Stat6 (res: %0.1f)", rTest))
    ggsave(sprintf("tree_alv_res%0.1f_%s.pdf", rTest, Sys.Date()));
}

## Redo without Stat6
c57.S6KO.alv.filt.clusterTest <- c57.S6KO.alv.noMx;
for(rTest in seq(0.4, 1.5, by=0.1)){
    #cat(rTest, "\n");
    idents.rt <- Idents(FindClusters(c57.S6KO.alv.filt.clusterTest, resolution=rTest, n.start=10, verbose=FALSE));
    res.b <- table(list(Idents(c57.S6KO.alv.noMx), Idents(c57.S6KO.alv.filt)));
    res.1 <- table(list(idents.rt, Idents(c57.S6KO.alv.filt)));
    res.2 <- table(list(idents.rt, Idents(c57.S6KO.alv.noMx)));
    res.cluster <- table(list(idents.rt, as.character(unlist(c57.S6KO.alv.filt[["orig.ident"]]))));
    write.csv(res.cluster, sprintf("cluster_counts_res%0.1f_noStat6_%s.csv", rTest, Sys.Date()));
    maxProp.1.1 <- max(apply(res.1, 1, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[1] / cmp.vals[2]}));
    maxProp.2.1 <- max(apply(res.2, 1, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[1] / cmp.vals[2]}));
    maxProp.1.2 <- max(apply(res.1, 2, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[1] / cmp.vals[2]}));
    maxProp.2.2 <- max(apply(res.2, 2, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[1] / cmp.vals[2]}));
    maxDiff.1.1 <- max(apply(res.1, 1, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[2] - cmp.vals[1]}));
    maxDiff.2.1 <- max(apply(res.2, 1, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[2] - cmp.vals[1]}));
    maxDiff.1.2 <- max(apply(res.1, 2, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[2] - cmp.vals[1]}));
    maxDiff.2.2 <- max(apply(res.2, 2, function(x){cmp.vals <- tail(sort(x), 2); cmp.vals[2] - cmp.vals[1]}));
    cat(sprintf("%s: %0.2f (%d) with Stat6; %0.2f (%d) without Stat6\n",
                rTest, maxProp.1.1, maxDiff.1.1, maxProp.2.1, maxDiff.2.1));
    cat(sprintf("%s: %0.2f (%d) with Stat6; %0.2f (%d) without Stat6\n",
                rTest, maxProp.1.2, maxDiff.1.2, maxProp.2.2, maxDiff.2.2));
    c57.S6KO.alv.filt[["cluster.new"]] <- idents.rt;
    c57.S6KO.alv.filt.clusterTest[["cluster.new"]] <- idents.rt;
    DimPlot(c57.S6KO.alv.filt.clusterTest, label=FALSE, group.by="cluster.new", pt.size=1.5, split.by="orig.ident") +
        ggtitle(sprintf("SCTransform without Stat6, cluster resolution %0.1f", rTest));
    ggsave(sprintf("M23_UMAP320_origValues_filt_res%0.1f_SCT_alv_noStat6_%s.png", rTest, Sys.Date()), width=11, height=8);
    ## DotPlot
    genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);
    geneIndex <- match(c(genes.OL13, "Zap70"), sub("-ENSM.*$", "", rownames(c57.S6KO.alv.filt)));
    listSub <- na.omit(rownames(c57.S6KO.alv.filt)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    subMat <- GetAssayData(c57.S6KO.alv.filt.clusterTest, slot="scale.data")[listSub,];
    cellOrder <- hclust(dist(t(subMat)))$order;
    cairo_pdf(sprintf("DotPlot_OL13_unscaled_unsorted_byCluster_res%0.1f_alv_NoStat6_%s.pdf", rTest, Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.filt.clusterTest, features=listSub,
                   group.by="cluster.new",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
            scale_colour_viridis() +
            ylab("Cluster") +
            ggtitle(sprintf("%s (res: %0.1f)", plotName, rTest)) +
            scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
    Idents(c57.S6KO.alv.filt.clusterTest) <- idents.rt;
    c57.S6KO.alv.filt.clusterTest <- BuildClusterTree(c57.S6KO.alv.filt.clusterTest, assay="SCT")
    myPhyTree <- Tool(object=c57.S6KO.alv.filt.clusterTest, slot = "BuildClusterTree")
    ggtree(myPhyTree) +
        geom_tiplab()+
        theme_tree() +
        ggtitle(sprintf("Dendrogram, without Stat6 (res: %0.1f)", rTest))
    ggsave(sprintf("tree_alv_res%0.1f_noStat6_%s.pdf", rTest, Sys.Date()));
}

## OL meeting 2021-Apr-01
## [use gene set OL13]

genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);

OL.labels <- read_xlsx("Cluster_annotation_2021-Apr-01.xlsx");
OL.names <- OL.labels$Name;
names(OL.names) <- OL.labels$`Cluster Number`;
OL.colours <- OL.labels$Colored;

clusterCols <- rep("", nrow(OL.labels));

clusterCols[OL.labels$Colored] <-
    brewer.pal(sum(OL.labels$Colored), "Dark2");
clusterCols[!OL.labels$Colored] <-
    grey_pal(end=0.85)(sum(!OL.labels$Colored));
names(clusterCols) <- OL.labels$Name;

clusterCols["C12"] <- "darkgreen";

c57.S6KO.alv.noMx[["OL.clusters"]] <-
    factor(OL.names[as.character(unlist(c57.S6KO.alv.noMx[["cluster"]]))],
           levels=OL.names);
c57.S6KO.alv.noMx[["OL.clustersREV"]] <-
    factor(OL.names[as.character(unlist(c57.S6KO.alv.noMx[["cluster"]]))],
           levels=rev(OL.names));
c57.S6KO.alv.noMx[["OL.clustersSplit"]] <-
    factor(paste(as.character(unlist(c57.S6KO.alv.noMx[["OL.clusters"]])),
                 as.character(unlist(c57.S6KO.alv.noMx[["orig.ident"]])),
                 sep="_"),
           levels=paste(rep(OL.names, each=2), c("WT", "Stat6KO"), sep="_"));
c57.S6KO.alv.noMx[["OL.clustersSplitREV"]] <-
    factor(paste(as.character(unlist(c57.S6KO.alv.noMx[["OL.clusters"]])),
                 as.character(unlist(c57.S6KO.alv.noMx[["orig.ident"]])),
                 sep="_"),
           levels=rev(paste(rep(OL.names, each=2), c("WT", "Stat6KO"), sep="_")));

## UMAP

DimPlot(c57.S6KO.alv.noMx, label=FALSE, group.by="OL.clusters", pt.size=1.5, split.by="orig.ident",
        cols=clusterCols) +
ggtitle("SCTransform without Stat6, cluster resolution 1.5");
ggsave(sprintf("M23_UMAP320_split_res1.5_SCT_alv_withoutStat6_%s.png", Sys.Date()), width=14, height=8);
ggsave(sprintf("M23_UMAP320_split_res1.5_SCT_alv_withoutStat6_%s.pdf", Sys.Date()), width=14, height=8);

DimPlot(c57.S6KO.alv.noMx, label=FALSE, group.by="OL.clusters", pt.size=1.5,
        cols=clusterCols) +
ggtitle("SCTransform without Stat6, cluster resolution 1.5");
ggsave(sprintf("M23_UMAP320_res1.5_SCT_alv_withoutStat6_%s.png", Sys.Date()), width=8, height=8);
ggsave(sprintf("M23_UMAP320_res1.5_SCT_alv_withoutStat6_%s.pdf", Sys.Date()), width=8, height=8);

## tSNE
c57.S6KO.alv.noMx <- RunTSNE(c57.S6KO.alv.noMx);

DimPlot(c57.S6KO.alv.noMx, label=FALSE, group.by="OL.clusters", pt.size=1.5, split.by="orig.ident",
        cols=clusterCols, reduction="tsne") +
ggtitle("SCTransform without Stat6, cluster resolution 1.5");
ggsave(sprintf("M23_tSNE_split_res1.5_SCT_alv_withoutStat6_%s.png", Sys.Date()), width=14, height=8);
ggsave(sprintf("M23_tSNE_split_res1.5_SCT_alv_withoutStat6_%s.pdf", Sys.Date()), width=14, height=8);

DimPlot(c57.S6KO.alv.noMx, label=FALSE, group.by="OL.clusters", pt.size=1.5,
        cols=clusterCols, reduction="tsne") +
ggtitle("SCTransform without Stat6, cluster resolution 1.5");
ggsave(sprintf("M23_tSNE_res1.5_SCT_alv_withoutStat6_%s.png", Sys.Date()), width=14, height=8);
ggsave(sprintf("M23_tSNE_res1.5_SCT_alv_withoutStat6_%s.pdf", Sys.Date()), width=14, height=8);

## DotPlot
{
    genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);
    geneIndex <- match(c(genes.OL13), sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    subMat <- GetAssayData(c57.S6KO.alv.noMx, slot="scale.data")[listSub,];
    cellOrder <- hclust(dist(t(subMat)))$order;
    cairo_pdf(sprintf("DotPlot_OL13_byLabel_res1.5_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
                   group.by="OL.clustersREV",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
        scale_colour_viridis() +
        ylab("Cluster") +
        ggtitle("Gene set OL13 (res: 1.5)") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
    cairo_pdf(sprintf("DotPlot_OL13_split_byLabel_res1.5_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DotPlot(c57.S6KO.alv.noMx, features=listSub,
                   group.by="OL.clustersSplitREV",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
        scale_colour_viridis() +
        ylab("Cluster") +
        ggtitle("Gene set OL13 (res: 1.5)") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
}

## HeatMap
{
    genes.OL13 <- read_xlsx('List of genes for figures.xlsx') %>% pull(`Gene Name`);
    geneIndex <- match(c(genes.OL13), sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    subMat <- GetAssayData(c57.S6KO.alv.noMx, slot="scale.data")[listSub,];
    cellOrder <- hclust(dist(t(subMat)))$order;
    cairo_pdf(sprintf("Heatmap_OL13_byLabel_res1.5_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub, cells=cellOrder,
                     group.by="OL.clusters", raster=FALSE, lines.width=20, group.colors=clusterCols) +
        scale_fill_viridis() +
        scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        guides(col=FALSE)
    print(res);
    invisible(dev.off());
    cairo_pdf(sprintf("Heatmap_OL13_split_byLabel_res1.5_alv_NoMx_%s.pdf", Sys.Date()), onefile=TRUE,
              width=10, height=8, pointsize=8);
    res <- DoHeatmap(c57.S6KO.alv.noMx, features=listSub, cells=cellOrder,
                     group.by="OL.clustersSplit", raster=FALSE, lines.width=20, group.colors=clusterCols) +
        scale_fill_viridis() +
        scale_y_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        guides(col=FALSE)
    print(res);
    invisible(dev.off());
}

## monocle
## https://satijalab.org/signac/articles/monocle.html#building-trajectories-with-monocle-3-1

##remotes::install_github('satijalab/seurat-wrappers');
library(SeuratWrappers);
library(monocle3);

##BiocManager::install("batchelor");
##BiocManager::install("fastmap");

"SummarizedExperiment, S4Vectors, limma, DelayedMatrixStats, DelayedArray, BiocGenerics, batchelor, SingleCellEx\
 periment, Biobase"

DefaultAssay(c57.S6KO.alv.noMx) <- "RNA";
c57.S6KO.cds <- as.cell_data_set(c57.S6KO.alv.noMx);
c57.S6KO.cds <- cluster_cells(cds = c57.S6KO.cds, reduction_method = "UMAP");
c57.S6KO.cds <- learn_graph(c57.S6KO.cds, use_partition = TRUE);
c57.S6KO.cds <- order_cells(c57.S6KO.cds, reduction_method = "UMAP");

cairo_pdf(sprintf("monocle_pseudotime_trajectory_%s.pdf", Sys.Date()), onefile=TRUE,
          width=11, height=8);
plot_cells(c57.S6KO.cds,
           cell_size=2,
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE);
invisible(dev.off());

## do mini-bulk with CD11b high vs low
c57.S6KO.alv.noMx[["CD11b_est"]] <-
    ifelse(unlist(c57.S6KO.alv.noMx[["OL.clusters"]]) %in% paste0("C", 1:8),
           "high", "low");

## CD11c high vs low
c57.S6KO.alv.noMx[["CD11c"]] <-
    ifelse(GetAssayData(c57.S6KO.alv.noMx, slot="counts")["Itgax-ENSMUSG00000030789.9",] == 0,
           "low", "high");


{
    res <- NULL;
    res.all <- NULL;
    for(marker in c("CD11b_est", "CD11c")){
        for(pop in c("high", "low")){
            cat(marker, pop,"\n");
            subCells <- which(unlist(c57.S6KO.alv.noMx[[marker]]) == pop);
            c57.S6KO.alv.noMx.distSub <- subset(c57.S6KO.alv.noMx, cells=subCells);
            test.markers <- FindMarkers(c57.S6KO.alv.noMx.distSub, group.by="orig.ident", test.use="DESeq2",
                                        ident.1 = "Stat6KO", ident.2 = "WT");
            res.all <- rbind(res.all, cbind(subCluster=paste0(marker, "_", pop), test="Stat6KO-WT",
                                            gene=rownames(test.markers), logSum=geneLogSums[rownames(test.markers)], test.markers));
            test.markers <- subset(test.markers, p_val_adj < 0.1);
            test.markers <- test.markers[order(-test.markers$avg_log2FC),];
            res <- rbind(res, cbind(subCluster=paste0(marker, "_", pop), test="Stat6KO-WT",
                                    gene=rownames(test.markers), logSum=geneLogSums[rownames(test.markers)], test.markers));
        }
    }
    write_csv(res, file=sprintf("ClusterDE_1.5_miniBulk_noMx_alv_%s.csv", Sys.Date()));
    write_csv(res.all, file=sprintf("ClusterDE_1.5_all_miniBulk_noMx_alv_%s.csv", Sys.Date()));
}

## Alternative monocle construction
## [https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/]
DefaultAssay(c57.S6KO.alv.noMx) <- "RNA";

c57.S6KO.cds <- as.cell_data_set(c57.S6KO.alv.noMx) %>%
    preprocess_cds();

##c57.S6KO.cds <- align_cds(c57.S6KO.cds, alignment_group="orig.ident");
c57.S6KO.cds <- reduce_dimension(c57.S6KO.cds, verbose=TRUE);

c57.S6KO.cds <- cluster_cells(cds = c57.S6KO.cds, num_iter=1, k=50, resolution=0.04, cluster_method="leiden");
length(levels(clusters(c57.S6KO.cds)));

c57.S6KO.cds <- learn_graph(c57.S6KO.cds);

c57.S6KO.cds <- order_cells(c57.S6KO.cds, reduction_method = "LSI");

str(colData(c57.S6KO.cds))


cairo_pdf(sprintf("monocle_UMAP_altConstruction_%s.pdf", Sys.Date()), onefile=TRUE,
          width=11, height=8);
plot_cells(c57.S6KO.cds,
           color_cells_by="cluster",
           label_groups_by_cluster=FALSE,
           group_label_size=4,
           graph_label_size=3,
           show_trajectory_graph=FALSE,
           cell_size=2) +
plot_cells(c57.S6KO.cds,
           color_cells_by="orig.ident",
           label_groups_by_cluster=FALSE,
           group_label_size=4,
           graph_label_size=3,
           show_trajectory_graph=FALSE,
           cell_size=2) +
    labs(title = "Monocle3 UMAP, cell tags")
plot_cells(c57.S6KO.cds,
           color_cells_by="OL.clusters",
           label_groups_by_cluster=FALSE,
           group_label_size=4,
           graph_label_size=3,
           show_trajectory_graph=FALSE,
           cell_size=2) +
    labs(title = "Monocle3 UMAP, OL clusters")
## plot_cells(c57.S6KO.cds,
##            color_cells_by="partition",
##            label_groups_by_cluster=FALSE,
##            group_label_size=4,
##            graph_label_size=3,
##            cell_size=2);
plot_cells(c57.S6KO.cds,
           color_cells_by="cluster",
           label_groups_by_cluster=FALSE,
           group_label_size=4,
           graph_label_size=3,
           show_trajectory_graph=FALSE,
           cell_size=2) +
    labs(title = "Monocle3 UMAP, Monocle3 clusters")
plot_cells(c57.S6KO.cds,
           color_cells_by="cluster",
           label_groups_by_cluster=FALSE,
           group_label_size=4,
           graph_label_size=3,
           show_trajectory_graph=FALSE,
           cell_size=2) +
    facet_wrap(~orig.ident, ncol=2) +
    labs(title = "Monocle3 UMAP, Monocle3 clusters (split by cell tag)")
## plot_cells(c57.S6KO.cds,
##            color_cells_by="pseudotime",
##            label_cell_groups=FALSE,
##            label_branch_points=FALSE,
##            graph_label_size=3,
##            cell_size=2);
invisible(dev.off());

library(monocle3);

cairo_pdf(sprintf("monocle_UMAP_FeaturePlots_%s.pdf", Sys.Date()), onefile=TRUE,
          width=11, height=8);
monocle::plot_cell_clusters(c57.S6KO.cds,
                   markers = listSub,
                   show_group_id = T, cell_size = 0.1);
invisible(dev.off());

marker_test_res <- top_markers(c57.S6KO.cds, group_cells_by="cluster",
                               cores=10);

c57.S6KO.alv.noMx.monocle3 <- AddMetaData(c57.S6KO.alv.noMx,
                                          metadata=c57.S6KO.cds@clusters@listData$UMAP$clusters,
                                          col.name="monocle50_0.04_cluster");

c57.S6KO.alv.noMx.monocle3 <- AddMetaData(c57.S6KO.alv.noMx.monocle3,
                                          metadata=as.character(c57.S6KO.cds@clusters@listData$UMAP$clusters) == "9",
                                          col.name="monocle_dump");

c57.S6KO.alv.noMx.monocle3[["monocleUMAP"]] <-
    CreateDimReducObject(
    embeddings = reducedDims(c57.S6KO.cds)$UMAP,
    key = "monocleUMAP_",
    assay = "RNA"
)

str(reducedDims(c57.S6KO.cds)$UMAP)

str(c57.S6KO.cds@clusters@listData$UMAP$cluster_result);

DefaultAssay(c57.S6KO.alv.noMx.monocle3) <- "SCT";
{
    cairo_pdf(sprintf("monocle_UMAP_asSeuratPlot_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    print(DimPlot(c57.S6KO.alv.noMx.monocle3, reduction="umap",
                  group.by="monocle50_0.04_cluster", pt.size=1.5));
    print(DimPlot(c57.S6KO.alv.noMx.monocle3, reduction="monocleUMAP",
                  group.by="monocle50_0.04_cluster", pt.size=1.5));
    print(DimPlot(c57.S6KO.alv.noMx.monocle3, reduction="monocleUMAP",
            group.by="monocle50_0.04_cluster", pt.size=1.5, split.by="orig.ident"));
    print(DimPlot(c57.S6KO.alv.noMx.monocle3, reduction="monocleUMAP",
            group.by="OL.clusters", cols=clusterCols, pt.size=1.5));
    print(DimPlot(c57.S6KO.alv.noMx.monocle3, reduction="monocleUMAP",
            group.by="OL.clusters", split.by="orig.ident", cols=clusterCols, pt.size=1.5));
    ## Create Gene Feature Plots
    geneList.OL14 <- c("Aldh1a2", "Ccl9", "Ccl17", "Ccrl2", "Fgl2", "Itgax",
                       "Mpeg1", "Slc7a11", "Crispld2", "Il9r", "Sned1",
                       "Zbtb10", "Il1b", "Ccl5", "Ccl22", "Cd74", "Cd1d1");
    geneIndex <- match(geneList.OL14, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    selectPoss <- seq(1, length(listSub), by=4);
    maxRes <- NULL;
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx.monocle3, reduction="monocleUMAP",
                                      col=c("lightgrey", "#e31836"),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=2, nrow=2)));
    }
    res <- DotPlot(c57.S6KO.alv.noMx.monocle3, features=listSub,
                   group.by="monocle50_0.04_cluster", assay="SCT",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
        scale_colour_viridis() +
        ylab("Cluster") +
        ggtitle("Gene set OL14") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    invisible(dev.off());
}


write_csv(marker_test_res, sprintf("top_markers_monocle50_0.04_clusters_%s.csv", Sys.Date()));

c57.S6KO.alv.noMx.monocle3.filt <-
    subset(c57.S6KO.alv.noMx.monocle3,
           subset=(monocle_dump == FALSE));

c57.S6KO.cds.filt <- as.cell_data_set(c57.S6KO.alv.noMx.monocle3.filt) %>%
    preprocess_cds();

## -> 0.04
c57.S6KO.cds.filt <- reduce_dimension(c57.S6KO.cds.filt, verbose=TRUE);

c57.S6KO.alv.noMx.monocle3.filt[["monocleUMAP"]] <-
    CreateDimReducObject(
        embeddings = reducedDims(c57.S6KO.cds.filt)$UMAP,
        key = "monocleUMAP_",
        assay = "RNA"
    )

{
    cairo_pdf(sprintf("monocle_clusterTest_UMAP_asSeuratPlot_%s.pdf", Sys.Date()), onefile=TRUE,
          width=11, height=8);
    for(tRes in round(seq(0.002, 0.04, length.out=10), 3)){
        c57.S6KO.cds.filt <-
            cluster_cells(cds = c57.S6KO.cds.filt, num_iter=1, k=50,
                          resolution=tRes, cluster_method="leiden");
        cat(tRes, "=", length(levels(clusters(c57.S6KO.cds.filt))), "\n");
        c57.S6KO.alv.noMx.monocle3.filt <- AddMetaData(c57.S6KO.alv.noMx.monocle3.filt,
                                                       metadata=c57.S6KO.cds.filt@clusters@listData$UMAP$clusters,
                                                       col.name="monocle_cluster");
        c57.S6KO.alv.noMx.monocle3.filt[["monocle_clusterSplit"]] <-
            factor(paste(as.character(unlist(c57.S6KO.alv.noMx.monocle3.filt[["monocle_cluster"]])),
                         as.character(unlist(c57.S6KO.alv.noMx.monocle3.filt[["orig.ident"]])),
                         sep="_"),
                   levels=rev(paste(rep(levels(unlist(c57.S6KO.alv.noMx.monocle3.filt[["monocle_cluster"]])), each=2),
                                    c("WT", "Stat6KO"), sep="_")));
        print(table(unlist(c57.S6KO.alv.noMx.monocle3.filt[["orig.ident"]]),
                    unlist(c57.S6KO.alv.noMx.monocle3.filt[["monocle_cluster"]])));
        print(DimPlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                      group.by="monocle_cluster", pt.size=1.5) +
              ggtitle(sprintf("k=50, res = %0.3f, Monocle UMAP", tRes)));
        print(DimPlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                      group.by="monocle_cluster", pt.size=1.5, split.by="orig.ident") +
              ggtitle(sprintf("k=50, res = %0.3f, Monocle UMAP (split)", tRes)));
        res <- DotPlot(c57.S6KO.alv.noMx.monocle3.filt, features=listSub,
                       group.by="monocle_clusterSplit", assay="SCT",
                       dot.min=0.0001, scale.by="size", scale=TRUE,
                       col.min=0, col.max=1) +
            scale_colour_viridis() +
            ylab("Cluster") +
            ggtitle(sprintf("k=50, res = %0.3f, Gene set OL14", tRes)) +
            scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            theme(axis.text.x=element_text(angle=45, hjust=1));
        print(res);
    }
    invisible(dev.off());
}


c57.S6KO.cds <- learn_graph(c57.S6KO.cds);

c57.S6KO.cds <- order_cells(c57.S6KO.cds, reduction_method = "LSI");

str(colData(c57.S6KO.cds))
