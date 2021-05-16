library(tximeta);
library(SingleCellExperiment);
library(org.Mm.eg.db);
library(tidyverse);
library(Seurat);
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
mtNames <- grep("(^LOC|[0-9][0-9]\\-|^LTO)", invert=TRUE, value=TRUE,
                grep("^[A-Z][A-Z]", rownames(counts.combined), value=TRUE));
rownames(counts.combined)[rownames(counts.combined) %in% mtNames] <-
    paste0("mt-", rownames(counts.combined)[rownames(counts.combined) %in% mtNames]);

## Okay, that's all set up; now onto Seurat
## [https://satijalab.org/seurat/articles/pbmc3k_tutorial.html]
## [https://satijalab.org/seurat/articles/sctransform_vignette.html]

c57.s6KO.all <- CreateSeuratObject(counts = counts.combined, project = "c57.s6KO",
                                   min.cells = 3, min.features = 200);

c57.s6KO.all[["orig.ident"]] <- Idents(c57.s6KO.all);


## Mouse mt genes start with lower case 'mt'
c57.s6KO.all[["percent.mt"]] <- PercentageFeatureSet(c57.s6KO.all, pattern = "^mt-");

VlnPlot(c57.s6KO.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3);
ggsave(sprintf("FeatureCount_alevin_all_%s.png", Sys.Date()), width=11, height=8);

plot1 <- FeatureScatter(c57.s6KO.all, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    geom_hline(yintercept=8) +
    geom_vline(xintercept=25000);
plot2 <- FeatureScatter(c57.s6KO.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(yintercept=c(200,3500)) +
    geom_vline(xintercept=25000);
plot3 <- FeatureScatter(c57.s6KO.all, feature1 = "nFeature_RNA", feature2 = "percent.mt");
plot1 | (plot2 / plot3);
ggsave(sprintf("FeatureScatter_alevin_all_%s.png", Sys.Date()), width=11, height=8);

c57.s6KO.filt <-
    subset(c57.s6KO.all,
           subset = nFeature_RNA > 200 &
               nFeature_RNA < 3500 &
               percent.mt < 8 &
           nCount_RNA < 25000);

c57.s6KO.filt <- SCTransform(c57.s6KO.filt, vars.to.regress="percent.mt");

## [iteration limit warnings]

c57.s6KO.filt <- RunPCA(c57.s6KO.filt, verbose = FALSE, npcs=100);

DimPlot(c57.s6KO.filt, reduction = "pca");
ggsave(sprintf("PCA_alevin_filt_SCT_%s.png", Sys.Date()), width=8, height=8);

png(sprintf("DimHeatmap_alevin_filt_SCT_%s.png", width=1280*2, height=720*2);
## Note: DimHeatmap is not a ggplot thing
DimHeatmap(c57.s6KO.filt, dims = 1:15, cells = 500, balanced = TRUE);
invisible(dev.off());

ElbowPlot(c57.s6KO.filt, ndims=100);
ggsave(sprintf("Elbowplot_alevin_filt_SCT_%s.png", Sys.Date()), width=8, height=8);

## Choose 50 dimensions for reductions
c57.s6KO.filt <- RunUMAP(c57.s6KO.filt, dims = 1:50, n.epochs=1000);

c57.s6KO.filt <- FindNeighbors(c57.s6KO.filt, dims = 1:50);
c57.s6KO.filt <- FindClusters(c57.s6KO.filt, resolution=0.8);
c57.s6KO.filt[["cluster"]] <- Idents(c57.s6KO.filt);

## Dims flipped to make the rabbit easier to see
DimPlot(c57.s6KO.filt, label=FALSE, group.by="orig.ident", dims=c(2,1),
        pt.size=1.5) +
    NoLegend();
ggsave(sprintf("UMAP_origValues_alevin_filt_SCT_noLabel_%s.png", Sys.Date()), width=12, height=6);

DimPlot(c57.s6KO.filt, label = TRUE, group.by="orig.ident", cols=DiscretePalette(4), dims=c(2,1),
        pt.size=1.5);
ggsave(sprintf("UMAP_origValues_alevin_filt_SCT_%s.png", Sys.Date()), width=12, height=6);
DimPlot(c57.s6KO.filt, label=FALSE, group.by="cluster", cols=DiscretePalette(8), dims=c(2,1),
        pt.size=1.5) +
   NoLegend();
ggsave(sprintf("UMAP_byCluster_alevin_filt_SCT_%s_noLabel.png", Sys.Date()), width=12, height=6);
DimPlot(c57.s6KO.filt, label = TRUE, group.by="cluster", cols=DiscretePalette(8), dims=c(2,1),
        pt.size=1.5);
ggsave(sprintf("UMAP_byCluster_alevin_filt_SCT_%s.png", Sys.Date()), width=12, height=6);
