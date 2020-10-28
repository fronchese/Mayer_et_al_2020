#!/usr/bin/env Rscript

library(dplyr);
library(Seurat);
library(patchwork);

dataDir <- "/mnt/BigBirdStornext/Scratch/public_data/GEO/GSE153760";

GSMDirs <- list.files(dataDir, pattern="GSM", full.names=TRUE);
names(GSMDirs) <- sub("^GSM.*_", "", list.files(dataDir, pattern="GSM"));

## Note: based on [https://satijalab.org/seurat/v3.0/integration.html]

##Takes ~5 mins
expressionMat <- Read10X(data.dir=GSMDirs);

## For secondary analysis, the R-package ‘Seurat’ was used (Seurat
## v3.1, Satija Lab, NYU, New York, USA).

sessionInfo()$otherPkgs[["Seurat"]]$Version;
## [1] "3.2.0"

## To remove unwanted variations in the scRNA-seq data, cells were
## first analyzed for their UMI and mitochondrial gene counts, and
## cells with low (<500) UMI counts were excluded from the data set.

seur.obj <- CreateSeuratObject(counts = expressionMat);

## To align data of all samples and remove batch effects, the data was
## integrated in a standardized workflow as recommended by the
## developers of the “Seurat”-package. This workflow (Hafemeister &
## Satija, preprint Bioarxiv, https://doi.org/10.1101/576827) includes
## data normalization

## as per [https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html]

## Determine mitochondrial expression percentage

seur.obj[["percent.mt"]] <-
    PercentageFeatureSet(seur.obj, pattern = "^MT-");

### based on Vera's email, 2020-Sep-03, Relabel HC1 and HC2 as HC1_2
## hc.list<-c(hc1_2, ad1, ad2, ad3, ad4, ad5, ad6, ad7, ad8, hc3, hc4, hc5)
seur.obj <- RenameIdents(seur.obj, "HC1" = "HC1_2", "HC2" = "HC1_2");

## remove HC6/HC7 [based on the same email]
samplesToFetch <- levels(Idents(seur.obj));

## remove AD samples
samplesToFetch <- samplesToFetch[!grepl("AD", samplesToFetch)];

hc.unfiltered <- subset(seur.obj, idents=samplesToFetch);

pdf(sprintf("percent_mt_HC_unfiltered_%s.pdf", Sys.Date()));
plot(density(unlist(hc.unfiltered[["percent.mt"]])),
     main="Percent Mitochondrial Transcripts (HC only)");
abline(v=c(0.5, 15), col="#80808050");
print(VlnPlot(hc.unfiltered, features=c("percent.mt"),
              pt.size=0.25, cols=DiscretePalette(14)) +
      geom_hline(yintercept=c(0.5, 15), col="#80808050"));
print(FeatureScatter(hc.unfiltered,
                     feature1="percent.mt", feature2="nCount_RNA",
                     cols=DiscretePalette(14), pt.size=0.25) +
      geom_vline(xintercept=c(0.5,15), col="#80808050") +
      geom_hline(yintercept=10000, col="#80808050"));
dev.off();

pdf(sprintf("Cluster_features_vs_counts_HC_unfiltered_RNA_%s.pdf", Sys.Date()), width=11, height=8);
print(FeatureScatter(hc.unfiltered, feature1="nCount_RNA", feature2="nFeature_RNA",
                     cols=DiscretePalette(14)) +
      geom_hline(yintercept=500, col="#80808050") +
      geom_vline(xintercept=10000, col="#80808050"));
invisible(dev.off());

### Questions:
## * Did you compare results with the SCTransform integration alternative
##   (i.e. using SCTransform instead of NormalizeData)?

hc.list <- list();
## Takes ~1hr
for(cellOrigin in samplesToFetch){
    ## filtering mostly based on Vera's email, 2020-Sep-03
    ## hc.list[[i]] <- subset(hc.list[[i]],
    ##   subset = nFeature_RNA > 500 & nCount_RNA < 6000 &
    ##      percent_mt<0.5)
    subset(seur.obj, idents=cellOrigin,
           subset=nFeature_RNA > 500 & nCount_RNA < 10000 &
               percent.mt > 0.5 & percent.mt < 15) %>%
        SCTransform(verbose=TRUE, return.only.var.genes=FALSE) ->
        hc.list[[cellOrigin]];
}

pdf("Cluster_features_vs_counts_RNA_2020-09-16.pdf", width=11, height=8);
VlnPlot(seur.obj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol=3, pt.size=0);

pdf(sprintf("Cluster_features_vs_counts_orig_RNA_%s.pdf", Sys.Date()), width=11, height=8);
print(FeatureScatter(seur.obj, feature1="nCount_RNA", feature2="nFeature_RNA",
                     cols=DiscretePalette(14)))
invisible(dev.off());

## takes ~3 mins
Sys.time();
hc.features <-
    SelectIntegrationFeatures(object.list=hc.list,
                              nfeatures=3000, fvf.nfeatures=3000);
Sys.time();

## fast
hc.list <-
    PrepSCTIntegration(hc.list,
                       anchor.features=hc.features,
                       verbose=TRUE);

## takes ~20 mins
Sys.time();
hc.anchors <-
    FindIntegrationAnchors(hc.list, normalization.method="SCT",
                           anchor.features=hc.features,
                           verbose=TRUE);
Sys.time();

## takes ~10 mins
Sys.time();
hc.all <- IntegrateData(anchorset=hc.anchors,
                        normalization.method="SCT", verbose=TRUE);
Sys.time();

## Set default searched matrix to the integrated data
DefaultAssay(hc.all) <- "integrated";

## Okay, now we're back in business
## takes ~2 mins
## Determine principal components
hc.all <- RunPCA(hc.all, npcs=100);

## Check PCA

png(sprintf("PCA_SCT_HCall_10k_%s.png", format(Sys.Date())));
DimPlot(hc.all, reduction="pca", cols=DiscretePalette(17), pt.size=1);
## triangular shape?
invisible(dev.off());

## Takes about 4h
Sys.time();
hc.all <- JackStraw(hc.all, dims=100);
Sys.time();

hc.all <- ScoreJackStraw(hc.all, dims=1:100);

png(sprintf("JackStraw_SCT_HCall_10k_%s.png", format(Sys.Date())),
    width=1600, height=800);
JackStrawPlot(hc.all, dims=1:100, ymax=0.5, xmax=0.2);
invisible(dev.off());

plot(qqnorm(JS(hc.all[["pca"]], slot = "empirical")[,1]), pch="",
     xlim=c(0,0.5), ylim=c(0,0.2));
for(di in 1:10){
    points(qqnorm(JS(hc.all[["pca"]], slot = "empirical")[,di],
                  plot.it=FALSE), type="p", pch=20);
}

JS(hc.all[["pca"]], slot = "empirical") %>%
    as.data.frame() %>%
    as_tibble(rownames="Contig") %>%
    tidyr::pivot_longer(cols=-Contig, names_to="PC", values_to="Value") %>%
    filter(PC == "PC30") %>%
    pull(Value) %>%
    qqunif(logscale=FALSE);

JS(hc.all[["pca"]], slot = "empirical") %>%
    as.data.frame() %>%
    as_tibble(rownames="Contig") %>%
    tidyr::pivot_longer(cols=-Contig, names_to="PC", values_to="Value") %>%
    filter(as.numeric(sub("^PC","", PC)) == 15) %>%
    ggplot() +
    aes(sample=Value, colour=factor(PC)) +
    stat_qq(distribution=qunif) +
    xlim(0, 0.5) +
    ylim(0, 0.25) +
    coord_flip() +
    geom_abline(intercept=0, slope=1, linetype="dashed", na.rm=TRUE)

plot(y=seq(0,1,length.out=3000),
     x=sort(JS(hc.all[["pca"]], slot = "empirical")[,15]),
     xlim=c(0,0.25), ylim=c(0,0.5), type="n");
abline(h=0.5);
for(pi in 1:15){
    pf <- sort(JS(hc.all[["pca"]], slot = "empirical")[,pi]);
    pf <- head(pf, 1500);
    points(type="l", y=seq(0,1,length.out=3000)[1:length(pf)],
           x=pf);
}

abline(v=apply(JS(hc.all[["pca"]], slot = "empirical"), 2, quantile,
probs=0.5), col=rainbow(100));


pdf(sprintf("Threshold_values_JackStraw_%s.pdf", format(Sys.Date())));
par(mfrow=c(2,1));
plot(apply(JS(hc.all[["pca"]], slot = "empirical"), 2, quantile, probs=0.5),
     type="s",
     xlab="PC", ylab = "Empirical Value",
     main="Threshold values (50% of data)");
plot(apply(JS(hc.all[["pca"]], slot = "empirical"), 2, quantile, probs=0.3),
     xlab="PC", ylab = "Empirical Value", type="s",
     main="Threshold values (30% of data)");
invisible(dev.off());

## Show also PC heatmap, Elbow plot
png(sprintf("Elbow_SCT_HCall_10k_%s.png", format(Sys.Date())),
    width=800, height=800);
ElbowPlot(hc.all, ndims=100);
invisible(dev.off());

pdf(sprintf("DimLoadings_SCT_HCall_10k_%s.pdf", format(Sys.Date())),
    height=8, width=8, onefile=TRUE);
for(di in 1:100){
    print(VizDimLoadings(hc.all, dims=di, reduction="pca", nfeatures=30));
}
invisible(dev.off());

pdf(sprintf("DimHeatmaps_SCT_HCall_10k_%s.pdf", format(Sys.Date())),
    height=8, width=8, onefile=TRUE);
for(di in 1:100){
    print(DimHeatmap(hc.all, dims=di, balanced=TRUE, cells=200));
}
invisible(dev.off());

## Generate cluster labels
hc.all %>%
    FindNeighbors(dims=1:50) %>%
    FindClusters(resolution=0.8) -> hc.all

## Generate cluster colouring
cluster.cols <- paste0(DiscretePalette(length(unique(Idents(hc.all)))), "50");
names(cluster.cols) <- 1:length(cluster.cols);

## takes ~3 mins
Sys.time();
hc.all <- RunUMAP(object=hc.all, dims=1:50, n.components=3,
                  n.epochs=200);
Sys.time();

library(ggplot2);

png(sprintf("UMAP_Clustered_SCT_HCall_10k_200epochs_bySample_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.all, cols=paste0(DiscretePalette(7, palette="polychrome"), "50"),
        group.by="orig.ident",
        pt.size=1, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());

png(sprintf("UMAP_Clustered_SCT_HCall_10k_200epochs_byCluster_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.all, cols=cluster.cols,
        pt.size=1, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());


## takes ~3 mins
Sys.time();
hc.all <- RunUMAP(object=hc.all, dims=1:50, n.components=3,
                  n.epochs=1000);
Sys.time();

png(sprintf("UMAP_Clustered_SCT_HCall_10k_1000epochs_bySample_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.all, cols=paste0(DiscretePalette(7, palette="polychrome"), "50"),
        group.by="orig.ident",
        pt.size=1, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());

png(sprintf("UMAP_Clustered_SCT_HCall_10k_1000epochs_byCluster_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.all, cols=cluster.cols,
        pt.size=1, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());


## takes ~3 mins
Sys.time();
hc.all <- RunUMAP(object=hc.all, dims=1:50, n.components=3,
                  n.epochs=500);
Sys.time();

png(sprintf("UMAP_Clustered_SCT_HCall_10k_500epochs_bySample_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.all, cols=paste0(DiscretePalette(7, palette="polychrome"), "50"),
        group.by="orig.ident",
        pt.size=1, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());

png(sprintf("UMAP_Clustered_SCT_HCall_10k_500epochs_byCluster_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.all, cols=cluster.cols,
        pt.size=1, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());


## takes ~3 mins
Sys.time();
hc.all <- RunUMAP(object=hc.all, dims=1:50, n.components=3,
                  n.epochs=400);
Sys.time();

png(sprintf("UMAP_Clustered_SCT_HCall_10k_400epochs_bySample_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.all, cols=paste0(DiscretePalette(7, palette="polychrome"), "50"),
        group.by="orig.ident",
        pt.size=1, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());

png(sprintf("UMAP_Clustered_SCT_HCall_10k_400epochs_byCluster_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.all, cols=cluster.cols,
        pt.size=1, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());


png(sprintf("UMAP_Clustered_SCT_HCall_10k_400epochs_byCluster_%s_split.png",
            format(Sys.Date())), width=800*24, height=800);
DimPlot(hc.all, cols=cluster.cols, split.by="ident",
        pt.size=1, label=TRUE, label.size=10, combine=TRUE) +
    theme(text=element_text(size=30));
invisible(dev.off());


## [400 epochs looks the best in terms of separation of DCs]


missing <- c("IL5", "IL17A");

preferred.genes <-
    c("AIF1", "CD3D", "CD4", "CD8A",
      "FOXP3", "CD14", "HLA-DRA", "LYVE1",
      "CD207", "CD1A", "CD1C", "KRT5", "KRT1",
      "PMEL", "TPSAB1", "COL1A1", "ACTA2",
      "PECAM1", "IL13", "IL13RA1", "CD141"="THBD", "SIRPA",
      "IRF4", "IRF8", "CLEC9A",
      "AXL", "BATF3", "EPCAM", "TROP2"="TACSTD2", "CD301a"="CLEC10A",
      "VEGFA", "CD32A"="FCGR2A", "CLEC9A", "IDO1", "CCL22",
      "CCL17", "CXCL8", "MMP9", "FSCN1", "TNFRSF1B", "ANXA1",
      "ICOS", "IL1RL1", "IL17B", "IL17D", "IL18R1");
preferred.genes %>% unique %>% sort -> preferred.genes;

## Example for validating if a gene is present
grep("CD141", rownames(GetAssayData(hc.all, slot="data")), ignore.case=TRUE,
     value=TRUE)

grep("MT-ND1", rownames(GetAssayData(hc.all, slot="data")), ignore.case=TRUE,
     value=TRUE)

DefaultAssay(hc.all) <- "SCT";

## MIMR red: #e31836; MIMR orange: #f37a1f

png(sprintf("Feature_preferred_UMAP_SCT_HCall_10k_%s_%%02d.png",
           format(Sys.Date())),
    width=1600, height=1200);
for(pi in seq(1, length(preferred.genes), by=12)){
    pn <- preferred.genes[pi:min(length(preferred.genes), (pi+11))];
    print(FeaturePlot(hc.all, slot="data", col=c("lightgrey", "#e31836"),
                      features=pn, order=TRUE, pt.size=1));
}
invisible(dev.off());

png(sprintf("Feature_preferred_scale_UMAP_SCT_HCall_10k_%s_%%02d.png",
           format(Sys.Date())),
    width=1600, height=1200);
for(pi in seq(1, length(preferred.genes), by=12)){
    pn <- preferred.genes[pi:min(length(preferred.genes), (pi+11))];
    print(FeaturePlot(hc.all, slot="scale.data", col=c("lightgrey", "#e31836"),
                      features=pn, order=TRUE, pt.size=1));
}
invisible(dev.off());

hc.all %>%
    BuildClusterTree(features=hc.features) -> hc.all;
phyTree <- Tool(hc.all, slot="BuildClusterTree");
png(sprintf("phyTree_integrationSet_UMAP_SCT_HCall_10k_%s.png",
            format(Sys.Date())), width=600, height=600);
plot(phyTree, tip.color=sub("..$", "", cluster.cols));
invisible(dev.off());


hc.all %>%
    BuildClusterTree(features=preferred.genes) -> hc.all;
phyTree <- Tool(hc.all, slot="BuildClusterTree");
png(sprintf("phyTree_preferred_UMAP_SCT_HCall_10k_%s.png",
            format(Sys.Date())), width=600, height=600);
plot(phyTree, tip.color=sub("..$", "", cluster.cols));
invisible(dev.off());


sub.hc.12345 <- subset(hc.all, subset=orig.ident %in%
                                   c("HC1", "HC2", "HC3", "HC4", "HC5"));
sub.hc.67 <- subset(hc.all, subset=orig.ident %in% c("HC6", "HC7"));

## Zoomed in plots for all clusters
pdf(sprintf("UMAP_Clustered_HCall_10k_zoom_byCluster_%s.pdf",
            format(Sys.Date())), width=10, height=10);
for(ci in levels(Idents(hc.all))){
    print(ci);
    sub.hc <- subset(hc.all, ident=ci);
    print(DimPlot(sub.hc,
                  cols=cluster.cols,
                  pt.size=2, label=TRUE, label.size=10) +
          theme(text=element_text(size=30)))
}
invisible(dev.off());

## Determine descriptive genes
comp.list.20 <- NULL;
comp.list.10 <- NULL;
comp.list.5 <- NULL;
library(gdata);
diffFileName <-sprintf("differentiatingGenes_10k_byCluster_%s.txt", Sys.Date());
cat("** Most-distinguishing cluster genes **\n", file=diffFileName);
for(cl in levels(Idents(hc.all))){
    hc.sub <- GetAssayData(hc.all, slot="scale.data")[,Idents(hc.all) == cl];
    ##hc.sub <- hc.sub[apply(hc.sub,1,function(x){(mad(x) < 0.05)}),];
    hc.sub <- hc.sub[order(-apply(hc.sub,1,
                                  function(x){abs(median(x, na.rm=TRUE))})),];
    cat(sprintf("Cluster %s: \n", cl));
    cat(sprintf("Cluster %s: \n", cl), file=diffFileName, append=TRUE);
    res <- head(apply(hc.sub,1,median), 20);
    write.fwf(as.data.frame(res), rownames=TRUE, colnames=FALSE, append=TRUE,
              file=diffFileName);
    comp.list.20 <- c(comp.list.20, head(rownames(hc.sub), 20));
    comp.list.10 <- c(comp.list.10, head(rownames(hc.sub), 10));
    comp.list.5 <- c(comp.list.5, head(rownames(hc.sub), 5));
}

comp.list.20 <- unique(comp.list.20);
comp.list.10 <- unique(comp.list.10);
comp.list.5 <- unique(comp.list.5);

selected.cells <-
    c(sapply(levels(Idents(hc.all)),
           function(x){sample(which(Idents(hc.all) == x), 20)}));

## Plot 10 most descriptive genes for each cluster
pdf(sprintf("Cluster_compositions_HCall_10k_sampled20_%s.pdf", format(Sys.Date())),
    height=20, width=8);
DoHeatmap(hc.all, features=comp.list.10, cells=selected.cells,
          group.colors=sub("..$", "", cluster.cols),
          size=4);
invisible(dev.off());

## ## Plot IL13-positive cells on UMAP dimensions
## plot(Embeddings(sub.hc.67, reduction="umap")[
##     GetAssayData(sub.hc.67, slot="data")["IL13",] > 0, c("UMAP_1","UMAP_2")],
##     xlim=c(-15,15), ylim=c(-15,15))

rojahn.genes <-
    c("CD3D", "CD8A", "CD4", "IL22", "IL13", "GSMB", "NKG7", "CCL5",
      "DUT", "TYMS", "FOXP3", "CTLA4", "TNFRSF18", "CD14", "F13A1",
      "CCL18", "CCL13", "TPSAB1", "TPSAB2", "CD207", "FCGBP", "CD1A",
      "CD1C", "HLA-DRA", "ITGAX", "AIF1", "IL1B", "IRF8", "LAMP3",
      "CCL17", "CCL22", "CD80", "CD86", "CCR7", "CD83", "LILRA4",
      "MRC1", "CLEC4C", "TCF4", "KRT5", "KRT1", "KRT14", "COL7A1",
      "COL17A1", "KRT15", "CXCL14", "CCL22", "LY6D", "KRT10", "KRT2",
      "PCNA", "TK1", "TOP2A", "MKI67", "DEFB1", "GJB6", "ATP1B1",
      "AQP5", "DCD", "COL1A1", "FBLN1", "ACTA2", "MYL9", "PECAM1",
      "SELE", "PMEL", "MLANA", "KIT", "SOX10", "TYR", "MITF", "LYVE1",
      "CD3E", "CD28", "IL12RB2", "ALOX5", "PPARG", "RXRG", "IL23R",
      "CHAD", "TNFRSF21", "RORC", "TNFRSF21", "AS3MT", "LNPP4B",
      "GZMK");

pdf(sprintf("Cluster_compositions_HCall_10k_rojahn_sampled20_%s.pdf", format(Sys.Date())),
    height=20, width=8);
DoHeatmap(hc.all, features=rojahn.genes, cells=selected.cells,
          group.colors=sub("..$", "", cluster.cols),
          size=4, slot="data");
invisible(dev.off());

pdf(sprintf("HeatMap_Cluster_compositions_HCall_10k_scaleData_rojahn_sampled20_%s.pdf",
            format(Sys.Date())),
    height=20, width=8);
DoHeatmap(hc.all, features=rojahn.genes, cells=selected.cells,
          group.colors=sub("..$", "", cluster.cols),
          size=4, slot="scale.data");
invisible(dev.off());

png(sprintf("Feature_rojahn_scale_UMAP_SCT_HCall_10k_%s_%%02d.png",
           format(Sys.Date())),
    width=1600, height=1200);
for(pi in seq(1, length(rojahn.genes), by=12)){
    pn <- rojahn.genes[pi:min(length(rojahn.genes), (pi+11))];
    print(FeaturePlot(hc.all, slot="data", col=c("lightgrey", "#e31836"),
                      features=pn, order=TRUE, pt.size=1));
}
invisible(dev.off());

## Not found: GSMB, TPSAB2, LNPP4B

hc.all %>%
    BuildClusterTree(features=rojahn.genes) -> hc.all;
phyTree <- Tool(hc.all, slot="BuildClusterTree");
png(sprintf("phyTree_rojahnSet_UMAP_SCT_HCall_10k_%s.png",
            format(Sys.Date())), width=600, height=600);
plot(phyTree, tip.color=sub("..$", "", cluster.cols));
invisible(dev.off());

## After discussion with OL, remove cluster 11
## [It has demonstrable batch effects; HC5/7 cells are separated from others]
hc.no11 <- subset(hc.all, idents="11", invert=TRUE);

hc.all[["OL_clusters"]] <-
    clusters.tbl$label[match(Idents(hc.all), clusters.tbl$clusterID)];

## Save hc.all dataset
saveRDS(hc.all, file = "/mnt/BigBirdStornext/R_data/RStudioProjects/Seurat/JACI/JACI_10k_400epochs.rds");
## Save hc.no11 dataset
saveRDS(hc.no11, file = "/mnt/BigBirdStornext/R_data/RStudioProjects/Seurat/JACI/JACI_10k_400epochs_no11.rds");

## create differential expression table
Idents(hc.no11) <- hc.no11[["OL_clusters"]];

map_dfr(as.character(unique(Idents(hc.no11))), function(x){
    cat(x,"\n");
    FindMarkers(hc.no11, ident.1=x) %>%
        rownames_to_column("gene") %>%
        mutate(cluster=x);
}) -> marker.DE.tbl;

write_csv(marker.DE.tbl, "/mnt/BigBirdStornext/R_data/RStudioProjects/Seurat/JACI/DEtable_10k_no11.csv.gz");

## Fetch new cluster assignments
library(readxl);

read_excel("Cluster_identification.xlsx") %>%
    mutate(clusterID=sub("^.", "", `Cluster Number`),
           label=paste(Group, Number, sep="_")) %>%
    mutate(label=sub("_NA$", "", label)) -> clusters.tbl;

hc.no11 <- subset(hc.all, idents="11", invert=TRUE);

hc.no11[["OL_clusters"]] <-
    clusters.tbl$label[match(as.character(Idents(hc.no11)), as.character(clusters.tbl$clusterID))];

OL.cols <- cluster.cols[-10];
names(OL.cols) <- sort(unlist(unique(hc.no11[["OL_clusters"]])));

OL.cols.sameColour <- OL.cols;
OL.cols.sameColour[grep("^KC_", names(OL.cols.sameColour))] <-
    "#9DCC0050";
OL.cols.sameColour[grep("^FB_", names(OL.cols.sameColour))] <-
    "#005C3150";
OL.cols.sameColour[grep("^T-cell_[12]", names(OL.cols.sameColour))] <-
    "#C2008850";

hc.no11[["shape"]] <- factor(16);

cairo_pdf(sprintf("UMAP_Clustered_SCT_HCall_no11_10k_400epochs_byCluster_%s.pdf",
            format(Sys.Date())), width=16, height=10);
DimPlot(hc.no11, cols=OL.cols,
        group.by="OL_clusters",
        pt.size=1.5, label=TRUE, label.size=5, shape.by="shape") +
    theme(text=element_text(size=20)) +
    guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)))
invisible(dev.off());

## With recoloured KC, T-cell 1/2, FB
cairo_pdf(sprintf("UMAP_Clustered_SCT_HCall_mergedCols_no11_10k_400epochs_byCluster_%s.pdf",
            format(Sys.Date())), width=16, height=10);
DimPlot(hc.no11, cols=OL.cols.sameColour,
        group.by="OL_clusters",
        pt.size=1.5, label=TRUE, label.size=5, shape.by="shape") +
    theme(text=element_text(size=20)) +
    guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)))
invisible(dev.off());

## JM gene set
DefaultAssay(hc.no11) <- "SCT";
cairo_pdf(sprintf("Feature_Rojahn_Clustered_SCT_HCall_no11_10k_400epochs_%s.pdf",
                  format(Sys.Date())), width=12, height=20);
rojahn.genes <- c("CD3D", "CD4", "CD8A", "FOXP3", "CD14", "CD207",
                  "CD1A", "CD1C", "KRT5", "KRT1", "PMEL", "TPSAB1",
                  "COL1A1", "ACTA2", "PECAM1");
res <- lapply(FeaturePlot(hc.no11, slot="data", col=c("lightgrey", "#e31836"),
                          features=rojahn.genes, order=TRUE, pt.size=1, shape.by="shape",
                          combine=FALSE),
              function(x){x + guides(shape=FALSE)});
print(wrap_plots(res, ncol=3));
invisible(dev.off());

## OL5 gene set
DefaultAssay(hc.no11) <- "SCT";
cairo_pdf(sprintf("Feature_OL5_Clustered_SCT_HCall_no11_10k_400epochs_%s.pdf",
                  format(Sys.Date())), width=12, height=24);
gene.set.5  <- c("HLA-DRA", "CD14", "CD207",
                 "XCR1", "CD1C", "CCR7",
                 "CD3D", "CD4", "CD8A",
                 "FOXP3", "NKG7", "TPSAB1",
                 "KRT5", "LYVE1", "COL1A1",
                 "ACTA2", "PECAM1", "PMEL");
res <- lapply(FeaturePlot(hc.no11, slot="data", col=c("lightgrey", "#e31836"),
                          features=gene.set.5, order=TRUE, pt.size=1, shape.by="shape",
                          combine=FALSE),
              function(x){x + guides(shape=FALSE)});
print(wrap_plots(res, ncol=3));
invisible(dev.off());


## OL gene set
gene.set.1 <- c("HLA-DRA", "ITGAX", "CD14", "F13A1", "CD207", "CD1C",
                "IRF4", "XCR1", "CLEC9A", "IRF8", "CCR7", "LAMP3",
                "CCL17", "CCL22");
gene.set.2 <- c("HLA-DRA", "ITGAX", "CD14", "F13A1", "CD207", "CD1C",
                "XCR1", "CLEC9A", "CCR7", "LAMP3", "IL4R", "IL13RA1",
                "CCL22", "RHOU", "CXCL8", "FOXO3", "BATF", "MAOA",
                "MMP9", "IL23A", "VEGFA", "JAK1", "LAMA5", "FSCN1",
                "IRF4", "MMP1", "SOCS1", "HMOX1", "FOXO1", "GATA3",
                "IL6", "MMP3", "BCL2", "IL1A", "BCL2L1", "TNFRSF1B",
                "MUC1", "IL10", "ANXA1");
gene.set.3 <- c("CCL2", "FOS", "SOCS3", "F13A1", "PIM1", "CEBPD",
                "TIMP1", "CD36", "IL4R", "POU2F1", "CDKN1A", "STAT1",
                "BATF", "TGFB1", "RHOU", "TYK2", "CXCL8", "HIF1A",
                "FOXO3", "JUNB", "MCL1", "STAT3", "ITGAM", "PIK3R1",
                "FN1", "CCL22");
gene.set.4 <- c("STAT6", "KLF4", "IL4R", "IL13RA1");
DefaultAssay(hc.no11) <- "SCT";
cairo_pdf(sprintf("Feature_OL_GS1_Clustered_SCT_HCall_no11_10k_400epochs_%s.pdf",
                  format(Sys.Date())), width=12, height=20);
res <- lapply(FeaturePlot(hc.no11, slot="data", col=c("lightgrey", "#e31836"),
                          features=gene.set.1, order=TRUE, pt.size=1, shape.by="shape",
                          combine=FALSE),
              function(x){x + guides(shape=FALSE)});
print(wrap_plots(res, ncol=3));
invisible(dev.off());


DefaultAssay(hc.no11) <- "SCT";
cairo_pdf(sprintf("Feature_OL_GS4_Clustered_SCT_HCall_no11_10k_400epochs_%s.pdf",
                  format(Sys.Date())), width=12, height=12, onefile=TRUE);
for(gene in gene.set.4){
    (FeaturePlot(hc.no11, slot="data", col=c("lightgrey", "#e31836"),
                features=gene, order=TRUE, pt.size=1, shape.by="shape") + guides(shape=FALSE)) %>% print();
}
invisible(dev.off());

DefaultAssay(hc.no11) <- "SCT";
plot.cells <- which(unlist(hc.no11[["OL_clusters"]]) %in%
                    c("Monocytes", "DC_LCs", "DC_cDC2", "DC_cDC1", "DC_Mature DCs"));
cairo_pdf(sprintf("Heatmap_OL_GS2_Clustered_SCT_HCall_no11_10k_400epochs_%s.pdf",
                  format(Sys.Date())), width=12, height=20);
print(DoHeatmap(hc.no11, slot="data", raster=TRUE, features=gene.set.2,
                group.colors=sub("..$", "", cluster.cols.TCs),
                group.by="OL_clusters", cells=plot.cells) +
      scale_fill_viridis());
invisible(dev.off());
png(sprintf("Heatmap_OL_GS2_Clustered_SCT_HCall_no11_10k_400epochs_%s.png",
                  format(Sys.Date())), width=1200, height=2000);
print(DoHeatmap(hc.no11, slot="data", raster=TRUE, features=gene.set.2,
                group.colors=sub("..$", "", cluster.cols.TCs),
                group.by="OL_clusters", cells=plot.cells) +
      scale_fill_viridis());
invisible(dev.off());

DefaultAssay(hc.no11) <- "SCT";
plot.cells <- which(unlist(hc.no11[["OL_clusters"]]) %in%
                    c("Monocytes", "DC_LCs", "DC_cDC2", "DC_cDC1", "DC_Mature DCs"));
cairo_pdf(sprintf("Heatmap_OL_GS3_Clustered_SCT_HCall_no11_10k_400epochs_%s.pdf",
                  format(Sys.Date())), width=12, height=20);
print(DoHeatmap(hc.no11, slot="data", raster=TRUE, features=gene.set.3,
                group.colors=sub("..$", "", cluster.cols.TCs),
                group.by="OL_clusters", cells=plot.cells) +
      scale_fill_viridis());
invisible(dev.off());
png(sprintf("Heatmap_OL_GS3_Clustered_SCT_HCall_no11_10k_400epochs_%s.png",
                  format(Sys.Date())), width=1200, height=2000);
print(DoHeatmap(hc.no11, slot="data", raster=TRUE, features=gene.set.3,
                group.colors=sub("..$", "", cluster.cols.TCs),
                group.by="OL_clusters", cells=plot.cells) +
      scale_fill_viridis());
invisible(dev.off());


unique(unlist(hc.no11[["OL_clusters"]]));

png(sprintf(paste0("UMAP_Clustered_SCT_HCall_no11_labelOnly_",
                   "10k_400epochs_byCluster_%s.png"),
            format(Sys.Date())), width=1600, height=1000);
DimPlot(hc.no11, cols=rep("white", 23),
        group.by="OL_clusters",
        pt.size=1, label=TRUE, label.size=7) +
    theme(text=element_text(size=30))
invisible(dev.off());

png(sprintf(paste0("UMAP_Clustered_SCT_HCall_no11_nolabel_",
                   "10k_400epochs_byCluster_%s.png"),
            format(Sys.Date())), width=1600, height=1000);
DimPlot(hc.no11, cols=OL.cols,
        group.by="OL_clusters",
        pt.size=1, label=FALSE) +
    theme(text=element_text(size=30))
invisible(dev.off());


il4.il13.genes <- c("IL4", "IL13", "IL4R", "IL13RA1");
cairo_pdf(sprintf("Feature_IL4_IL13_UMAP_SCT_HCall_no11_10k_%s.pdf",
           format(Sys.Date())),
    width=12, height=12);
print(FeaturePlot(hc.no11, slot="data", col=c("lightgrey", "#e31836"),
                  features=il4.il13.genes, order=TRUE, pt.size=1, shape.by="shape") +
          guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16))));
invisible(dev.off());

png(sprintf("Ridge_IL4_IL13_UMAP_SCT_HCall_no11_10k_%s.png",
           format(Sys.Date())),
    width=1600, height=1200);
print(RidgePlot(hc.no11, slot="data", group.by="OL_clusters",
                col=OL.cols,
                features=il4.il13.genes));
invisible(dev.off());

png(sprintf("Violin_IL4_IL13_UMAP_SCT_HCall_no11_10k_%s.png",
           format(Sys.Date())),
    width=1600, height=1200);
print(VlnPlot(hc.no11, slot="data", group.by="OL_clusters",
                col=OL.cols,
                features=il4.il13.genes));
invisible(dev.off());

### Process DC subsets

hc.DCs <- subset(hc.all, idents=c("4", "9", "19", "22", "23"));
hc.DCs[["hc.cluster"]] <- Idents(hc.DCs);
DefaultAssay(hc.DCs) <- "integrated";
## Determine principal components
hc.DCs <- RunPCA(hc.DCs, npcs=100);
## Takes 8 mins
hc.DCs <- JackStraw(hc.DCs, dims=30);
## Add JackStraw Scores
hc.DCs <- ScoreJackStraw(hc.DCs, dims=1:30);

png(sprintf("PCA_SCT_HCall_DCs_10k_%s.png", format(Sys.Date())));
DimPlot(hc.DCs, reduction="pca", cols=DiscretePalette(17), pt.size=1);
## Yay, that's a good PCA plot
invisible(dev.off());

png(sprintf("JackStraw_SCT_HCall_DCs_10k_%s.png", format(Sys.Date())),
    width=1600, height=800);
JackStrawPlot(hc.DCs, dims=1:30, ymax=0.3, xmax=0.1);
invisible(dev.off());

## Run UMAP for 1000 epochs
Sys.time();
hc.DCs <- RunUMAP(object=hc.DCs, dims=1:15, n.components=3,
                  n.epochs=1000);
Sys.time();

## Generate cluster labels
hc.DCs %>%
    FindNeighbors(dims=1:15) %>%
    FindClusters(resolution=0.8) -> hc.DCs

## Tack on original clusters
hc.DCs[["comb.cluster"]] <-
    paste0("CL", unlist(hc.DCs[["hc.cluster"]]),
           "_", Idents(hc.DCs));
hcc <- unlist(hc.DCs[["comb.cluster"]]);

highCount.cells <-
    which(table(hcc)[hcc] > 30);

## Generate cluster colouring
cluster.cols.DCs <- paste0(DiscretePalette(length(unique(hcc))), "50");
names(cluster.cols.DCs) <- unique(hcc);
png(sprintf("UMAP_Clustered_SCT_HCall_DCs_10k_1000epochs_byCluster_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.DCs, cols=cluster.cols.DCs,
        cells=highCount.cells,
        group.by="comb.cluster",
        pt.size=3, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());

selected.cells <-
    c(sapply(unique(names(highCount.cells)),
             function(x){sample(which(hc.DCs[["comb.cluster"]] == x),
                                25)}));
hc.DCs.sub <- subset(hc.DCs, cells=selected.cells);
hccs <- unlist(hc.DCs.sub[["comb.cluster"]]);
DC.assayData <- GetAssayData(hc.DCs.sub, slot="scale.data");

comp.list.DCs.10 <- NULL;
comp.list.DCs.20 <- NULL;
diffFileName <-
    sprintf("differentiatingGenes_10k_DCs_byCluster_%s.txt", Sys.Date());
cat("** Most-distinguishing cluster genes **\n", file=diffFileName);
for(cl in unique(hccs)){
    print(cl);
    hc.sub <- apply(DC.assayData[,hccs == cl], 1, median, na.rm=TRUE);
    hc.other <- apply(DC.assayData[,hccs != cl], 1, median, na.rm=TRUE);
    gene.diffRanked <- hc.sub - hc.other;
    cat(sprintf("Cluster %s: \n", cl), file=diffFileName, append=TRUE);
    names(gene.diffRanked) <- names(hc.sub);
    gene.diffRanked <- gene.diffRanked[order(-abs(gene.diffRanked))];
    comp.list.DCs.10 <- c(comp.list.DCs.10, names(head(gene.diffRanked, 10)));
    res <- head(gene.diffRanked, 20);
    comp.list.DCs.20 <- c(comp.list.DCs.20, names(res));
    write.fwf(as.data.frame(res), rownames=TRUE, colnames=FALSE, append=TRUE,
              file=diffFileName);
}
comp.list.DCs.10 <- unique(comp.list.DCs.10);
comp.list.DCs.20 <- unique(comp.list.DCs.20);

features.DCs <- tail(readLines("OL_list_DCs.txt"), -1);

for(pi in seq(1, length(preferred.genes), by=12)){
    pn <- preferred.genes[pi:min(length(preferred.genes), (pi+11))];
    print(FeaturePlot(hc.all, slot="data", col=c("lightgrey", "#e31836"),
                      features=pn, order=TRUE, pt.size=1));
}
invisible(dev.off());

png(sprintf("Feature_OL_DClist_UMAP_SCT_HCall_DCs_10k_%s_%%02d.png",
           format(Sys.Date())),
    width=1600, height=1200);
DefaultAssay(hc.DCs) <- "SCT";
DC.assayData <- GetAssayData(hc.DCs.sub, slot="data");
for(pi in seq(1, length(features.DCs), by=12)){
    pn <- features.DCs[pi:min(length(features.DCs), (pi+11))];
    print(FeaturePlot(hc.DCs, slot="data", col=c("lightgrey", "#e31836"),
                      features=pn, order=TRUE, pt.size=1));
}
invisible(dev.off());

### Process T-Cell subsets

## plot IL13

pdf(sprintf("FeaturePlot_T-cells_IL13_%s.pdf", Sys.Date()));
DefaultAssay(hc.TCs) <- "SCT";
print(FeaturePlot(hc.TCs, slot="data", col=c("lightgrey", "#e31836"),
                features=c("IL13"), order=TRUE));
dev.off();

pdf(sprintf("RidgePlot_T-cells_%s.pdf", Sys.Date()));
DefaultAssay(hc.TCs) <- "SCT";
print(RidgePlot(hc.TCs, slot="data", col=c("lightgrey", "#e31836"),
                features=c("IL13"), order=TRUE));
dev.off();

hc.TCs <- subset(hc.all, idents=c("14","1","13","17"));
hc.TCs[["hc.cluster"]] <- Idents(hc.TCs);
DefaultAssay(hc.TCs) <- "integrated";
## Determine principal components
hc.TCs <- RunPCA(hc.TCs, npcs=100);
## Takes 8 mins
hc.TCs <- JackStraw(hc.TCs, dims=30);
## Add JackStraw Scores
hc.TCs <- ScoreJackStraw(hc.TCs, dims=1:30);

png(sprintf("PCA_SCT_HCall_TCs_10k_%s.png", format(Sys.Date())));
DimPlot(hc.TCs, reduction="pca", cols=DiscretePalette(17), pt.size=1);
## Yay, that's a good PCA plot
invisible(dev.off());

png(sprintf("JackStraw_SCT_HCall_TCs_10k_%s.png", format(Sys.Date())),
    width=1600, height=800);
JackStrawPlot(hc.TCs, dims=1:30, ymax=0.3, xmax=0.1);
invisible(dev.off());

## Run UMAP for 1000 epochs
Sys.time();
hc.TCs <- RunUMAP(object=hc.TCs, dims=1:15, n.components=3,
                  n.epochs=1000);
Sys.time();

## Generate cluster labels
hc.TCs %>%
    FindNeighbors(dims=1:15) %>%
    FindClusters(resolution=0.8) -> hc.TCs

## Tack on original clusters
hc.TCs[["comb.cluster"]] <-
    paste0("CL", unlist(hc.TCs[["hc.cluster"]]),
           "_", Idents(hc.TCs));
hcc <- unlist(hc.TCs[["comb.cluster"]]);

highCount.cells <-
    which(table(hcc)[hcc] > 30);

## Generate cluster colouring
cluster.cols.TCs <- paste0(DiscretePalette(length(unique(hcc))), "50");
names(cluster.cols.TCs) <- unique(hcc);
png(sprintf("UMAP_Clustered_SCT_HCall_TCs_10k_1000epochs_byCluster_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.TCs, cols=cluster.cols.TCs,
        cells=highCount.cells,
        group.by="comb.cluster",
        pt.size=3, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());

selected.cells <-
    c(sapply(unique(names(highCount.cells)),
             function(x){sample(which(hc.TCs[["comb.cluster"]] == x),
                                25)}));
hc.TCs.sub <- subset(hc.TCs, cells=selected.cells);
hccs <- unlist(hc.TCs.sub[["comb.cluster"]]);
TC.assayData <- GetAssayData(hc.TCs.sub, slot="scale.data");

comp.list.TCs.10 <- NULL;
comp.list.TCs.20 <- NULL;
diffFileName <-
    sprintf("differentiatingGenes_10k_TCs_byCluster_%s.txt", Sys.Date());
cat("** Most-distinguishing cluster genes **\n", file=diffFileName);
for(cl in unique(hccs)){
    print(cl);
    hc.sub <- apply(TC.assayData[,hccs == cl], 1, median, na.rm=TRUE);
    hc.other <- apply(TC.assayData[,hccs != cl], 1, median, na.rm=TRUE);
    gene.diffRanked <- hc.sub - hc.other;
    cat(sprintf("Cluster %s: \n", cl), file=diffFileName, append=TRUE);
    names(gene.diffRanked) <- names(hc.sub);
    gene.diffRanked <- gene.diffRanked[order(-abs(gene.diffRanked))];
    comp.list.TCs.10 <- c(comp.list.TCs.10, names(head(gene.diffRanked, 10)));
    res <- head(gene.diffRanked, 20);
    comp.list.TCs.20 <- c(comp.list.TCs.20, names(res));
    write.fwf(as.data.frame(res), rownames=TRUE, colnames=FALSE, append=TRUE,
              file=diffFileName);
}
comp.list.TCs.10 <- unique(comp.list.TCs.10);
comp.list.TCs.20 <- unique(comp.list.TCs.20);

pdf(sprintf("Cluster_compositions_HCall_TCs_10k_sampled25_%s.pdf",
            format(Sys.Date())),
    height=20, width=8);
DoHeatmap(hc.TCs, cells=highCount.cells,
          features=comp.list.TCs.10,
          group.colors=sub("..$", "", cluster.cols.TCs),
          group.by="comb.cluster",
          size=3, raster=FALSE);
invisible(dev.off());

features.TCs <- tail(readLines("OL_list_TCs.txt"), -1);

png(sprintf("Feature_OL_TClist_UMAP_SCT_HCall_TCs_10k_%s_%%02d.png",
           format(Sys.Date())),
    width=1600, height=1200);
DefaultAssay(hc.TCs) <- "SCT";
for(pi in seq(1, length(features.TCs), by=12)){
    pn <- features.TCs[pi:min(length(features.TCs), (pi+11))];
    print(FeaturePlot(hc.TCs, slot="data", col=c("lightgrey", "#e31836"),
                      features=pn, order=TRUE, pt.size=1));
}
invisible(dev.off());


## Additional experimentation with cluster 11

hc.cl11s <- subset(hc.all, idents=c("11"));
DefaultAssay(hc.cl11s) <- "integrated";
## Determine principal components
hc.cl11s <- RunPCA(hc.cl11s, npcs=100);
## Takes 8 mins
hc.cl11s <- JackStraw(hc.cl11s, dims=30);
## Add JackStraw Scores
hc.cl11s <- ScoreJackStraw(hc.cl11s, dims=1:30);

png(sprintf("PCA_SCT_HCall_cl11s_10k_%s.png", format(Sys.Date())));
DimPlot(hc.cl11s, reduction="pca", cols=DiscretePalette(17), pt.size=1);
## Also looks reasonable
invisible(dev.off());

png(sprintf("JackStraw_SCT_HCall_cl11s_10k_%s.png", format(Sys.Date())),
    width=1600, height=800);
JackStrawPlot(hc.cl11s, dims=1:30, ymax=0.3, xmax=0.1);
invisible(dev.off());

## Generate cluster labels
hc.cl11s %>%
    FindNeighbors(dims=1:15) %>%
    FindClusters(resolution=0.8) -> hc.cl11s
## Generate cluster colouring
cluster.cols.cl11s <- paste0(DiscretePalette(length(unique(Idents(hc.cl11s)))), "50");
names(cluster.cols.cl11s) <- 1:length(cluster.cols.cl11s);

## Run UMAP for 4000 epochs
Sys.time();
hc.cl11s <- RunUMAP(object=hc.cl11s, dims=1:20, n.components=3,
                  n.epochs=2000);
Sys.time();

png(sprintf("UMAP_Clustered_SCT_HCall_cl11s_10k_2000epochs_byCluster_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.cl11s, cols=cluster.cols.cl11,
        pt.size=3, label=TRUE, label.size=10) +
    theme(text=element_text(size=30))
invisible(dev.off());

png(sprintf("UMAP_Clustered_SCT_HCall_cl11s_10k_2000epochs_byIdent_%s.png",
            format(Sys.Date())), width=1024, height=1024);
DimPlot(hc.cl11s,
        pt.size=2, label=TRUE, label.size=10, group.by="orig.ident") +
    theme(text=element_text(size=30))
invisible(dev.off());

cl11.cells <- which(Idents(hc.all) == "11");

hc.cl11 <- subset(hc.all, idents="11");
hc.cl14 <- subset(hc.all, idents="14");
hc.cl9 <- subset(hc.all, idents="9");
cl11.mat <-  GetAssayData(subset(hc.all, idents="11"), slot="scale.data");
cl14.mat <-  GetAssayData(subset(hc.all, idents="14"), slot="scale.data");
cl7.mat <-  GetAssayData(subset(hc.all, idents="7"), slot="scale.data");

pdf(sprintf("plot_cluster11_10k_%s.pdf",
            format(Sys.Date())), width=12, height=6);
par(mfrow=c(2,3));
plot(density(cl11.mat["KRT1",]));
plot(density(cl14.mat["KRT1",]));
plot(density(cl7.mat["KRT1",]));
plot(density(cl11.mat["KRT10",]));
plot(density(cl14.mat["KRT10",]));
plot(density(cl7.mat["KRT10",]));
invisible(dev.off());

png(sprintf("Feature_preferred_UMAP_SCT_cluster11_10k_%s.png",
            format(Sys.Date())), width=1024, height=1024);

FeaturePlot(hc.all, features="MT-CO3", cells=cl11.cells,
            pt.size=4, slot="data", order=TRUE);
invisible(dev.off());

cl11.krt <- cl11.mat[,cl11.mat["KRT1",] > 10];

diff.val <- (apply(cl7.mat,1,median) - apply(cl11.krt,1,median));

diff.val <- diff.val[order(abs(diff.val))];

tail(diff.val, 20);

## IL13-expressing cells - barplot of clusters
cells.IL13 <- which(GetAssayData(hc.no11)["IL13",] > 0);
hc.no11.il13 <- subset(hc.no11, cells=cells.IL13);


tapply(GetAssayData(hc.no11.il13)["IL13",], hc.no11.il13[["OL_clusters"]], median);

DefaultAssay(hc.no11.il13) <- "RNA";
bind_cols(
    tibble(cluster=tapply(unlist(hc.no11.il13[["OL_clusters"]]), hc.no11.il13[["OL_clusters"]], first)),
    count=tapply(GetAssayData(hc.no11.il13)["IL13",], hc.no11.il13[["OL_clusters"]], length),
    median=tapply(GetAssayData(hc.no11.il13)["IL13",], hc.no11.il13[["OL_clusters"]], median),
    min=tapply(GetAssayData(hc.no11.il13)["IL13",], hc.no11.il13[["OL_clusters"]], min),
    max=tapply(GetAssayData(hc.no11.il13)["IL13",], hc.no11.il13[["OL_clusters"]], max)
) -> metadata.il13;

## pie chart

## Compute the position of labels
metadata.il13 %>%
    arrange(-count) %>%
    mutate(pct=count / sum(count) * 100) %>%
    mutate(ypos=cumsum(pct) - 0.5 * pct) %>%
    mutate(bar.col=OL.cols[cluster])-> bar.data

bar.data$cluster <- factor(bar.data$cluster, levels=bar.data$cluster)

## Basic barplot
bar.data %>%
    ggplot() +
    aes(x=cluster, y=count, fill=cluster) +
    scale_fill_manual(values=sub("..$","",bar.data$bar.col)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf("IL13_counts_barPlot_SCT_HCall_no11_10k_%s.pdf",
               Sys.Date()));

## Pie chart
## see https://www.r-graph-gallery.com/piechart-ggplot2.html
bar.data %>%
    ggplot() +
    aes(x="", y=pct, fill=cluster) +
    scale_fill_manual(values=sub("..$","",bar.data$bar.col)) +
    geom_col(width=1, color="grey90") +
    coord_polar("y", start=0) +
    theme_minimal() +
    theme(legend.position="bottom",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.x=element_blank(),
          axis.title=element_blank())
ggsave(sprintf("IL13_counts_pieChart_SCT_HCall_no11_10k_%s.pdf",
               Sys.Date()));

pdf(sprintf("UMAP_Clustered_SCT_HCall_no11_10k_400epochs_byCluster_%s.pdf",
            format(Sys.Date())), width=16, height=10);
DimPlot(hc.no11, cols=OL.cols,
        group.by="OL_clusters",
        pt.size=1, label=TRUE, label.size=5) +
    theme(text=element_text(size=20))
invisible(dev.off());

hc.no11[["shape"]] <- factor(16);
hc.TCs.noCluster <- subset(hc.no11, idents=c("14","1","13","17"));
hc.DCs.noCluster <- subset(hc.no11, idents=c("4", "9", "19", "22", "23"));

cairo_pdf(sprintf("UMAP_Clustered_SCT_HCall_no11_10k_400epochs_DCzoom_byCluster_%s.pdf",
            format(Sys.Date())), width=12, height=8);
DimPlot(hc.DCs.noCluster, cols=OL.cols,
        group.by="OL_clusters",
        pt.size=2, label=TRUE, label.size=5, shape.by="shape") +
    xlim(-12, -6) +
    ylim(-3, 2.5) +
    theme(text=element_text(size=20)) +
    guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)))
invisible(dev.off());

### Feature plot with IL13RA1, IL13RA2, and IL4RA
cairo_pdf(sprintf("Feature_IL4_IL13_SCT_HCall_no11_10k_400epochs_DCzoom_%s.pdf",
            format(Sys.Date())), width=10, height=8);
hc.DCs.zoom <- subset(hc.no11,
                      subset=((UMAP_1 > -12) & (UMAP_1 < -6) &
                              (UMAP_2 > -3) & (UMAP_2 < 2.5)));
DefaultAssay(hc.DCs.zoom) <- "SCT";
res <- lapply(FeaturePlot(hc.DCs.zoom, slot="data", features=il4.il13.genes,
                          col=c("lightgrey", "#e31836"), shape.by="shape",
                          pt.size=1, order=TRUE, label.size=5, combine=FALSE),
              function(x){x + guides(shape=FALSE)});
print(res[[1]] + res[[2]] + res[[3]] + res[[4]]);
invisible(dev.off());

cairo_pdf(sprintf("Feature_IL4_IL13_SCT_HCall_no11_10k_400epochs_largeUMAP_%s.pdf",
            format(Sys.Date())), width=10, height=8, onefile=T);
DefaultAssay(hc.no11) <- "SCT";
res <- lapply(FeaturePlot(hc.no11, slot="data", features=il4.il13.genes,
                  col=c("lightgrey", "#e31836"), shape.by="shape",
                  pt.size=1, order=TRUE, label.size=5, combine=FALSE),
              function(x){x + guides(shape=FALSE)});
print(res[[1]] + res[[2]] + res[[3]] + res[[4]]);
invisible(dev.off());

library(tidyverse);

IL4.IL13.counts <-
    GetAssayData(hc.no11)[il4.il13.genes,] %>%
    as.matrix() %>%
    t() %>%
    as_tibble() %>%
    bind_cols(hc.no11[["OL_clusters"]]) %>%
    pivot_longer(names_to="Gene", values_to="expression", cols=-starts_with("OL_clusters")) %>%
    mutate(expression=expression > 0) %>%
    group_by(OL_clusters, Gene) %>%
    summarise(count=sum(expression)) %>%
    mutate(bar.col=OL.cols[OL_clusters])

IL4.IL13.counts %>%
    select(-bar.col) %>%
    pivot_wider(names_from=Gene, values_from=count) %>%
    write_csv(sprintf("counts_IL4_IL13_SCT_HCall_no11_10k_%s.csv",
                      Sys.Date()));

pdf(sprintf("IL4_IL13_counts_barPlot_SCT_HCall_no11_10k_%s.pdf", Sys.Date()),
    width=16, height=8);
IL4.IL13.counts %>%
    ggplot() +
    aes(x = OL_clusters, y=count, fill=OL_clusters) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values=sub("..$","",unique(IL4.IL13.counts$bar.col))) +
    facet_wrap(~ Gene)
invisible(dev.off());

hc.TCs.noCluster[["shape"]] <- factor(19);

cairo_pdf(sprintf("UMAP_Clustered_SCT_HCall_no11_10k_400epochs_TCzoom_byCluster_%s.pdf",
            format(Sys.Date())), width=10, height=10);
DimPlot(hc.TCs.noCluster, cols=OL.cols,
        group.by="OL_clusters",
        pt.size=2, label=TRUE, label.size=5, shape.by="shape") +
    guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16))) +
    xlim(-6.5, -1) +
    ylim(-3, 7) +
    theme(text=element_text(size=20))
invisible(dev.off());

### Attempt clustering with IL13 from cluster 1

Cluster.TC1 <- subset(hc.no11, subset=(OL_clusters=="T-cell_1"));
DefaultAssay(Cluster.TC1) <- "SCT";
cells.TC1.IL13 <- which(GetAssayData(Cluster.TC1)["IL13",] > 0);
cells.TC1.noIL13 <- which(GetAssayData(Cluster.TC1)["IL13",] == 0);
Idents(Cluster.TC1) <- ifelse(GetAssayData(Cluster.TC1)["IL13",] > 0, "IL13pos", "IL13neg");

res <- FindMarkers(Cluster.TC1, ident.1="IL13pos") %>%
    filter(p_val_adj < 0.01) %>%
    arrange(-abs(avg_logFC))
write_csv(res,
          sprintf("FindMarkers_Cluster1_IL13_wilcox_SCT_HCall_no11_10k_400epochs_%s.csv",
                  Sys.Date()));

res <- FindMarkers(Cluster.TC1, ident.1="IL13pos", test.use="DESeq2") %>%
    filter(p_val_adj < 0.1) %>%
    arrange(-abs(avg_logFC))

write_csv(res,
          sprintf("FindMarkers_Cluster1_SCT_HCall_no11_10k_400epochs_%s.csv",
                  Sys.Date()));


Cluster.TC1.IL13only <- subset(hc.no11, subset=(OL_clusters=="T-cell_1"),
                               cells=which(GetAssayData(hc.no11)["IL13",] > 0));
ncol(Cluster.TC1.IL13only);

Cluster.TC1.IL13only %>%
    FindVariableFeatures(selection.method="mvp") %>%
    RunPCA(npcs=20) %>%
    RunUMAP(dims=1:20, n.components=3, n.epochs=400) %>%
    FindNeighbors(dims=1:20) %>%
    FindClusters(resolution=0.8) -> Cluster.TC1.IL13only;

pdf(sprintf("DimPlot_HeatMap_Cluster_mvp_UMAP_400epochs_TC1.IL13only_%s.pdf", Sys.Date()));
library(RColorBrewer);
DimPlot(Cluster.TC1.IL13only,
        cols=brewer.pal(7, "Dark2"),
        group.by="orig.ident", reduction="pca",
        pt.size=1, label=TRUE, label.size=5) +
    theme(text=element_text(size=15))
DimPlot(Cluster.TC1.IL13only,
        cols=brewer.pal(7, "Dark2"),
        reduction="pca",
        pt.size=1, label=TRUE, label.size=5) +
    theme(text=element_text(size=15))
DimPlot(Cluster.TC1.IL13only,
        cols=brewer.pal(7, "Dark2"),
        group.by="orig.ident", reduction="umap",
        pt.size=1, label=TRUE, label.size=5) +
    theme(text=element_text(size=15))
DimPlot(Cluster.TC1.IL13only,
        cols=brewer.pal(7, "Dark2"),
        reduction="umap",
        pt.size=1, label=TRUE, label.size=5) +
    theme(text=element_text(size=15))
res <- lapply(FeaturePlot(Cluster.TC1.IL13only, slot="data", features=c("TRAC", "CD3D", "TMSB4X", "RPL21"),
                          col=c("lightgrey", "#e31836"), shape.by="shape",
                          reduction="pca",
                          pt.size=1, order=TRUE, label.size=5, combine=FALSE),
              function(x){x + guides(shape=FALSE)});
print(wrap_plots(res, ncol=2));
res <- lapply(FeaturePlot(Cluster.TC1.IL13only, slot="data", features=c("TRAC", "CD3D", "TMSB4X", "RPL21"),
                          col=c("lightgrey", "#e31836"), shape.by="shape",
                          reduction="umap",
                          pt.size=1, order=TRUE, label.size=5, combine=FALSE),
              function(x){x + guides(shape=FALSE)});
print(wrap_plots(res, ncol=2));
DoHeatmap(Cluster.TC1.IL13only, slot="data",
        features=c("TRAC", "CD3D", "TMSB4X", "RPL21"));
invisible(dev.off());

IL13.assayData <- GetAssayData(Cluster.TC1);
comp.list.IL13.100 <- NULL;
diffFileName <-
    sprintf("differentiatingGenes_10k_Cluster1_IL13_byCluster_%s.txt", Sys.Date());
cat("** Most-distinguishing IL13 / cluster 1 genes **\n", file=diffFileName);
for(bootStrap in (1:100)){
    print(bootStrap);
    cl1.IL13 <- apply(IL13.assayData[, sample(cells.TC1.IL13, 50, replace=TRUE)], 1, mean, na.rm=TRUE);
    cl1.other <- apply(IL13.assayData[, sample(cells.TC1.noIL13, 50, replace=TRUE)], 1, mean, na.rm=TRUE);
    gene.diffRanked <- cl1.IL13 - cl1.other;
    cat(sprintf("Bootstrap %d: \n", bootStrap), file=diffFileName, append=TRUE);
    names(gene.diffRanked) <- names(cl1.IL13);
    gene.diffRanked <- gene.diffRanked[order(-abs(gene.diffRanked))];
    comp.list.IL13.100 <- c(comp.list.IL13.100, names(head(gene.diffRanked, 100)));
    res <- head(gene.diffRanked, 100);
    write.fwf(as.data.frame(res), rownames=TRUE, colnames=FALSE, append=TRUE,
              file=diffFileName);
    cat(head(names(res), 10), fill=TRUE);
}

pdf(sprintf("consistent_set_IL13_differentiating_cluster1_%s.pdf", Sys.Date()), width=16, height=8);
par(mar=c(8,5,1,1));
barplot(rev(sort(table(comp.list.IL13.100)))[1:50], las=2,
        ylab="Number of bootstrap sub-samples", xlab="");
dev.off();

IL13.80pct <- rev(sort(table(comp.list.IL13.100)[table(comp.list.IL13.100) > 80]));
IL13.60pct <- rev(sort(table(comp.list.IL13.100)[table(comp.list.IL13.100) > 60]));
IL13.20pct <- rev(sort(table(comp.list.IL13.100)[table(comp.list.IL13.100) > 20]));
IL13.any <- rev(sort(table(comp.list.IL13.100)));
cat("Consistent set (> 60% of bootstraps): \n", file=diffFileName, append=TRUE);
cat(sort(names(IL13.60pct)), file=diffFileName, fill=TRUE, append=TRUE);
cat("Consistent set (> 80% of bootstraps): \n", file=diffFileName, append=TRUE);
cat(sort(names(IL13.80pct)), file=diffFileName, fill=TRUE, append=TRUE);

## Plot 60% consistent set
DefaultAssay(Cluster.TC1) <- "SCT";
library(viridis);
selected.cells <- c(which(GetAssayData(Cluster.TC1)["IL13",] > 0),
                    sample(which(GetAssayData(Cluster.TC1)["IL13",] == 0), 200));
selected.cells <- selected.cells[order(GetAssayData(Cluster.TC1)["IL13",selected.cells])];
pdf(sprintf("Heatmap_IL13_differentiating_genes_10k_bootstrap60_%s.pdf", format(Sys.Date())),
    height=8, width=8);
DoHeatmap(Cluster.TC1, features=names(IL13.60pct), cells=selected.cells,
          group.by="IL13.pos", size=4, raster=FALSE, slot="data") +
    scale_fill_viridis();
invisible(dev.off());

pdf(sprintf("Feature_IL13_differentiating_genes_10k_bootstrap80_%s.pdf",
           format(Sys.Date())),
    width=16, height=12);
    print(FeaturePlot(Cluster.TC1, col=c("lightgrey", "#e31836"),
                      features=names(IL13.80pct), order=TRUE, pt.size=1));
invisible(dev.off());

pdf(sprintf("Violin_IL13_differentiating_genes_10k_bootstrap80_%s.pdf",
           format(Sys.Date())),
    width=16, height=12);
print(VlnPlot(Cluster.TC1, features=names(IL13.80pct),
              pt.size=1, split.by="IL13.pos"));
invisible(dev.off());

Cluster.TC1 <- RunPCA(Cluster.TC1, features=setdiff(names(IL13.60pct), "IL13"));
Cluster.TC1 <- RunUMAP(object=Cluster.TC1, features=setdiff(names(IL13.60pct), "IL13"),
                       n.components=3, n.epochs=1000);

pdf(sprintf("PCA_UMAP_IL13_differentiating_genes_10k_bootstrap60_%s.pdf",
            format(Sys.Date())));
DimPlot(Cluster.TC1, reduction="pca", group.by="IL13.pos", pt.size=1, label=TRUE);
DimPlot(Cluster.TC1, reduction="umap", group.by="IL13.pos", pt.size=1, label=TRUE);
invisible(dev.off());


Cluster.TC1 <- RunPCA(Cluster.TC1, features=setdiff(names(IL13.80pct), "IL13"));
Cluster.TC1 <- RunUMAP(object=Cluster.TC1, features=setdiff(names(IL13.80pct), "IL13"),
                       n.components=3, n.epochs=1000);

DefaultAssay(Cluster.TC1) <- "integrated";
cells.IL13subset <- c(which(!!Cluster.TC1[["IL13.pos"]]), sample(which(!Cluster.TC1[["IL13.pos"]]), 200));
Cluster.TC1.any <- subset(Cluster.TC1, cells=cells.IL13subset);
Cluster.TC1.any <- ScaleData(Cluster.TC1.any, features=names(IL13.any));
Cluster.TC1.any <- RunPCA(Cluster.TC1.any, features=setdiff(names(IL13.any), "IL13"));
Cluster.TC1.any <- RunUMAP(object=Cluster.TC1.any, features=setdiff(names(IL13.any), "IL13"),
                           n.components=3, n.epochs=1000);
Cluster.TC1.any[["enriched.cluster"]] <- Embeddings(Reductions(Cluster.TC1.any, "umap"))[,"UMAP_2"] > 2;
Cluster.TC1.assay <- GetAssayData(Cluster.TC1.any);

diffFileName <-
    sprintf("differentiatingGenes_10k_Cluster1_IL13enriched_%s.txt", Sys.Date());
comp.list.enriched.100 <- NULL;
res <- NULL;
for(bootStrap in (1:100)){
    print(bootStrap);
    cl1.enriched <- apply(Cluster.TC1.assay[,sample(which(!!Cluster.TC1.any[["enriched.cluster"]]), 40)] , 1, mean, na.rm=TRUE);
    cl1.other <- apply(Cluster.TC1.assay[,sample(which(!Cluster.TC1.any[["enriched.cluster"]]), 40)] , 1, mean, na.rm=TRUE);
    gene.diffRanked <- cl1.enriched - cl1.other;
    cat(sprintf("Bootstrap %d: \n", bootStrap), file=diffFileName, append=TRUE);
    names(gene.diffRanked) <- names(cl1.enriched);
    gene.diffRanked <- gene.diffRanked[order(-abs(gene.diffRanked))];
    comp.list.enriched.100 <- c(comp.list.enriched.100, names(head(gene.diffRanked, 100)));
    res <- head(gene.diffRanked, 100);
    write.fwf(as.data.frame(res), rownames=TRUE, colnames=FALSE, append=TRUE,
              file=diffFileName);
    cat(head(names(res), 10), fill=TRUE);
}

enriched.80pct <- rev(sort(table(comp.list.enriched.100)[table(comp.list.enriched.100) > 80]));
enriched.60pct <- rev(sort(table(comp.list.enriched.100)[table(comp.list.enriched.100) > 60]));
enriched.20pct <- rev(sort(table(comp.list.enriched.100)[table(comp.list.enriched.100) > 20]));
enriched.any <- rev(sort(table(comp.list.enriched.100)));
cat("Consistent set (> 60% of bootstraps): \n", file=diffFileName, append=TRUE);
cat(sort(names(enriched.60pct)), file=diffFileName, fill=TRUE, append=TRUE);
cat("Consistent set (> 80% of bootstraps): \n", file=diffFileName, append=TRUE);
cat(sort(names(enriched.80pct)), file=diffFileName, fill=TRUE, append=TRUE);

head(rev(sort(table(comp.list.enriched.100))));

pdf(sprintf("consistent_set_IL13_enriched_cluster1_%s.pdf", Sys.Date()), width=16, height=8);
par(mar=c(8,5,1,1));
barplot(rev(sort(table(comp.list.enriched.100)))[1:50], las=2,
        ylab="Number of bootstrap sub-samples", xlab="");
dev.off();

DefaultAssay(Cluster.TC1) <- "SCT";
pdf(sprintf("Feature_IL13_enriched_differentiating_genes_10k_bootstrap80_%s.pdf",
           format(Sys.Date())),
    width=16, height=12);
    print(FeaturePlot(Cluster.TC1, col=c("lightgrey", "#e31836"),
                      features=c("IL13", names(enriched.80pct)), order=TRUE, pt.size=1));
invisible(dev.off());


pdf("out.pdf");
DimPlot(Cluster.TC1.any, reduction="pca", group.by="IL13.pos", pt.size=1, label=TRUE);
DimPlot(Cluster.TC1.any, reduction="umap", group.by="IL13.pos", pt.size=1, label=TRUE);
invisible(dev.off());

pdf(sprintf("PCA_UMAP_IL13_differentiating_genes_10k_bootstrap80_%s.pdf",
            format(Sys.Date())));
DimPlot(Cluster.TC1, reduction="pca", group.by="IL13.pos", pt.size=1, label=TRUE);
DimPlot(Cluster.TC1, reduction="umap", group.by="IL13.pos", pt.size=1, label=TRUE);
invisible(dev.off());


### Attempt clustering with SRGN from Cluster 9 (Monocytes)

Cluster.Mon <- subset(hc.no11, subset=(OL_clusters=="Monocytes" &
                                       UMAP_1 < -6 & UMAP_1 > -12 &
                                       UMAP_2 < 2.5 & UMAP_2 > -3));
DefaultAssay(Cluster.Mon) <- "SCT";

cells.Mon.SRGN <- which(GetAssayData(Cluster.Mon)["SRGN",] > 0);
cells.Mon.noSRGN <- which(GetAssayData(Cluster.Mon)["SRGN",] == 0);
Cluster.Mon[["SRGN.pos"]] <- GetAssayData(Cluster.Mon)["SRGN",] > 0;

SRGN.assayData <- GetAssayData(Cluster.Mon);
comp.list.SRGN.100 <- NULL;
diffFileName <-
    sprintf("differentiatingGenes_10k_Cluster1_SRGN_byCluster_%s.txt", Sys.Date());
cat("** Most-distinguishing SRGN / cluster 1 genes **\n", file=diffFileName);
for(bootStrap in (1:100)){
    print(bootStrap);
    cl1.SRGN <- apply(SRGN.assayData[, sample(cells.Mon.SRGN, 50, replace=TRUE)], 1, mean, na.rm=TRUE);
    cl1.other <- apply(SRGN.assayData[, sample(cells.Mon.noSRGN, 50, replace=TRUE)], 1, mean, na.rm=TRUE);
    gene.diffRanked <- cl1.SRGN - cl1.other;
    cat(sprintf("Bootstrap %d: \n", bootStrap), file=diffFileName, append=TRUE);
    names(gene.diffRanked) <- names(cl1.SRGN);
    gene.diffRanked <- gene.diffRanked[order(-abs(gene.diffRanked))];
    comp.list.SRGN.100 <- c(comp.list.SRGN.100, names(head(gene.diffRanked, 100)));
    res <- head(gene.diffRanked, 100);
    write.fwf(as.data.frame(res), rownames=TRUE, colnames=FALSE, append=TRUE,
              file=diffFileName);
    cat(head(names(res), 10), fill=TRUE);
}

pdf(sprintf("consistent_set_SRGN_differentiating_cluster1_%s.pdf", Sys.Date()), width=16, height=8);
par(mar=c(8,5,1,1));
barplot(rev(sort(table(comp.list.SRGN.100)))[1:50], las=2,
        ylab="Number of bootstrap sub-samples", xlab="");
dev.off();

SRGN.100pct <- rev(sort(table(comp.list.SRGN.100)[table(comp.list.SRGN.100) >= 100]));
SRGN.80pct <- rev(sort(table(comp.list.SRGN.100)[table(comp.list.SRGN.100) > 80]));
SRGN.60pct <- rev(sort(table(comp.list.SRGN.100)[table(comp.list.SRGN.100) > 60]));
SRGN.20pct <- rev(sort(table(comp.list.SRGN.100)[table(comp.list.SRGN.100) > 20]));
SRGN.1pct <- rev(sort(table(comp.list.SRGN.100)[table(comp.list.SRGN.100) >= 1]));
cat("Consistent set (> 60% of bootstraps): \n", file=diffFileName, append=TRUE);
cat(sort(names(SRGN.60pct)), file=diffFileName, fill=TRUE, append=TRUE);
cat("Consistent set (> 80% of bootstraps): \n", file=diffFileName, append=TRUE);
cat(sort(names(SRGN.80pct)), file=diffFileName, fill=TRUE, append=TRUE);

## Plot 60% consistent set
DefaultAssay(Cluster.Mon) <- "SCT";
library(viridis);
selected.cells <- c(which(GetAssayData(Cluster.Mon)["SRGN",] > 0),
                    which(GetAssayData(Cluster.Mon)["SRGN",] == 0));
selected.cells <- selected.cells[order(GetAssayData(Cluster.Mon)["SRGN",selected.cells])];
pdf(sprintf("Heatmap_SRGN_differentiating_genes_10k_bootstrap60_%s.pdf", format(Sys.Date())),
    height=8, width=8);
DoHeatmap(Cluster.Mon, features=names(SRGN.60pct), cells=selected.cells,
          group.by="SRGN.pos", size=4, raster=FALSE, slot="data") +
    scale_fill_viridis();
invisible(dev.off());

pdf(sprintf("Violin_SRGN_differentiating_genes_10k_bootstrap80_%s.pdf",
           format(Sys.Date())),
    width=16, height=12);
print(VlnPlot(Cluster.Mon, col=c("lightgrey", "#e31836"), group.by="SRGN.pos",
                  features=head(names(SRGN.80pct), 12), pt.size=1));
print(VlnPlot(Cluster.Mon, col=c("lightgrey", "#e31836"), group.by="SRGN.pos",
                      features=head(tail(names(SRGN.80pct), -12), 12), pt.size=1));
print(VlnPlot(Cluster.Mon, col=c("lightgrey", "#e31836"), group.by="SRGN.pos",
                  features=head(tail(names(SRGN.80pct), -24), 12), pt.size=1));
invisible(dev.off());

Cluster.Mon <- RunPCA(Cluster.Mon, features=setdiff(names(SRGN.60pct), "SRGN"));
Cluster.Mon <- RunUMAP(object=Cluster.Mon, features=setdiff(names(SRGN.60pct), "SRGN"),
                       n.components=3, n.epochs=1000);

pdf(sprintf("PCA_UMAP_SRGN_differentiating_genes_10k_bootstrap60_%s.pdf",
            format(Sys.Date())));
DimPlot(Cluster.Mon, reduction="pca", group.by="SRGN.pos", pt.size=1, label=TRUE);
DimPlot(Cluster.Mon, reduction="umap", group.by="SRGN.pos", pt.size=1, label=TRUE);
invisible(dev.off());

pdf(sprintf("Feature_SRGN_differentiating_genes_10k_bootstrap80_%s.pdf",
           format(Sys.Date())),
    width=16, height=12);
print(FeaturePlot(Cluster.Mon, col=c("lightgrey", "#e31836"), reduction="umap",
                  features=head(names(SRGN.80pct), 12), order=TRUE, pt.size=1));
print(FeaturePlot(Cluster.Mon, col=c("lightgrey", "#e31836"), reduction="umap",
                      features=head(tail(names(SRGN.80pct), -12), 12), order=TRUE, pt.size=1));
print(FeaturePlot(Cluster.Mon, col=c("lightgrey", "#e31836"), reduction="umap",
                  features=head(tail(names(SRGN.80pct), -24), 12), order=TRUE, pt.size=1));
invisible(dev.off());

