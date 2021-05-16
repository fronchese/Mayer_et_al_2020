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
library(SeuratWrappers);
library(monocle3);

## Loading base Seurat structure from source file
c57.S6KO.alv.filt.noStat6 <- readRDS(file="FR_2.3k_320epochs_noStat6_Alevin.rds");

c57.S6KO.alv.noMx <- c57.S6KO.alv.filt.noStat6;

## OL meeting 2021-Apr-01
## [use gene set OL13]

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

##BiocManager::install("batchelor");
##BiocManager::install("fastmap");

## Alternative monocle construction
## [https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/]
DefaultAssay(c57.S6KO.alv.noMx) <- "RNA";

grep("Stat6", rownames(c57.S6KO.alv.noMx), value = TRUE);

c57.S6KO.cds <- as.cell_data_set(c57.S6KO.alv.noMx) %>%
    preprocess_cds(use_genes = setdiff(rownames(c57.S6KO.alv.noMx), "Stat6-ENSMUSG00000002147.18"));

##c57.S6KO.cds <- align_cds(c57.S6KO.cds, alignment_group="orig.ident");
c57.S6KO.cds <- reduce_dimension(c57.S6KO.cds, verbose=TRUE);

c57.S6KO.cds <- cluster_cells(cds = c57.S6KO.cds, num_iter=1, k=50, resolution=0.04, cluster_method="leiden");
length(levels(clusters(c57.S6KO.cds)));

cairo_pdf(sprintf("monocle_UMAP_altConstruction_%s.pdf", Sys.Date()), onefile=TRUE,
          width=11, height=8);
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
invisible(dev.off());

cairo_pdf(sprintf("monocle_UMAP_FeaturePlots_%s.pdf", Sys.Date()), onefile=TRUE,
          width=11, height=8);
monocle::plot_cell_clusters(c57.S6KO.cds,
                   markers = listSub,
                   show_group_id = T, cell_size = 0.1);
invisible(dev.off());


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

marker_test_res <- top_markers(c57.S6KO.cds, group_cells_by="cluster",
                               cores=10);
write_csv(marker_test_res, sprintf("top_markers_monocle50_0.04_clusters_%s.csv", Sys.Date()));

c57.S6KO.alv.noMx.monocle3.filt <-
    subset(c57.S6KO.alv.noMx.monocle3,
           subset=(monocle_dump == FALSE));

DimPlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
        group.by="monocle50_0.04_cluster", pt.size=1.5)

c57.S6KO.cds.filt <- as.cell_data_set(c57.S6KO.alv.noMx.monocle3.filt) %>%
    preprocess_cds(use_genes = setdiff(rownames(c57.S6KO.alv.noMx.monocle3.filt), "Stat6-ENSMUSG00000002147.18"));

## -> 0.04
c57.S6KO.cds.filt <- reduce_dimension(c57.S6KO.cds.filt, verbose=TRUE);

c57.S6KO.alv.noMx.monocle3.filt[["monocleUMAP"]] <-
    CreateDimReducObject(
        embeddings = reducedDims(c57.S6KO.cds.filt)$UMAP,
        key = "monocleUMAP_",
        assay = "RNA"
    )

{
    cairo_pdf(sprintf("monocle_clusterTest_default_UMAP_reFilt_asSeuratPlot_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    c57.S6KO.cds.filt <-
        cluster_cells(cds = c57.S6KO.cds.filt);
    cat("Default =", length(levels(clusters(c57.S6KO.cds.filt))), "\n");
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
              ggtitle(sprintf("k=[default], res = [default], Monocle UMAP")));
    print(DimPlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                  group.by="orig.ident", pt.size=1.5) +
              ggtitle(sprintf("k=[default], res = [default], Monocle UMAP")));
    print(DimPlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                  group.by="monocle_cluster", pt.size=1.5, split.by="orig.ident") +
              ggtitle(sprintf("k=[default], res = [default], Monocle UMAP (split)")));
    cat("Finding top markers... ");
    marker_test_res <- top_markers(c57.S6KO.cds.filt, group_cells_by="cluster",
                                   cores=10);
    write_csv(marker_test_res, sprintf("top_markers_reFilt_monocle50_default_clusters_%s.csv", Sys.Date()));
    cat("done\n");
    invisible(dev.off());
}


{
    cairo_pdf(sprintf("monocle_clusterTest_0.002-0.4_UMAP_asSeuratPlot_%s.pdf", Sys.Date()), onefile=TRUE,
          width=11, height=8);
    for(tRes in round(seq(0.002, 0.04, length.out=10), 4)){
        c57.S6KO.cds.filt <-
            cluster_cells(cds = c57.S6KO.cds.filt, num_iter=1, random_seed=42, k=50,
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
            ggtitle(sprintf("k=50, res = %0.4f, Gene set OL14", tRes)) +
            scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            theme(axis.text.x=element_text(angle=45, hjust=1));
        print(res);
        cat("Finding top markers... ");
        marker_test_res <- top_markers(c57.S6KO.cds.filt, group_cells_by="cluster",
                                       cores=10);
        write_csv(marker_test_res, sprintf("top_markers_reFilt_monocle50_%0.4f_clusters_%s.csv", tRes, Sys.Date()));
        cat("done\n");
    }
    invisible(dev.off());
}

{ ## Including FindMarkers / DESeq2 testing
    cairo_pdf(sprintf("monocle_clusterTest_0.019_UMAP_asSeuratPlot_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    for(tRes in round(seq(0.019, 0.019, length.out=1), 4)){
        c57.S6KO.cds.filt <-
            cluster_cells(cds = c57.S6KO.cds.filt, num_iter=1, random_seed=42, k=50,
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
        DefaultAssay(c57.S6KO.alv.noMx.monocle3.filt) <- "SCT";
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
            ggtitle(sprintf("k=50, res = %0.4f, Gene set OL14", tRes)) +
            scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
            theme(axis.text.x=element_text(angle=45, hjust=1));
        print(res);
        cat("Finding top markers (Monocle)... ");
        marker_test_res <- top_markers(c57.S6KO.cds.filt, group_cells_by="cluster",
                                       cores=10);
        write_csv(marker_test_res, sprintf("top_markers_reFilt_monocle50_%0.4f_clusters_%s.csv", tRes, Sys.Date()));
        cat("done\n");
        cat("Finding top markers (Seurat)... ");
        res.all <- NULL;
        for(cl in levels(clusters(c57.S6KO.cds.filt))){
            cells.1 <- which(clusters(c57.S6KO.cds.filt) == cl);
            cells.2 <- which(clusters(c57.S6KO.cds.filt) != cl);
            c57.S6KO.alv.noMx.monocle3.filt[["clusterTest"]] <- ifelse(clusters(c57.S6KO.cds.filt) == cl, cl, "other");
            res.cl <- cbind(cluster=cl, FindMarkers(c57.S6KO.alv.noMx.monocle3.filt, group.by="clusterTest", ident.1=cl,
                                  test.use="DESeq2"));
            res.all <- rbind(res.all, res.cl);
        }
        write_csv(res.all %>% rownames_to_column("gene"),
                  sprintf("FindMarkers_all_reFilt_monocle50_%0.4f_clusters_%s.csv", tRes, Sys.Date()));
        top10 <- res.all %>% rownames_to_column("gene") %>% filter(p_val_adj < 0.25) %>%
            arrange(avg_log2FC) %>% group_by(cluster) %>% slice_head(n=10);
        bot10 <- res.all %>% rownames_to_column("gene") %>% filter(p_val_adj < 0.25) %>%
            arrange(avg_log2FC) %>% group_by(cluster) %>% slice_tail(n=10);
        write_csv(bind_rows(top10, bot10) %>% arrange(cluster, -avg_log2FC) %>% distinct(),
                  sprintf("FindMarkers_MaxMin10_reFilt_monocle50_%0.4f_clusters_%s.csv", tRes, Sys.Date()));
        cat("done\n");
    }
    invisible(dev.off());
}

## Generate features plots for the top 20 differentiating markers find by Seurat and Moncle for each clusters.

Seurat.genes <- read_csv("FindMarkers_MaxMin10_reFilt_monocle50_0.0190_clusters_2021-04-20.csv") %>% pull(gene);
monocle.genes <- read_csv("top_markers_reFilt_monocle50_0.0190_clusters_2021-04-20.csv") %>% pull(gene_id);

{
    cairo_pdf(sprintf("FeaturePlot_top_markers_monocle50_0.0190_UMAP_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    listSub <- monocle.genes;
    selectPoss <- seq(1, length(listSub), by=12);
    maxRes <- NULL;
    cat("To print:",length(listSub),"\n");
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+11)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                                      col=c("lightgrey", "#e31836"),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=4, nrow=3)));
    }
    invisible(dev.off());
}

{
    cairo_pdf(sprintf("FeaturePlot_FindMarkers_10MaxMin_monocle50_0.0190_UMAP_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    listSub <- Seurat.genes;
    geneIndex <- match(sub("-ENSM.*$", "", listSub),
                       sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    
    selectPoss <- seq(1, length(listSub), by=12);
    maxRes <- NULL;
    cat("To print:",length(listSub),"\n");
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+11)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                                      col=c("lightgrey", "#e31836"),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=4, nrow=3)));
    }
    invisible(dev.off());
}

{ ## DESeq2 testing, Stat6 vs WT
    for(tRes in round(seq(0.019, 0.019, length.out=1), 4)){
        c57.S6KO.cds.filt <-
            cluster_cells(cds = c57.S6KO.cds.filt, num_iter=1, random_seed=42, k=50,
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
        DefaultAssay(c57.S6KO.alv.noMx.monocle3.filt) <- "SCT";
        print(table(unlist(c57.S6KO.alv.noMx.monocle3.filt[["orig.ident"]]),
                    unlist(c57.S6KO.alv.noMx.monocle3.filt[["monocle_cluster"]])));
        cat("looking for differentially-expressed markers (WT vs Stat6)... ");
        res.all <- NULL;
        for(cl in levels(clusters(c57.S6KO.cds.filt))){
            cat(cl, "\n");
            cells.1 <- which((clusters(c57.S6KO.cds.filt) == cl) & (c57.S6KO.alv.noMx.monocle3.filt[["orig.ident"]] == "WT"));
            cells.2 <- which((clusters(c57.S6KO.cds.filt) == cl) & (c57.S6KO.alv.noMx.monocle3.filt[["orig.ident"]] == "Stat6KO"));
            data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, cells=c(cells.1, cells.2));
            res.cl <- cbind(cluster=cl, FindMarkers(data.sub, group.by="orig.ident", ident.1="Stat6KO",
                                                    test.use="DESeq2") %>% rownames_to_column("gene"));
            res.all <- rbind(res.all, res.cl);
        }
        write_csv(res.all,
                  sprintf("WholeGroupDE_all_STat6KO-WT_reFilt_monocle50_%0.4f_clusters_%s.csv", tRes, Sys.Date()));
        selected.genes <- res.all %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > log2(1.5)) %>%
            arrange(avg_log2FC) %>% group_by(cluster) %>% arrange(cluster, -avg_log2FC) %>% distinct();
        write_csv(selected.genes,
                  sprintf("WholeGroupDE_0.58_0.05_Stat6KO-WT_reFilt_monocle50_%0.4f_clusters_%s.csv", tRes, Sys.Date()));
        cat("done\n");
    }
    invisible(dev.off());
}

{
    Seurat.genes <- read_csv("WholeGroupDE_0.58_0.05_Stat6KO-WT_reFilt_monocle50_0.0190_clusters_2021-04-21.csv") %>% pull(gene);
    cairo_pdf(sprintf("FeaturePlot_DESeq2_0.58_0.05_monocle50_0.0190_UMAP_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    listSub <- Seurat.genes;
    geneIndex <- match(sub("-ENSM.*$", "", listSub),
                       sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    
    selectPoss <- seq(1, length(listSub), by=12);
    maxRes <- NULL;
    cat("To print:",length(listSub),"\n");
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+11)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                                      col=c("lightgrey", "#e31836"),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=4, nrow=3)));
    }
    invisible(dev.off());
}

{
    ## Create Gene Feature Plots
    geneList.OL15 <- c("Aldh1a2", "Ccl9", "Ccl17", "Ccrl2", "Fgl2", "Itgax",
                       "Mpeg1", "Slc7a11", "Crispld2", "Il9r", "Sned1",
                       "Zbtb10");
    geneIndex <- match(geneList.OL15, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx)[geneIndex]);
    cairo_pdf(sprintf("FeaturePlot_OL15_monocle_UMAP_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    selectPoss <- seq(1, length(listSub), by=6);
    maxRes <- NULL;
    cat("To print:",length(listSub),"\n");
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+5)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                                      col=c("lightgrey", "#e31836"),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=3, nrow=2)));
    }
    invisible(dev.off());
}


## Mini-bulk DESeq2 testing
{ 
    c57.S6KO.alv.noMx.monocle3.filt[["CD11b"]] <-
        ifelse(unlist(c57.S6KO.alv.noMx.monocle3.filt[["monocle_cluster"]]) %in% c("1", "2", "3"),
               "high", "low");
    res.all <- NULL;
    for(pop11b in c("high", "low")){
        cat(sprintf("CD11b[%s]: looking for differentially-expressed markers (WT vs Stat6)...\n", pop11b));
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, subset=(CD11b == pop11b));
        res.cl <- cbind(cluster=paste0("CD11b", pop11b), test="Stat6KO-WT",
                        FindMarkers(data.sub, group.by="orig.ident", ident.1="Stat6KO",
                                    test.use="DESeq2") %>% rownames_to_column("gene"));
        res.all <- rbind(res.all, res.cl);
    }
    for(tag in c("Stat6KO", "WT")){
        cat(sprintf("%s: looking for differentially-expressed markers (CD11bhigh vs CD11blow)...\n", tag));
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, subset=(orig.ident == tag));
        res.cl <- cbind(cluster=tag, test="CD11bhigh-CD11blow",
                        FindMarkers(data.sub, group.by="CD11b", ident.1="high",
                                    test.use="DESeq2") %>% rownames_to_column("gene"));
        res.all <- rbind(res.all, res.cl);
    }
    write_csv(res.all,
              sprintf("MiniBulkDE_all_reFilt_monocle50_%0.4f_clusters_%s.csv", tRes, Sys.Date()));
    selected.genes <- res.all %>% filter(p_val_adj < 0.05, abs(avg_log2FC) > log2(1.5)) %>%
        arrange(avg_log2FC) %>% group_by(cluster) %>% arrange(cluster, -avg_log2FC) %>% distinct();
    write_csv(selected.genes,
              sprintf("MiniBulkDE_0.58_0.05_reFilt_monocle50_%0.4f_clusters_%s.csv", tRes, Sys.Date()));
    cat("done\n");
}


## MA plots
marker.list.withGeneNames$logSum <-geneLogSums[marker.list.withGeneNames$gene];
marker.list.withGeneNames$shortName <- sub("-ENSM.*$", "", marker.list.withGeneNames$gene);

## Mini-bulk MA plots
{ 
    c57.S6KO.alv.noMx.monocle3.filt[["CD11b"]] <-
        ifelse(unlist(c57.S6KO.alv.noMx.monocle3.filt[["monocle_cluster"]]) %in% c("1", "2", "3"),
               "high", "low");
    cairo_pdf(sprintf("MAplot_MiniBulkDE_all_reFilt_monocle50_0.0190_clusters_%s.pdf", Sys.Date()), 
              width=11, height=8, onefile = TRUE);
    res.all <- read_csv("MiniBulkDE_all_reFilt_monocle50_0.0190_clusters_2021-04-22.csv");
    for(pop11b in c("high", "low")){
        cat(sprintf("CD11b[%s]: looking for differentially-expressed markers (WT vs Stat6)...\n", pop11b));
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, subset=(CD11b == pop11b));
        count.mat <- as.matrix(GetAssayData(data.sub, "counts", "RNA"));
        geneLogSums <- log2(rowSums(count.mat));
        res.sub <- res.all %>% filter(cluster==paste0("CD11b", pop11b), test=="Stat6KO-WT") %>%
            mutate(logSum = geneLogSums[gene], shortName = sub("-ENSM.*$", "", gene));
        res.sub %>%
            ggplot() +
            aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
            geom_point() +
            ggtitle(sprintf("CD11b%s: Stat6KO-WT", pop11b)) +
            geom_text_repel(data=res.sub %>% filter((p_val_adj < 0.05) & abs(avg_log2FC) > log2(1.5))) -> myPlot;
            print(myPlot);
    }
    for(tag in c("Stat6KO", "WT")){
        cat(sprintf("%s: looking for differentially-expressed markers (CD11bhigh vs CD11blow)...\n", tag));
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, subset=(orig.ident == tag));
        count.mat <- as.matrix(GetAssayData(data.sub, "counts", "RNA"));
        geneLogSums <- log2(rowSums(count.mat));
        res.sub <- res.all %>% filter(cluster==tag, test=="CD11bhigh-CD11blow") %>%
            mutate(logSum = geneLogSums[gene], shortName = sub("-ENSM.*$", "", gene));
        res.sub %>%
            ggplot() +
            aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
            geom_point() +
            ggtitle(sprintf("%s: CD11bhigh-CD11blow", tag)) +
            geom_text_repel(data=res.sub %>% filter((p_val_adj < 0.05) & abs(avg_log2FC) > log2(1.5))) -> myPlot;
        print(myPlot);
    }
    invisible(dev.off());
    cairo_pdf(sprintf("Volcano_MiniBulkDE_all_reFilt_monocle50_0.0190_clusters_%s.pdf", Sys.Date()), 
              width=11, height=8, onefile = TRUE);
    res.all <- read_csv("MiniBulkDE_all_reFilt_monocle50_0.0190_clusters_2021-04-22.csv");
    for(pop11b in c("high", "low")){
        cat(sprintf("CD11b[%s]: looking for differentially-expressed markers (WT vs Stat6)...\n", pop11b));
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, subset=(CD11b == pop11b));
        count.mat <- as.matrix(GetAssayData(data.sub, "counts", "RNA"));
        geneLogSums <- log2(rowSums(count.mat));
        res.sub <- res.all %>% filter(cluster==paste0("CD11b", pop11b), test=="Stat6KO-WT") %>%
            mutate(logSum = geneLogSums[gene], shortName = sub("-ENSM.*$", "", gene));
        res.sub %>%
            ggplot() +
            aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
            scale_y_log10() +
            geom_point() +
            ggtitle(sprintf("CD11b%s: Stat6KO-WT", pop11b)) +
            geom_text_repel(data=res.sub %>% filter((p_val_adj < 0.05) & abs(avg_log2FC) > log2(1.5))) -> myPlot;
        print(myPlot);
    }
    for(tag in c("Stat6KO", "WT")){
        cat(sprintf("%s: looking for differentially-expressed markers (CD11bhigh vs CD11blow)...\n", tag));
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, subset=(orig.ident == tag));
        count.mat <- as.matrix(GetAssayData(data.sub, "counts", "RNA"));
        geneLogSums <- log2(rowSums(count.mat));
        res.sub <- res.all %>% filter(cluster==tag, test=="CD11bhigh-CD11blow") %>%
            mutate(logSum = geneLogSums[gene], shortName = sub("-ENSM.*$", "", gene));
        res.sub %>%
            ggplot() +
            aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
            scale_y_log10() +
            geom_point() +
            ggtitle(sprintf("%s: CD11bhigh-CD11blow", tag)) +
            geom_text_repel(data=res.sub %>% filter((p_val_adj < 0.05) & abs(avg_log2FC) > log2(1.5))) -> myPlot;
        print(myPlot);
    }
    invisible(dev.off());
    cat("done\n");
}

## Cluster MA plots
{ 
    cairo_pdf(sprintf("MAplot_FindMarkers_all_reFilt_monocle50_0.0190_clusters_%s.pdf", Sys.Date()), 
              width=11, height=8, onefile = TRUE);
    res.all <- read_csv("FindMarkers_all_reFilt_monocle50_0.0190_clusters_2021-04-21.csv");
    for(cl in levels(unlist(c57.S6KO.alv.noMx.monocle3.filt[["monocle_cluster"]]))){
        cat(sprintf("Cluster %s: looking for differentiating markers...\n", cl));
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, subset=(monocle_cluster == cl));
        count.mat <- GetAssayData(data.sub, "counts", "RNA");
        geneLogSums <- log2(rowSums(count.mat));
        names(geneLogSums) <- sub("-ENSM.*$", "", names(geneLogSums));
        res.sub <- res.all %>% filter(cluster==cl) %>%
            mutate(shortName = sub("-ENSM.*$", "", gene), logSum = geneLogSums[shortName]);
        res.sub %>%
            ggplot() +
            aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
            geom_point() +
            ggtitle(sprintf("Cluster %s: differentiating markers", cl)) +
            geom_text_repel(data=res.sub %>% filter((p_val_adj < 0.05) & abs(avg_log2FC) > log2(1.5))) -> myPlot;
        print(myPlot);
    }
    invisible(dev.off());
    cairo_pdf(sprintf("Volcano_FindMarkers_all_reFilt_monocle50_0.0190_clusters_%s.pdf", Sys.Date()), 
              width=11, height=8, onefile = TRUE);
    res.all <- read_csv("FindMarkers_all_reFilt_monocle50_0.0190_clusters_2021-04-21.csv");
    for(cl in levels(unlist(c57.S6KO.alv.noMx.monocle3.filt[["monocle_cluster"]]))){
        cat(sprintf("Cluster %s: looking for differentiating markers...\n", cl));
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, subset=(monocle_cluster == cl));
        count.mat <- GetAssayData(data.sub, "counts", "RNA");
        geneLogSums <- log2(rowSums(count.mat));
        names(geneLogSums) <- sub("-ENSM.*$", "", names(geneLogSums));
        res.sub <- res.all %>% filter(cluster==cl) %>%
            mutate(shortName = sub("-ENSM.*$", "", gene), logSum = geneLogSums[shortName]);
        res.sub %>%
            ggplot() +
            aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
            scale_y_log10() +
            geom_point() +
            ggtitle(sprintf("Cluster %s: differentiating markers", cl)) +
            geom_text_repel(data=res.sub %>% filter((p_val_adj < 0.05) & abs(avg_log2FC) > log2(1.5))) -> myPlot;
        print(myPlot);
    }
    invisible(dev.off());
    cat("done\n");
}

## Save cluster information as a new file
###saveRDS(c57.S6KO.alv.noMx.monocle3.filt, "FR_Monocle_Filt_2.3k_320epochs_noStat6_Alevin.rds");


## Trajectory analysis
{
    c57.S6KO.cds.filt.test <-
        cluster_cells(cds = c57.S6KO.cds.filt);
    #c57.S6KO.cds.filt.test <-
    #    cluster_cells(cds = c57.S6KO.cds.filt, num_iter=1, random_seed=42, k=50,
    #                  resolution=tRes, cluster_method="leiden");
    c57.S6KO.cds.filt.test <- learn_graph(c57.S6KO.cds.filt.test, use_partition = TRUE);
    c57.S6KO.cds.filt.test <- order_cells(c57.S6KO.cds.filt.test, reduction_method = "UMAP");
    cairo_pdf(sprintf("monocle_reFilt_default_UMAP_test_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    res <- plot_cells(c57.S6KO.cds.filt.test,
               color_cells_by="pseudotime",
               label_groups_by_cluster=FALSE,
               group_label_size=4,
               graph_label_size=3,
               show_trajectory_graph=TRUE,
               cell_size=2) +
        labs(title = "Monocle3 UMAP, Monocle3 pseudotime trajectories")
    print(res);
    res <- plot_cells(c57.S6KO.cds.filt.test,
               color_cells_by="pseudotime",
               label_groups_by_cluster=FALSE,
               group_label_size=4,
               graph_label_size=3,
               show_trajectory_graph=TRUE,
               cell_size=2) +
        facet_wrap(~orig.ident, ncol=2) +
        labs(title = "Monocle3 UMAP, Monocle3 pseudotime trajectories (split by cell tag)")
    print(res);
    invisible(dev.off());
}

## Trajectory analysis
{
    #c57.S6KO.cds.filt.test <-
    #    cluster_cells(cds = c57.S6KO.cds.filt);
    c57.S6KO.cds.filt.test <-
        cluster_cells(cds = c57.S6KO.cds.filt, num_iter=1, random_seed=42, k=50,
                      resolution=0.019, cluster_method="leiden");
    c57.S6KO.cds.filt.test <- learn_graph(c57.S6KO.cds.filt.test, use_partition = TRUE);
    c57.S6KO.cds.filt.test <- order_cells(c57.S6KO.cds.filt.test, reduction_method = "UMAP");
    cairo_pdf(sprintf("monocle_reFilt_0.019_UMAP_test_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    res <- plot_cells(c57.S6KO.cds.filt.test,
                      color_cells_by="pseudotime",
                      label_groups_by_cluster=FALSE,
                      group_label_size=4,
                      graph_label_size=3,
                      show_trajectory_graph=TRUE,
                      cell_size=2) +
        labs(title = "Monocle3 UMAP, Monocle3 pseudotime trajectories")
    print(res);
    res <- plot_cells(c57.S6KO.cds.filt.test,
                      color_cells_by="pseudotime",
                      label_groups_by_cluster=FALSE,
                      group_label_size=4,
                      graph_label_size=3,
                      show_trajectory_graph=TRUE,
                      cell_size=2) +
        facet_wrap(~orig.ident, ncol=2) +
        labs(title = "Monocle3 UMAP, Monocle3 pseudotime trajectories (split by cell tag)")
    print(res);
    invisible(dev.off());
}


## Module test (0.019)
module_res <- graph_test(c57.S6KO.cds.filt.test, neighbor_graph="principal_graph", cores=10)
pr_deg_ids <- row.names(subset(module_res, morans_test_statistic > 70))

{
    cairo_pdf(sprintf("monocle_graphTest_reFilt_0.019_UMAP_test_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    listSub <- pr_deg_ids;
    selectPoss <- seq(1, length(listSub), by=4);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+3)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                                      col=viridis(100),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=2, nrow=2)));
    }
    invisible(dev.off());
}
    
## Module test (default)
module_res <- graph_test(c57.S6KO.cds.filt.test, neighbor_graph="principal_graph", cores=10);
write.csv(module_res, sprintf("graph_test_default_reFilt_default_UMAP_test_%s.csv", Sys.Date()));
module_res %<>% arrange(-morans_test_statistic);
pr_deg_ids <- rownames(head(module_res, 100));

{
    cairo_pdf(sprintf("monocle_graphTest_reFilt_default_UMAP_test_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    listSub <- pr_deg_ids;
    selectPoss <- seq(1, length(listSub), by=6);
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+5)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                                      col=viridis(100),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=3, nrow=2)));
    }
    invisible(dev.off());
}


## [doesn't work; genes not found]
#plot_cells(c57.S6KO.cds.filt.test, genes=c(10746),
#           show_trajectory_graph=FALSE,
#           label_cell_groups=FALSE,
#           label_leaves=FALSE)

lineage_cds <- c57.S6KO.cds.filt.test[rownames(c57.S6KO.cds.filt.test) %in% rownames(module_res)[module_res$q_value < 0.05],];
gene_module_df <- find_gene_modules(lineage_cds, verbose=TRUE, cores=10);
plot(gene_module_df$dim_1, gene_module_df$dim_2, col=gene_module_df$module)

write.csv(gene_module_df, sprintf("gene_modules_default_reFilt_default_UMAP_test_%s.csv", Sys.Date()));

## Email 2021-May-05

## Loading clustered data from RDS
c57.S6KO.alv.noMx.monocle3.filt <- readRDS(file="FR_Monocle_Filt_2.3k_320epochs_noStat6_Alevin.rds");

## Combine and rename clusters
c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]] <-
c("1" = "CD11bhigh DC2-1",
  "6" = "CD11bhigh DC2-1",
  "2" = "CD11bhigh DC2-2",
  "3" = "CD11bhigh DC2-2",
  "4" = "CD11blow DC2",
  "5" = "CD11blow DC2")[as.character(unlist(c57.S6KO.alv.noMx.monocle3.filt[["monocle_cluster"]]))];

## UMAP WT and STAT6-KO cells combined and split. Please, could you choose a colour palette that is colour blind friendly.
ol.colours <- brewer.pal(3, "Dark2");
names(ol.colours) <- sort(unique(unlist(c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]])));
{
    DefaultAssay(c57.S6KO.alv.noMx.monocle3.filt) <- "SCT";
    cairo_pdf(sprintf("monocle_UMAP_asSeuratPlot_OLlabels_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    print(DimPlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                  group.by="monocle_cluster", pt.size=1.5, label=TRUE));
    print(DimPlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                  group.by="OL_mon_cluster", pt.size=1.5, cols = ol.colours));
    print(DimPlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                  group.by="OL_mon_cluster", pt.size=1.5, cols = ol.colours,
                  split.by="orig.ident"));
    invisible(dev.off());
}

##Number of cells in each cluster in WT and STAT6-KO.
table(unlist(c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]]),
      unlist(c57.S6KO.alv.noMx.monocle3.filt[["orig.ident"]]))

#                Stat6KO  WT
#CD11bhigh DC2-1     435 205
#CD11bhigh DC2-2     575 286
#CD11blow DC2        110 492

DefaultAssay(c57.S6KO.alv.noMx) <- "SCT";

{
    ## Create Gene Feature Plots
    geneList.OL16 <- c("Aldh1a2", "Ccl9", "Ccl17", "Fgl2", "Itgax", "Mpeg1", "Slc7a11", "Cdhr1", 
                       "Clmn", "Crispld2", "Il9r", "Khdc1a", "Sned1", "Zbtb10");
    geneIndex <- match(geneList.OL16, sub("-ENSM.*$", "", rownames(c57.S6KO.alv.noMx.monocle3.filt)));
    listSub <- na.omit(rownames(c57.S6KO.alv.noMx.monocle3.filt)[geneIndex]);
    cairo_pdf(sprintf("DotPlot_FeaturePlot_OL16_monocle_UMAP_OLlabels_%s.pdf", Sys.Date()), onefile=TRUE,
              width=11, height=8);
    res <- DotPlot(c57.S6KO.alv.noMx.monocle3.filt, features=listSub,
                   group.by="OL_mon_cluster",
                   dot.min=0.0001, scale.by="size", scale=TRUE,
                   col.min=0, col.max=1) +
        scale_colour_viridis() +
        ylab("Cluster") +
        ggtitle("Gene set OL16 (res: 0.019)") +
        scale_x_discrete(labels=function(x){sub("-ENSM.*", "", x)}) +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    print(res);
    selectPoss <- seq(1, length(listSub), by=6);
    maxRes <- NULL;
    cat("To print:",length(listSub),"\n");
    for(listStart in selectPoss){
        cat(listStart, "\n");
        geneSubset <- na.omit(listSub[listStart:(listStart+5)]);
        res <- list();
        for(gene in geneSubset){
            res <- c(res, FeaturePlot(c57.S6KO.alv.noMx.monocle3.filt, reduction="monocleUMAP",
                                      col=c("lightgrey", "#e31836"),
                                      features=gene, order=TRUE, pt.size=1, shape.by="shape",
                                      combine=FALSE));
        }
        res <- lapply(res, function(x){
            x$labels$title <- sub("-ENSM.*$", "", x$labels$title);
            x + guides(shape=FALSE, colour=guide_legend(override.aes=list(shape=16)));
        });
        suppressWarnings(print(wrap_plots(res, ncol=3, nrow=2)));
    }
    invisible(dev.off());
}

## Monocle 25 top markers genes for the 3 clusters.

c57.S6KO.cds.OLnames <- as.cell_data_set(c57.S6KO.alv.noMx.monocle3.filt);
monocle_marker_test_res <- top_markers(c57.S6KO.cds.OLnames, group_cells_by="OL_mon_cluster",
                               cores=10);
write_csv(monocle_marker_test_res, sprintf("top_markers_monocle50_0.0190_OLlabels_clusters_%s.csv", Sys.Date()));

## Seurat Differentiating Markers for the 3 clusters
cat("Finding top markers (Seurat)... ");
res.all <- NULL;
for(cl in unique(unlist(c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]]))){
    ##cells.1 <- which(unlist(c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]]) == cl);
    ##cells.2 <- which(unlist(c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]]) != cl);
    res.cl <- cbind(cluster=cl, FindMarkers(c57.S6KO.alv.noMx.monocle3.filt, group.by="OL_mon_cluster", ident.1=cl,
                                            test.use="DESeq2"));
    res.all <- rbind(res.all, res.cl);
}
write_csv(res.all %>% rownames_to_column("gene"),
          sprintf("FindMarkers_all_reFilt_monocle50_0.019_OLlabels_clusters_%s.csv", Sys.Date()));
top10 <- res.all %>% rownames_to_column("gene") %>% dplyr::filter(p_val_adj < 0.25) %>%
    arrange(avg_log2FC) %>% group_by(cluster) %>% slice_head(n=10);
bot10 <- res.all %>% rownames_to_column("gene") %>% dplyr::filter(p_val_adj < 0.25) %>%
    arrange(avg_log2FC) %>% group_by(cluster) %>% slice_tail(n=10);
write_csv(bind_rows(top10, bot10) %>% arrange(cluster, -avg_log2FC) %>% distinct(),
          sprintf("FindMarkers_MaxMin10_reFilt_monocle50_0.0190_OLlabels_clusters_%s.csv", Sys.Date()));
cat("done\n");

## DESeq2 comparison for STAT6-KO vs WT for the 3 clusters with DEGs list, MA and volcano plots.
{
    cat("looking for differentially-expressed markers (WT vs Stat6)... ");
    res.all <- NULL;
    for(cl in unique(unlist(c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]]))){
        cat(cl, "\n");
        cells.1 <- which((unlist(c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]]) == cl) & 
                             (c57.S6KO.alv.noMx.monocle3.filt[["orig.ident"]] == "WT"));
        cells.2 <- which((unlist(c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]]) == cl) & 
                             (c57.S6KO.alv.noMx.monocle3.filt[["orig.ident"]] == "Stat6KO"));
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, cells=c(cells.1, cells.2));
        res.cl <- cbind(cluster=cl, FindMarkers(data.sub, group.by="orig.ident", ident.1="Stat6KO",
                                                test.use="DESeq2") %>% rownames_to_column("gene"));
        res.all <- rbind(res.all, res.cl);
    }
    write_csv(res.all,
              sprintf("WholeGroupDE_all_STat6KO-WT_reFilt_monocle50_0.0190_OLlabels_clusters_%s.csv", Sys.Date()));
    selected.genes <- res.all %>% dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > log2(1.5)) %>%
        arrange(avg_log2FC) %>% group_by(cluster) %>% arrange(cluster, -avg_log2FC) %>% distinct();
    write_csv(selected.genes,
              sprintf("WholeGroupDE_0.58_0.05_Stat6KO-WT_reFilt_monocle50_0.0190_OLlabels_clusters_%s.csv", Sys.Date()));
    cat("done\n");
}
## MA and volcano plots
{
    res.all <- read_csv(sprintf("WholeGroupDE_all_STat6KO-WT_reFilt_monocle50_0.0190_OLlabels_clusters_2021-05-10.csv"));
    cairo_pdf(sprintf("Volcano_MA_all_reFilt_monocle50_0.0190_OLlabels_clusters_%s.pdf", Sys.Date()), 
              width=11, height=8, onefile = TRUE);
    for(cl in unique(unlist(c57.S6KO.alv.noMx.monocle3.filt[["OL_mon_cluster"]]))){
        cat("Cluster =", cl,"\n");
        data.sub <- subset(c57.S6KO.alv.noMx.monocle3.filt, subset=(OL_mon_cluster == cl));
        count.mat <- as.matrix(GetAssayData(data.sub, "counts", "RNA"));
        geneLogSums <- log2(rowSums(count.mat));
        res.sub <- res.all %>% dplyr::filter(cluster==cl) %>%
            mutate(logSum = geneLogSums[gene], shortName = sub("-ENSM.*$", "", gene));
        res.sub %>%
            ggplot() +
            aes(x=logSum, y=avg_log2FC, col=p_val_adj, label=shortName) +
            geom_point() +
            ggtitle(sprintf("%s: Stat6KO-WT", cl)) +
            geom_text_repel(data=res.sub %>% dplyr::filter((p_val_adj < 0.05) & abs(avg_log2FC) > log2(1.5))) -> myPlot;
        print(myPlot);
        res.sub %>%
            ggplot() +
            aes(x=avg_log2FC, y=1/p_val_adj, col=logSum, label=shortName) +
            scale_y_log10() +
            geom_point() +
            ggtitle(sprintf("%s: Stat6KO-WT", cl)) +
            geom_text_repel(data=res.sub %>% dplyr::filter((p_val_adj < 0.05) & abs(avg_log2FC) > log2(1.5))) -> myPlot;
        print(myPlot);
    }
    invisible(dev.off());
}

saveRDS(c57.S6KO.alv.noMx.monocle3.filt, "FR_Monocle_reFilt_2.3k_OLlabels.rds");
