#!/usr/bin/env Rscript

library(Seurat);
library(tidyverse);
library(sctransform);
library(patchwork);

dateStr <- format(Sys.Date(), "%Y-%b-%d");

## TOKEEP: information here regarding SCT normalising data, up to SCTransform
## [used for Figure 7b]

## In the blood and skin data, clusters were called at a resolution of
## 1.8 using the first six principal components, generated using
## highly variable genes.

bloodSkin.tbl <-
    read_csv("data/private_download/CSV/BloodSkin_rawExpr.csv");

bloodSkin.meta <-
    read_csv("data/private_download/CSV/BloodSkin_meta.csv");

metadata.tbl <- read_delim("data/E-MTAB-8498.sdrf.txt", delim="\t") %>%
    select(-`Comment[FASTQ_URI]`, -`Scan Name`,
           -`Comment[SUBMITTED_FILE_NAME]`) %>%
    distinct %>%
    mutate(cellID = sub("^", "cell", sub("_#", "_", `Source Name`)),
           cellType = `Characteristics[inferred cell type]`,
           cellPart = `Characteristics[organism part]`,
           stimulus = `Characteristics[stimulus]`) %>%
    right_join(bloodSkin.meta, by=c("cellID" = "X1"))

write_csv(metadata.tbl, "data/JEM_metadata_combined.csv");

bloodSkin.mat <- as.matrix(bloodSkin.tbl[,-1]);
rownames(bloodSkin.mat) <- gsub("_",".",bloodSkin.tbl$X1);
colnames(bloodSkin.mat) <- gsub("_",".",colnames(bloodSkin.mat));

## sctransform vignette from
## https://satijalab.org/seurat/v3.1/sctransform_vignette.html

bs <- CreateSeuratObject(counts = bloodSkin.mat);

## store mitochondrial percentage in object meta data
bs <- PercentageFeatureSet(bs, pattern = "\\.MT-", col.name = "percent.mt");

## run sctransform
bs <- SCTransform(bs, vars.to.regress = "percent.mt", verbose = FALSE);

## ^^^^ STOP HERE ^^^^

## These are now standard steps in the Seurat workflow for
## visualization and clustering
bs <- RunPCA(bs, verbose = FALSE);
bs <- RunTSNE(bs, dims = 1:30, verbose = FALSE)
bs <- RunUMAP(bs, dims = 1:30, verbose = FALSE)

## Note: resolution = 2.4 is needed to distinguish top left cells
## Using consistent workflow as per saline/HDM:
## - resolution 2.4
## - 6 dimensions
## - all genes
bs <- FindNeighbors(bs, dims = 1:6, verbose = FALSE);
bs <- FindClusters(bs, verbose = FALSE, resolution=2.4);

bs[["cellOrigin"]] <- metadata.tbl$cellPart[match(colnames(bs),
                                                  metadata.tbl$cellID)];
bs[["cellType"]] <- metadata.tbl$cellType[match(colnames(bs),
                                                metadata.tbl$cellID)];

(TSNEPlot(bs) /
 TSNEPlot(bs, group.by = "cellOrigin") /
 TSNEPlot(bs, group.by = "cellType")) |
(UMAPPlot(bs) /
 UMAPPlot(bs, group.by = "cellOrigin") /
 UMAPPlot(bs, group.by = "cellType"))

ggsave(sprintf("Cluster_plot_BloodSkin_%s.pdf", dateStr),
       width=13, height=11);

bs.markers <- FindAllMarkers(bs, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25);

bs.top10 <- bs.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC);
DoHeatmap(bs, features=bs.top10$gene);

ggsave(sprintf("Heatmap_BloodSkin_%s.pdf", dateStr),
       width=13, height=11);

## In the saline and HDM data, clusters were called at a resolution of
## 2.4 using the first six principal components, generated using all
## genes.

## [ I'm using res=0.8 ]

salineHDM.tbl <-
    read_csv("data/private_download/CSV/AllSkin_rawExpr.csv");

salineHDM.mat <- as.matrix(salineHDM.tbl[,-1]);
rownames(salineHDM.mat) <- gsub("_",".",salineHDM.tbl$X1);
colnames(salineHDM.mat) <- gsub("_",".",colnames(salineHDM.mat));

## follow through the same procedure as BloodSkin
sh <- CreateSeuratObject(counts = salineHDM.mat);
sh <- PercentageFeatureSet(sh, pattern = "\\.MT-", col.name = "percent.mt");
sh <- SCTransform(sh, vars.to.regress = "percent.mt", verbose = FALSE);
## dimension reduction
sh <- RunPCA(sh, verbose = FALSE);
sh <- RunTSNE(sh, dims = 1:30, verbose = FALSE)
sh <- RunUMAP(sh, dims = 1:30, verbose = FALSE)

## clustering
sh <- FindNeighbors(sh, dims = 1:6, verbose = FALSE)
sh <- FindClusters(sh, verbose = FALSE, resolution=0.8)
## additional labels
sh[["cellOrigin"]] <- metadata.tbl$cellPart[match(colnames(sh),
                                                  metadata.tbl$cellID)];
sh[["stimulus"]] <- metadata.tbl$stimulus[match(colnames(sh),
                                                  metadata.tbl$cellID)];
sh[["cellType"]] <- metadata.tbl$cellType[match(colnames(sh),
                                                metadata.tbl$cellID)];

## plotting
TSNEPlot(sh) +
    UMAPPlot(sh) +
    TSNEPlot(sh, group.by = "stimulus") +
    UMAPPlot(sh, group.by = "stimulus") +
    TSNEPlot(sh, group.by = "cellType") +
    UMAPPlot(sh, group.by = "cellType") +
    plot_layout(ncol=2);

ggsave(sprintf("Cluster_plot_AllBlood_%s.pdf", dateStr),
       width=13, height=11);

## Find distinctive markers for identified clusters
sh.markers <- FindAllMarkers(sh, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25);

## draw plot
sh.top10 <- sh.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC);
DoHeatmap(sh, features=sh.top10$gene)

ggsave(sprintf("Heatmap_AllBlood_%s.pdf", dateStr),
       width=13, height=11);



#### GSEA for skin vs blood (BloodSkin dataset) ####

## TOKEEP: information here regarding generating the normalised SCT count matrix
## [used for Figure 7b]

## Create Gene comparison scores
sct.gene.mat <- sct.count.mat <- bs[["SCT"]]@data;

## ^^^^ STOP HERE ^^^^
rownames(sct.gene.mat) <- sub("^ENSG.*\\.", "", rownames(sct.count.mat));
sct.gene.mat <- sct.gene.mat[rownames(sct.gene.mat) != "NA",];
sct.gene.mat <- sct.gene.mat[!grepl("^[0-9]+$", rownames(sct.gene.mat)),];
sct.gene.mat <-
    sct.gene.mat[match(unique(rownames(sct.gene.mat)),
                       rownames(sct.gene.mat)),];

write.csv(round(sct.gene.mat,4),
          gzfile("data/JEM_normalised_SCT_countMat.csv.gz"));

geneList.tbl <- read_csv("Aggregated_OL.csv");

## Carry out GSEA on populations vs all others

pdf(sprintf("GSEA_absDiff_OLsets_vs_JEM-clusters-all-vs-all_%s.pdf", dateStr),
    width=12, height=6);
for(cell.pop in sort(unique(Idents(bs)))){
    print(cell.pop);
    pop1.cols <- which(Idents(bs) == cell.pop);
    pop2.cols <- which(Idents(bs) != cell.pop);
    pop.diff <- sort(abs(rowMeans(as.matrix(sct.count.mat[,pop1.cols])) -
                         rowMeans(as.matrix(sct.count.mat[,pop2.cols]))),
                     decreasing=TRUE);
    names(pop.diff) <- sub("^.*\\.", "", names(pop.diff));
    for(listName in unique(geneList.tbl$ListName)){
        (geneList.tbl %>%
         filter(ListName == listName))$HumanOrth -> listGeneNames;
        gsea(names(pop.diff), listGeneNames,
             listName=listName,
             rankStatName=sprintf("JEM %s vs others", cell.pop));
    }
}
invisible(dev.off());


pdf(sprintf("GSEA_absDiff_OLSets_vs_cellType_%s.pdf", dateStr),
    width=12, height=6);
ncp <- bs[["cellType"]]$cellType;
names(ncp) <- rownames(bs[["cellType"]]);
for(cell.pop in sort(unique(ncp))){
    print(cell.pop);
    pop1.cols <- names(which(ncp == cell.pop));
    pop2.cols <- names(which(ncp != cell.pop));
    pop.diff <- sort(abs(rowMeans(as.matrix(sct.count.mat[,pop1.cols])) -
                         rowMeans(as.matrix(sct.count.mat[,pop2.cols]))),
                     decreasing=TRUE);
    names(pop.diff) <- sub("^.*\\.", "", names(pop.diff));
    for(listName in unique(geneList.tbl$ListName)){
        (geneList.tbl %>%
         filter(ListName == listName))$HumanOrth -> listGeneNames;
        gsea(names(pop.diff), listGeneNames,
             listName=listName,
             rankStatName=sprintf("JEM %s vs others", cell.pop));
    }
}
invisible(dev.off());


pdf(sprintf("GSEA_absDiff_OLsets_vs_skinBlood_%s.pdf", dateStr),
    width=12, height=6);
ncp <- bs[["cellOrigin"]]$cellType;
names(ncp) <- rownames(bs[["cellOrigin"]]);
pop1.cols <- names(which(ncp == cell.pop));
pop2.cols <- names(which(ncp != cell.pop));
pop.diff <- sort(abs(rowMeans(as.matrix(sct.count.mat[,pop1.cols])) -
                     rowMeans(as.matrix(sct.count.mat[,pop2.cols]))),
                 decreasing=TRUE);
names(pop.diff) <- sub("^.*\\.", "", names(pop.diff));
for(listName in unique(geneList.tbl$ListName)){
    print(listName);
    (geneList.tbl %>%
     filter(ListName == listName))$HumanOrth -> listGeneNames;
    gsea(names(pop.diff), listGeneNames,
         listName=listName,
         rankStatName="JEM skin vs blood");
}
invisible(dev.off());

## re-attempt with Olivier's classifications 2019-Feb-19 (for Blood/Skin)

pdf(sprintf("GSEA_absDiff_OLSets_vs_OLClassification_%s.pdf", dateStr),
    width=6, height=6);
sct.count.mat <- bs[["SCT"]]@data;
colSets <- list(
    ## "JEM activated DC blood vs skin" =
    ##     list("pop1.cols" = which((bs[["cellType"]]$cellType ==
    ##                               "activated dendritic cell") &
    ##                              (bs[["cellOrigin"]]$cellOrigin == "blood")),
    ##          "pop2.cols" = which((bs[["cellType"]]$cellType ==
    ##                               "activated dendritic cell") &
    ##                              (bs[["cellOrigin"]]$cellOrigin == "skin"))),
    ## "JEM cDC/monocyte blood vs activated DC skin" =
    ##     list("pop1.cols" = which((bs[["cellType"]]$cellType == "monocyte") &
    ##                              (bs[["cellOrigin"]]$cellOrigin == "blood")),
    ##          "pop2.cols" = which((bs[["cellType"]]$cellType ==
    ##                               "activated dendritic cell") &
    ##                              (bs[["cellOrigin"]]$cellOrigin == "skin"))),
    "JEM cDC/monocyte; skin vs blood" =
        list("pop1.cols" = which((bs[["cellType"]]$cellType == "monocyte") &
                                 (bs[["cellOrigin"]]$cellOrigin == "skin")),
             "pop2.cols" = which((bs[["cellType"]]$cellType == "monocyte") &
                                 (bs[["cellOrigin"]]$cellOrigin == "blood")))
    ## "JEM activated DC vs pDC" =
    ##     list("pop1.cols" = which((bs[["cellType"]]$cellType ==
    ##                               "activated dendritic cell")),
    ##          "pop2.cols" = which((bs[["cellType"]]$cellType ==
    ##                               "plasmacytoid dendritic cell"))),
    ## "JEM cDC/monocyte vs pDC" =
    ##     list("pop1.cols" = which((bs[["cellType"]]$cellType ==
    ##                               "monocyte")),
    ##          "pop2.cols" = which((bs[["cellType"]]$cellType ==
    ##                               "plasmacytoid dendritic cell")))
    );
for(li in names(colSets)){
    cat(sprintf("** %s **\n", li));
    pop1.cols <- colSets[[li]][["pop1.cols"]];
    pop2.cols <- colSets[[li]][["pop2.cols"]];
    pop.diff <- rowMeans(as.matrix(sct.count.mat[,pop1.cols])) -
        rowMeans(as.matrix(sct.count.mat[,pop2.cols]));
    pop.diff <- pop.diff[(rowMeans(as.matrix(sct.count.mat[,pop1.cols])) > 0) |
                         (rowMeans(as.matrix(sct.count.mat[,pop2.cols]) > 0))];
    pop.diff <- pop.diff[order(-abs(pop.diff))];
    ## filter out duplicates
    names(pop.diff) <- sub("\\.[0-9]+$","",
                           sub("^.*?\\.", "", names(pop.diff)));
    pop.diff <- pop.diff[names(pop.diff) != "NA"];
    ## retain the most significant of each named gene
    pop.diff <- pop.diff[match(unique(names(pop.diff)), names(pop.diff))];
    ## reorder, just in case unique fiddled with the order
    pop.diff <- pop.diff[order(-(pop.diff))];
    write_csv(tibble(
        "Dataset" = li,
        "Gene" = names(pop.diff),
        "Difference" = round(pop.diff, 3)),
        path=sprintf("geneRank_%s_%s.csv", gsub("[/ ;]+",".",li), dateStr));
    for(listName in unique(geneList.tbl$ListName)){
        print(listName);
        (geneList.tbl %>%
         filter(ListName == listName))$HumanOrth -> listGeneNames;
        gsea(names(pop.diff), listGeneNames,
             listName=listName,
             rankStatName=li, zeroPoint=sum(pop.diff > 0));
    }
}
invisible(dev.off());

geneLookup <- rownames(bs);
names(geneLookup) <- sub("\\.[0-9]+$","",
                         sub("^.*?\\.", "", geneLookup));
geneLookup <- geneLookup[names(geneLookup) != "NA"];

pwNames <- sub(".genes.txt$", "",
               list.files("results/olivier/my_analysis.GseaPreranked.1582500542678",
                          pattern=".genes.txt$"));
geneList <- sapply(pwNames, function(x){
    readLines(sprintf("results/olivier/my_analysis.GseaPreranked.1582500542678/%s.genes.txt",x));
});


## Set up GSEA function

gsea <- function(rankedGenes, targetList, listName="<listName>",
                 rankStatName="<rankStatName>", zeroPoint=length(rankedGenes)){
    targetList <- unique(targetList);
    targetList <- targetList[targetList %in% rankedGenes];
    ## Determine adjustment for positively-associated genes
    posRankedGenes <- rankedGenes[seq_along(rankedGenes) < zeroPoint];
    posListGenes <- targetList[targetList %in% posRankedGenes];
    increasePerPosGene <- 1 - (seq_along(posRankedGenes) / length(posRankedGenes));
    increasePosSum <- sum(increasePerPosGene[posRankedGenes %in% posListGenes]);
    increasePerPosGene <- increasePerPosGene * 1 / (increasePosSum);
    ## Determine adjustment for negatively-associated genes
    negRankedGenes <- rankedGenes[seq_along(rankedGenes) >= zeroPoint];
    negListGenes <- targetList[targetList %in% negRankedGenes];
    increasePerNegGene <- rev(1 - (seq_along(negRankedGenes) / length(negRankedGenes)));
    increaseNegSum <- sum(increasePerNegGene[negRankedGenes %in% negListGenes]);
    increasePerNegGene <- increasePerNegGene * 1 / (increaseNegSum);
    increasePerListGene <- c(increasePerPosGene, increasePerNegGene);
    increasePerNonList <-
        c(rep(-(1 / (length(posRankedGenes) - length(posListGenes))),
              length(posRankedGenes)),
          rep(-(1 / (length(negRankedGenes) - length(negListGenes))),
              length(negRankedGenes)));
    tibble(Gene=rankedGenes) %>%
        mutate(inList = Gene %in% targetList,
               inListVal = ifelse(inList, increasePerListGene, 0),
               nonListVal = ifelse(inList, 0, increasePerNonList)) %>%
        mutate(cModList = cumsum(inListVal),
               cModNon = cumsum(nonListVal),
               cMod = cModList + cModNon) -> geneOrder;
    area.traps <-
        (head(geneOrder$cMod, -1) + tail(geneOrder$cMod, -1)) / 2;
    integrated.area.pos <- sum(area.traps[head(seq_along(rankedGenes) < zeroPoint, -1)]);
    integrated.area.neg <- sum(area.traps[head(seq_along(rankedGenes) >= zeroPoint, -1)]);
    print(c(length(posRankedGenes), length(negRankedGenes), integrated.area.pos, integrated.area.neg));
    max.area.pos <- 0.5 * (length(posRankedGenes));
    max.area.neg <- 0.5 * (length(negRankedGenes));
    par(mar=c(8,4.5,4,0.5));
    plot(NA, xlim=c(0,length(rankedGenes)), type="l", axes=FALSE, xlab="",
         ylab = "Enrichment",
         ylim=c(-1,1),
         main=sprintf(paste0("Gene Set: %s\n",
                             "Rank Statistic: %s\n",
                             "AUC: %0.3f; %0.3f"),
                      listName,
                      rankStatName,
                      integrated.area.pos / max.area.pos,
                      integrated.area.neg / max.area.neg));
    points(geneOrder$cMod, type="l", col="black");
    ##points(geneOrder$cModList-1, type="l", col="darkgreen");
    ##points(geneOrder$cModNon+1, type="l", col="magenta");
    posCorePoint <- head(which(geneOrder$cMod == max(geneOrder$cMod)), 1);
    posCoreSize <- sum(rankedGenes[1:posCorePoint] %in% targetList);
    negCorePoint <- tail(which(geneOrder$cMod == min(geneOrder$cMod)), 1);
    negCoreSize <- sum(rankedGenes[negCorePoint:length(rankedGenes)] %in% targetList);
    abline(h=c(max(geneOrder$cMod), min(geneOrder$cMod)), col="#A0A0A080");
    text(x=posCorePoint, y=max(geneOrder$cMod),
         labels=sprintf("%0.3f (%d genes)", max(geneOrder$cMod), posCoreSize), adj=c(0.1, -0.5));
    text(x=negCorePoint, y=min(geneOrder$cMod),
         labels=sprintf("(%d genes) %0.3f", negCoreSize, min(geneOrder$cMod)), adj=c(0.9, 1.5));
    abline(h=0);
    if(zeroPoint < length(rankedGenes)){
        segments(x0=zeroPoint, y0=-0.1, y1=0.1, col="#208020A0");
    }
    axis(2);
    axis(1, at=which(geneOrder$inList), las=3,
         labels=geneOrder$Gene[geneOrder$inList], cex.axis=0.5);
    return(max(integrated.area.pos / max.area.pos, integrated.area.neg / max.area.neg));
}

## TOKEEP: information here regarding calculating rowMeans and creating output
## [used for Figure 7b]

pdf(sprintf("GSEA_relDiff_GSEASets_vs_OLClassification_%s.pdf", dateStr),
    width=6, height=6);
sct.count.mat <- bs[["SCT"]]@data;
colSets <- list(
    "JEM cDC/monocyte; skin vs blood" =
        list("pop1.cols" = which((bs[["cellType"]]$cellType == "monocyte") &
                                 (bs[["cellOrigin"]]$cellOrigin == "skin")),
             "pop2.cols" = which((bs[["cellType"]]$cellType == "monocyte") &
                                 (bs[["cellOrigin"]]$cellOrigin == "blood")))
    );
for(li in names(colSets)){
    cat(sprintf("** %s **\n", li));
    pop1.cols <- colSets[[li]][["pop1.cols"]];
    pop2.cols <- colSets[[li]][["pop2.cols"]];
    pop.diff <- rowMeans(as.matrix(sct.count.mat[,pop1.cols])) -
        rowMeans(as.matrix(sct.count.mat[,pop2.cols]));
    pop.diff <- pop.diff[(rowMeans(as.matrix(sct.count.mat[,pop1.cols])) > 0) |
                         (rowMeans(as.matrix(sct.count.mat[,pop2.cols]) > 0))];
    pop.diff <- pop.diff[order(-abs(pop.diff))];
    ## filter out duplicates
    names(pop.diff) <- sub("\\.[0-9]+$","",
                           sub("^.*?\\.", "", names(pop.diff)));
    pop.diff <- pop.diff[names(pop.diff) != "NA"];
    ## retain the most significant of each named gene
    pop.diff <- pop.diff[match(unique(names(pop.diff)), names(pop.diff))];
    ## reorder, just in case unique fiddled with the order
    pop.diff <- pop.diff[order(-(pop.diff))];
    write_csv(tibble(
        "Dataset" = li,
        "Gene" = names(pop.diff),
        "Difference" = round(pop.diff, 3)),
        path=sprintf("geneRank_%s_%s.csv", gsub("[/ ;]+",".",li), dateStr));
    ## ^^^^ STOP HERE ^^^^

    for(listName in names(geneList)){
        print(listName);
        gsea(names(pop.diff), geneList[[listName]],
             listName=listName,
             rankStatName=li, zeroPoint=sum(pop.diff > 0));
    }
}
invisible(dev.off());


coreGenes <-
    list("IL-4 AND IL-13 PATHWAY (REACTOME)" =
             c("CCL2", "FOS", "SOCS3", "F13A1", "IL13RA1", "PIM1",
               "CEBPD", "TIMP1", "CD36", "IL4R", "POU2F1", "CDKN1A",
               "STAT1", "BATF", "TGFB1", "RHOU", "TYK2", "HIF1A",
               "FOXO3", "JUNB", "MCL1", "STAT3", "ITGAM", "PIK3R1",
               "FN1", "CCL22"),
         "T-HELPER 2 CELL" =
             c("BCL3", "IL4R", "RSAD2", "HLX", "BATF", "IL17RA",
               "LGALS9"),
         "TYPE 2 IMMUNE RESPONSE" =
             c("BCL3", "IL4R", "CCR2", "RSAD2", "HLX", "BATF"),
         "IL4 SIGNALING (WIKI PATHWAY)" =
             c("FOS", "SOCS3", "RELA", "IL4R", "NFKBIA", "STAT1",
               "MAPK1", "TYK2", "STAT3", "PIK3CD", "ATF2", "PIK3R1",
               "FES", "STAT5A", "MAPK3"),
         "TN SIGNATURE 2016" =
             c("CALM3", "CLMN", "SHISA9", "TVP23A", "CLCN5", "SCD",
               "NETO2", "SPECC1", "TCAF2", "TMEM25", "EYA1", "CNNM2",
               "DNAJC22", "MGLL", "CD109", "MSRB3", "GLCCI1",
               "DENND5B", "STAT4", "CEACAM1", "GDPD1", "FAM107B",
               "LRMP", "MYB", "FAM234B", "SRL", "HLA-DOA", "CCDC120"),
         "TN SIGNATURE 2019" =
             c("HLA-DMA", "HLA-DMB", "PRNP", "CLMN", "SPATS2",
               "STOML1", "GAA", "MILR1", "CLCN5", "SCD", "NETO2",
               "ZNF532", "TCAF2", "CIITA", "EYA1", "CNNM2", "ZHX3",
               "SEMA6A", "NREP", "MSRB3", "GLCCI1", "DENND5B",
               "GPRC5A", "FOXRED2", "ENO3", "MYB", "GCNT2", "RAB29",
               "FAM234B", "SRL", "HLA-DOA", "PHC1", "CRISPLD2",
               "CAMK2B", "TSHZ2", "PTPRO", "IL18R1"),
         "CD11B SIGNATURE 2016" =
             c("ENG", "ISG15", "THBS1", "SLC43A2", "DRAM1", "IFITM3",
               "IFIT1", "ITPR2", "SLC4A7", "CMPK2", "CMKLR1", "IL1RN",
               "PLXDC2", "BLVRB", "SIRPA", "UNC93B1", "LRP1",
               "CLEC4E", "PDCD1LG2", "IFI27", "MX2", "FCGR2B", "MSR1",
               "IFIT3", "BCAM", "TBC1D2B", "OAS2", "TIMP2",
               "LGALS3BP", "F11R", "SMAD3", "ACSL1", "TMEM159", "VDR",
               "AXL", "CSF3R", "TFEC", "CCL5", "S100A10", "PHF11",
               "ITGAM", "IRF7", "ITGB8", "CLU", "ATP13A2", "FN1",
               "ADAP2", "SLC7A11", "FES", "EMB", "FAM83F", "PGD",
               "AHRR", "DTX4", "ARHGAP26", "VRK2", "KLRK1",
               "MAPKAPK3", "PLXNB2", "ABCC3", "CSF1R", "GBGT1",
               "MRPL38", "CASS4", "CD53", "P2RY14"),
         "CD11B SIGNATURE 2019" =
             c("LILRB4", "CD36", "IL1RN", "PLXDC2", "BLVRB", "CCR2",
               "CLEC4E", "FCGR2B", "BCAM", "LGALS3BP", "F11R",
               "KCTD12", "TFEC", "CCL5", "PARVB", "S100A10", "IRF7"));

coreGenes.tbl <- as.tbl(data.frame(Reduce(rbind, sapply(names(coreGenes), function(x){cbind(Pathway=x, Gene=coreGenes[[x]])}))));
write_csv(coreGenes.tbl, "data/JEM_coreGenes.csv");

dhm <- function(listName){
    DoHeatmap(bs, slot = "data", assay="SCT",
              cells=names(Idents(bs))[(bs[["cellType"]]$cellType == "monocyte")],
              group.by = "cellOrigin",
              features=geneLookup[geneList[[listName]]]) +
        ggtitle(listName) + theme(plot.title = element_text(hjust = 0.5, size=20));
}

ggsave("out.pdf", width=20, height=30,
       plot=dhm("TYPE 2 IMMUNE RESPONSE") +
       dhm("T-HELPER 2 CELL") +
       dhm("IL-4 AND IL-13 PATHWAY (REACTOME)") +
       dhm("IL4 SIGNALING (WIKI PATHWAY)") +
       dhm("CD11B SIGNATURE 2016") +
       dhm("CD11B SIGNATURE 2019") +
       dhm("TN SIGNATURE 2019") +
       dhm("TN SIGNATURE 2016"));

dhm2 <- function(myGeneList){
    DoHeatmap(bs, slot = "data",## assay="SCT",
              cells = which(Idents(bs) != "5"),
              group.by = "cellType",
              features=geneLookup[myGeneList]) +
        ggtitle(listName) + theme(plot.title = element_text(hjust = 0.5, size=20));
}

mgl <- c("C12orf75", "ALOX5AP", "GZMB", "ITM2C", "PLD4", "JCHAIN",
       "MZB1", "IRF7", "IRF8", "PLAC8", "TCF4", "SERPINF1", "UGCG",
       "PTGDS", "TCL1A", "S100A9", "S100A8", "CD14", "TIMP1",
       "SERPINA1", " FCGR3A", "S100A4", "COTL1", "CTSS", "FCN1",
       "FGR", "GLUL", " CTSL", "FOS", "CD163", "BIRC3", "IDO1",
       "CD83", "CCR7", "LAMP3 ", "CST7", "RGS1", "WFDC21P", "TBC1D4",
       "DAPP1", "POGLUT1", " CCL19", "RASSF4", "CCL17", "IL7R");
dhm2(mgl)
