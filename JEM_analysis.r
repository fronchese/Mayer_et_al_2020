#!/usr/bin/env Rscript

library(Seurat);
library(tidyverse);
library(sctransform);
library(patchwork);

dateStr <- format(Sys.Date(), "%Y-%b-%d");

## In the blood and skin data, clusters were called at a resolution of
## 1.8 using the first six principal components, generated using
## highly variable genes.
## [Note: these were obtained from the paper authors, and are not
##  included in this git repository]

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

#### GSEA for skin vs blood (BloodSkin dataset) ####

## Create Gene comparison scores
sct.gene.mat <- sct.count.mat <- bs[["SCT"]]@data;

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
}


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
