#!/usr/bin/env Rscript

library(tidyverse);

dirNames <- list.dirs("microarray_2020-Aug-19", recursive=FALSE);
names(dirNames) <- list.dirs("microarray_2020-Aug-19", recursive=FALSE,
                             full.names=FALSE);

results.all <- list();

dateStr <- format(Sys.Date(), "%Y-%b-%d");

pdf(file=sprintf("GSEA_plots_combined_%s.pdf", dateStr), width=11, height=8);
layout(matrix(1:4,2,2));
for(dsName in names(dirNames)){
    dn <- dirNames[dsName];
    allGenesFile <- list.files(dn, pattern="^ranked");
    geneListFile <- grep("^ranked", list.files(dn, pattern="\\.xls$"),
                         invert=TRUE, value=TRUE);
    imageFile <- list.files(dn, pattern="^enplot.*\\.png");
    ## Collect input files
    all.genes <- read_delim(paste(dn, allGenesFile, sep="/"), delim="\t");
    list.genes <- read_delim(paste(dn, geneListFile, sep="/"), delim="\t");
    ## Fetch known scores and fill in rank score changes for all genes
    rank.scores <- list.genes$`RANK METRIC SCORE`;
    rank.val <- round(1 / sum(abs(rank.scores)) * abs(rank.scores), 5);
    es.drop <- 1 / (nrow(all.genes) - nrow(list.genes));
    add.value <- rep(0, length.out=nrow(all.genes));
    add.adj.value <- rep(0, length.out=nrow(all.genes));
    all.genes$inSet <- FALSE;
    all.genes$inSet[(list.genes$`RANK IN GENE LIST`+1)] <- TRUE;
    zero.point <- which.min(all.genes$SCORE > 0);
    add.value[(list.genes$`RANK IN GENE LIST`+1)] <- rank.val;
    vals.pos <- sum(head(add.value, zero.point));
    vals.neg <- sum(tail(add.value, -zero.point));
    drop.pos <- vals.pos / sum(head(add.value == 0, zero.point));
    drop.neg <- vals.neg / sum(tail(add.value == 0, zero.point));
    add.value[-(list.genes$`RANK IN GENE LIST`+1)] <- -es.drop;
    ## Calculate adjusted values (zero point at zero crossover)
    add.adj.value <-
        ifelse(seq_along(add.adj.value) < zero.point, -drop.pos, -drop.neg);
    add.adj.value[(list.genes$`RANK IN GENE LIST`+1)] <- rank.val;
    ## Calculate running scores
    options(scipen=15);
    all.genes$addValue <- signif(add.value, 4);
    all.genes$addAdjValue <- signif(add.adj.value, 4);
    all.genes$runningScore <- signif(cumsum(add.value), 4);
    all.genes$runningAdjScore <- signif(cumsum(add.adj.value), 4);
    ## Determine PNG location
    all.genes$ypos <- cumsum(add.value);
    all.genes$ypos <- (all.genes$ypos - min(all.genes$ypos)) /
        diff(range(all.genes$ypos));
    all.genes$ypos <- round((1 - all.genes$ypos) * (237-80) + (80), 1);
    all.genes$xpos <- round(seq(62, 471, length.out=nrow(all.genes)), 1);
    all.genes$listID <- names(dirNames)[dirNames == dn];
    write_csv(all.genes,
              gzfile(paste0(dn, "/withES_", sub("\\.[^\\.]+$", "", allGenesFile),
                            ".csv.gz")));
    ## Add to results array
    results.all[[names(dirNames)[dirNames == dn]]] <- all.genes;
    ## Create plot of Running ES
    plot(y=c(0,all.genes$runningScore), x=0:nrow(all.genes),
         xlab = "Gene Rank", ylab="Running Enrichment Score",
         type="l", lwd=2, main=dsName, ylim=c(-0.2,0.7));
    ## Add gene set markers
    segments(x0=which(all.genes$inSet), y0=-0.2, y1=-0.16,
             col=c("darkred", "dodgerblue")[(all.genes$SCORE > 0)[all.genes$inSet] + 1]);
    ## Add Zero point mark
    points(x=which.min(all.genes$SCORE > 0),
           y=all.genes$runningScore[which.min(all.genes$SCORE > 0)]);
}
for(dsName in names(dirNames)){
    dn <- dirNames[dsName];
    allGenesFile <- list.files(dn, pattern="^ranked");
    geneListFile <- grep("^ranked", list.files(dn, pattern="\\.xls$"),
                         invert=TRUE, value=TRUE);
    imageFile <- list.files(dn, pattern="^enplot.*\\.png");
    ## Collect input files
    all.genes <- read_delim(paste(dn, allGenesFile, sep="/"), delim="\t");
    list.genes <- read_delim(paste(dn, geneListFile, sep="/"), delim="\t");
    ## Fetch known scores and fill in rank score changes for all genes
    rank.scores <- list.genes$`RANK METRIC SCORE`;
    rank.val <- round(1 / sum(abs(rank.scores)) * abs(rank.scores), 5);
    es.drop <- 1 / (nrow(all.genes) - nrow(list.genes));
    add.value <- rep(0, length.out=nrow(all.genes));
    add.adj.value <- rep(0, length.out=nrow(all.genes));
    all.genes$inSet <- FALSE;
    all.genes$inSet[(list.genes$`RANK IN GENE LIST`+1)] <- TRUE;
    zero.point <- which.min(all.genes$SCORE > 0);
    add.value[(list.genes$`RANK IN GENE LIST`+1)] <- rank.val;
    vals.pos <- sum(head(add.value, zero.point));
    vals.neg <- sum(tail(add.value, -zero.point));
    drop.pos <- vals.pos / sum(head(add.value == 0, zero.point));
    drop.neg <- vals.neg / sum(tail(add.value == 0, -zero.point));
    add.value[-(list.genes$`RANK IN GENE LIST`+1)] <- -es.drop;
    ## Calculate adjusted values (zero point at zero crossover)
    add.adj.value <-
        ifelse(seq_along(add.adj.value) < zero.point, -drop.pos, -drop.neg);
    add.adj.value[(list.genes$`RANK IN GENE LIST`+1)] <- rank.val;
    ## Calculate running scores
    options(scipen=15);
    all.genes$addValue <- signif(add.value, 4);
    all.genes$addAdjValue <- signif(add.adj.value, 4);
    all.genes$runningScore <- signif(cumsum(add.value), 4);
    all.genes$runningAdjScore <- signif(cumsum(add.adj.value), 4);
    ## Determine PNG location
    all.genes$ypos <- cumsum(add.value);
    all.genes$ypos <- (all.genes$ypos - min(all.genes$ypos)) /
        diff(range(all.genes$ypos));
    all.genes$ypos <- round((1 - all.genes$ypos) * (237-80) + (80), 1);
    all.genes$xpos <- round(seq(62, 471, length.out=nrow(all.genes)), 1);
    all.genes$listID <- names(dirNames)[dirNames == dn];
    write_csv(all.genes,
              gzfile(paste0(dn, "/withES_", sub("\\.[^\\.]+$", "", allGenesFile),
                            ".csv.gz")));
    ## Add to results array
    results.all[[names(dirNames)[dirNames == dn]]] <- all.genes;
    ## Create adjusted plot of Running ES
    plot(y=c(0,all.genes$runningAdjScore), x=0:nrow(all.genes),
         xlab = "Gene Rank", ylab="Running Enrichment Score [Adj]",
         type="l", lwd=2, main=dsName, ylim=c(-0.3,0.7));
    ## Add gene set markers
    segments(x0=which(all.genes$inSet), y0=-0.3, y1=-0.26,
             col=c("darkred", "dodgerblue")[(all.genes$SCORE > 0)[all.genes$inSet] + 1]);
    ## Add Zero point Mark
    points(x=which.min(all.genes$SCORE > 0),
           y=all.genes$runningAdjScore[which.min(all.genes$SCORE > 0)]);
}
invisible(dev.off());
