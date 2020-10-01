library(tidyverse) # v1.3.0
library(UpSetR)    # v1.4.0


setwd("/Users/siold/Proj/FR_RNA-Seq/shiny_august2020/FR19/data/")

all_normalised <- read_csv("allnormalised.csv")
all_pval       <- read_csv("all_pval.csv")
all_l2fc       <- read_csv("all_l2fc.csv")
gene_info      <- read_csv("gene_info.csv")


# Recreate UpsetR plot for Johannes ---------------------------------------

dataSelected  <- c("KOm24pvsB6m24p","KOm24nvsB6m24n","KOs103pvsB6s103p","KOs103nvsB6s103n","KOe11bpvsB6e11bp","KOeTNvsB6eTN" )
l2fc_cutoff   <- 1
pval_cutoff   <- 0.05
biotypeFilter <- "protein_coding"

genes_of_biotype <- gene_info %>%
  dplyr::filter(gene_biotype %in% biotypeFilter)

pval_data <- all_pval %>%
  dplyr::filter(gene_name %in% genes_of_biotype$external_gene_name) %>%
  dplyr::select(all_of(dataSelected)) %>% 
  dplyr::rename(`KO vs B6 m24+` = KOm24pvsB6m24p,
                `KO vs B6 m24-` = KOm24nvsB6m24n,
                `KO vs B6 s103+` = KOs103pvsB6s103p,
                `KO vs B6 s103-` = KOs103nvsB6s103n,
                `KO vs B6 e11bHi` = KOe11bpvsB6e11bp,
                `KO vs B6 e11bLo` = KOeTNvsB6eTN)

l2fc_data <- all_l2fc %>%
  dplyr::filter(gene_name %in% genes_of_biotype$external_gene_name) %>%
  dplyr::select(all_of(dataSelected)) %>% 
  dplyr::rename(`KO vs B6 m24+` = KOm24pvsB6m24p,
                `KO vs B6 m24-` = KOm24nvsB6m24n,
                `KO vs B6 s103+` = KOs103pvsB6s103p,
                `KO vs B6 s103-` = KOs103nvsB6s103n,
                `KO vs B6 e11bHi` = KOe11bpvsB6e11bp,
                `KO vs B6 e11bLo` = KOeTNvsB6eTN)

# Prepare UpSetR input data frame using selected cutoffs

pval_data_cut <- pval_data < pval_cutoff
l2fc_data_neg <- l2fc_data < -l2fc_cutoff
l2fc_data_pos <- l2fc_data >  l2fc_cutoff
  
l2fc_data_com <- l2fc_data_neg + l2fc_data_pos

final_data    <- l2fc_data_com + pval_data_cut
final_data    <- final_data == 2
final_data    <- final_data[rowSums(final_data) > 0,]
final_data    <- as.data.frame(final_data + 1 - 1)

# Create final image

pdf(file = "/Users/siold/Desktop/test_upset.pdf", height = 8, width = 9)

upset(data = final_data,
      nintersects = NA,
      nsets = length(dataSelected),
      sets = colnames(l2fc_data),
      point.size = 5,
      group.by = "degree",
      matrix.color = "black",
      main.bar.color = "black",
      sets.bar.color = "black",
      shade.color = "black",
      text.scale = 2)
      
dev.off() 
