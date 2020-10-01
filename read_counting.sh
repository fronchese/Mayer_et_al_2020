## Gene counting using Rsubread version 1.28.1 on R version 3.4.4 (2018-03-15) "Someone to Lean On"
docker run -v /stornext/BioinfoSAN/Bioinformatics/Illumina/OG4953/OG4953_fastq/combined/trimmed/STAR_mapped:/data/ -it rsubread/afterstar:latest R
library(Rsubread)
sampleVector = c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
samFileNames = NULL
for(sampleName in sampleVector){
  samFileNames = c(samFileNames, paste0("/data/output/bamfiles/", sampleName, ".bam"))
}
FCFile = featureCounts(nthreads = 12,
                       files = samFileNames,
                       
                       primaryOnly = FALSE,
                       isPairedEnd = TRUE,
                       checkFragLength = TRUE,
                       allowMultiOverlap = TRUE,
                       countMultiMappingReads = TRUE,
                       
                       GTF.featureType = "gene",
                       GTF.attrType = "gene_name",
                       isGTFAnnotationFile = TRUE,
                       annot.ext = "/data/star_genome/gencode.vM16.primary_assembly.annotation.gff3")
write.csv(FCFile$counts, file="/data/output/RNASeq_raw_counts_gene_STAR_FR_2019-11-01.csv")
write.csv(FCFile$stat, file="/data/output/FCFile$stat_gene_STAR_FR_2019-11-01.csv")
