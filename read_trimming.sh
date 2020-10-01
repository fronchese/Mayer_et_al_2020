# Concatenate all RNA-Seq exeperiment HiSeq fastq files prior to read trimming
for x in {1..9};
do
  echo merging $x
  cat CDT11ANXX-4953-0$x-50-1_S$x_L00*_R1_001.fastq.gz > combined/4953_S0${x}_R1.fastq.gz
  cat CDT11ANXX-4953-0$x-50-1_S$x_L00*_R2_001.fastq.gz > combined/4953_S0${x}_R2.fastq.gz
done

for x in {10..36};
do
  echo merging $x;
  cat CDT11ANXX-4953-$x-50-1_S$x_L00*_R1_001.fastq.gz > combined/4953_S${x}_R1.fastq.gz
  cat CDT11ANXX-4953-$x-50-1_S$x_L00*_R2_001.fastq.gz > combined/4953_S${x}_R2.fastq.gz
done

# Trimmomatic 0.36 docker and code for read trimming
# docker pull comics/trimmomatic:0.36
docker run -v /stornext/BioinfoSAN/Bioinformatics/Illumina/OG4953/OG4953_fastq/combined:/data/ -it comics/trimmomatic:0.36 bash
ADAPTfile=/software/applications/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa;
for x in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36;
echo trimming sample $x;
java -jar $TRIMMOMATIC PE -threads 12 4953_S${x}_R1.fastq.gz 4953_S${x}_R2.fastq.gz -baseout trimmed/4953_S${x}.fq.gz ILLUMINACLIP:${ADAPTfile}:2:15:8:5:true  LEADING:3  TRAILING:3  SLIDINGWINDOW:10:20  MINLEN:20;
done

# EOF
