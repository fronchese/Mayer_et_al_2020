# STAR 2.7.1a https://github.com/alexdobin/STAR/releases/tag/2.7.1a docker and code for alignment to mouse genome mm10 M16 downloaded from https://www.gencodegenes.org/
docker run -v /stornext/BioinfoSAN/Bioinformatics/Illumina/OG4953/OG4953_fastq/combined/trimmed:/data/ -it samold_staraligner/star2.7.1a:latest bash
STAR --runThreadN 12 -- runMode genomeGenerate --genomeDir /data/STAR_Mapped/star_genome_directory --genomeFastaFiles /data/STAR_Mapped/star_genome/GRCm38.primary_assembly.genome.fa --sjdbGTFfile /data/STAR_Mapped/star_genome/gencode.vM16.primary_assembly.annotation.gff3
for x in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36;
do echo mapping $x;
STAR --runThreadN 12 --genomeDir /data/STAR_mapped/star_genome_directory --readFilesIn /data/4953_S${x}_1P.fq.gz /data/4953_S${x}_2P.fq.gz --readFilesCommand gunzip -c --outFileNamePrefix /data/STAR_mapped/output/${x}_STAR_9.Oct.19_
done

# Convert to sorted bam
for x in 01 02 03 04 05 06 07 08 09 {10..36};
do echo converting $x;
samtools view -u /data/${x}_STAR_9.Oct.19_Aligned.out.sam | samtools sort -@ 4 -o ${x}.bam
done

# EOF
