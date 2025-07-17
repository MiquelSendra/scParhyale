#!/bin/bash

for i in 1 2 3 4 5 6; do
  echo "Processing library $i..."

  ~/workspace/scRNAseq/tools/STAR-2.7.11a/source/STAR \
    --runThreadN 60 \
    --genomeDir ~/workspace/scRNAseq/genome_index \
    --readFilesIn \
      ~/workspace/scRNAseq/Data/fastq/PIPseq_lib0${i}_S${i}_R2_001.fastq.gz \
      ~/workspace/scRNAseq/Data/fastq/PIPseq_lib0${i}_S${i}_R1_001.fastq.gz \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --soloCBwhitelist None \
    --soloFeatures Gene \
    --soloBarcodeReadLength 0 \
    --outFileNamePrefix ~/workspace/scRNAseq/results/lib0${i}_
done
