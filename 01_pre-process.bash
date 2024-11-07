#!/bin/bash
########################################################
# 1.Preprocessing of raw sequencing data for Graft-seq #
########################################################
SEQLIBS=(Sample_CON1 Sample_CON2 Sample_JUMP1 Sample_JUMP2) # The file name
branch="AGTCACAC" # Custom branch sequence
for seqlib in ${SEQLIBS[@]};
do
len=${#branch}
gunzip ${seqlib}_R1.fastq.gz # Read1 data was processed
echo ${seqlib#*_} >> /data/nohup.out
bash /code/delNoBranch.bash $branch ${seqlib} >> /data/nohup.out
nohup fastp -f ${len} -a CCCATTCACTCTGCGTTGATACCACTGCTT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -x -q 25 -w 4 -l 25 -i ${seqlib}_R1_${branch}.fastq -o ${seqlib}_R1_clean.fastq >> /data/nohup.out
nohup hisat2 -p 24 -k 1 -x /reference/hg38 -U ${seqlib}_R1_clean.fastq -S ${seqlib}_clean.sam >> /data/nohup.out
samtools view -@ 16 -b -S ${seqlib}_clean.sam -o ${seqlib}.bam
samtools sort -@ 16 ${seqlib}.bam > ${seqlib}.sort.bam
samtools index ${seqlib}.sort.bam ${seqlib}.sort.bam.bai
gzip ${seqlib}_R1.fastq
rm ${seqlib}_R1_${branch}.fastq
rm ${seqlib}_clean.sam
done
