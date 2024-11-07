#!/bin/bash
####################################
# 1.1 Filter reads without branche #
####################################
branch=$1
file=$2
echo -e "\e[32m Start Filtering non-branch reads. \e[0m"
line=`cat ${file}_R1.fastq | wc -l`
let line=line/4
echo -e "\e[32m Orign reads number: $line \e[0m"
grep -A 2 -B 1 "^$branch" ${file}_R1.fastq > ${file}_R1_${branch}.fastq
line=`cat ${file}_R1_${branch}.fastq | wc -l`
let line=line/4
echo -e "\e[32m Filtering reads number: $line \e[0m"
