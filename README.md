# Graft-seq

## 01_pre-process.bash

   Code for pre-processing raw sequencing data for Graft-seq:
   (1) Filter reads starting with branch using delNoBranch.bash;
   (2) Trim the sequence of TSO and adapters using fastp;
   (3) Align the sequences using Hisat2;
   (4) Convert and sort alignment files.

```
bash 01_pre-process.bash
```

## 02_FindPos.R

   Code for identification A sites in transcriptome. 
   (1) Identify the position of the first base of reads;
   (2) Identify the position of nearest A.
   
```
Rscript 02_FindPos.R
```

## 03_SiteSelect.R

   Code for Filter high-confidence m6A/m6mA sites.
   (1) Filter sites by calculating Exp/Ctrl RPM fold change;
   (2) Filter sites by comparing the results of two reps.

```
Rscript 03_SiteSelect.R
```
