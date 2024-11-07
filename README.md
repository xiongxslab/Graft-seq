# Graft-seq

## 01_pre-process.bash

   Code for pre-processing raw sequencing data for Graft-seq.
- Filter reads starting with branch using delNoBranch.bash;
- Trim the sequence of TSO and adapters using fastp;
- Align the sequences using Hisat2;
- Convert and sort alignment files.

```
bash 01_pre-process.bash
```

## 02_FindPos.R

   Code for identification A sites in transcriptome. 
- Identify the position of the first base of reads;
- Identify the position of nearest A.
   
```
Rscript 02_FindPos.R
```

## 03_SiteSelect.R

   Code for Filter high-confidence m6A/m6mA sites.
- Filter sites by calculating Exp/Ctrl RPM fold change;
- Filter sites by comparing the results of two reps.

```
Rscript 03_SiteSelect.R
```
