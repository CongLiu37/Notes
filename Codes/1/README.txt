Differential expression analysis of RNA-seq by Hisat2, StringTie and DESeq2 (reference genome available)

Dependencies:
Softwares: SRA-toolkit, FastQC, Trimmomatic, Cufflinks, Hisat2, SAMtools, StringTie, eggNOG-mapper. 
R packages: DESeq2, ggplot2, ComplexHeatmap, circlize, annotate, stringr, dplyr

main.R
Major modules.

Run.R
Need to run on HPC.

Hpc.sh
Submitted to HPC.

prepDE.py3
From https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual. It computes expression matrix from output of StringTie.

metadata.txt
An example of metadata. 