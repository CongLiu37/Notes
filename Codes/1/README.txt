Differential expression analysis of RNA-seq with reference genome

Dependencies:
Softwares: SRA-toolkit, FastQC, Trimmomatic, Gffread, Hisat2, SAMtools, StringTie, eggNOG-mapper. 
R packages: DESeq2, ggplot2, ComplexHeatmap, circlize, annotate, stringr, dplyr, pathview.

Functions:
(1) DownloadSRA
Download fastq from NCBI SRA by SRA-toolkit.
Dependent on SRA-toolkit.

(2) QualityCheck
Assess quality of fastq by FastQC.
Dependent on FastQC.

(3) QualityFilter
Quality control of sequencing data by Trimmomatic.
Dependent on Trimmomatic and fasta of adaptor used for adaptor trimming.

(4) gffread
Extract cds from gff.
Dependent on Gffread.

(5) Hisat2Build
Build hisat2 index of reference genome by Hisat2.
Dependent on hisat2.

(6) Hisat
Map fastq to reference genome index by Hisat2. 
Compress SAM to BAM and sort BAM by SAMtools.
Dependent on hisat2.

(7) StringTie
Assemble transcripts by STringTie.
Dependent on StringTie.

(8) run.prepDE.py
Run prepDE.py3 to extract gene/transcript read count matrix from output of stringtie.
Dependent on prepDE.py3 from https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual.

(9) DESeq
Conduct differential expression analysis by DESeq2.
Dependent on DESeq2.

(10) Volcano
Plot a volcano plot in pdf.
In volcano plot, log2FoldChange is plotted against -log10(padj).
Dependent on ggplot2.

(11) PCA
Plot principle component analysis.
Dependent on ggplot2.

(12) Heatmap
Plot gene expression heatmap.
Dependent on ComplexHeatmap and circlize.

(13) eggNOG
CDS functional annotation by eggNOG-mapper.
Dependent on eggNOG-mapper.

(14) enrichGO
GO enrichment via Fisher's exact test, adjust p-value to false discovery rate (Benjamini & Hochberg), and trim GO by removing enriched parent if its child is enriched.
Dependent on annotate, stringr, dplyr.

(15) enrichGOplot
Bar plot of 10 most significantly enriched BP/MF/CC GO terms.
Dependent on ggplot2.

(16) PathInf
Infer KEGG pathway of whole genome by Minpath.
Dependent on Minpath, stringr, dplyr.

(17) enrichKEGG
KEGG pathway enrichment via Fisher's exact test, and adjust p-value to false discovery rate (Benjamini & Hochberg).
Dependent on stringr, dplyr.

(18) enrichKEGGplot
Bar plot for 30 most significantly enriched KEGG pathways.
Dependent on ggplot2.

(19) PathwayView
Visualize pathways enriched by DEGs.
Dependent on stringr, dplyr, pathview.

(20) enrichKO
KO enrichment via Fisher's exact test, adjust p-value to false discovery rate (Benjamini & Hochberg).
Dependent on stringr, dplyr.

(21) enrichKOplot
Bar plot of 30 most significantly enriched KOs.
Dependent on ggplot2.
