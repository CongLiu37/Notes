Genome/Metagenome assembly, binning and assembly quality assessment

Dependencies:
   Softwares: PhyloFlash, FastQC, Trimmomatic, BUSCO
   R packages: ape, Biostrings

(1) PhyloFlash
Taxon profiling of reads/assembly via reconstructing the SSU rRNAs, using PhyloFlash.
Dependent on PhyloFlash. 
PhyloFlash is dependent on Perl >= 5.13.2, EMIRGE and its dependencies, BBmap, Vsearch, SPAdes, Bedtools, Mafft, Barrnap.

(2) QualityCheck
Check quality of sequencing data by FastQC.
Dependent on FastQC. 

(3) QualityFilter
Quality control of sequencing data by Trimmomatic.
Dependent on Trimmomatic.

(4) SPAdes
Genome/metagenome assembly by SPAdes.
Dependent on SPAdes.

(5) BSG
Compute statistics of assembly, including N50, min/max/median contig size, gap percent, GC content. 
Dependent on R packages: ape, Biostrings.

(6) BUSCO
Assess genome completeness via searching universal single-copy orthologue genes by BUSCO.
Dependent on BUSCO.

(7) checkm
Assess completeness and contamination of genomic bins via using collocated sets of genes that are ubiquitous and single-copy within a phylogenetic lineage by CheckM.
Dependent on CheckM. 
CheckM is dependent on python3, HMMER (>=3.1b1), prodigal (2.60 or >=2.6.1), pplacer (>=1.1). 

(8) Bowtie2Build
Build bowtie2 index of genome by Bowtie2.
Dependent on Bowtie2.

(9) Bowtie2
Map reads to genome by Bowtie2, compress SAM to BAM and sort BAM by SAMtools.
Dependent on Bowtie2 and SAMtools.

(10) SprayNPray
Binning based on tetranucleotide frequency, GC-content, codon usage bias, and read coverage, by SprayNPray.
Dependent on SprayNPray.
SprayNPray is dependent on DIAMOND, Prodigal, Metabat, Python3, Biopython3, Joblib.
