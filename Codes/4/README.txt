Genome/Metagenome assembly, binning and assembly quality assessment

(1) PhyloFlash
Taxon profiling of reads/assembly via reconstructing the SSU rRNAs by PhyloFlash.
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

(5) LengthFilter
Extract contigs longer than threshold.
Dependent on R packages: ape.

(6) BSG
Compute statistics of assembly, including N50, min/max/median contig size, gap percent, GC content. 
Dependent on R packages: ape, Biostrings.

(7) BUSCO
Assess genome completeness via searching universal single-copy orthologue genes by BUSCO.
Dependent on BUSCO.

(8) checkm
Assess completeness and contamination of genomic bins via using collocated sets of genes that are ubiquitous and single-copy within a phylogenetic lineage by CheckM.
Dependent on CheckM. 
CheckM is dependent on python3, HMMER (>=3.1b1), prodigal (2.60 or >=2.6.1), pplacer (>=1.1). 

(9) Bowtie2Build
Build bowtie2 index of genome by Bowtie2.
Dependent on Bowtie2.

(10) Bowtie2
Map reads to genome by Bowtie2, compress SAM to BAM, sort and index BAM by SAMtools.
Dependent on Bowtie2, SAMtools.

(11) Metabat
Unsupervised binning by Metabat2.
Dependent on Metabat2.

(12) blastn_blobtools
DNA-DNA search by blastn.
Compute contig length, coverage, GC-content, taxonomy at major ranks.
Dependent on blast, blobtools.

(13) Diamond_Megan
DNA-protein search by diamond.
Compute taxonomy at major ranks by Megan.
Diamond and Megan in long-read mode.
Dependent on Diamond and Megan.

(14) SprayNPray
Predict genes by Prodigal, and compute contig length, coding density, coverage, GC-content, average AAI, domain of closest DIAMOND hits.
Dependent on SprayNPray.
SprayNPray is dependent on DIAMOND, Prodigal, Metabat, Python3, Biopython3, Joblib.

(15) TAGC
Plot taxonomy-annotated GC content-coverage plot.
Dependent on R package ggplot2, stringr.

(16) ContaminationPlot
Evaluate contamination level of draft genome by GC-coverage plot, coverage distribution and GC distribution.
Dependent on R package ggplot2, ggExtra.

(17) Extract_fq
Extract short reads mapped to an assembly.
Dependent on Bowtie2, SAMtools.

(18) kmc
Compute kmer spectrum.
Dependent on KMC.

(19) Genomescope
Estimate homozygous, heterozygous and haploid size.
Dependent on GenomeScope.

(20) Smudgeplot
Estimate ploidity.
Dependent on Smudgeplot.
