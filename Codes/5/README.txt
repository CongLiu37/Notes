Genome annotation

(1) GenomeMask
Identify and mask repeat fragments.
Dependent on RepeatModeler, RepeatMasker.
RepeatModeler and RepeatMasker can be installed by singularity container.

(2) GenomeThreader
Homology-based gene prediction.
Dependent on GenomeThreader.

(3) Hisat2Build
Build Hisat2 index.
Dependent on Hisat2.

(4) Hisat
Map short reads to assembly, sort BAM and index BAM.
Dependent on Hisat2, SAMtools.

(5) MergeBAM
Merge BAM files.
Dependent on SAMtools.

(6) StringTie
Assemble transcripts for gene prediction.
Dependent on StringTie.

(6) TransDecoder
Gene prediction from transcripts.
Dependent on TransDecoder.

(7) MergeGff3
Merge two gff3 files.
Dependent on gff3toolkit.

(8) gtf2read
Extract transcripts and protein seuqences from gtf.
Dependent on Cufflinks.

(9) Augustus
De novo gene prediction.
Dependent on AUgustus.
