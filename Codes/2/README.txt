Search for specific pattern in promoters of a genome
Promoter: 1000 bp upstream to 500 bp downstream of the transcription start site

Dependencies:
  Softwares: BEDtools
  R packages: Biostrings

(1) gff2bed_promoters
Extract promoter locations (bed) from gff.

(2) bed2fasta_promoter
Extract promoter fasta from promoter locations (bed) by BEDtools.
Dependent on BEDtools.

(3) PatternInPromoters
Search for specific pattern in promoter fasta.
Dependent on Biostrings.
