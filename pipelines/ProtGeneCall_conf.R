# Genome assembly
genome="/bucket/BourguignonU/Cong/public_db/model_insects/Drosophila_melanogaster_GCF_000001215.4.fna" # <file>
RNA1=paste(#"/bucket/BourguignonU/Cong/public_db/test_db/Drosophila/RNA_seq/SRR23047821_1.trimmed.fastq.gz",
           "/bucket/BourguignonU/Cong/public_db/test_db/Drosophila/RNA_seq/SRR23047822_1.trimmed.fastq.gz",
           "/bucket/BourguignonU/Cong/public_db/test_db/Drosophila/RNA_seq/SRR23047823_1.trimmed.fastq.gz",
           sep=",") # comma-list of files
RNA2=paste(#"/bucket/BourguignonU/Cong/public_db/test_db/Drosophila/RNA_seq/SRR23047821_2.trimmed.fastq.gz",
           "/bucket/BourguignonU/Cong/public_db/test_db/Drosophila/RNA_seq/SRR23047822_2.trimmed.fastq.gz",
           "/bucket/BourguignonU/Cong/public_db/test_db/Drosophila/RNA_seq/SRR23047823_2.trimmed.fastq.gz",
           sep=",") # comma-list of fq
ref_protein="/bucket/BourguignonU/Cong/public_db/model_insects/Apis_mellifera_GCF_003254395.2.faa" # comma-list of faa
# Output:
out_dir="/flash/BourguignonU/Cong/pipeline_eval/"
save_dir="/bucket/BourguignonU/Cong/pipeline_eval/"
label="Drosophila"
# External database
UniRef_dmnd="/apps/unit/BioinfoUgrp/DB/diamondDB/uniref90/2020_06/uniref90.dmnd"
PfamA="/bucket/BourguignonU/Cong/public_db/pfam/Pfam-A.hmm"
# External configure files
pasa.alignAssembly.conf="/home/c/c-liu/Softwares/PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt" # Path to PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt
pasa.annotationCompare.conf="/home/c/c-liu/Softwares/PASApipeline/pasa_conf/pasa.annotationCompare.Template.txt" # Path to PASApipeline/pasa_conf/pasa.annotationCompare.Template.txt
# Others
threads=64
step=8
