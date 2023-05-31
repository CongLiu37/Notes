tab="/bucket/BourguignonU/Cong/gut_microbiome/summaryGutTrans.tsv" # tsv, no header
        # Example:
        # #tag  fq1 fq2 hostGenome  bam
        # Aaac  Aaac1.1.fq.gz,Aaac2.1.fq.gz Aaac1.2.fq.gz,Aaac2.2.fq.gz Aaac.fna Aaac1.bam,Aaac2.bam
        # Aban  Aban1.1.fq.gz,Aban2.1.fq.gz Aban1.2.fq.gz,Aban2.2.fq.gz Aban.fna none
# External database
sortmerna.db="/bucket/BourguignonU/Cong/public_db/sortmerna/smr_v4.3_default_db.fasta"
uniprot.dmdb="/apps/unit/BioinfoUgrp/DB/diamondDB/uniref90/2020_06/uniref90.dmnd" # DIAMOND protein db, uniprot
nr.dmdb="/apps/unit/BioinfoUgrp/DB/diamondDB/ncbi/2022-07/nr.dmnd"
pfam="/bucket/BourguignonU/Cong/public_db/pfam/Pfam-A.hmm" # Pfam-A.hmm
megan.db="/bucket/BourguignonU/Cong/public_db/megan.db/megan-map-Feb2022.db"
out_dir="/flash/BourguignonU/Cong/gut_microbiome/"
save_dir="/bucket/BourguignonU/Cong/gut_microbiome/"
threads=64
Step=1
