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


# genome='/bucket/BourguignonU/Cong/termite_genome_annotation/Cpun_mark/cryptocercus_punctulatus.scaffolded.fasta'
# RNA1=paste(system("ls /bucket/BourguignonU/Cong/termite_genome_annotation/Cpun_mark/DRR05870*_1*",intern=TRUE),collapse = ",")
# RNA2=paste(system("ls /bucket/BourguignonU/Cong/termite_genome_annotation/Cpun_mark/DRR05870*_2*",intern=TRUE),collapse = ",")
# ref_protein1='/bucket/BourguignonU/Cong/public_db/model_insects/Acyrthosiphon_pisum_GCF_005508785.2.faa,/bucket/BourguignonU/Cong/public_db/model_insects/Drosophila_melanogaster_GCF_000001215.4.faa,/bucket/BourguignonU/Cong/public_db/model_insects/Aeder_aegypti_GCF_002204515.2.faa,/bucket/BourguignonU/Cong/public_db/model_insects/Manduca_sexta_GCF_014839805.1.faa,/bucket/BourguignonU/Cong/public_db/model_insects/Anopheles_gambiae_GCF_000005575.2.faa,/bucket/BourguignonU/Cong/public_db/model_insects/Nasonia_vitripennis_GCF_009193385.2.faa,/bucket/BourguignonU/Cong/public_db/model_insects/Apis_mellifera_GCF_003254395.2.faa,/bucket/BourguignonU/Cong/public_db/model_insects/Solenopsis_invicta_GCF_016802725.1.faa,/bucket/BourguignonU/Cong/public_db/model_insects/Bombyx_mori_GCF_014905235.1.faa,/bucket/BourguignonU/Cong/public_db/model_insects/Tribolium_castaneum_GCF_000002335.3.faa,/bucket/BourguignonU/Cong/public_db/Blattodea/Blattella_germanica.anno.pep.faa,/bucket/BourguignonU/Cong/public_db/Blattodea/Cryptotermes_secundus.anno.pep.faa,/bucket/BourguignonU/Cong/public_db/Blattodea/Coptotermes_formosanus.anno.pep.faa,/bucket/BourguignonU/Cong/public_db/Blattodea/Zootermopsis_nevadensis.anno.pep.faa'
# ref_protein2=paste(system("ls /bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/*/*.faa",intern=TRUE),collapse = ",")
# ref_protein=paste(ref_protein1,ref_protein2,sep=",")
# galba_protein='/bucket/BourguignonU/Cong/public_db/Blattodea/Blattella_germanica.anno.pep.faa,/bucket/BourguignonU/Cong/public_db/Blattodea/Cryptotermes_secundus.anno.pep.faa,/bucket/BourguignonU/Cong/public_db/Blattodea/Coptotermes_formosanus.anno.pep.faa,/bucket/BourguignonU/Cong/public_db/Blattodea/Zootermopsis_nevadensis.anno.pep.faa'
# out_dir='/flash/BourguignonU/Cong/termite_genome_annotation/Cpun_mark'
# save_dir='/bucket/BourguignonU/Cong/termite_genome_annotation/Cpun_mark'
# label='Cpun'
# BUSCO_db='/flash/BourguignonU/Cong/insecta_odb10/'
# UniRef_dmnd='/apps/unit/BioinfoUgrp/DB/diamondDB/uniref90/2020_06/uniref90.dmnd'
# PfamA='/bucket/BourguignonU/Cong/public_db/pfam/Pfam-A.hmm'
# miRNAture.db='/bucket/BourguignonU/Cong/public_db/miRNAture/Dataset_mirnature_Sept21_2022/Data'
# Rfam.clanin='/bucket/BourguignonU/Cong/public_db/Rfam/Rfam.clanin'
# Rfam.cm='/bucket/BourguignonU/Cong/public_db/Rfam/Rfam.cm'
# pasa.alignAssembly.conf='/home/c/c-liu/Softwares/PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt'
# pasa.annotationCompare.conf='/home/c/c-liu/Softwares/PASApipeline/pasa_conf/pasa.annotationCompare.Template.txt'
# threads=64
# step= 1