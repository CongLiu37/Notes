# Input
genome="/bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/Aaca/Aaca_maskedGenome.fna"

# Database
edta.sif="/home/c/c-liu/Softwares/EDTA.sif"
DeepTE.model="/bucket/BourguignonU/Cong/public_db/DeepTE/Metazoans_model"
DeepTE.sp="M"
TEsorted.db="rexdb-metazoa"

# Output:
out_dir="/flash/BourguignonU/Cong/termite_genome_annotation/RepCall"
save_dir="/bucket/BourguignonU/Cong/termite_genome_annotation/RepCall"
label="Aaca"

# Others
threads=64
step=1

# d=read.table("~/quality_genomePeptide.tsv",sep="\t",header=TRUE,quote="")
# sp.lst=d[,"Label"]
# dire="/flash/BourguignonU/Cong/termite_genome_annotation/shell"
# for (sp in sp.lst){
#   conf=paste(dire,"/conf/",sp,"_conf.R",sep="")
#   write(paste("genome='/bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/",sp,"/",sp,"_maskedGenome.fna'",sep=""),
#         conf,append=FALSE)
#   write("edta.sif='/home/c/c-liu/Softwares/EDTA.sif'",conf,append=TRUE)
#   write("DeepTE.model='/bucket/BourguignonU/Cong/public_db/DeepTE/Metazoans_model'",conf,append = TRUE)
#   write("DeepTE.sp='M'",conf,append=TRUE)
#   write("TEsorted.db='rexdb-metazoa'",conf,append=TRUE)
#   write("out_dir='/flash/BourguignonU/Cong/termite_genome_annotation/RepCall'",
#         conf,append=TRUE)
#   write("save_dir='/bucket/BourguignonU/Cong/termite_genome_annotation/RepCall'",
#         conf,append=TRUE)
#   write(paste("label='",sp,"'",sep=""),
#         conf,append=TRUE)
#   write("threads=64",
#         conf,append=TRUE)
#   write("step=4",
#         conf,append=TRUE)
# 
#   shell=paste(dire,"/hpc/",sp,".sh",sep="")
#   write("#!/bin/bash",shell,append=FALSE)
#   write("#SBATCH --partition=compute",shell,append=TRUE)
#   write(paste("#SBATCH --job-name=",sp,sep=""),shell,append=TRUE)
#   write(paste("#SBATCH --output=",sp,".o%j",sep=""),shell,append=TRUE)
#   write(paste("#SBATCH --error=",sp,".e%j",sep=""),shell,append=TRUE)
#   write("#SBATCH --time=95:00:00",shell,append=TRUE)
#   write("#SBATCH --mem=256G",shell,append=TRUE)
#   write("#SBATCH --ntasks=1",shell,append=TRUE)
#   write("#SBATCH --cpus-per-task=64",shell,append=TRUE)
#   write("#SBATCH --mail-type=END,FAIL",shell,append=TRUE)
#   write("#SBATCH --mail-user=c.liu@oist.jp",shell,append=TRUE)
#   write("module load bioinfo-ugrp-modules",shell,append=TRUE)
#   write("module load DebianMed/11.2",shell,append=TRUE)
#   write("module load singularity",shell,append=TRUE)
#   write(paste("Rscript ~/pipelines/RepCall_pipeline.R ~/pipelines/RepCall_main.R ",conf,sep=""),shell,append=TRUE)
# }


