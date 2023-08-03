arg=commandArgs(trailingOnly = TRUE)

main=arg[1]
conf=arg[2]

source(main)
source(conf)

out_dir=sub("/$","",out_dir)
save_dir=sub("/$","",save_dir)
threads=as.character(threads)
if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}

############################################
# Repeat masking
if (1 %in% step){
  Stamp1=paste(out_dir,"/masking/",label,"/",label,"_masking.finished",sep="")
  Stamp2=paste(save_dir,"/masking/",label,"/",label,"_masking.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 1: Repeat masking FINISHED")
  }else{
    print("Step 1: Repeat masking START")
    if (!file.exists(paste(out_dir,"/masking/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/masking/",sep=""),wait=TRUE)
    }
    repeatmodeler(fna=genome,
                  out_dir=paste(out_dir,"/masking/",label,sep=""),
                  out_prefix=label,
                  Threads=threads)
    repeatmasker(fna=genome,# Fasta file of genome.
                 out_dir=paste(out_dir,"/masking/",label,sep=""),
                 RepeatLib.fa=paste(out_dir,"/masking/",label,"/",label,"_RepeatModeler.db-families.fa",sep=""),
                 Threads=threads)
    system(paste("mv"," ",
                 out_dir,"/masking/",label,"/",genome_file,".masked"," ",
                 out_dir,"/masking/",label,"/",label,".fasta.masked",sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# braker
if (2 %in% step){
  Stamp1=paste(out_dir,"/braker/",label,"/",label,"_braker.finished",sep="")
  Stamp2=paste(save_dir,"/braker/",label,"/",label,"_braker.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 2: braker FINISHED")
  }else{
    print("Step 2: braker START")
    if (!file.exists(paste(out_dir,"/braker/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/braker/",sep=""))
    }
    braker(genome=paste(save_dir,"/masking/",label,"/",basename(genome_file),".fasta.masked",sep=""),
           bam="none", # comma-list
           ref_proteins=odb.faa,
           species=label,
           tsebra.conf="none", # TSEBRA/config/default.cfg 
           out_dir=paste(out_dir,"/braker/",label,sep=""),
           threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}


