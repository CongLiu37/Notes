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
# tRNAscan-SE
if (1 %in% step){
  Stamp1=paste(out_dir,"/tRNAscanSE/",label,"/",label,"_tRNAscanSE.finished",sep="")
  Stamp2=paste(save_dir,"/tRNAscanSE/",label,"/",label,"_tRNAscanSE.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 1: tRNAscanSE FINISHED")
  }else{
    print("Step 1: tRNAscanSE START")
  
    if (!file.exists(paste(out_dir,"/tRNAscanSE/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/tRNAscanSE/",sep=""))
    }
    tRNAscan(genome=genome,
             mode="eukaryotic", # eukaryotic/bacterial/archaeal/mitochondrial_mammal/mitochondrial_vert/other_organellar/general
             out_dir=paste(out_dir,"/tRNAscanSE/",label,sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# miRNAture
if (2 %in% step){
  Stamp1=paste(out_dir,"/miRNAture/",label,"/",label,"_miRNAture.finished",sep="")
  Stamp2=paste(save_dir,"/miRNAture/",label,"/",label,"_miRNAture.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 2: miRNAture FINISHED")
  }else{
    print("Step 2: miRNAture START")
    
    if (!file.exists(paste(out_dir,"/miRNAture/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/miRNAture/",sep=""))
    }
    miRNAture(genome=genome,
              dataF=miRNAture.db, # pre-calculated data directory from miRNAture
              species=paste(label,"_mirna",sep=""), # Genus_species
              out_dir=paste(out_dir,"/miRNAture/",label,sep=""),
              threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# miRanda
if (3 %in% step){
  Stamp1=paste(out_dir,"/miranda/",label,"/",label,"_miranda.finished",sep="")
  Stamp2=paste(save_dir,"/miranda/",label,"/",label,"_miranda.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 3: miranda FINISHED")
  }else{
    print("Step 3: miranda START")
  
    if (!file.exists(paste(out_dir,"/miranda/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/miranda/",sep=""))
    }
    miranda(gff=paste(save_dir,"/protein/",label,"/",label,"_genes.gff3",sep=""), # three_prime_UTR feature required
            genome=genome,
            miRNA.fna=paste(save_dir,"/miRNAture/",label,"/miRNA.fna",sep=""),
            threads=threads,
            out_dir=paste(out_dir,"/miranda/",label,sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}  
############################################
# Infernal
if (4 %in% step){
  Stamp1=paste(out_dir,"/infernal/",label,"/",label,"_infernal.finished",sep="")
  Stamp2=paste(save_dir,"/infernal/",label,"/",label,"_infernal.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 4: infernal FINISHED")
  }else{
    print("Step 4: infernal START")
    
    if (!file.exists(paste(out_dir,"/infernal/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/infernal/",sep=""))
    }
    infernal(genome=genome,
             domain="eukarya", # eukarya/archaea/bacteria Which domain shall be retained?
             Rfam.clanin=Rfam.clanin,
             Rfam.cm=Rfam.cm,
             out_dir=paste(out_dir,"/infernal/",label,sep=""),
            threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# Final
if (5 %in% step){
  print("Step 5: ncRNA")
  if (!file.exists(paste(out_dir,"/ncRNA/",sep=""))){
    system(paste("mkdir"," ",out_dir,"/ncRNA/",sep=""))
  }
  ncRNA(tRNA.gff3=paste(save_dir,"/tRNAscanSE/",label,"/tRNA.gff3",sep=""),
        miRNAture.gff3=paste(save_dir,"/miRNAture/",label,"/miRNAture.gff3",sep=""),
        miRNA.fa=paste(save_dir,"/miRNAture/",label,"/miRNA.fna",sep=""),
        infernal.gff3=paste(save_dir,"/infernal/",label,"/Infernal.gff3",sep=""),
        miRanda.tsv=paste(save_dir,"/miranda/",label,"/miRanda.tsv",sep=""),
        species=label,
        out_dir=paste(out_dir,"/ncRNA/",label,sep=""))
}