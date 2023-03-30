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
# mite_hunter
if (1 %in% step){
  Stamp1=paste(out_dir,"/mite_hunter/",label,"/",label,"_mite_hunter.finished",sep="")
  Stamp2=paste(save_dir,"/mite_hunter/",label,"/",label,"_mite_hunter.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 1: mite_hunter FINISHED")
  }else{
    print("Step 1: mite_hunter START")
    if (!file.exists(paste(out_dir,"/mite_hunter/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/mite_hunter/",sep=""))
    }
    mite_hunter(genome=genome,
                out_dir=paste(out_dir,"/mite_hunter/",label,sep=""),
                out_basename=label,
                threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# ltr
if (2 %in% step){
  Stamp1=paste(out_dir,"/ltr/",label,"/",label,"_ltr.finished",sep="")
  Stamp2=paste(save_dir,"/ltr/",label,"/",label,"_ltr.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 2: ltr FINISHED")
  }else{
    print("Step 2: ltr START")
    if (!file.exists(paste(out_dir,"/ltr/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/ltr/",sep=""))
    }
    ltr(fna=genome,
        out_dir=paste(out_dir,"/ltr/",label,sep=""),
        out_prefix=label,
        threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# RepeatMasker: Mask MITE and LTR
if (3 %in% step){
  Stamp1=paste(out_dir,"/repeatmasker/",label,"/",label,"_repeatmasker.finished",sep="")
  Stamp2=paste(save_dir,"/repeatmasker/",label,"/",label,"_repeatmasker.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 3: repeatmasker FINISHED")
  }else{
    print("Step 3: repeatmasker START")
    if (!file.exists(paste(out_dir,"/repeatmasker/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/repeatmasker/",sep=""))
    }
    repeatmasker(fna=genome,# Fasta file of genome.
                 out_dir=paste(out_dir,"/repeatmasker/",label,sep=""),
                 RepeatLib.fa=paste(save_dir,"/mite_hunter/",label,"/",label,"_MITE.fa"," ",
                                    save_dir,"/ltr/",label,"/LTR.fa",
                                    sep=""), # space-list
                 Threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# RepeatModeler no LTR
if (4 %in% step){
  Stamp1=paste(out_dir,"/repeatmodeler_noLTR/",label,"/",label,"_repeatmodeler_noLTR.finished",sep="")
  Stamp2=paste(save_dir,"/repeatmodeler_noLTR/",label,"/",label,"_repeatmodeler_noLTR.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 4: repeatmodeler_noLTR FINISHED")
  }else{
    print("Step 4: repeatmodeler_noLTR START")
    if (!file.exists(paste(out_dir,"/repeatmodeler_noLTR/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/repeatmodeler_noLTR/",sep=""))
    }
    repeatmodeler_noLTR(fna=paste(save_dir,"/repeatmasker/",label,"/genome.fa.masked",sep=""),
                        out_dir=paste(out_dir,"/repeatmodeler_noLTR/",label,sep=""),
                        out_prefix=label,
                        Threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# Annotate repeats with repeatmasker
if (5 %in% step){
  Stamp1=paste(out_dir,"/repeat_elements/",label,"/",label,"_repeat_elements.finished",sep="")
  Stamp2=paste(save_dir,"/repeat_elements/",label,"/",label,"_repeat_elements.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 5: repeat_elements FINISHED")
  }else{
    print("Step 5: repeat_elements START")
    if (!file.exists(paste(out_dir,"/repeat_elements/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/repeat_elements/",sep=""))
    }
    repeatmasker(fna=genome,# Fasta file of genome.
                 out_dir=paste(out_dir,"/repeat_elements/",label,sep=""),
                 RepeatLib.fa=paste(save_dir,"/mite_hunter/",label,"/",label,"_MITE.fa"," ",
                                    save_dir,"/ltr/",label,"/LTR.fa"," ",
                                    save_dir,"/repeatmodeler_noLTR/",label,"/noLTR_rep.lib.fa",
                                    sep=""), # space-list
                 Threads=threads)
    system(paste("rm ",out_dir,"/repeat_elements/",label,"/genome.fa.cat.gz"))
    system(paste("rm ",out_dir,"/repeat_elements/",label,"/genome.fa.masked"))
    system(paste("rm ",out_dir,"/repeat_elements/",label,"/genome.fa.out.gff"))
  
    rep_elements(repeatmasker.out=paste(out_dir,"/repeat_elements/",label,"/genome.fa.out",sep=""),
                 # Clarify source of repeat elements
                 MITE_Hunter.replib=paste(save_dir,"/mite_hunter/",label,"/",label,"_MITE.fa",sep=""), # RepLib from MITE-Hunter
                 LTR_retriever.replib=paste(save_dir,"/ltr/",label,"/LTR.fa",sep=""), # RepLib from LTR_retriever
                 RepeatModeler.replib=paste(save_dir,"/repeatmodeler_noLTR/",label,"/noLTR_rep.lib.fa",sep=""), # RepLib from RepeatModeler
                 all.replib=paste(out_dir,"/repeat_elements/",label,"/",label,"_replib.fa",sep=""),
                 out.gff3=paste(out_dir,"/repeat_elements/",label,"/",label,"_repeat_elements.gff3",sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}