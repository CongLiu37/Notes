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
# trash
if (1 %in% step){
  Stamp1=paste(out_dir,"/trash/",label,"/",label,"_trash.finished",sep="")
  Stamp2=paste(save_dir,"/trash/",label,"/",label,"_trash.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 1: trash FINISHED")
  }else{
    print("Step 1: trash START")
    if (!file.exists(paste(out_dir,"/trash/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/trash/",sep=""))
    }
    trash(genome=genome,
          out_dir=paste(out_dir,"/trash/",label,sep=""),
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
# RepeatModeler
if (3 %in% step){
  Stamp1=paste(out_dir,"/repeatmodeler/",label,"/",label,"_repeatmodeler.finished",sep="")
  Stamp2=paste(save_dir,"/repeatmodeler/",label,"/",label,"_repeatmodeler.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 3: repeatmodeler FINISHED")
  }else{
    print("Step 3: repeatmodeler START")
    if (!file.exists(paste(out_dir,"/repeatmodeler/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/repeatmodeler/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/repeatmodeler/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/repeatmodeler/",label,sep=""))
    }
    RMlib=paste(out_dir,"/repeatmodeler/",label,"/",label,"_RepeatModeler.db-families.fa",sep="")
    if (!file.exists(RMlib)){
      cmd=paste("seqkit","seq","-u",
                "-j",threads,
                genome,">",
                paste(out_dir,"/repeatmodeler/",label,"/genome.fa",sep=""),
                sep=" ")
      print(cmd);system(cmd,wait=TRUE)
      repeatmodeler(fna=paste(out_dir,"/repeatmodeler/",label,"/genome.fa",sep=""),
                    out_dir=paste(out_dir,"/repeatmodeler/",label,sep=""),
                    out_prefix=label,
                    Threads=threads)
    }
    # cmd=paste("seqkit grep -vnrp '#LTR' ",paste(out_dir,"/repeatmodeler/",label,"/",label,"_RepeatModeler.db-families.fa",sep=""),
    #           " | seqkit grep -vnrp '#Satellite'",
    #           " > ",paste(out_dir,"/repeatmodeler/",label,"/filtered_replib.fa",sep=""))
    # print(cmd);system(cmd,wait=TRUE)
    system(paste("rm",paste(out_dir,"/repeatmodeler/",label,"/genome.fa",sep=""),sep=" "))
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# EDTA
if (4 %in% step){
  Stamp1=paste(out_dir,"/edta/",label,"/",label,"_edta.finished",sep="")
  Stamp2=paste(save_dir,"/edta/",label,"/",label,"_edta.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 4: edta FINISHED")
  }else{
    print("Step 4: edta START")
    if (!file.exists(paste(out_dir,"/edta/",sep=""))){
      system(paste("mkdir ",out_dir,"/edta/",sep=""))
    }
    edta(genome=genome,
         edta.sif=edta.sif,
         threads=threads,
         out_dir=paste(out_dir,"/edta/",label,sep=""))
  }
}

# MITE
# if (4 %in% step){
#   Stamp1=paste(out_dir,"/miteFinder/",label,"/",label,"_miteFinder.finished",sep="")
#   Stamp2=paste(save_dir,"/miteFinder/",label,"/",label,"_miteFinder.finished",sep="")
#   if (file.exists(Stamp1) | file.exists(Stamp2)){
#     print("Step 4: miteFinder FINISHED")
#   }else{
#     print("Step 4: miteFinder START")
#     if (!file.exists(paste(out_dir,"/miteFinder/",sep=""))){
#       system(paste("mkdir"," ",out_dir,"/miteFinder/",sep=""))
#     }
#     miteFinder(genome=genome,
#                pattern_scoring=pattern_scoring,
#                out_dir=paste(out_dir,"/miteFinder/",label,sep=""))
#   }
# }
# if (4 %in% step){
#   Stamp1=paste(out_dir,"/mite_hunter/",label,"/",label,"_mite_hunter.finished",sep="")
#   Stamp2=paste(save_dir,"/mite_hunter/",label,"/",label,"_mite_hunter.finished",sep="")
#   if (file.exists(Stamp1) | file.exists(Stamp2)){
#     print("Step 4: mite_hunter FINISHED")
#   }else{
#     print("Step 4: mite_hunter START")
#     if (!file.exists(paste(out_dir,"/mite_hunter/",sep=""))){
#       system(paste("mkdir ",out_dir,"/mite_hunter/",sep=""))
#     }
#     mite_hunter(genome=genome,
#                 out_dir=paste(out_dir,"/mite_hunter/",label,sep=""),
#                 out_basename=label,
#                 threads=threads)
#   }
# }
############################################
# Classify transposons and wrangle TE library of different origin
if (5 %in% step){
  Stamp1=paste(out_dir,"/DeepTE/",label,"/",label,"_DeepTE.finished",sep="")
  Stamp2=paste(save_dir,"/DeepTE/",label,"/",label,"_DeepTE.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 5: DeepTE FINISHED")
  }else{
    print("Step 5: DeepTE START")
    if (!file.exists(paste(out_dir,"/DeepTE/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/DeepTE/",sep=""))
    }
    
    if (!file.exists(paste(out_dir,"/DeepTE/mite",sep=""))){
      system(paste("mkdir"," ",out_dir,"/DeepTE/mite",sep=""))
    }
    DeepTE(in.fna=paste(save_dir,"/miteFinder/",label,"/MITE_sorted.fa",sep=""),
           out_dir=paste(out_dir,"/DeepTE/mite/",label,sep=""),
           TE.fam="ClassII", # none
           DeepTE.sp=DeepTE.sp,# M:Metazoans, F:Fungi, and O: Others.
           DeepTE.model=DeepTE.model)
    
    if (!file.exists(paste(out_dir,"/DeepTE/ltr",sep=""))){
      system(paste("mkdir"," ",out_dir,"/DeepTE/ltr",sep=""))
    }
    DeepTE(in.fna=paste(save_dir,"/ltr/",label,"/LTR.fa",sep=""),
           out_dir=paste(out_dir,"/DeepTE/ltr/",label,sep=""),
           TE.fam="LTR", # none
           DeepTE.sp=DeepTE.sp,# M:Metazoans, F:Fungi, and O: Others.
           DeepTE.model=DeepTE.model)
    
    if (!file.exists(paste(out_dir,"/DeepTE/rmodeler",sep=""))){
      system(paste("mkdir"," ",out_dir,"/DeepTE/rmodeler",sep=""))
    }
    DeepTE(in.fna=paste(save_dir,"/repeatmodeler/",label,"/",label,"_RepeatModeler.db-families.fa",sep=""),
           out_dir=paste(out_dir,"/DeepTE/rmodeler/",label,sep=""),
           TE.fam="none", # none
           DeepTE.sp=DeepTE.sp,# M:Metazoans, F:Fungi, and O: Others.
           DeepTE.model=DeepTE.model)
  }
}
tandem.gff3=paste(save_dir,"/trash/",label,"/TRASH.gff3",sep="")


# >repeatname#class/subclass
# In this format, the data will be processed (overlapping repeats are
# merged etc), alternative output (.ace or .gff) can be created and an
# overview .tbl file will be created. Classes that will be displayed in
# the .tbl file are 'SINE', 'LINE', 'LTR', 'DNA', 'Satellite', anything
# with 'RNA' in it, 'Simple_repeat', and 'Other' or 'Unknown' (the
#                                                              latter defaults when class is missing). Subclasses are plentiful. They
# are not all tabulated in the .tbl file or necessarily spelled
# identically as in the repeat files, so check the RepeatMasker.embl
# file for names that can be parsed into the .tbl file.

############################################
# Annotate repeats with repeatmasker
if (6 %in% step){
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
                                    save_dir,"/repeatmodeler/",label,"/filtered_replib.fa",
                                    sep=""), # space-list
                 Threads=threads)
    #system(paste("rm ",out_dir,"/repeat_elements/",label,"/genome.fa.cat.gz"))
    #system(paste("rm ",out_dir,"/repeat_elements/",label,"/genome.fa.masked"))
    #system(paste("rm ",out_dir,"/repeat_elements/",label,"/genome.fa.out.gff"))
  
    rep_elements(repeatmasker.out=paste(out_dir,"/repeat_elements/",label,"/genome.fa.out",sep=""),
                 # Clarify source of repeat elements
                 TRASH.gff3=paste(save_dir,"/trash/",label,"/TRASH.gff3",sep=""),
                 MITE_Hunter.replib=paste(save_dir,"/mite_hunter/",label,"/",label,"_MITE.fa",sep=""), # RepLib from MITE-Hunter
                 LTR_retriever.replib=paste(save_dir,"/ltr/",label,"/LTR.fa",sep=""), # RepLib from LTR_retriever
                 RepeatModeler.replib=paste(save_dir,"/repeatmodeler/",label,"/filtered_replib.fa",sep=""), # RepLib from RepeatModeler
                 all.replib=paste(out_dir,"/repeat_elements/",label,"/",label,"_replib.fa",sep=""),
                 out.gff3=paste(out_dir,"/repeat_elements/",label,"/",label,"_repeat_elements.gff3",sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}