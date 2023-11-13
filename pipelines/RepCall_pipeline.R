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
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# nr TE library
if (5 %in% step){
  Stamp1=paste(out_dir,"/seqNR/",label,"/",label,"_seqNR.finished",sep="")
  Stamp2=paste(save_dir,"/seqNR/",label,"/",label,"_seqNR.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 5: seqNR FINISHED")
  }else{
    print("Step 5: seqNR START")
    if (!file.exists(paste(out_dir,"/seqNR/",sep=""))){
      system(paste("mkdir ",out_dir,"/seqNR/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/seqNR/",label,sep=""))){
      system(paste("mkdir ",out_dir,"/seqNR/",label,sep=""))
    }
    
    repeatmodeler.lib=paste(save_dir,"/repeatmodeler/",label,"/",label,"_RepeatModeler.db-families.fa",sep="")
    ltr.lib=paste(save_dir,"/ltr/",label,"/LTR.fa",sep="")
    edta.lib=paste(save_dir,"/edta/",label,"/genome.fa.mod.EDTA.TElib.fa",sep="")
    
    repeatmodeler.oldIDs=system(paste("grep '>' ",repeatmodeler.lib," | sed 's/ .*$//' | sed 's/>//'",sep=""),intern=TRUE)
    repeatmodeler.oldClass=sub("^.*#","",repeatmodeler.oldIDs)
    repeatmodeler.newIDs=paste("RM",as.character(1:length(repeatmodeler.oldIDs)),"#",repeatmodeler.oldClass,sep="")
    write.table(data.frame(repeatmodeler.oldIDs,repeatmodeler.newIDs),
                paste(out_dir,"/seqNR/",label,"/alias.txt",sep=""),
                sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
    cmd=paste("seqkit replace -p '^(\\S+)' -r '{kv}$2'",
              "-k",paste(out_dir,"/seqNR/",label,"/alias.txt",sep=""),
              repeatmodeler.lib,
              ">",
              paste(out_dir,"/seqNR/",label,"/lib.fna",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    
    ltr.oldIDs=system(paste("grep '>' ",ltr.lib," | sed 's/>//'",sep=""),intern=TRUE)
    ltr.oldClass=sub("^.*#","",ltr.oldIDs)
    ltr.newIDs=paste("ltr",as.character(1:length(ltr.oldIDs)),"#",ltr.oldClass,sep="")
    write.table(data.frame(ltr.oldIDs,ltr.newIDs),
                paste(out_dir,"/seqNR/",label,"/alias.txt",sep=""),
                sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
    cmd=paste("seqkit replace -p '^(\\S+)' -r '{kv}$2'",
              "-k",paste(out_dir,"/seqNR/",label,"/alias.txt",sep=""),
              ltr.lib,
              ">>",
              paste(out_dir,"/seqNR/",label,"/lib.fna",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    
    edta.oldIDs=system(paste("grep '>' ",edta.lib," | sed 's/>//'",sep=""),intern=TRUE)
    edta.oldClass=sub("^.*#","",edta.oldIDs)
    edta.newIDs=paste("edta",as.character(1:length(edta.oldIDs)),"#",edta.oldClass,sep="")
    write.table(data.frame(edta.oldIDs,edta.newIDs),
                paste(out_dir,"/seqNR/",label,"/alias.txt",sep=""),
                sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
    cmd=paste("seqkit replace -p '^(\\S+)' -r '{kv}$2'",
              "-k",paste(out_dir,"/seqNR/",label,"/alias.txt",sep=""),
              edta.lib,
              ">>",
              paste(out_dir,"/seqNR/",label,"/lib.fna",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    
    system(paste("rm",paste(out_dir,"/seqNR/",label,"/alias.txt",sep=""),sep=" "))
    sum.lib.fna=data.frame(oldID=c(repeatmodeler.oldIDs,ltr.oldIDs,edta.oldIDs),
                           newID=c(repeatmodeler.newIDs,ltr.newIDs,edta.newIDs),
                           oldClass=c(repeatmodeler.oldClass,ltr.oldClass,edta.oldClass),
                           source=c(rep("repeatmodeler",length(repeatmodeler.oldIDs)),
                                    rep("ltr",length(ltr.oldIDs)),
                                    rep("edta",length(edta.oldIDs))))
    write.table(sum.lib.fna,paste(out_dir,"/seqNR/",label,"/lib.fna_info.tsv",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
    
    seqNR(in.fasta=paste(out_dir,"/seqNR/",label,"/lib.fna",sep=""),
          out_dir=paste(out_dir,"/seqNR/",label,sep=""),
          Identity=0.8, # [0.0,1.0]
          cov_mode=0, # 0: alignment covers ${coverage} of target and of query
                   # 1: alignment covers ${coverage} of target
                   # 2: alignment covers ${coverage} of query
                   # 3: target is of ${coverage} query length
          coverage=0.8, # [0.0,1.0]
          threads=threads)
    system(paste("mv",
                 paste(out_dir,"/seqNR/",label,"/rep.fasta",sep=""),
                 paste(out_dir,"/seqNR/",label,"/nrlib.fna",sep=""),
                 sep=" "))
    repIDs=system(paste("grep '>' ",out_dir,"/seqNR/",label,"/nrlib.fna | sed 's/>//' | sed 's/ .*$//'",sep=""),intern=TRUE)
    sum.nrlib.fna=sum.lib.fna[sum.lib.fna[,"newID"] %in% repIDs,]
    write.table(sum.nrlib.fna,paste(out_dir,"/seqNR/",label,"/nrlib.fna_info.tsv",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
    
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# TE classification
if (6 %in% step){
  Stamp1=paste(out_dir,"/classify/",label,"/",label,"_classify.finished",sep="")
  Stamp2=paste(save_dir,"/classify/",label,"/",label,"_classify.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 6: classify FINISHED")
  }else{
    print("Step 6: classify START")
    rfsb=function(in.fna=paste(out_dir,"/seqNR/",label,"/nrlib.fna",sep=""),
                  DeepTE.sp=DeepTE.sp, # M:Metazoans, F:Fungi, and O: Others.
                  DeepTE.model=DeepTE.model,
                  out_dir=paste(out_dir,"/classify/",label,sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# TE annotation
if (7 %in% step){
  Stamp1=paste(out_dir,"/annotation/",label,"/",label,"_annotation.finished",sep="")
  Stamp2=paste(save_dir,"/annotation/",label,"/",label,"_annotation.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 7: annotation FINISHED")
  }else{
    print("Step 7: annotation START")
    if (!file.exists(paste(out_dir,"/annotation/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/annotation/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/annotation/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/annotation/",label,sep=""))
    }
    # source("/flash/BourguignonU/Cong/termite_genome_annotation/shell/conf/Aaca_conf.R")
    nrlib=paste(save_dir,"/seqNR/",label,"/nrlib.fna",sep="")
    repeatmasker(fna=genome,# Fasta file of genome.
                 out_dir=paste(out_dir,"/annotation/",label,sep=""),
                 RepeatLib.fa=paste(save_dir,"/seqNR/",label,"/nrlib.fna",sep=""), # space-list
                 Threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# Final results
# if (8 %in% step){
#   Stamp1=paste(out_dir,"/repeat_element/",label,"/",label,"_repeat_element.finished",sep="")
#   Stamp2=paste(save_dir,"/repeat_element/",label,"/",label,"_repeat_element.finished",sep="")
#   if (file.exists(Stamp1) | file.exists(Stamp2)){
#     print("Step 8: repeat_element FINISHED")
#   }else{
#     print("Step 8: repeat_element START")
#     if (!file.exists(paste(out_dir,"/repeat_element/",sep=""))){
#       system(paste("mkdir"," ",out_dir,"/repeat_element/",sep=""))
#     }
#     if (!file.exists(paste(out_dir,"/repeat_element/",label,sep=""))){
#       system(paste("mkdir"," ",out_dir,"/repeat_element/",label,sep=""))
#     }
#     
#     tandem_repeat=paste(save_dir,"/trash/",label,"/TRASH.gff3",sep="")
#     system("cp",tandem_repeat,
#            paste(out_dir,"/repeat_element/",label,"/",label,"_tandem.gff3",sep=""),
#            sep=" ")
#     replib=paste(save_dir,"/seqNR/",label,"/nrlib.fna",sep="")
#     system("cp",replib,
#            paste(out_dir,"/repeat_element/",label,"/",label,"_interspersedLib.fna",sep=""),
#            sep=" ")
#     
#     classes=read.table(paste(save_dir,"/annotation/",label,"/",label,"_classification.tsv",sep=""),
#                        sep="\t",header=TRUE,quote="",comment.char="")
#     rownames(classes)=sub("#.*$","",classes[,"newID"])
#     interspersed=read.table(paste(save_dir,"/annotation/",label,"/genome.fa.out.gff",sep=""),
#                             sep="\t",header=FALSE,quote="")
#     interspersed[,9]=sapply(1:nrow(interspersed),
#                             function(i){
#                               target=sub("Target \"Motif:","",interspersed[i,9])
#                               target=sub("\".*$","",target)
#                               source=classes[target,"source"]
#                               class=classes[target,"Class"]
#                               order=classes[target,"Order"]
#                               res=NA
#                               if (!is.na(source)){
#                                 res=paste("Target=",target,";",
#                                           "Source=",source,";",
#                                           "Class=",class,";",
#                                           "Order=",order,";",
#                                           sep="")
#                               }
#                               return(res)
#                             })
#     write.table(interspersed[!is.na(interspersed[,9]),],
#                 paste(out_dir,"/repeat_element/",label,"/",label,"_interspersed.gff3",sep=""),
#                 sep="\t",col.names = FALSE,row.names = FALSE,quote=FALSE)
#   }
# }