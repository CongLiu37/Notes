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
  Stamp1=paste(out_dir,"/TEsorter/",label,"/",label,"_TEsorter.finished",sep="")
  Stamp2=paste(save_dir,"/TEsorter/",label,"/",label,"_TEsorter.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 6: TEsorter FINISHED")
  }else{
    print("Step 6: TEsorter START")
    if (!file.exists(paste(out_dir,"/TEsorter/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/TEsorter/",sep=""))
    }
    TEsorter(in.fna=paste(save_dir,"/seqNR/",label,"/nrlib.fna",sep=""),
             out_dir=paste(out_dir,"/TEsorter/",label,"/",sep=""),
             db=TEsorted.db, # gydb,rexdb,rexdb-plant,rexdb-metazoa,rexdb-pnas,rexdb-line,sine
             threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
  
  Stamp1=paste(out_dir,"/DeepTE/",label,"/",label,"_DeepTE.finished",sep="")
  Stamp2=paste(save_dir,"/DeepTE/",label,"/",label,"_DeepTE.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 6: DeepTE FINISHED")
  }else{
    print("Step 6: DeepTE START")
    if (!file.exists(paste(out_dir,"/DeepTE/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/DeepTE/",sep=""))
    }
    DeepTE(in.fna=paste(save_dir,"/seqNR/",label,"/nrlib.fna",sep=""),
           out_dir=paste(out_dir,"/DeepTE/",label,"/",sep=""),
           TE.fam="none", # none
           # ClassI: the input sequence is ClassI TEs
           # ClassII: the input sequence is ClassII subclass1 TEs
           # LTR: the input sequence is LTR TEs
           # nLTR: the input sequence is nLTR TEs
           # LINE: the input sequence is LINE TEs
           # SINE: the input sequence is SINE TEs
           # Domain: the input sequence is Class II subclass1 TEs with specified super families
           sp=DeepTE.sp,# M:Metazoans, F:Fungi, and O: Others.
           DeepTE.model=DeepTE.model)
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
    
    sum.nrlib=read.table(paste(save_dir,"/seqNR/",label,"/nrlib.fna_info.tsv",sep=""),
                         sep="\t",header=TRUE,quote="",comment.char="")
    TEsorter.res=read.table(paste(save_dir,"/TEsorter/",label,"/nrlib.fna.",TEsorted.db,".cls.tsv",sep=""),
                            sep="\t",header=TRUE,quote="",comment.char="")
    colnames(TEsorter.res)=c("newID","Order","Superfamily","Clade","Complete","Strand","Domains")
    TEsorter.res=TEsorter.res[,c("newID","Order","Superfamily","Clade")]
    d=merge(sum.nrlib,TEsorter.res,by="newID",all=TRUE)
    d[is.na(d)]="unknown"
    d[,"Class"]=sapply(d[,"Order"],
                       function(order){
                         if (order %in% c("LTR","DIRS","Penelope","LINE","SINE")){
                           return("ClassI")
                         }else if (order %in% c("Maverick","Helitron")){
                           return("ClassII_subclass2")
                         }else if (order %in% c("TIR","Crypton")){
                           return("ClassII_subclass1")
                         }else{
                           return("unknown")
                         }
                       })
    d=d[,c("newID","oldID","oldClass","source","Class","Order","Superfamily","Clade")]
    DeepTE=read.table(paste(save_dir,"/DeepTE/",label,"/opt_DeepTE.txt",sep=""),
                      sep="\t",header=FALSE,quote="",comment.char="")
    colnames(DeepTE)=c("newID","classification")
    rownames(DeepTE)=DeepTE[,"newID"]
    classOrder=t(sapply(1:nrow(d),
                        function(i){
                          res=d[i,c("Class","Order")]
                          deepte=unlist(strsplit(DeepTE[d[i,"newID"],"classification"],"_"))
                          deepte=deepte[!grepl("MITE",deepte)]
                          if (deepte=="unknown"){
                            return(res)
                          }else{
                            if (res[1]=="unknown" & res[2]=="unknown"){
                              # Class
                              res[1]=deepte[1]
                              if (res[1]=="ClassII"){res[1]="ClassII_sub1"}
                              if (res[1]=="ClassIII"){res[1]="ClassII_sub2"}
                              # Order
                              res[2]=deepte[2]
                              if (is.na(res[2])){res[2]="unknown"}
                              if (res[2]=="LTR"){res[2]="LTR"}
                              if (res[2]=="nLTR"){res[2]=deepte[3]}
                              if (is.na(res[2])){res[2]="unknown"}
                              if (res[2]=="PLE"){res[2]="Penelope"}
                              if (res[2]=="DNA"){res[2]="TIR"}
                              if (res[2]=="Helitron"){res[2]="Helitron"}
                              if (is.na(res[2])){res[2]="unknown"}
                            }
                            return(res)
                          }
                        }))
    d[,"Class"]=unlist(classOrder[,"Class"])
    d[,"Order"]=unlist(classOrder[,"Order"])
    write.table(d,paste(out_dir,"/annotation/",label,"/label_classification.tsv",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
    repeatmasker(fna=genome,# Fasta file of genome.
                 out_dir=paste(out_dir,"/annotation/",label,sep=""),
                 RepeatLib.fa=paste(save_dir,"/seqNR/",label,"/nrlib.fna",sep=""), # space-list
                 Threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}

