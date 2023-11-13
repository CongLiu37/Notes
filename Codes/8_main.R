# More genome annotation (repeat elements and non-coding RNA)

######
# Repeat elements
######
# TRASH: tandem repeats
trash=function(genome=genome,
               out_dir=out_dir,
               threads=threads){
  threads=as.character(threads)
  wd=getwd()
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  setwd(out_dir)
  
  if (!file.exists("all.repeats.from.genome.fa.csv")){
    cmd=paste("seqkit","seq","-u",
              "-j",threads,
              genome,">",
              paste(out_dir,"/genome.fa",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    
    cmd=paste("TRASH_run.sh",
              "genome.fa",
              "--o",out_dir,
              "--par",threads,
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  system("rm -r genome.fa_out")
  d=read.csv("all.repeats.from.genome.fa.csv",header=TRUE)
  library(stringr)
  d[,"block"]=rep(NA,nrow(d))
  d[1,"block"]=1
  for (i in 2:nrow(d)){
    if (d[i,"seq.name"]==d[i-1,"seq.name"] &
        d[i,"strand"]==d[i-1,"strand"] &
        d[i,"start"]==d[i-1,"end"]+1){
      d[i,"block"]=d[i-1,"block"]
    }else{
      d[i,"block"]=d[i-1,"block"]+1
    }
  }
  d[,"block"]=str_pad(d[,"block"],8,side="left",pad="0")
  d[,"ID"]=rep(NA,nrow(d))
  d[1,"ID"]=1
  for (i in 2:nrow(d)){
    if (d[i,"block"]!=d[i-1,"block"]){
      d[i,"ID"]=1
    }else{
      d[i,"ID"]=1+d[i-1,"ID"]
    }
  }
  type=tapply(d[,"width"],d[,"block"],
         function(x){
           m=mean(x);n=length(x)
           if (m>60){return("Satellite")}
           if (m>10 & m<=60){return("miniSatellite")}
           if (m<=10){return("microSatellite")}
         })
  d[,"type"]=type[d[,"block"]]
  d[,"block"]=as.character(d[,"block"])
  d[,"ID"]=as.character(d[,"ID"])
  
  gff=data.frame(d[,"seq.name"],
                 rep("TRASH",nrow(d)),
                 rep("Tandem_repeat",nrow(d)),
                 d[,"start"],
                 d[,"end"],
                 rep(".",nrow(d)),
                 d[,"strand"],
                 rep(".",nrow(d)),
                 paste("ID=",d[,"type"],d[,"block"],".",d[,"ID"],";",
                       "Tandem_block=",d[,"type"],d[,"block"],";",
                       "Type=",d[,"type"],
                       sep=""))
  write.table(gff,"TRASH.gff3",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

  system("rm genome.fa")
  setwd(wd)
}

# MITE-Hunter
# Dependencies: MITE-hunter, seqkit
mite_hunter=function(genome=genome,
                     out_dir=out_dir,
                     out_basename=out_basename,
                     threads=threads){
  threads=as.character(threads)
  wd=getwd()
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  setwd(out_dir)
  
  cmd=paste("seqkit","seq","-u",
            "-j",threads,
            genome,">",
            paste(out_dir,"/genome.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("MITE_Hunter_manager.pl",
            "-i",paste(out_dir,"/genome.fa",sep=""),
            "-g",out_basename,
            "-n",threads,
            "-S 12345678",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("cat",
            paste(out_basename,"_Step8*",sep=""),
            ">",
            paste(out_basename,"_MITE.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm"," ",out_basename,".*",sep=""))
  system(paste("rm"," ",out_basename,"_raw*",sep=""))
  system(paste("rm"," ",out_basename,"_[0-9]*",sep=""))
  system(paste("rm"," ",out_basename,"_Step*",sep=""))
  system("rm genome.*")
  system("rm formatdb.log")
  setwd(wd)
}

# LTR: long terminal repeats
# Dependencies: LTR_finder, LTR_Harvest (gt), LTR_retriever, seqkit
ltr=function(fna=fna,
             out_dir=out_dir,
             out_prefix=out_prefix,
             threads=threads){
  threads=as.character(threads)
  wd=getwd()
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  setwd(out_dir)
  
  cmd=paste("seqkit","seq","-u",
            "-j",threads,
            genome,">",
            paste(out_dir,"/genome.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  fna=paste(out_dir,"/genome.fa",sep="")
  
  if (!file.exists(paste(fna,".finder.combine.scn",sep=""))){
    cmd=paste("LTR_FINDER_parallel",
              "-seq",fna,
              "-harvest_out",
              "-threads",threads,
              "-size 1000000 -time 300",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  if (!file.exists(paste(out_prefix,".LTRharvest",sep=""))){
    cmd=paste("gt suffixerator",
              "-db",fna,
              "-indexname",paste(out_prefix,".suffixerator",sep=""),
              "-tis -suf -lcp -des -ssp -sds -dna",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("gt ltrharvest",
              "-index",paste(out_prefix,".suffixerator",sep=""),
              "-minlenltr 100 -maxlenltr 7000 -mintsd 4",
              "-maxtsd 6 -motif TGCA -motifmis 1 -similar 85",
              "-vic 10 -seed 20 -seqids yes",
              ">",paste(out_prefix,".LTRharvest",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("cat",
            paste(fna,".finder.combine.scn",sep=""),
            paste(out_prefix,".LTRharvest",sep=""),
            "> rawLTR.scn")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("LTR_retriever",
            "-genome",fna,
            "-inharvest","rawLTR.scn",
            "-threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("mv genome.fa.LTR.gff3 LTR.gff3")
  system("mv genome.fa.LTRlib.fa LTR.fa")
  
  system(paste("rm",fna,sep=" "),wait=TRUE)
  system(paste("rm","*.suffixerator.*",sep=" "))
  system("rm genome.fa")
  system("rm rawLTR.scn")
  setwd(wd)
}

# RepeatModeler
# Dependencies: RepeatModeler, Singularity
repeatmodeler=function(fna=fna,
                       out_dir=out_dir,
                       out_prefix=out_prefix,
                       Threads=Threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  pwd_begin=getwd();setwd(out_dir)
  Threads=as.character(Threads)
  
  system(paste("cp",fna,out_dir,sep=" "),wait=TRUE)
  fna_name=unlist(strsplit(fna,"/"));fna_name=fna_name[length(fna_name)]
  fna=paste(out_dir,"/",fna_name,sep="")
  #####################################################################
  # RepeatModeler & RepeatMasker installed in Singularity container
  path="singularity run /home/c/c-liu/Softwares/dfam-tetools-latest.sif"
  #####################################################################
  
  # RepeatModeler: de novo repeat library
  cmd=paste(path,
            "BuildDatabase",
            "-name",paste(out_prefix,"_RepeatModeler.db",sep=""),
            "-engine","ncbi",
            fna,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste(path,
            "RepeatModeler",
            "-database",paste(out_prefix,"_RepeatModeler.db",sep=""),
            "-engine","ncbi",
            "-pa",Threads,
            "-LTRStruct",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",fna,sep=" "),wait=TRUE)
  setwd(pwd_begin)
  return(paste(out_dir,"/",out_prefix,"_RepeatModeler.db-families.fa",sep=""))
}

# EDTA
# Dependencies: EDTA, singularity
edta=function(genome=genome,
              edta.sif="/home/c/c-liu/Softwares/EDTA.sif",
              threads=threads,
              out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("seqkit","seq","-u",
            "-j",threads,
            genome,">",
            paste(out_dir,"/genome.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("singularity exec",
            edta.sif,"EDTA.pl",
            "--genome genome.fa",
            "--species others",
            "--step all",
            "--overwrite 0", # enable restart from previous run
            "--sensitive 0", # do not invoke repeatmodeler
            "--anno 0",
            "--threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("rm genome.fa")
  setwd(wd)
}

# MMseq2: remove redundancy
# Dependencies: mmseq2
seqNR=function(in.fasta=in.fasta,
               out_dir=out_dir,
               Identity=Identity, # [0.0,1.0]
               cov_mode=cov_mode, # 0: alignment covers ${coverage} of target and of query
               # 1: alignment covers ${coverage} of target
               # 2: alignment covers ${coverage} of query
               # 3: target is of ${coverage} query length
               # quert as representative
               coverage=coverage, # [0.0,1.0]
               threads=threads){
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("mmseqs createdb",in.fasta,"sequenceDB",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cov_mode=as.character(cov_mode)
  coverage=as.character(coverage)
  Identity=as.character(Identity)
  cmd=paste("mmseqs cluster sequenceDB clusterDB", out_dir,
            "--cov-mode",cov_mode,
            "-c",coverage,
            "--min-seq-id",Identity,
            "--threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="mmseqs createsubdb clusterDB sequenceDB clusterDB_rep"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="mmseqs convert2fasta clusterDB_rep rep.fasta"
  print(cmd);system(cmd,wait=TRUE)
  
  # cmd="mmseqs createtsv sequenceDB sequenceDB clusterDB_rep clusters.tsv"
  setwd(wd)
}

# TE classification by transposon_classifier_RFSB & DeepTE 
# Dependencies: transposon_classifier_RFSB, DeepTE, seqkit
rfsb=function(in.fna=in.fna,
              DeepTE.sp="M", # M:Metazoans, F:Fungi, and O: Others.
              DeepTE.model="/bucket/BourguignonU/Cong/public_db/DeepTE/Metazoans_model",
              out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("transposon_classifier_RFSB",
            "-mode classify",
            "-fastaFile",in.fna,
            "-outputPredictionFile classification.txt",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  res=readLines("classification.txt")
  res=res[res!="" & !grepl("^#",res)]
  res=sub("^>","",res)
  TEs=res[seq(1,length(res),2)]
  class=res[seq(2,length(res),2)]
  class=sapply(class,function(i){ return(unlist(strsplit(i," "))[2]) })
  class=unname(class)
  
  res=data.frame(TE=TEs,rfsbClass=class)
  res[,c("class","order","superfamily")]=
    t(sapply(res[,"rfsbClass"],
             function(i){
               a=c(NA,NA,NA)
               if (i=="1"){a=c("Retrotransposons",NA,NA)}
               if (i=="1/1"){a=c("Retrotransposons","LTR",NA)}
               if (i=="1/1/1"){a=c("Retrotransposons","LTR","Copia")}
               if (i=="1/1/2"){a=c("Retrotransposons","LTR","Gypsy")}
               if (i=="1/1/3"){a=c("Retrotransposons","LTR","ERV")}
               if (i=="1/2"){a=c("Retrotransposons","nLTR",NA)}
               if (i=="1/2/1"){a=c("Retrotransposons","LINE",NA)}
               if (i=="1/2/2"){a=c("Retrotransposons","SINE",NA)}
               if (i=="2"){a=c("DNAtransposons",NA,NA)}
               if (i=="2/1"){a=c("DNAtransposons","TIR",NA)}
               if (i=="2/1/1"){a=c("DNAtransposons","TIR","Tc1-Mariner")}
               if (i=="2/1/2"){a=c("DNAtransposons","TIR","hAT")}
               if (i=="2/1/3"){a=c("DNAtransposons","TIR","CMC")}
               if (i=="2/1/4"){a=c("DNAtransposons","TIR","Sola")}
               if (i=="2/1/5"){a=c("DNAtransposons","TIR","Zator")}
               if (i=="2/1/6"){a=c("DNAtransposons","TIR","Novosib")}
               if (i=="2/2"){a=c("DNAtransposons","Helitron","Helitron")}
               if (i=="2/3"){a=c("DNAtransposons","MITE","MITE")}
               return(a)
             }))
  rownames(res)=res[,"TE"]
  
  system("mkdir LINE")
  system("mkdir SINE")
  LINE=res[res[,"order"]=="LINE",];SINE=res[res[,"order"]=="SINE",]
  writeLines(paste(LINE[,"TE"],sep=""),"LINE/lst")
  writeLines(paste(SINE[,"TE"],sep=""),"SINE/lst")
  cmd=paste("seqkit grep","-f LINE/lst",in.fna,"> LINE/LINE.fna",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("seqkit grep","-f SINE/lst",in.fna,"> SINE/SINE.fna",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("DeepTE.py",
            "-d LINE",
            "-o LINE",
            "-i LINE/LINE.fna",
            "-sp",DeepTE.sp,
            "-m_dir",DeepTE.model,
            "-fam LINE",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("DeepTE.py",
            "-d SINE",
            "-o SINE",
            "-i SINE/SINE.fna",
            "-sp",DeepTE.sp,
            "-m_dir",DeepTE.model,
            "-fam SINE",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  LINE=read.table("LINE/opt_DeepTE.txt",sep="\t",header=FALSE,quote="")
  LINE[,2]=sub("^ClassI_nLTR_LINE_","",LINE[,2])
  LINE[,2][LINE[,2]=="ClassI_nLTR_LINE"]=NA
  rownames(LINE)=LINE[,1]
  SINE=read.table("SINE/opt_DeepTE.txt",sep="\t",header=FALSE,quote="")
  SINE[,2]=sub("^ClassI_nLTR_SINE_","",SINE[,2])
  SINE[,2][SINE[,2]=="ClassI_nLTR_SINE"]=NA
  rownames(SINE)=SINE[,1]
  
  tmp1=res[res[,"order"]=="LINE",]
  tmp1[,"superfamily"]=sapply(tmp1[,"TE"],
                              function(te){return(LINE[te,2])})
  tmp2=res[res[,"order"]=="SINE",]
  tmp2[,"superfamily"]=sapply(tmp2[,"TE"],
                              function(te){return(SINE[te,2])})
  tmp3=res[res[,"order"]!="LINE" & res[,"order"]!="SINE",]
  
  d=rbind(tmp1,tmp2,tmp3)
  d=d[order(d[,"rfsbClass"]),]
  write.table(d,"TEclassification.tsv",sep="\t",row.names=FALSE,quote=FALSE)
  
  oldID=d[,"TE"]
  newID=sapply(1:nrow(d),
               function(i){
                 te=d[i,"TE"]
                 class=d[i,"class"];if (is.na(class)){class="Unknown"}
                 order=d[i,"order"];if (is.na(order)){order="Unknown"}
                 superfamily=d[i,"superfamily"];if (is.na(superfamily)){superfamily="Unknown"}
                 newName=NA
                 if (class=="Retrotransposons"){
                   newName=paste(te,"#",order,"/",superfamily,sep="")
                 }
                 if (class=="DNAtransposons"){
                   newName=paste(te,"#","DNA","/",superfamily,sep="")
                 }
                 return(newName)
               })
  write.table(data.frame(oldID,newID),"alias.tsv",
              sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  cmd=paste("seqkit replace -p '^(\\S+)' -r '{kv}'",
            "-k alias.tsv",in.fna,"> TEclassification.fna",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# RepeatMasker: rush jobs
repeatmaskerQQ=function(fna=fna,# Fasta file of genome.
                      out_dir=out_dir,
                      RepeatLib.fa=RepeatLib.fa, # space-list
                      Threads=Threads){
  Threads=as.character(Threads)
  wd=getwd()
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  setwd(out_dir)
  
  cmd=paste("seqkit","seq","-u",
            "-j",Threads,
            fna,">",
            paste(out_dir,"/genome.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  fna=paste(out_dir,"/genome.fa",sep="")
  
  cmd=paste("cat",RepeatLib.fa,"> repeat.fa",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  RepeatLib.fa=paste(out_dir,"/repeat.fa",sep="")
  
  #####################################################################
  # RepeatModeler & RepeatMasker installed in Singularity container
  path="singularity run /home/c/c-liu/Softwares/dfam-tetools-latest.sif"
  #####################################################################
  
  # RepeatMasker: Mask repeats found by RepeatModeler
  cmd=paste(path,
            "RepeatMasker",
            "-xsmall", # soft masking
            "-lib",RepeatLib.fa,
            "-pa",Threads,
            "-gff",
            "-nolow", # Does not mask low_complexity DNA or simple repeats
            "-norna", # Does not mask small RNA (pseudo) genes
            "-no_is", # Skips bacterial insertion element check
            "-qq", #Rush job; about 10% less sensitive, 4->10 times faster than default (quick searches are fine under most circumstances) repeat options
            fna,
            "-dir",out_dir)
  print(cmd);system(cmd,wait=TRUE)
  
  # # RepeatMasker: Mask repeats in RepBase
  # cmd=paste(path,
  #           "RepeatMasker",
  #           "-xsmall", # soft masking
  #           "-pa",Threads,
  #           "-gff",
  #           paste(fna,".masked",sep=""),
  #           "-dir",out_dir,
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  
  # Remove temporaries
  system(paste("rm",fna,sep=" "),wait=TRUE)
  system(paste("rm",RepeatLib.fa,sep=" "),wait=TRUE)
  setwd(wd)
  
  # masked genome
  return(paste(out_dir,"/",fna,".masked.masked",sep=""))
}

# Final repeat annotation: repeatmasker out to gff3
# Dependencies: parallel (R)
rep_elements=function(repeatmasker.gff2=repeatmasker.gff2, # interspersed repeats
                      trash.gff3=trash.gff3, # tandem repeats
                      TEclassification.tsv=TEclassification.tsv, # from rfsb
                      out_prefix=out_prefix,
                      threads=threads){
  cmd=paste("cp",trash.gff3,
            paste(out_prefix,"_tandem.gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  TEclassification=read.table(TEclassification.tsv,header=TRUE,sep="\t",quote="")
  rownames(TEclassification)=TEclassification[,"TE"]
  
  interspersed=read.table(repeatmasker.gff2,sep="\t",header=FALSE,quote="")
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  clusterExport(clus,list("interspersed","TEclassification"),envir=environment())
  interspersed[,9]=parSapply(clus,
                             1:nrow(interspersed),
                              function(i){
                                target=sub("Target \"Motif:","",interspersed[i,9])
                                target=sub("\".*$","",target)
                                class=TEclassification[target,"class"]
                                order=TEclassification[target,"order"]
                                superfamily=TEclassification[target,"superfamily"]
                                res=paste("Target=",target,";",
                                          "Class=",class,";",
                                          "Order=",order,";",
                                          "Superfamily=",superfamily,";",
                                          sep="")
                                return(res)
                              })
  write.table(interspersed,
              paste(out_prefix,"_interspersed.gff3",sep=""),
              sep="\t",col.names = FALSE,row.names = FALSE,quote=FALSE)
}

######
# non-coding RNA
######
# tRNAscan-SE
# Dependencies: tRNAscan-SE, biocode, seqkit
tRNAscan=function(genome=genome,
                  mode="eukaryotic", # eukaryotic/bacterial/archaeal/mitochondrial_mammal/mitochondrial_vert/other_organellar/general
                  out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("seqkit","seq","-u",
            genome,">",
            paste(out_dir,"/genome.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  genome=paste(out_dir,"/genome.fa",sep="")
  
  if (mode=="eukaryotic"){mode="-E"}
  if (mode=="bacterial"){mode="-B"}
  if (mode=="archaeal"){mode="-A"}
  if (mode=="mitochondrial_mammal"){mode="-M mammal"}
  if (mode=="mitochondrial_vert"){mode="-M vert"}
  if (mode=="other_organellar"){mode="-O"}
  if (mode=="general"){mode="G"}
  
  cmd=paste("tRNAscan-SE",
            mode,
            "-o ./final_results",
            "-f ./secondary_structures",
            "-m ./summary",
            genome,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="convert_tRNAScanSE_to_gff3.pl -g --input=final_results > tRNA.gff3" # from biocode
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",genome,sep=" "))
  
  setwd(out_dir)
}

# miRNAture: Identify microRNA
miRNAture=function(genome=genome,
                   dataF=dataF, # pre-calculated data directory from miRNAture
                   species=species, # Genus_species
                   out_dir=out_dir,
                   threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  cmd=paste("cp -r",dataF,"./dataF",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",genome,"./genome.fa",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("miRNAture",
            "-stage complete",
            "-dataF","dataF",
            "-speG","./genome.fa",
            "-speN",species,
            "-speT mrna",
            "-w",out_dir,
            "-m hmm,rfam,mirbase,infernal,final",
            "-pe 0",
            "-nbitscore_cut 1.0",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("mv ./Final_miRNA_evaluation/miRNA_annotation_mrna_accepted_conf.gff3 miRNAture.gff3")
  system("cat ./Final_miRNA_evaluation/Fasta/Drosophila_mirna_miRNAs_high_confidence.fasta Final_miRNA_evaluation/Fasta/Drosophila_mirna_miRNAs_medium_confidence.fasta > miRNA.fna")
  system("rm -r ./dataF")
  system("rm ./genome.fa")
  system("rm -r LOGS")
  system("rm -r miRNA_prediction")
  system("rm miRNAture_configuration_mrna.yaml")
  system("rm -r miRNA_validation")
  system("rm -r TemporalFiles")
  system("rm -r Final_miRNA_evaluation")
  setwd(wd)
}

# miRanda: microRNA target gene identification
# Dependencies: miRanda, bedtools
miranda=function(gff=gff, # three_prime_UTR feature required
                 genome=genome,
                 miRNA.fna=miRNA.fna,
                 threads=threads,
                 out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  cmd=paste("grep 'three_prime_UTR' ",gff,
            " | awk -F '\t' -v OFS='\t' '{print $1,$4-1,$5,$9,200,$7}' | sed 's/ID=.*;Parent=//' | sed 's/;//' > 3_UTR.bed",
            sep="")
  print(cmd);system(cmd,wait=TRUE)
  system(paste("cp",genome,"./genome.fa",sep=" "))
  cmd=paste("bedtools getfasta -s -name",
            "-fi","./genome.fa",
            "-bed 3_UTR.bed",
            " | sed 's/::.*//' > 3_UTR.fa")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("seqkit split2",
            "3_UTR.fa",
            "--by-part",as.character(as.numeric(threads)-1),
            "--out-dir","tmp",
            "--threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("cp",miRNA.fna,"./miRNA.fna",sep=" "))
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  parSapply(clus,
            1:(as.numeric(threads)-1),
            function(i){
              utr=system("ls tmp/*.fa",intern=TRUE)[i]
              cmd=paste("miranda","miRNA.fna",utr,
                        "-sc 150","-en -20",
                        "| grep '>>' | sed 's/>>//' | sort -k 5 -n -r",
                        "| awk -F '\t' -v OFS='\t' '{print $1,$2}' | sort | uniq -u >",
                        paste("tmp/",as.character(i),".tsv",sep=""),
                        sep=" ")
              system(cmd,wait=TRUE)
            })
  stopCluster(clus)
  
  o=system("ls tmp/*.tsv",intern=TRUE)
  cmd=paste("cat",paste(o,collapse=" "),"> miRanda.tsv")
  system(cmd)
  
  system("rm ./miRNA.fna")
  system("rm ./3_UTR.fa")
  system("rm ./3_UTR.bed")
  system("rm ./genome.fa")
  system("rm ./genome.fa.fai")
  system("rm -r tmp")
  setwd(wd)
}

# Infernal: non-coding RNA gene identification
# cmscan treats uppercase and lowercase nucleotides identically
# Dependencies: Infernal, Rfam database, stringr (R)
infernal=function(genome=genome,
                  domain=domain, # eukarya/archaea/bacteria Which domain shall be retained?
                  Rfam.clanin=Rfam.clanin,
                  Rfam.cm=Rfam.cm,
                  out_dir=out_dir,
                  threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)

  cmd=paste("cmscan",
            "--cut_ga","--rfam","--nohmmonly",
            "--tblout Infernal.tblout",
            "--fmt 2",
            "--clanin",Rfam.clanin,
            "--cpu",threads,
            Rfam.cm,genome,
            "> Infernal.cmscan")
  print(cmd);system(cmd,wait=TRUE)
  
  domains=c("eukarya","archaea","bacteria")
  for (i in 1:3){if (grepl(domain,domains[i])){domains=domains[-i]}}
  cmd=paste("grep -v '=' Infernal.tblout | ",
            "grep -v 'tRNA' | ",
            "grep -v '",domains[1],"' | ",
            "grep -v '",domains[2],"' | ",
            "grep -v 'microRNA'",
            " > final.tblout",
            sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="awk 'NR>2 {if ($1!=\"#\") print $4\"\tInfernal\t\"$2\"\t\"$10\"\t\"$11\"\t\"$17\"\t\"$12\"\t\\.\tID=\"$2\"_\"$1\";Class=\"$2\";Rfam=\"$3\";Description=\"}' final.tblout  | awk -F '\t' -v OFS='\t' '{if ($7==\"-\") print $1,$2,$3,$5,$4,$6,$7,$8,$9; else print $0}' > Infernal.gff3"
  print(cmd);system(cmd,wait=TRUE)
  Des=system("sed 's/.*- //' final.tblout",intern=TRUE)
  Des=sapply(1:length(Des),
             function(i){
               if (grepl("#",Des[i])){return(NA)}else{return(Des[i])}
             })
  Des=Des[!is.na(Des)]
  Des=Des[-1]
  gff=readLines("Infernal.gff3")
  
  out=sapply(1:length(Des),
             function(i){
               return(paste(gff[i],Des[i],sep=""))
             })
  writeLines(out,"Infernal.gff3")
  setwd(wd)
}

# Integrate ncRNA results from tRNAscan-SE, miRNAture, infernal and miRanda
ncRNA=function(tRNA.gff3=tRNA.gff3,
               miRNAture.gff3=miRNAture.gff3,
               miRNA.fa=miRNA.fa,
               infernal.gff3=infernal.gff3,
               miRanda.tsv=miRanda.tsv,
               species=species,
               out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",tRNA.gff3,
            paste(species,"_tRNA.gff3",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",miRNAture.gff3,
            paste(species,"_miRNA.gff3",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",miRNA.fa,
            paste(species,"_miRNA.fa",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",infernal.gff3,
            paste(species,"_Rfam.gff3",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",miRanda.tsv,
            paste(species,"_miRNA_target.tsv",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}


# 6-frame translation (DNA to protein)
# hmmer build (protein)
# hmmer search (protein-protein)




# # mirdeep2: Identify microRNA from microRNA-seq
# # Dependencies: seqkit,mirdeep2
# mirdeep=function(genome=genome,
#                  fq=fq, # comma-list of RNA-seq fastq
#                  threads=threads,
#                  out_dir=out_dir){
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   wd=getwd();setwd(out_dir)
#   
#   fqs=unlist(strsplit(fq,","))
#   for (j in 1:length(fqs)){
#     i=fqs[j]
#     cmd=paste("seqkit fq2fa",
#               i,
#               "-j",threads,
#               "-o",paste(as.character(j),".fa",sep=""),
#               sep=" ")
#     print(cmd);system(cmd,wait=TRUE)
#     
#     system(paste("cat ",as.character(j),".fa"," >> reads.fa",sep=""))
#     system(paste("rm"," ",as.character(j),".fa",sep=""))
#   }
#   
#   cmd="collapse_reads_md.pl reads.fa seq > reads.col.fa"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("sed","'s/ .$//'",genome,">","genome.fa",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("bowtie-build","genome.fa","./bowtie_index",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("mapper.pl reads.col.fa -c -p bowtie_index -t map.arf -o",threads,sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("miRDeep2.pl reads.col.fa",genome,"map.arf none none none ‑v 2>report.log",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   setwd(wd)
# }

# # RepeatModeler (no LTR)
# # Dependencies: RepeatModeler, Singularity
# repeatmodeler_noLTR=function(fna=fna,
#                              out_dir=out_dir,
#                              out_prefix=out_prefix,
#                              Threads=Threads){
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   pwd_begin=getwd();setwd(out_dir)
#   Threads=as.character(Threads)
#   
#   system(paste("cp",fna,out_dir,sep=" "),wait=TRUE)
#   fna_name=unlist(strsplit(fna,"/"));fna_name=fna_name[length(fna_name)]
#   fna=paste(out_dir,"/",fna_name,sep="")
#   #####################################################################
#   # RepeatModeler & RepeatMasker installed in Singularity container
#   path="singularity run /home/c/c-liu/Softwares/dfam-tetools-latest.sif"
#   #####################################################################
#   
#   # RepeatModeler: de novo repeat library
#   cmd=paste(path,
#             "BuildDatabase",
#             "-name",paste(out_prefix,"_RepeatModeler.db",sep=""),
#             "-engine","ncbi",
#             fna,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste(path,
#             "RepeatModeler",
#             "-database",paste(out_prefix,"_RepeatModeler.db",sep=""),
#             "-engine","ncbi",
#             "-pa",Threads,
#             # "-LTRStruct",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   system(paste("rm",fna,sep=" "),wait=TRUE)
#   system("rm -r RM_*")
#   
#   cmd=paste("seqkit grep -vnrp '#LTR' ",label,"_RepeatModeler.db-families.fa | ",
#             "seqkit grep -vnrp '#[tr]RNA'",
#             " > noLTR_rep.lib.fa")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   setwd(pwd_begin)
#   return(paste(out_dir,"/noLTR_rep.lib.fa",sep=""))
# }

# # miteFinder
# # Dependencies: miteFinder, seqkit
# miteFinder=function(genome=genome,
#                     pattern_scoring="~/Softwares/miteFinder/profile/pattern_scoring.txt",
#                     out_dir=out_dir){
#   threads=as.character(threads)
#   wd=getwd()
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   setwd(out_dir)
#   
#   cmd=paste("seqkit","seq","-u",
#             "-j",threads,
#             genome,">",
#             paste(out_dir,"/genome.fa",sep=""),
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("miteFinder_linux_x64",
#             "-input genome.fa",
#             "-output MITE.fa",
#             "-pattern_scoring",pattern_scoring,
#             "-threshold 0.5",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="seqkit replace -p .+ -r 'Mite{nr}' MITE.fa > MITE_sorted.fa"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   system("rm genome.fa")
#   
#   setwd(wd)
# }

# # Generic Repeat Finder: TIR, TDRs, interspersed repeats, MITEs, and LTR
# grf=function(genome.fna=genome.fna,
#              out_dir=out_dir,
#              threads=threads){
#   # grf-main -i tests/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -o tir -c 0 --min_tr 10 -t 8
#   # grf-main -i tests/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -o tdr -c 2 --min_tr 10 -t 8
#   # grf-main -i tests/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -o mite -c 1 --min_tr 10 -t 8
#   # grf-intersperse -i tests/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa -o interspersed -t 8
# }

# # DeepTE: classify transposons
# DeepTE=function(in.fna=in.fna,
#                 out_dir=out_dir,
#                 TE.fam=TE.fam, # none
#                 # ClassI: the input sequence is ClassI TEs
#                 # ClassII: the input sequence is ClassII subclass1 TEs
#                 # LTR: the input sequence is LTR TEs
#                 # nLTR: the input sequence is nLTR TEs
#                 # LINE: the input sequence is LINE TEs
#                 # SINE: the input sequence is SINE TEs
#                 # Domain: the input sequence is Class II subclass1 TEs with specified super families
#                 sp="M",# M:Metazoans, F:Fungi, and O: Others.
#                 DeepTE.model="/bucket/BourguignonU/Cong/public_db/DeepTE/Metazoans_model"){
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   pwd_begin=getwd();setwd(out_dir)
#   
#   cmd=paste("DeepTE.py",
#             "-d",out_dir,
#             "-o",out_dir,
#             "-i",in.fna,
#             "-sp",sp,
#             "-m_dir",DeepTE.model,
#             sep=" ")
#   if (TE.fam!="none"){
#     cmd=paste(cmd,
#               "-fam",TE.fam,
#               sep=" ")
#   }
#   print(cmd);system(cmd,wait=TRUE)
#   
#   setwd(pwd_begin)
# }

# # TEsorter: classify transposons
# # Dependencies: TEsorter
# TEsorter=function(in.fna=in.fna,
#                   out_dir=out_dir,
#                   db=db, # gydb,rexdb,rexdb-plant,rexdb-metazoa,rexdb-pnas,rexdb-line,sine
#                   threads=threads){
#   threads=as.character(threads)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   pwd_begin=getwd();setwd(out_dir)
#   
#   cmd=paste("TEsorter",
#             in.fna,
#             "-db",db,
#             "-p",threads,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   setwd(pwd_begin)
# }
# 
# DeepTE=function(in.fna=in.fna,
#                 out_dir=out_dir,
#                 TE.fam=TE.fam, # none
#                 # ClassI: the input sequence is ClassI TEs
#                 # ClassII: the input sequence is ClassII subclass1 TEs
#                 # LTR: the input sequence is LTR TEs
#                 # nLTR: the input sequence is nLTR TEs
#                 # LINE: the input sequence is LINE TEs
#                 # SINE: the input sequence is SINE TEs
#                 # Domain: the input sequence is Class II subclass1 TEs with specified super families
#                 sp="M",# M:Metazoans, F:Fungi, and O: Others.
#                 DeepTE.model="/bucket/BourguignonU/Cong/public_db/DeepTE/Metazoans_model"){
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   pwd_begin=getwd();setwd(out_dir)
#   
#   cmd=paste("DeepTE.py",
#             "-d",out_dir,
#             "-o",out_dir,
#             "-i",in.fna,
#             "-sp",sp,
#             "-m_dir",DeepTE.model,
#             sep=" ")
#   if (TE.fam!="none"){
#     cmd=paste(cmd,
#               "-fam",TE.fam,
#               sep=" ")
#   }
#   print(cmd);system(cmd,wait=TRUE)
#   
#   setwd(pwd_begin)
# }