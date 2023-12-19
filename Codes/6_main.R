# Genomic analysis of protein-coding genes: 
# genome collinearity, HGT, phylogeny, gene content, selective pressure, copy number evolution

# For peptide sets from NCBI
# Retain the longest protein isoform of each gene by header lines of faa.
# Proteins are considered as isoforms only when the header line declares isoform.
# e.g. these are the same gene as they declare isoform and share same name "pyruvate decarboxylase 2-like"
# >XP_024356342.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
# >XP_024356343.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
# these are different genes:
# >XP_024390255.1 40S ribosomal protein S12-like [Physcomitrella patens]
# >XP_024399722.1 40S ribosomal protein S12-like [Physcomitrella patens]
# Dependencies: seqkit
Isoform_filter=function(faa_in=faa_in,
                        faa_out=faa_out,
                        threads=threads){
  threads=as.character(threads)
  tmp=paste(faa_out,"_tmp",sep="")
  system(paste("mkdir",tmp,sep=" "))
  
  # Sequence header (everything except >) to length (tabular)
  cmd=paste("seqkit","fx2tab",
            "--threads",threads,
            "--length --name --header-line",
            faa_in,
            ">",
            paste(tmp,"/Header2Length",sep=""))
  system(cmd,wait=TRUE)
  
  Header2Length=read.table(paste(tmp,"/Header2Length",sep=""),sep="\t",header=FALSE,quote="",
                           col.names=c("Header","Length"))
  Header2Length[,"ID"]=sapply(Header2Length[,"Header"],
                              function(Header){
                                Header=unlist(strsplit(Header," "))
                                return(Header[1])
                              })

  # Select longest isoform
  Isoforms=Header2Length[grepl("[Ii]soform",Header2Length[,"Header"]),] # Proteins that declare isoform
  Isoforms[,"Header"]=sub("[Ii]soform [0-9, A-Z]*","",Isoforms[,"Header"])
  if (nrow(Isoforms)!=0){
  Isoforms[,"Header"]=sapply(1:nrow(Isoforms),
                             function(i){ # convert header to name
                               return(sub(Isoforms[i,"ID"],"",Isoforms[i,"Header"]))
                             })
  IDs_LongestIso=sapply(names(table(Isoforms[,"Header"])),
                        function(Header){
                          d=Isoforms[Isoforms[,"Header"]==Header,]
                          d=d[d[,"Length"]==max(d[,"Length"]),]
                          ID=sample(d[,"ID"],1)
                          return(ID)
                        }) # IDs of longest isoforms
  IDs_LongestIso=unname(IDs_LongestIso)
  }else{IDs_LongestIso=c()}
  IDs_noIso=Header2Length[!grepl("[Ii]soform",Header2Length[,"Header"]),][,"ID"] # Protein IDs that do not declare isoform
  
  target_IDs=c(IDs_noIso,IDs_LongestIso)
  target_IDs=data.frame(target_IDs)
  write.table(target_IDs,paste(tmp,"/target_IDs",sep=""),
              sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  # Extract fasta by IDs
  cmd=paste("seqkit","grep",
            "--threads",threads,
            "-f",paste(tmp,"/target_IDs",sep=""),
            faa_in,
            ">",
            faa_out,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm","-r",tmp,sep=" "),wait=TRUE)
}

#####
# genome collinearity
#####
# MCScanX: classify duplicated gene pairs
# Dependencies: blastp, MCScanX, stringr (R)
MCScanX_dupliClass=function(gff=gff,
                            blast=blast, # self2self blast of proteins; IDs same with gene ID in gff
                            out_dir=out_dir,
                            out_basename=out_basename){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  system(paste("cp",blast,".",sep=" "))
  library(stringr)
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3==\"gene\") print $1,$9,$4,$5}' ",gff," > ",out_basename,".gff",sep="")
  print(cmd);system(cmd,wait=TRUE)
  df=read.table(paste(out_basename,".gff",sep=""),sep="\t",header=FALSE,quote="")
  oldChr=names(table(df[,1]));newChr=1:length(oldChr);names(newChr)=oldChr
  IDs=sapply(1:nrow(df),function(i){return( paste("sp",as.character(newChr[df[i,1]]),sep="") )})
  df[,1]=IDs
  geneIDs=str_extract(df[,2],"ID=[^;]*;");geneIDs=sub(";$","",geneIDs);geneIDs=sub("ID=","",geneIDs)
  df[,2]=geneIDs
  write.table(df,paste(out_basename,".gff",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  
  cmd=paste("MCScanX"," ",out_dir,"/",out_basename,sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("duplicate_gene_classifier"," ",out_dir,"/",out_basename,sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  df=read.table(paste(out_basename,".gene_type",sep=""),sep="\t",header=FALSE,quote="")
  types=sapply(1:nrow(df),function(i){
                            code=df[i,2]
                            if (code==0){return("singleton")}
                            if (code==1){return("dispersed")}
                            if (code==2){return("proximal")}
                            if (code==3){return("tandem")}
                            if (code==4){return("segmental")}})
  df[,3]=types
  write.table(df,paste(out_basename,".gene_type",sep=""),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  
  setwd(wd)
}

# WGDI input files
# Dependencies: blastp
wgdi_input=function(genome.fna1=genome.fna1,gff1=gff1,proteins.faa1=proteins.faa1,cds.fna1=cds.fna1, 
                    genome.fna2=genome.fna2,gff2=gff2,proteins.faa2=proteins.faa2,cds.fna2=cds.fna2,
                    # No protein isoforms on proteins.faa[12] and cds.fna[12]
                    # Corresponding proetin & cds has same ID
                    threads=threads,
                    out_dir=out_dir,
                    out_basename=out_basename){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  self2self=genome.fna1==genome.fna2
  
  # Prepare fake gff and len file required by wgdi
  # Dependencies: https://github.com/xuzhougeng/myscripts/blob/master/comparative/generate_conf.py
  wgdi_input=function(genome.fna=genome.fna,
                      gff=gff,
                      out_prefix=out_prefix){
    cmd=paste("generate_conf.py",
              "-p",out_prefix,
              genome.fna,gff,
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  wgdi_input(genome.fna=genome.fna1,
             gff=gff1,
             out_prefix=out_basename)
  fake_gff1=paste(out_dir,"/",out_basename,".gff",sep="")
  len1=paste(out_dir,"/",out_basename,".len",sep="")
  # Remove contigs without genes
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3!=0) print $0}'",len1,"> len1",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system(paste("rm",len1,sep=" "));system(paste("mv","len1",len1,sep=" "))
  # Format gff
  cmd=paste("awk -F '\t' -v OFS='\t' '{print $1,$7,$3,$4,$5,$6,$2}'",fake_gff1,"> fake_gff1",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system(paste("rm",fake_gff1,sep=" "));system(paste("mv","fake_gff1",fake_gff1,sep=" "))
  if (!self2self){
    wgdi_input(genome.fna=genome.fna2,
               gff=gff2,
               out_prefix=out_basename)
    fake_gff2=paste(out_dir,"/",out_basename,".gff",sep="")
    len2=paste(out_dir,"/",out_basename,".len",sep="")
    # Remove contigs without genes
    cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3!=0) print $0}'",len2,"> len2",sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system(paste("rm",len2,sep=" "));system(paste("mv","len2",len2,sep=" "))
    # Format gff
    cmd=paste("awk -F '\t' -v OFS='\t' '{print $1,$7,$3,$4,$5,$6,$2}'",fake_gff2,"> fake_gff2",sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system(paste("rm",fake_gff2,sep=" "));system(paste("mv","fake_gff2",fake_gff2,sep=" "))
  }else{
    fake_gff2=paste(out_dir,"/",out_basename,".gff",sep="")
    len2=paste(out_dir,"/",out_basename,".len",sep="")
  }
  
  # Run blast
  system("mkdir tmp",wait=TRUE)
  setwd("tmp/")
  cmd=paste("cp"," ",proteins.faa1," ",".",sep="")
  print(cmd);system(cmd,wait=TRUE)
  proteins.faa1=unlist(strsplit(proteins.faa1,"/"));proteins.faa1=proteins.faa1[length(proteins.faa1)]
  if (!self2self){
    cmd=paste("cp"," ",proteins.faa2," ",".",sep="")
    print(cmd);system(cmd,wait=TRUE)
  }
  proteins.faa2=unlist(strsplit(proteins.faa2,"/"));proteins.faa2=proteins.faa2[length(proteins.faa2)]
  if (!file.exists(paste("../",out_basename,".blast",sep=""))){
    cmd=paste("makeblastdb",
              "-in",proteins.faa1,
              "-dbtype prot",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("blastp",
              "-num_threads",threads,
              "-db",proteins.faa1,
              "-query",proteins.faa2,
              "-outfmt 6",
              "-evalue 1e-5",
              "-num_alignments 20",
              "-out",paste(out_basename,".blast",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    blast1=paste(out_basename,".blast",sep="")
    system(paste("mv",blast1,"../",sep=" "),wait=TRUE)
  }
  if (!self2self){
    if (!file.exists(paste("../",out_basename,"_reverse.blast",sep=""))){
      cmd=paste("makeblastdb",
                "-in",proteins.faa2,
                "-dbtype prot",
                sep=" ")
      print(cmd);system(cmd,wait=TRUE)
      cmd=paste("blastp",
                "-num_threads",threads,
                "-db",proteins.faa2,
                "-query",proteins.faa1,
                "-outfmt 6",
                "-evalue 1e-5",
                "-num_alignments 20",
                "-out",paste(out_basename,"_reverse.blast",sep=""),
                sep=" ")
      print(cmd);system(cmd,wait=TRUE)
      blast2=paste(out_basename,"_reverse.blast",sep="")
      system(paste("mv",blast2,"../",sep=" "),wait=TRUE)
    }
  }else{
    #blast2=paste(out_dir,"/",out_basename,".blast",sep="")
    blast2="false"
  }
  system(paste("mv",proteins.faa1,"../",sep=" "),wait=TRUE)
  if (!self2self){system(paste("mv",proteins.faa2,"../",sep=" "),wait=TRUE)}
  setwd("../")
  system("rm -r tmp",wait=TRUE)
  system(paste("cp"," ",cds.fna1," ",".",sep=""),wait=TRUE)
  system(paste("cp"," ",cds.fna2," ",".",sep=""),wait=TRUE)
  
  setwd(wd)
}

# Genome-genome dotplot
# Dependencies: wgdi
wgdi_dot=function(blast=blast,blast_reverse="false",
                  fake.gff1=fake.gff1,fake.gff2=fake.gff2,
                  lens1=lens1,lens2=lens2,
                  genome1_name=genome1_name,genome2_name=genome2_name,
                  out_basename=out_basename,
                  out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  getFile=function(i){i=unlist(strsplit(i,"/"));return(i[length(i)])}
  
  cmd="wgdi -d \\?"
  conf=system(cmd,wait=TRUE,intern=TRUE)
  conf[2]=paste("blast = ",blast,sep="")
  conf[3]=paste("gff1 = ",fake.gff1,sep="")
  conf[4]=paste("gff2 = ",fake.gff2,sep="")
  conf[5]=paste("lens1 = ",lens1,sep="")
  conf[6]=paste("lens2 = ",lens2,sep="")
  conf[7]=paste("genome1_name = ",genome1_name,sep="")
  conf[8]=paste("genome2_name = ",genome2_name,sep="")
  conf[9]="multiple  = 1"
  conf[10]="score = 100"
  conf[11]="evalue = 1e-5"
  conf[12]="repeat_number = 10"
  conf[13]="position = order"
  conf[14]=paste("blast_reverse = ",blast_reverse,sep="")
  conf[15]="ancestor_left = none"
  conf[16]="ancestor_top = none"
  conf[17]="markersize = 20"
  conf[18]="figsize = 100,100"
  conf[19]=paste("savefig = ",out_basename,"_dotplot.pdf",sep="")
  writeLines(conf,paste(out_basename,"_dotplot.conf",sep=""))
  
  cmd=paste("wgdi -d",
            paste(out_basename,"_dotplot.conf",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# colinearity
# Dependencies: wgdi
wgdi_col=function(blast=blast,blast_reverse="false",
                  fake.gff1=fake.gff1,fake.gff2=fake.gff2,
                  lens1=lens1,lens2=lens2,
                  genome1_name=genome1_name,genome2_name=genome2_name,
                  pvalue=1,
                  out_basename=out_basename,
                  out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  pvalue=as.character(pvalue)
  
  cmd="wgdi -icl \\?"
  conf=system(cmd,wait=TRUE,intern=TRUE)
  conf[2]=paste("gff1 = ",fake.gff1,sep="")
  conf[3]=paste("gff2 = ",fake.gff2,sep="")
  conf[4]=paste("lens1 = ",lens1,sep="")
  conf[5]=paste("lens2 = ",lens2,sep="")
  conf[6]=paste("blast = ",blast,sep="")
  conf[7]=paste("blast_reverse = ",blast_reverse,sep="")
  conf[8]="multiple  = 1"
  conf[9]="process = 8"
  conf[10]="evalue = 1e-5"
  conf[11]="score = 100"
  conf[12]="grading = 50,40,25"
  conf[13]="mg = 40,40"
  conf[14]=paste("pvalue = ",pvalue,sep="")
  conf[15]="repeat_number = 10"
  conf[16]="positon = order"
  conf[17]=paste("savefile = ",out_basename,"_colinearity.txt",sep="")
  writeLines(conf,paste(out_basename,"_colinearity.conf",sep=""))
  
  cmd=paste("wgdi -icl",
            paste(out_basename,"_colinearity.conf",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# ka, ks
# Dependencies: wgdi
wgdi_kaks=function(wgdi_colinearity=wgdi_colinearity,
                   cds.fna=cds.fna,
                   protein.faa=protein.faa,
                   out_basename=out_basename,
                   out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd="wgdi -ks \\?"
  conf=system(cmd,wait=TRUE,intern=TRUE)
  conf[2]=paste("cds_file = ",cds.fna,sep="")
  conf[4]=paste("pep_file = ",protein.faa,sep="")
  conf[7]=paste("pairs_file = ",wgdi_colinearity,sep="")
  conf[8]=paste("ks_file = ",out_basename,"_kaks.txt",sep="")
  writeLines(conf,paste(out_basename,"_kaks.conf",sep=""))
  
  cmd=paste("wgdi -ks",
            paste(out_basename,"_kaks.conf",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# wgdi blockinfo
# Dependencies: wgdi
wgdi_bi=function(blast=blast,
                 fake.gff1=fake.gff1,fake.gff2=fake.gff2,
                 lens1=lens1,lens2=lens2,
                 wgdi_colinearity=wgdi_colinearity,
                 wgdi_kaks=wgdi_kaks,
                 ks_col="ks_NG86", # ks_NG86/ks_YN00
                 out_basename=out_basename,
                 out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd="wgdi -bi \\?"
  conf=system(cmd,wait=TRUE,intern=TRUE)
  conf[2]=paste("blast = ",blast,sep="")
  conf[3]=paste("gff1 = ",fake.gff1,sep="")
  conf[4]=paste("gff2 = ",fake.gff2,sep="")
  conf[5]=paste("lens1 = ",lens1,sep="")
  conf[6]=paste("lens2 = ",lens2,sep="")
  conf[7]=paste("collinearity = ",wgdi_colinearity,sep="")
  conf[8]="score = 100" 
  conf[9]="evalue = 1e-5" 
  conf[10]="repeat_number = 10"
  conf[11]="position = order"  
  conf[12]=paste("ks = ",wgdi_kaks,sep="")
  conf[13]=paste("ks_col =",ks_col,sep="")
  conf[14]=paste("savefile =",out_basename,"_BlockInfo.csv",sep="")
  writeLines(conf,paste(out_basename,"_BlockInfo.conf",sep=""))
  
  cmd=paste("wgdi -bi",
            paste(out_basename,"_BlockInfo.conf",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}

# ks+collinearity
# Dependencies: wgdi
wgdi_bk=function(wgdi_BlockInfo=wgdi_BlockInfo,
                 lens1=lens1,lens2=lens2,
                 genome1_name=genome1_name,genome2_name=genome2_name,
                 p_block=1, # 0-1,
                 out_basename=out_basename,
                 out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd="wgdi -bk \\?"
  conf=system(cmd,wait=TRUE,intern=TRUE)
  conf[2]=paste("lens1 = ",lens1,sep="")                    
  conf[3]=paste("lens2 = ",lens2,sep="")                    
  conf[4]=paste("genome1_name = ",genome1_name,sep="")          
  conf[5]=paste("genome2_name = ",genome2_name,sep="")          
  conf[6]=paste("blockinfo = ",wgdi_BlockInfo,sep="")
  conf[7]=paste("pvalue = ",as.character(p_block),sep="")
  conf[8]="tandem = true"
  conf[9]="tandem_length = 200"
  conf[10]="markersize = 20"  
  conf[11]="area = 0,10"
  conf[12]="block_length = 0"
  conf[13]="figsize = 100,100"   
  conf[14]=paste("savefig = ",out_basename,"_BlockKs.pdf",sep="")
  writeLines(conf,paste(out_basename,"_BlockKs.conf",sep=""))
  
  cmd=paste("wgdi -bk",
            paste(out_basename,"_BlockKs.conf",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# ks peaks of collinearity block
# Dependencies: wgdi
wgdi_kspeak=function(wgdi_BlockInfo=wgdi_BlockInfo,
                     p_block=1, # 0-1,
                     out_basename=out_basename,
                     out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd="wgdi -kp \\?"
  conf=system(cmd,wait=TRUE,intern=TRUE)
  conf[2]=paste("blockinfo = ",wgdi_BlockInfo,sep="")
  conf[3]=paste("pvalue = ",as.character(p_block),sep="")
  conf[4]="tandem = true"  
  conf[5]="block_length = 5"
  conf[6]="ks_area = 0,10"
  conf[7]="multiple  = 1" 
  conf[8]="homo = 0,1"
  conf[9]="fontsize = 9"  
  conf[10]="area = 0,10"
  conf[11]="figsize = 10,6.18" 
  conf[12]=paste("savefig = ",out_basename,"_KsPeaks.pdf",sep="")
  conf[13]=paste("savefile = ",out_basename,"_KsPeaks.csv",sep="")
  writeLines(conf,paste(out_basename,"_KsPeaks.conf",sep=""))
  
  cmd=paste("wgdi -kp",
            paste(out_basename,"_KsPeaks.conf",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# # Ks peaks
# cmd="wgdi -kp \\?"
# conf=system(cmd,wait=TRUE,intern=TRUE)
# conf[2]=paste("blockinfo = ",wgdi_BlockInfo,sep="")
# conf[3]="pvalue = 0.05"
# conf[5]="block_length = 5"
# conf[6]="ßks_area = 0,10"
# conf[11]="figsize = 10,6.18"
# conf[12]=paste("savefig = ",out_basename,"_KsPeak.pdf",sep="")
# conf[13]=paste("savefile = ",out_basename,"_KsPeak.csv",sep="")
# writeLines(conf,paste(out_basename,"_KsPeak.conf",sep=""))
# cmd=paste("wgdi -kp",
#           paste(out_basename,"_KsPeak.conf",sep=""),
#           sep=" ")
# print(cmd);system(cmd,wait=TRUE)
# 
# # Ks peak fit
# cmd="wgdi -pf \\?"
# conf=system(cmd,wait=TRUE,intern=TRUE)
# conf[2]=paste("blockinfo = ",wgdi_BlockInfo,sep="")
# conf[10]=paste("savefig = ",out_basename,"_PeakFit.pdf",sep="")
# writeLines(conf,paste(out_basename,"_PeakFit.conf",sep=""))
# cmd=paste("wgdi -pf",
#           paste(out_basename,"_PeakFit.conf",sep=""),
#           sep=" ")
# print(cmd);system(cmd,wait=TRUE)
# 
# # Ks figure
# cmd="wgdi -kf \\?"
# conf=system(cmd,wait=TRUE,intern=TRUE)
# conf[2]=paste("ksfit = ",out_basename,"_KsPeak.csv",sep="")
# conf[11]=paste("savefig = ",out_basename,"_KsFig.pdf",sep="")
# writeLines(conf,paste(out_basename,"_KsFig.conf",sep=""))
# cmd=paste("wgdi -kf",
#           paste(out_basename,"_KsFig.conf",sep=""),
#           sep=" ")
# print(cmd);system(cmd,wait=TRUE)

#####
# HGT
#####
# Alien index
# Dependencies: DIAMOND, parallel (R)
AI=function(pep.faa=pep.faa,
            pep.kingdom="Metazoa",
            pep.phylum="Arthropoda",
            nr.dmdb="/apps/unit/BioinfoUgrp/DB/diamondDB/ncbi/2022-07/nr.dmnd",
            ai=45, # threshold for alien index
            out_pct=0.9, # threshold for out_pct
            out_dir=out_dir,out_basename=out_basename,
            threads=threads){
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  out_dir=sub("/$","",out_dir)
  wd=getwd();setwd(out_dir)
  
  # blast search
  if (!file.exists("blast.tsv")){
    cmd=paste("diamond blastp",
              "--threads",threads,
              "--db",nr.dmdb,
              "--query",pep.faa,
              "--out blast.tsv",
              "--max-target-seqs 500",
              "--min-score 50",
              "--outfmt 6 qseqid sseqid evalue bitscore length pident skingdoms sphylums",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  if (!file.exists(paste(out_basename,"_AI.tsv",sep=""))){
    pepID.lst=system("awk -F '\t' -v OFS='\t' '{print $1}' blast.tsv",intern=TRUE)
    pepID.lst=pepID.lst[!duplicated(pepID.lst)]
    res=data.frame(ID=pepID.lst,
                   AI=rep(NA,length(pepID.lst)),
                   out_pct=rep(NA,length(pepID.lst)),
                   kingdom=rep(NA,length(pepID.lst)),
                   phylum=rep(NA,length(pepID.lst)))
    rownames(res)=res$ID
    blast=read.table("blast.tsv",sep="\t",header=FALSE,quote="",comment.char="")
    colnames(blast)=c("qseqid","sseqid","evalue","bitscore","length","pident","skingdoms","sphylums")
    
    library(parallel)
    clus=makeCluster(as.numeric(threads))
    clusterExport(clus,list("pep.kingdom","pep.phylum","res","blast"),envir=environment())
    res[,c("AI","out_pct","kingdom","phylum")]=
      t(parSapply(clus,
                pepID.lst,
                function(ID){
                 d=blast[blast[,"qseqid"]==ID,]
                 if (nrow(d)<25){
                   return(c(NA,NA,NA,NA))
                 }else{
                   ingroup.eval=d[grepl(pep.kingdom,d[,"skingdoms"]) & !grepl(pep.phylum,d[,"sphylums"]),
                                  "evalue"]
                   ingroup.eval=as.numeric(ingroup.eval)
                   if (length(ingroup.eval)==0){ingroup.eval=Inf}
        
                   outgroup.eval=d[!grepl(pep.kingdom,d[,"skingdoms"]) & !grepl(pep.phylum,d[,"sphylums"]),
                                   "evalue"]
                   outgroup.eval=as.numeric(outgroup.eval)
                   if (length(outgroup.eval)==0){outgroup.eval=Inf}
                   
                   ai=log(min(ingroup.eval)+1e-1000)-log(min(outgroup.eval)+1e-1000)
                   if (is.na(ai)){ai=-Inf}
                    
                   out_pct=length(outgroup.eval)/nrow(d)
                   if (is.na(out_pct)){out_pct=0}
                                 
                   kingdom=system(paste("awk -F '\t' -v OFS='\t'",
                                        paste("'{if ($1==\"",ID,"\") print $7}'",sep=""),
                                        "blast.tsv",sep=" "),
                                  intern=TRUE) 
                   kingdom=paste(kingdom,collapse=";")
                   kingdom=unlist(strsplit(kingdom,";"))
                   kingdom=kingdom[!duplicated(kingdom)]
                                
                   phyla=system(paste("awk -F '\t' -v OFS='\t'",
                                paste("'{if ($1==\"",ID,"\") print $8}'",sep=""),
                                "blast.tsv",sep=" "),
                                intern=TRUE) 
                   phyla=paste(phyla,collapse=";")
                   phyla=unlist(strsplit(phyla,";"))
                   phyla=phyla[!duplicated(phyla)]
                   return(c(ai,out_pct,paste(kingdom,collapse = ";"),paste(phyla,collapse = ";")))
                 }
    }))
    stopCluster(clus)
    write.table(res,paste(out_basename,"_AI.tsv",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
  }
  
  res=read.table(paste(out_basename,"_AI.tsv",sep=""),header=TRUE,sep="\t",quote="")
  HGT=res[res[,"AI"]>ai,] # Might be too strict
  HGT=HGT[HGT[,"out_pct"]>out_pct,]
  write.table(HGT,paste(out_basename,"_HGT.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
  
  setwd(wd)
}
# for (sp in sp.lst){
#   out_basename=sp
#   res=read.table(paste(sp,"/",out_basename,"_AI.tsv",sep=""),header=TRUE,sep="\t",quote="")
#   HGT=res[res[,"AI"]>ai,] # Might be too strict
#   HGT=HGT[HGT[,"out_pct"]>out_pct,]
#   write.table(HGT,paste(sp,"/",out_basename,"_HGT.tsv",sep=""),
#               sep="\t",row.names=FALSE,quote=FALSE)
# }

# Find homology sequences from nr by DIAMOND
# for HGT
# Dependencies: DIAMOND+nr, seqkit, stringr (R)
findHomo=function(pep.faa=pep.faa,
                  nr.dmdb="/apps/unit/BioinfoUgrp/DB/diamondDB/ncbi/2022-07/nr.dmnd",
                  threads=threads,
                  out_prefix=out_prefix){
  threads=as.character(threads)
  # /flash/BourguignonU/Cong/termite_pca/AI/hgtTree/allPep.faa
  blast=paste(out_prefix,".blast",sep="")
  if (!file.exists(blast)){
    cmd=paste("diamond blastp",
              "--threads",threads,
              "--db",nr.dmdb,
              "--query",pep.faa,
              "--out",blast,
              "--max-target-seqs 100",
              "--min-score 50",
              "--outfmt 6 qseqid sseqid evalue bitscore length pident qcovhsp scovhsp skingdoms sphylums sscinames staxids stitle full_sseq",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  if (file.size(blast)==0){
    print("No hits found in");print(blast)
  }else{
    blast=read.table(blast,sep="\t",header=FALSE,quote="",comment.char="")
    colnames(blast)=c("qseqid","sseqid","evalue","bitscore","length","pident","qcovhsp","scovhsp","skingdoms","sphylums","sscinames","staxids","stitle","full_sseq")
    blast=blast[!duplicated(blast[,"sseqid"]),]
    blast=blast[!duplicated(blast[,"staxids"]),]
    blast=blast[!grepl("[BZ]",blast[,"full_sseq"]),]
    blast=blast[!grepl("partial",blast[,"stitle"]),]
    library(stringr)
    blast[,"species"]=str_extract(blast[,"stitle"],"\\[.*\\]$")
    blast=blast[!duplicated(blast[,"species"]),]
    
    for (i in 1:nrow(blast)){
      sp=sub("\\[","",blast[i,"species"]);sp=sub("\\]","",sp)
      sp=paste(unlist(strsplit(sp," ")),collapse=".")
      header=paste(blast[i,"sseqid"],"_",sp,"_",blast[i,"sphylums"],"_",blast[i,"skingdoms"],sep="")
      write(paste(">",header,sep=""),
            paste(out_prefix,".faa",sep=""),
            append=TRUE)
      write(blast[i,"full_sseq"],
            paste(out_prefix,".faa",sep=""),
            append=TRUE)
    }
    
    system(paste("mv",
                 paste(out_prefix,".faa",sep=""),
                 paste(out_prefix,"_nr.faa",sep=""),
                 sep=" "))
    # system(paste("mkdir ",out_prefix,".tmp",sep=""))
    # system(paste("cp ",out_prefix,".faa"," ",out_prefix,".tmp/seq.faa",sep=""))
    # wd=getwd();setwd(paste(out_prefix,".tmp",sep=""))
    # cmd="mmseqs createdb seq.faa sequenceDB"
    # print(cmd);system(cmd,wait=TRUE)
    # cmd=paste("mmseqs cluster sequenceDB clusterDB .",
    #           "--cov-mode 0",
    #           "-c 0.9",
    #           "--min-seq-id 0.9",
    #           "--threads",threads,
    #           sep=" ")
    # print(cmd);system(cmd,wait=TRUE)
    # cmd="mmseqs createsubdb clusterDB sequenceDB clusterDB_rep"
    # print(cmd);system(cmd,wait=TRUE)
    # cmd="mmseqs convert2fasta clusterDB_rep rep.fasta"
    # print(cmd);system(cmd,wait=TRUE)
    # system(paste("cp rep.fasta ../",basename(out_prefix),"_nr.faa",sep=""))
    # setwd(wd)
    # system(paste("rm -r ",out_prefix,".tmp",sep=""))
    
    cmd=paste("seqkit seq -w0",
              pep.faa,">>",paste(out_prefix,"_nr.faa",sep=""),
              sep=" ")
    system(cmd)
  }
}

# With NCBI protein ID, get corresponding CDS and taxonomic lineage
# For each protein ID, ONE CDS is RANDOMLY picked
# Dependencies: Entrez Direct installed by $ sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
#               Biostrings (R), taxonomizr (R), parallel (R)
protID2cds=function(protID.lst=protID.lst, # comma list
                    taxonomizr.sql="/bucket/BourguignonU/Cong/public_db/ncbiTaxonomy/Taxonomy.sql",
                    threads=threads,
                    out_prefix=out_prefix){
  protID.lst=unlist(strsplit(protID.lst,","))
  # protID.lst="WP_255322253.1";out_prefix="protID2cds_cds"
  for (prot in protID.lst){
    cmd=paste("elink",
              "-target nuccore",
              "-db protein",
              "-name protein_nuccore",
              "-id",prot,"|",
              "efetch",
              "-format fasta_cds_na",">>", 
              paste(out_prefix,"_cds.fna",sep=""),
              sep=" ")
    print(cmd)
    system(cmd,wait=TRUE)
    Sys.sleep(1)
  }
  
  library(Biostrings)
  cds.fna=readDNAStringSet(paste(out_prefix,"_cds.fna",sep=""))
  cds.IDs=names(cds.fna)
  protein.IDs=gsub(pattern=".*\\[protein_id=(.*?)\\].*",replacement = "\\1",cds.IDs)
  genome.IDs=gsub(pattern="^lcl|(.*?)_cds_.*",replacement = "\\1",cds.IDs)
  genome.IDs=sub("^lcl\\|","",genome.IDs)
  df=data.frame(cds.ID=cds.IDs,protein.ID=protein.IDs,genome.ID=genome.IDs)
  df=df[df$protein.ID %in% protID.lst,]
  df=df[!duplicated(df$protein.ID),]
  print("No CDS for these proteins:")
  print(protID.lst[!(protID.lst %in% df$protein.ID)])
  rownames(df)=df$cds.ID
  cds.fna=cds.fna[names(cds.fna) %in% df$cds.ID]
  names(cds.fna)=paste(df[names(cds.fna),"protein.ID"],df[names(cds.fna),"genome.ID"],sep="-")
  writeXStringSet(cds.fna,paste(out_prefix,".fna",sep=""))
  
  genome.ID=df$genome.ID
  genome.ID=genome.ID[!duplicated(genome.ID)]
  library(taxonomizr)
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  clusterExport(cl=clus,varlist = list("genome.ID",
                                       "accessionToTaxa","getRawTaxonomy","taxonomizr.sql"))
  txid=parSapply(clus,genome.ID,function(i){return(accessionToTaxa(i,taxonomizr.sql))}) 
  names(txid)=genome.ID
  taxaPath=parSapply(clus,txid,function(i){i=getRawTaxonomy(i,taxonomizr.sql);return(paste(i[[1]],collapse=","))})
  names(taxaPath)=genome.ID
  
  df$txid=txid[df$genome.ID]
  df$taxaPath=taxaPath[df$genome.ID]
  write.table(df,paste(out_prefix,"_taxon.tsv",sep=""),
              sep="\t",row.names = FALSE,quote = FALSE)
}

#####
# Phylogeny
#####
# Collect single-copy busco sequences from busco runs
# Dependencies:
getBuscoSeq=function(tab=tab, # tsv. 
                              # First column is species label and second column is run_/busco_sequences/single_copy_busco_sequences/
                     out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  
  tab=read.table(tab,sep="\t",header=TRUE,quote="")
  system(paste("mkdir"," ",out_dir,"/seqs",sep=""))
  for (i in 1:nrow(tab)){
    sp=tab[i,1]
    seqs=system(paste("ls",tab[i,2],sep=" "),intern=TRUE)
    seqs=paste(tab[i,2],"/",seqs,sep="")
    for (j in seqs){
      out=unlist(strsplit(j,"/"));out=paste(out_dir,"/seqs/",out[length(out)],sep="")
      sed=paste("'s/>.*/>",sp,"/'",sep="")
      cmd=paste("sed",sed,j,">>",out,sep=" ")
      system(cmd,wait=TRUE)
    }
  }
  system("mv seqs/* ./")
}

# Wrangle orthofinder results to be ready for building gene trees
# Dependencies: seqkit, parallel (R)
getOrthofinderSeq=function(N0.tsv=N0.tsv,
                           allPep.faa=allPep.faa,
                           threads=threads,
                           out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  d=read.table(N0.tsv,sep="\t",header=TRUE,quote="")
  HOG=sub("^N0.","",d[,"HOG"])
  d=d[,4:ncol(d)]
  rownames(d)=HOG
  
  library(parallel)
  clus=makeCluster(as.numeric(threads)-1)
  clusterExport(clus,list("d"),envir=environment())
  parSapply(clus,
            rownames(d),
            function(hog){
              IDs=paste(d[hog,],collapse=",")
              IDs=unlist(strsplit(IDs,","))
              IDs=IDs[which(IDs!="")]
              IDs=sub(" ","",IDs)
              writeLines(IDs,paste(hog,".lst",sep=""))
              cmd=paste("seqkit grep",
                        "-f",paste(hog,".lst",sep=""),
                        allPep.faa,
                        ">",paste(hog,".faa",sep=""),
                        sep=" ")
              system(cmd,wait=TRUE)
              system(paste("rm ",hog,".lst",sep=""))
            })
  stopCluster(clus)
  setwd(wd)
}

# MAFFT: multiple sequence alignment
# Dependencies: MAFFT
mafft=function(in.fa=in.faa,
               align.fa=align.fa,
               threads=threads){
  threads=as.character(threads)
  cmd=paste("mafft",
            "--auto",
            "--thread",threads,
            in.fa,">",
            align.fa,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Convert protein alignments to CDS alignments
# Dependencies: PAL2NAL, seqkit
pal2nal=function(protAlign.fa=protAlign.fa,
                 cds.fna=cds.fna, # Including IDs in protAlign.fa
                 geneticCode=1,
                 # 1  Universal code (default)
                 # 2  Vertebrate mitochondrial code
                 # 3  Yeast mitochondrial code
                 # 4  Mold, Protozoan, and Coelenterate Mitochondrial code
                 # and Mycoplasma/Spiroplasma code
                 # 5  Invertebrate mitochondrial
                 # 6  Ciliate, Dasycladacean and Hexamita nuclear code
                 # 9  Echinoderm and Flatworm mitochondrial code
                 # 10  Euplotid nuclear code
                 # 11  Bacterial, archaeal and plant plastid code
                 # 12  Alternative yeast nuclear code
                 # 13  Ascidian mitochondrial code
                 # 14  Alternative flatworm mitochondrial code
                 # 15  Blepharisma nuclear code
                 # 16  Chlorophycean mitochondrial code
                 # 21  Trematode mitochondrial code
                 # 22  Scenedesmus obliquus mitochondrial code
                 # 23  Thraustochytrium mitochondrial code
                 outAlign.fa=outAlign.fa){
  cmd=paste("grep '>' ",protAlign.fa," | sed 's/>//' > ",outAlign.fa,".lst",sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("seqkit faidx ",cds.fna," --infile-list ",outAlign.fa,".lst > ",outAlign.fa,"_CDS.fa",sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("pal2nal.pl",
            protAlign.fa,
            paste(outAlign.fa,"_CDS.fa",sep=""),
            "-output fasta",
            "-codontable",as.character(geneticCode),
            ">",
            outAlign.fa,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm ",outAlign.fa,".lst",sep=""))
  system(paste("rm ",outAlign.fa,"_CDS.fa",sep=""))
}

# concatenate sequences with same ID from multiple files
# Dependencies: seqkit
catMSA=function(align.fa=align.fa, # comma list
                cat.fa=cat.fa){
  align.fa=unlist(strsplit(align.fa,","))
  align.fa=paste(align.fa,collapse=" ")
  cmd=paste("seqkit concat",align.fa,">",cat.fa,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# trimAL: curate multiple sequence alignments for phylogenetic analysis
# Dependencies: trimAL, seqkit
trimAL=function(inMSA.fa=inMSA.fa,
                outMSA.fa=outMSA.fa){
  cmd=paste("trimal",
            "-in",inMSA.fa,
            "-automated1",
            #"-keepheader",
            "| seqkit seq -u",
            ">",outMSA.fa,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Phylogenetic tree
# Dependencies: iqtree
iqtree=function(msa.fa=msa.fa, 
                type="dna", # dna/protein
                out_prefix=out_prefix,
                threads=threads){
  threads=as.character(threads)
  cmd=paste("iqtree",
            "-s",msa.fa,
            "--prefix",out_prefix,
            "-T",threads,
            # "-m TEST", # Standard model selection followed by tree inference
            # "-m MFP", # Extended model selection followed by tree inference
            # "-mset GTR", # test GTR+... models (DNA) only
            # "--msub nuclear", # Amino-acid model source (nuclear, mitochondrial, chloroplast or viral)
            "-B 1000",
            sep=" ")
  if (type=="dna"){cmd=paste(cmd,"-mset GTR",sep=" ")}
  if (type=="protein"){cmd=paste(cmd,"--msub nuclear",sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
}

# ASTRAL: supertree
# Dependencies: ASTRAL
astral=function(trees=trees, # comma-list, nwk trees
                path2astral="/home/c/c-liu/Softwares/ASTRAL/astral.5.7.8.jar",
                out_prefix=out_prefix){
  trees=unlist(strsplit(trees,","))
  for (tree in trees){
    cmd=paste("cat ",tree," >> ",out_prefix,"_inTrees.nwk",sep="")
    system(cmd,wait=TRUE)
  }
  cmd=paste("java -jar",path2astral,
            "-i",paste(out_prefix,"_inTrees.nwk",sep=""),
            "-o",paste(out_prefix,"_astralTree.nwk",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# STAG: supertree
# Dependencies: STAG
stag=function(gene2species.txt=gene2species.txt,# space-separated
              path2Stag="/home/c/c-liu/Softwares/STAG-1.0.0/stag/stag.py",
              tree_dir=tree_dir,
              out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  cmd=paste("python2",path2Stag,gene2species.txt,tree_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}

#####
# Gene content
#####
# Notung: root gene trees
# Dependencies: Notung, parallel (R), phytools (R)
rootNotung=function(Notung.jar="/home/c/c-liu/Softwares/Notung/Notung-2.9.1.5.jar",
                    speciesTree.nwk=speciesTree.nwk, # rooted, absolute path
                    geneTrees.nwk=geneTrees.nwk, # unrooted, comma-list, absolute path
                    # Link species with genes:
                    # musgene1, bov.gene1, dros-gene1, etc.
                    inferTransfers="false", # true/false. Consider gene transfers or not
                    out_dir=out_dir,
                    threads=threads){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  out_dir=sub("/$","",out_dir)
  wd=getwd();setwd(out_dir)
  
  system("unset DISPLAY")
  geneTrees.nwk=unlist(strsplit(geneTrees.nwk,","))
  
  library(parallel)
  clus=makeCluster(as.numeric(threads)-1)
  clusterExport(clus,list("speciesTree.nwk"),envir=environment())
  parSapply(clus,
            geneTrees.nwk,
            function(geneTree){
              cmd=paste("java -jar",
                        Notung.jar,
                        "-g",geneTree,
                        "-s",speciesTree.nwk,
                        "--root",
                        "--outputdir",out_dir,
                        "--infertransfers",inferTransfers,
                        "--absfilenames",
                        "--nolosses",
                        "--speciestag prefix",
                        "--treeoutput newick",
                        sep=" ")
              system(cmd,wait=TRUE)
            })
  stopCluster(clus)
  
  setwd(wd)
}

# Notung: gene content evolution
phyloNotung=function(Notung.jar="/home/c/c-liu/Softwares/Notung/Notung-2.9.1.5.jar",
                     speciesTree.nwk=speciesTree.nwk, # rooted, absolute path
                     geneTrees.nwk=geneTrees.nwk, # rooted, comma-list, absolute path
                     # Link species with genes:
                     # musgene1, bov.gene1, dros-gene1, etc.
                     inferTransfers="false", # true/false. Consider gene transfers or not
                     out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  system("unset DISPLAY")
  
  trees=paste(speciesTree.nwk,geneTrees.nwk,sep=",")
  trees=unlist(strsplit(trees,","))
  writeLines(trees,"inputTrees.lst")
  
  cmd=paste("java -jar",
            Notung.jar,
            "-b inputTrees.lst",
            "--reconcile",
            "--outputdir",out_dir,
            "--infertransfers",inferTransfers,
            "--absfilenames",
            "--speciestag prefix",
            "--treeoutput newick",
            "--phylogenomics",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

#####
# Selective pressure
#####
# absrel of hyphy
# Dependencies: hyphy
absrel=function(codon.align=codon.align,
                geneTree=geneTree,
                out_prefix=out_prefix){
  cmd=paste("hyphy absrel",
            "--alignment",codon.align,
            "--tree",geneTree,
            "--output",paste(out_prefix,".json",sep=""))
  print(cmd);system(cmd,wait=TRUE)
}
# hyphy absrel --alignment hiv1_transmission.fna --tree tree.nwk

# relax of hyphy
# hyphy relax --alignment pb2.fna --tree tree.nwk --test test

#####
# gain, birth, death and innovation gene (or DNA element) family rates
#####
# Badirate: Birth, death and innovation model
# Dependencies: Badirate, badirater (R), parallel (R), 
# readr (R), dplyr (R), stringr (R), ggtree (R), treeio (R)
badirate=function(perl.5.16="~/Softwares/perl-5.16.3/perl",
                  badirate.pl="~/Softwares/badirate/BadiRate.pl",
                  tree.nwk=tree.nwk,
                  famSize.tsv=famSize.tsv, # Fields: FAM_ID sp1 sp2
                  out_dir=out_dir,
                  threads=threads){
  library(badirater)
  library(parallel)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggtree)
  library(treeio)
  # perl.5.16="~/Softwares/perl-5.16.3/perl"
  # badirate.pl="~/Softwares/badirate/BadiRate.pl"
  # tree.nwk="/flash/BourguignonU/Cong/termite_pca/badirate/test/droso.12sp.tamura.nwk"
  # famSize.tsv="/flash/BourguignonU/Cong/termite_pca/badirate/test/og.tsv"
  # out_dir="/flash/BourguignonU/Cong/termite_pca/badirate/test/"
  # threads=2
  
  #turnOverRates=c("BDI","GD")
  #branchModels=c("FR","GR")
  #startValues=c("0","1","2","3")
  
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  out_dir=sub("/$","",out_dir)
  
  # split gene families
  famSize=read.table(famSize.tsv,header=TRUE,sep="\t",quote="")
  if (!file.exists(paste(out_dir,"/og.split",sep=""))){
    system(paste("mkdir ",out_dir,"/og.split",sep=""))
  }
  og.lst=system(paste("ls ",out_dir,"/og.split/*.tsv",sep=""),intern=TRUE)
  if (nrow(famSize)!=length(og.lst)){
    sapply(1:nrow(famSize), 
           function(i){
             if (!file.exists(paste(out_dir,"/og.split/",famSize[i,1],".tsv",sep=""))){
               write.table(famSize[i,],
                           paste(out_dir,"/og.split/",famSize[i,1],".tsv",sep=""),
                           sep="\t",row.names = FALSE,quote=FALSE)
             }})
  }
  
  # branch ID
  if (!file.exists(paste(out_dir,"/branchID.nwk",sep=""))){
    cmd=paste(perl.5.16,
              badirate.pl,
              "-print_ids",
              "-treefile",tree.nwk,
              "-sizefile",famSize.tsv,
              ">",paste(out_dir,"/branchID.nwk",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }

  # Setup badirate
  library(badirater)
  setupTab.0=prepare_badirate(og_path=paste(out_dir,"/og.split/",sep=""),
                              tree=tree.nwk,
                              branch_models=c(gr="GR",fr="FR"),
                              rate_model="BDI",
                              estimation="ML",
                              out_dir=paste(out_dir,"/rawOutput.0",sep=""),
                              script_dir=paste(out_dir,"/scripts.0",sep=""),
                              replicates=3,
                              ancestral=TRUE,
                              outlier=TRUE,
                              seed=20231212,
                              start_value=0,
                              pbs_q="smps",
                              badirate_path=badirate.pl,
                              create_scripts="pbs")
  if (!file.exists(paste(out_dir,"/setupTab.0.tsv",sep=""))){
    write.table(setupTab.0,paste(out_dir,"/setupTab.0.tsv",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
  }
  setupTab.1=prepare_badirate(og_path=paste(out_dir,"/og.split/",sep=""), 
                              tree=tree.nwk,
                              branch_models=c(gr="GR",fr="FR"),
                              rate_model="BDI", 
                              estimation="ML",
                              out_dir=paste(out_dir,"/rawOutput.1",sep=""),
                              script_dir=paste(out_dir,"/scripts.1",sep=""),
                              replicates=3,
                              ancestral=TRUE, 
                              outlier=TRUE,
                              seed=20231212,
                              start_value=1,
                              pbs_q="smps",
                              badirate_path=badirate.pl,
                              create_scripts="pbs")
  if (!file.exists(paste(out_dir,"/setupTab.1.tsv",sep=""))){
    write.table(setupTab.1,paste(out_dir,"/setupTab.1.tsv",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
  }
  
  # badirate cmd
  og.lst=system(paste("ls ",out_dir,"/og.split/*.tsv",sep=""),intern=TRUE)
  if (!file.exists(paste(out_dir,"/rawOutput.0/gr/",sep=""))){
    system(paste("mkdir ",out_dir,"/rawOutput.0/gr/",sep=""))
    system(paste("mkdir ",out_dir,"/rawOutput.0/fr/",sep=""))
    system(paste("mkdir ",out_dir,"/rawOutput.1/gr/",sep=""))
    system(paste("mkdir ",out_dir,"/rawOutput.1/fr/",sep=""))
  }
  
  library(stringr)
  cmd.0.GR.1=paste(perl.5.16,badirate.pl,
                 "-seed 20231212",
                 "-start_val 0",
                 "-rmodel BDI",
                 "-bmodel GR",
                 "-ep ML",
                 "-treefile",tree.nwk,
                 "-sizefile",og.lst,
                 "-anc -outlier",
                 "-out",paste(out_dir,"/rawOutput.0/gr/",
                              basename(og.lst),".gr01.bd",sep=""),
                  sep=" ")
  cmd.0.GR.2=str_replace(cmd.0.GR.1,"gr01","gr02")
  cmd.0.GR.3=str_replace(cmd.0.GR.1,"gr01","gr03")
  cmd.0.FR.1=paste(perl.5.16,badirate.pl,
                   "-seed 20231212",
                   "-start_val 0",
                   "-rmodel BDI",
                   "-bmodel FR",
                   "-ep ML",
                   "-treefile",tree.nwk,
                   "-sizefile",og.lst,
                   "-anc -outlier",
                   "-out",paste(out_dir,"/rawOutput.0/fr/",
                                basename(og.lst),".fr01.bd",sep=""),
                  sep=" ")
  cmd.0.FR.2=str_replace(cmd.0.FR.1,"fr01","fr02")
  cmd.0.FR.3=str_replace(cmd.0.FR.1,"fr01","fr03")
  cmd.1.GR.1=paste(perl.5.16,badirate.pl,
                   "-seed 20231212",
                   "-start_val 1",
                   "-rmodel BDI",
                   "-bmodel GR",
                   "-ep ML",
                   "-treefile",tree.nwk,
                   "-sizefile",og.lst,
                   "-anc -outlier",
                   "-out",paste(out_dir,"/rawOutput.1/gr/",
                               basename(og.lst),".gr01.bd",sep=""),
                   sep=" ")
  cmd.1.GR.2=str_replace(cmd.1.GR.1,"gr01","gr02")
  cmd.1.GR.3=str_replace(cmd.1.GR.1,"gr01","gr03")
  cmd.1.FR.1=paste(perl.5.16,badirate.pl,
                   "-seed 20231212",
                   "-start_val 1",
                   "-rmodel BDI",
                   "-bmodel FR",
                   "-ep ML",
                   "-treefile",tree.nwk,
                   "-sizefile",og.lst,
                   "-anc -outlier",
                   "-out",paste(out_dir,"/rawOutput.1/fr/",
                                basename(og.lst),".fr01.bd",sep=""),
                  sep=" ")
  cmd.1.FR.2=str_replace(cmd.1.FR.1,"fr01","fr02")
  cmd.1.FR.3=str_replace(cmd.1.FR.1,"fr01","fr03")
  cmd=c(cmd.0.GR.1,cmd.0.GR.2,cmd.0.GR.3,
        cmd.0.FR.1,cmd.0.FR.2,cmd.0.FR.3,
        cmd.1.GR.1,cmd.1.GR.2,cmd.1.GR.3,
        cmd.1.FR.1,cmd.1.FR.2,cmd.1.FR.3)
  out.file=sapply(cmd,
                  function(i){
                    out.file=unlist(strsplit(i," "))
                    out.file=out.file[length(out.file)]
                    return(out.file)
                  })
  out.file=unname(out.file)
  err.file=paste(out.file,".stderr",sep="")
  cmd=paste(cmd,">",err.file,sep=" ")
  end.output=sapply(out.file,
                    function(i){
                      res=TRUE
                      if (!file.exists(i)){res=FALSE}
                      if (file.exists(i)){
                        if (file.size(i)==0){res=FALSE}
                      }
                      return(res)
                    })
  cmd=cmd[which(!end.output)]
  
  # run badirate
  print("Start badirate jobs...")
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  parSapply(clus,cmd,
            function(i){system(i,wait=TRUE)})
  print("Collect badirate results...")
  bd_collect(setupTab.0,out_dir=paste(out_dir,"/results.0",sep=""))
  bd_collect(setupTab.1,out_dir=paste(out_dir,"/results.1",sep=""))
  
  # model selection
  bd_model_select <- function(setup_table, results_dir){
    dplyr::tibble(file_path = list.files(path = results_dir, pattern =  ".likelihood.txt")) %>%
      dplyr::mutate(model = file_path %>% stringr::str_sub(0,2),
                    replicates = file_path %>% stringr::str_sub(3,4) %>% as.numeric()) %>%
      dplyr::left_join(setup_table %>% select(-replicates), by = "model") %>%
      dplyr::mutate(data = file_path %>% purrr::map(function(x){paste(results_dir, "/", x, sep = "") %>%
          readr::read_delim(col_names=c("filename", "nann", "likelihood"), delim = " ", col_types = cols(col_character(), col_character(), col_double()))})) %>%
      tidyr::unnest() %>%
      dplyr::mutate(og = filename) %>%
      dplyr::select(og, model, parameters, replicates, likelihood) %>%
      dplyr::mutate(aic = 2*(parameters - likelihood)) %>%
      dplyr::group_by(og, model) %>%
      dplyr::slice(which.min(aic)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(og) %>%
      dplyr::mutate(min_aic = min(aic),
                    waic_num = exp((min_aic - aic)/2),
                    waic = waic_num / sum(waic_num),
                    best_waic = max(waic),
                    second_waic = sort(waic)[-2],
                    waic_ratio = best_waic / second_waic) %>%
      dplyr::slice(which.max(waic)) %>%
      dplyr::select(og, model, replicates, aic, waic, waic_ratio) %>%
      dplyr::mutate(signif = waic_ratio > 2.7) %>%
      return()
  }
  bestModels.0=bd_model_select(setup_table=setupTab.0, 
                               results_dir=paste(out_dir,"/results.0/",sep=""))
  bestModels.0=as.data.frame(bestModels.0)
  bestModels.0[,"fam"]=sub("\\.tsv.*$","",basename(bestModels.0[,"og"]))
  #rownames(bestModels.0)=paste(bestModels.0[,"fam"],"-",as.numeric(bestModels.0[,"replicates"]),sep="")
  write.table(bestModels.0,paste(out_dir,"/bestModels.0.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
  bestModels.1=bd_model_select(setup_table=setupTab.1, 
                               results_dir=paste(out_dir,"/results.1",sep=""))
  bestModels.1=as.data.frame(bestModels.1)
  bestModels.1[,"fam"]=sub("\\.tsv.*$","",basename(bestModels.1[,"og"]))
  #rownames(bestModels.1)=paste(bestModels.1[,"fam"],"-",bestModels.1[,"replicates"],sep="")
  write.table(bestModels.1,paste(out_dir,"/bestModels.1.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
  
  fam=c(bestModels.1$fam,bestModels.0$fam)
  fam=fam[!duplicated(fam)]
  bestModels=data.frame(fam=fam)
  bestModels[,"bd"]=sapply(bestModels[,"fam"],
                            function(fam){
                              aic.0=bestModels.0[bestModels.0[,"fam"]==fam,c("aic","og")]
                              aic.1=bestModels.1[bestModels.1[,"fam"]==fam,c("aic","og")]
                              aic=rbind(aic.1,aic.0)
                              bd=aic[aic$aic==max(aic$aic),"og"][1]
                              return(bd)
                            })
#  ~/Softwares/perl-5.16.3/perl ~/Softwares/badirate/BadiRate.pl -seed 20231212 -start_val 0 -rmodel BDI bmodel GR -ep ML -treefile /bucket/BourguignonU/Cong/termite_pca/timeTree_label.nwk -sizefile /flash/BourguignonU/Cong/termite_pca/AI/badirate.true_geneANDpgeneMat_copyMat/og.split/HOG0019764.1.tsv -anc -outlier -out test.bd
  bestModels[,"rateModel"]=rep("BDI",nrow(bestModels))
  bestModels[,"branchModel"]=toupper(basename(dirname(bestModels[,"bd"])))
  bestModels[,"startValue"]=sub("rawOutput\\.","",str_extract(bestModels[,"bd"],"rawOutput\\.[01]*"))
  
  library(ggtree)
  library(treeio)
  branchID.tree=read.newick(paste(out_dir,"/branchID.nwk",sep=""))
  branchID.tree=as.data.frame(ggtree(branchID.tree)$data)
  nodeID.bd=sapply(branchID.tree$label,
                   function(i){
                     if (grepl("_",i)){i=unlist(strsplit(i,"_"));i=i[length(i)]}
                     return(as.numeric(i))
                   })
  node2nodeID.bd=nodeID.bd;names(node2nodeID.bd)=branchID.tree$node
  parentID.bd=node2nodeID.bd[branchID.tree$parent]
  branchID.tree=cbind(branchID.tree,nodeID.bd=nodeID.bd,parentID.bd=parentID.bd)
  ref.tree=as.data.frame(ggtree(read.newick(tree.nwk))$data)
  branchID.tree$label=ref.tree$label
  branchID.tree$branchID.bd=paste(as.character(branchID.tree$parentID.bd),
                                               "->",as.character(branchID.tree$nodeID.bd),
                                               sep="")
  clus=makeCluster(as.numeric(threads))
  clusterExport(cl=clus,varlist = list("branchID.tree","bestModels","read.tree","ggtree","out_dir"))
  res=parSapply(clus,1:nrow(bestModels),
                function(i){
                  #i=2
                  df=branchID.tree
                  fam=rep(bestModels[i,"fam"],nrow(df))
                  branchModel.bd=rep(bestModels[i,"branchModel"],nrow(df))
                  rateModel.bd=rep(bestModels[i,"rateModel"],nrow(df))
                  startValue.bd=rep(bestModels[i,"startValue"],nrow(df))
                  
                  bd=bestModels[i,"bd"]
                  repli=unlist(strsplit(basename(bd),"\\."));repli=repli[length(repli)-1]
                  bd=readLines(bd)
                  bd=sub("^\t*","",bd)
                  
                  anc=bd[grepl(paste(bestModels[i,"fam"],"\t",sep=""),bd)]
                  anc=unlist(strsplit(anc,"\t"))[2]
                  anc=read.tree(text=anc)
                  anc=as.data.frame(ggtree(anc)$data)
                  anc.bd=sapply(anc$label,
                                function(i){
                                  if (grepl("_",i)){i=unlist(strsplit(i,"_"));i=i[length(i)]}
                                  return(as.numeric(i))
                                })
                  branchID.bd2branch_code.bd=paste("awk -F '\t' -v OFS='\t' '{if ($1==\"",
                                                   bestModels[i,"bd"]," \") print $3,$4}' ",
                                                   out_dir,"/results.",as.character(startValue.bd[1]),
                                                   "/",repli,".branch_code.txt",
                                                   sep="")
                  branchID.bd2branch_code.bd=system(branchID.bd2branch_code.bd,intern=TRUE)
                  branch_code.bd=sapply(df$branchID.bd,
                                        function(j){
                                          a=branchID.bd2branch_code.bd[grepl(j,branchID.bd2branch_code.bd)]
                                          if (length(a)==0){return(NA)}else{return(unlist(strsplit(a,"\t"))[2])}
                                        })
                  # branch_code.bd=unlist(unname(branch_code.bd))
                  rates=bd[which(bd=="#Branch_Group	Birth	Death	Innovation"):which(bd=="##Ancestral Family Size")]
                  rates=rates[c(-1,-length(rates),-length(rates)+1)]
                  Birth.bd=sapply(branch_code.bd,
                                  function(j){
                                    if (is.na(j)){return(NA)}else{
                                      a=unlist(strsplit(rates[grepl(paste("^",j,"\t",sep=""),rates)],"\t"))[2]
                                      return(as.numeric(a))
                                    }
                                  })
                  Death.bd=sapply(branch_code.bd,
                                  function(j){
                                    if (is.na(j)){return(NA)}else{
                                      a=unlist(strsplit(rates[grepl(paste("^",j,"\t",sep=""),rates)],"\t"))[3]
                                      return(as.numeric(a))
                                    }
                                  })
                  Innovation.bd=sapply(branch_code.bd,
                                       function(j){
                                         if (is.na(j)){return(NA)}else{
                                           a=unlist(strsplit(rates[grepl(paste("^",j,"\t",sep=""),rates)],"\t"))[4]
                                           return(as.numeric(a))
                                         }
                                       })
                  #Birth.bd=unname(Birth.bd)
                  #Death.bd=unname(Death.bd)
                  #Innovation.bd=unname(Innovation.bd)
                  
                  df=cbind(df,
                           fam=fam,branch_code.bd=branch_code.bd,
                           Birth.bd=Birth.bd,Death.bd=Death.bd,Innovation.bd=Innovation.bd,
                           anc.bd=anc.bd,branchModel.bd=branchModel.bd,rateModel.bd=rateModel.bd,
                           startValue.bd=startValue.bd)
                  return(df)
                })
  res=data.frame(parent=unlist(res["parent",]),
                 node=unlist(res["node",]),
                 branch.length=unlist(res["branch.length",]),
                 label=unlist(res["label",]),
                 isTip=unlist(res["isTip",]),
                 x=unlist(res["x",]),
                 y=unlist(res["y",]),
                 branch=unlist(res["branch",]),
                 angle=unlist(res["angle",]),
                 nodeID.bd=unlist(res["nodeID.bd",]),
                 parentID.bd=unlist(res["parentID.bd",]),
                 branchID.bd=unlist(res["branchID.bd",]),
                 fam=unlist(res["fam",]),
                 branch_code.bd=unlist(res["branch_code.bd",]),
                 Birth.bd=unlist(res["Birth.bd",]),
                 Death.bd=unlist(res["Death.bd",]),
                 Innovation.bd=unlist(res["Innovation.bd",]),
                 anc.bd=unlist(res["anc.bd",]),
                 branchModel.bd=unlist(res["branchModel.bd",]),
                 rateModel.bd=unlist(res["rateModel.bd",]),
                 startValue.bd=unlist(res["startValue.bd",]))
  write.table(res,paste(out_dir,"/badirate_results.tsv",sep=""),
              sep="\t",row.names = FALSE,quote = FALSE)
}




















# # HGTphyloDetect: detect HGTs
# # Super slow
# # Dependencies: HGTphyloDetect (modified), DIAMOND, parallel (R), seqkit
# HGTphyloDetect=function(pep.faa=pep.faa,
#                         pep.kingdom="33208",# txid. 33208: Metazoa
#                         pep.subphylum="6960", # txid. 6960: Hexapoda
#                         nr.dmdb="/apps/unit/BioinfoUgrp/DB/diamondDB/ncbi/2022-07/nr.dmnd",
#                         HGT_workflow_distant.py="/home/c/c-liu/Softwares/HGTphyloDetect/main/HGT_workflow_distant.py",
#                         #HGT_workflow_close.py="/home/c/c-liu/Softwares/HGTphyloDetect/main/HGT_workflow_close.py",
#                         out_dir=out_dir,
#                         out_basename=out_basename,
#                         threads=threads){
#   threads=as.character(threads)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
#   out_dir=sub("/$","",out_dir)
#   wd=getwd();setwd(out_dir)
#   
#   # blast search
#   if (!file.exists("blast.tsv")){
#     cmd=paste("diamond blastp",
#               "--threads",threads,
#               "--db",nr.dmdb,
#               "--query",pep.faa,
#               "--out blast.tsv",
#               "--max-target-seqs 400",
#               "--outfmt 6 qseqid sseqid evalue bitscore length pident",
#               sep=" ")
#     print(cmd);system(cmd,wait=TRUE)
#   }
#   
#   system("mkdir parallel/")
#   pepID.lst=system(paste("grep '>' ",pep.faa," | sed 's/>//'",sep=""),intern=TRUE)
#   library(parallel)
#   clus=makeCluster(as.numeric(threads))
#   parSapply(clus,pepID.lst,
#             function(ID){ # split blast results & run HGTphyloDetect
#               cmd0=paste("mkdir",
#                          paste("parallel/",ID,sep=""),
#                          sep=" ")
#               system(cmd0)
#               cmd1=paste("mkdir",
#                          paste("parallel/",ID,"/blastp_files",sep=""),
#                          sep=" ")
#               system(cmd1)
#               cmd2=paste("seqkit grep",
#                          "-n","-p",ID,
#                          pep.faa,
#                          ">",paste("parallel/",ID,"/",ID,".faa",sep=""),
#                          sep=" ")
#               system(cmd2)
#               cmd3=paste("echo '# Fields: qseqid sseqid evalue bitscore length pident' > ",
#                          "parallel/",ID,"/blastp_files/",ID,".txt",
#                          sep="")
#               system(cmd3)
#               cmd4=paste("awk -F '\t' -v OFS='\t' '{if ($1==\"",ID,"\") print $0}' ",
#                          "blast.tsv >> ",
#                          "parallel/",ID,"/blastp_files/",ID,".txt",
#                          sep="")
#               system(cmd4)
#               # cmd5=paste("cd",
#               #            paste("parallel/",ID,"/",sep=""),
#               #            sep=" ")
#               # system(cmd5)
#               setwd(paste("parallel/",ID,"/",sep=""))
#               cmd6=paste("python",
#                          HGT_workflow_distant.py,
#                          paste(ID,".faa",sep=""),
#                          "AI=45 out_pct=0.9",
#                          paste("gene_kingdom=",as.character(pep.kingdom),sep=""),
#                          paste("gene_subphylum=",as.character(pep.subphylum),sep=""),
#                          sep=" ")
#               system(cmd6)
#               # cmd7=paste("python",
#               #            HGT_workflow_close.py,
#               #            paste(ID,".faa",sep=""),
#               #            "bitscore_parameter=100 HGTIndex=0.5 out_pct=0.8",
#               #            paste("gene_kingdom=",as.character(pep.kingdom),sep=""),
#               #            paste("gene_subphylum=",as.character(pep.subphylum),sep=""),
#               #            sep=" ")
#               # cmd8="cd ../../"
#               # system(cmd8)
#               setwd("../../")
#               #cmd=paste(cmd0,cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,cmd7,cmd8,sep="; ")
#               #return(cmd)
#             })
#   #jobs=unname(jobs)
#   #writeLines(jobs,"jobs")
#   
#   # cmd=paste("split",
#   #           "-d",
#   #           "-n",paste("l/",as.character(as.numeric(threads)-1),sep=""),
#   #           "jobs",
#   #           sep=" ")
#   # print(cmd);system(cmd,wait=TRUE)
#   # library(parallel)
#   # clus=makeCluster(as.numeric(threads))
#   # parSapply(clus,
#   #           1:(as.numeric(threads)-1),
#   #           function(i){
#   #             scr=system("ls x*",wait=TRUE,intern=TRUE)[i]
#   #             cmd=paste("bash"," ",scr,sep="")
#   #             print(cmd);system(cmd,wait=TRUE)})
#   # stopCluster(clus)
#   
#   # output_close=paste("parallel/",pepID.lst,"/output_close_HGT.tsv",sep="")
#   # system(paste("head -n1 ",output_close[1]," > ",out_basename,"_closeHGT.tsv",sep=""))
#   # for (i in output_close){
#   #   wcl=as.numeric(unlist(strsplit(system(paste("wc -l ",i,sep=""),intern=TRUE)," "))[1])
#   #   if (wcl>1){
#   #     system(paste("tail -n1 ",i," >> ",out_basename,"_closeHGT.tsv",sep=""))
#   #   }
#   # }
#   
#   output_distant=paste("parallel/",pepID.lst,"/output_distant_HGT.tsv",sep="")
#   system(paste("head -n1 ",output_distant[1]," > ",out_basename,"_distantHGT.tsv",sep=""))
#   for (i in output_distant){
#     wcl=as.numeric(unlist(strsplit(system(paste("wc -l ",i,sep=""),intern=TRUE)," "))[1])
#     if (wcl>1){
#       system(paste("tail -n1 ",i," >> ",out_basename,"_distantHGT.tsv",sep=""))
#     }
#   }
#   
#   setwd(wd)
# }
# # Incomplete
# # AvP
# # Detect horizontal gene transfer (HGT)
# # Dependencies: AvP
# Avp=function(query.faa=query.faa,
#              blast_filename=blast_filename,
#              dmbd=dmbd, # not used if blast_filename exists
#              group.yaml=group.yaml,
#              config.yaml=config.yaml,
#              out_dir=out_dir,
#              threads=threads
# ){
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   threads=as.character(threads)
#   wd=getwd()
#   setwd(out_dir)
#   
#   if (!file.exists(blast_filename)){
#     cmd=paste("diamond blastp",
#               "--threads",threads,
#               "-d",dmbd,
#               "--max-target-seqs 500",
#               "-q",query.faa,
#               "--evalue 1e-5",
#               "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids",
#               ">",
#               blast_filename)
#     print(cmd);system(cmd,wait=TRUE)
#   }
#   
#   cmd=paste("calculate_ai.py",
#             "-i",blast_filename,
#             "-x",group.yaml,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("avp prepare",
#             "-a",paste(blast_filename,"_ai.out",sep=""),
#             "-o",".",
#             "-f",query.faa,
#             "-b",blast_filename,
#             "-x",group.yaml,
#             "-c",config.yaml,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("avp detect",
#             "-i","./mafftgroups/",
#             "-o",".",
#             "-g","./groups.tsv",
#             "-t","./tmp/taxonomy_nexus.txt",
#             "-c",config.yaml,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   setwd(wd)
# }

# # HGT inference by DIAMOND-nr-MEGAN
# hgt=function(proteins.faa=proteins.faa,
#              ref_diamond=ref_diamond,
#              ref_megan=ref_megan,
#              taxonExclude=taxonExclude, # exclude list of taxon ids (comma-separated)
#              out_prefix=out_prefix,
#              threads=threads){
#   threads=as.character(threads)
#   
#   cmd=paste("diamond blastp",
#             "--outfmt 100",
#             "-p",threads,
#             "-d",ref_diamond,
#             "-q",proteins.faa,
#             "--taxon-exclude",taxonExclude,
#             "--out",paste(out_prefix,".daa",sep=""),
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("daa-meganizer",
#             "-i",paste(out_prefix,".daa",sep=""),
#             "-mdb",ref_megan,
#             "-t",threads,
#             "-ram readCount",
#             "-supp 0",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("daa2info",
#             "-i",paste(out_prefix,".daa",sep=""),
#             "-o",paste(out_prefix,"_taxon.tsv",sep=""),
#             "-r2c Taxonomy",
#             "-n true",
#             "-p true",
#             "-r true",
#             "-mro true",
#             "-u false",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
# }