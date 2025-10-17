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
# self-self blastp
self2self_blastp=function(protein.faa=protein.faa,
                          threads=threads,
                          out_prefix=out_prefix){
  cmd=paste("makeblastdb",
            "-dbtype prot",
            "-in",protein.faa,
            "-out",paste(out_prefix,".blastp.db",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("blastp",
            "-query",protein.faa,
            "-db",paste(out_prefix,".blastp.db",sep=""),
            "-out",paste(out_prefix,".self2self.blast",sep=""),
            "-outfmt 6",
            "-evalue 1e-5",
            "-num_alignments 20",
            "-num_threads",as.character(threads),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

self2self_diamond=function(protein.faa=protein.faa,
                          threads=threads,
                          out_prefix=out_prefix){
  cmd=paste("diamond makedb",
            "--db",paste(out_prefix,".dmdb",sep=""),
            "--in",protein.faa,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("diamond blastp",
            "--query",protein.faa,
            "--db",paste(out_prefix,".dmdb",sep=""),
            "--threads",as.character(threads),
            "--out",paste(out_prefix,".self2self.blast",sep=""),
            "--evalue 1e-5",
            "--outfmt 6",
            "--max-target-seqs 20",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# MCScanX: classify duplicated gene pairs
# Dependencies: blastp, MCScanX, stringr (R)
MCScanX_dupliClass=function(gff=gff,
                            blast=blast, # self2self blast of proteins; IDs same with gene ID in gff
                            out_dir=out_dir,
                            out_basename=out_basename){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp ",blast," ./",out_basename,".blast",sep="")
  print(cmd);system(cmd,wait=TRUE)
  
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
                    out_basename1=out_basename1,
                    out_basename2=out_basename2 # not used for intragenomic comparisons
                    ){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  self2self=genome.fna1==genome.fna2
  
  
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
  
  # Prepare fake gff and len file required by wgdi
  wgdi_input(genome.fna=genome.fna1,
             gff=gff1,
             out_prefix=out_basename1)
  fake_gff1=paste(out_dir,"/",out_basename1,".gff",sep="")
  len1=paste(out_dir,"/",out_basename1,".len",sep="")
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
               out_prefix=out_basename2)
    fake_gff2=paste(out_dir,"/",out_basename2,".gff",sep="")
    len2=paste(out_dir,"/",out_basename2,".len",sep="")
    # Remove contigs without genes
    cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3!=0) print $0}'",len2,"> len2",sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system(paste("rm",len2,sep=" "));system(paste("mv","len2",len2,sep=" "))
    # Format gff
    cmd=paste("awk -F '\t' -v OFS='\t' '{print $1,$7,$3,$4,$5,$6,$2}'",fake_gff2,"> fake_gff2",sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system(paste("rm",fake_gff2,sep=" "));system(paste("mv","fake_gff2",fake_gff2,sep=" "))
  }else{
    fake_gff2=paste(out_dir,"/",out_basename1,".gff",sep="")
    len2=paste(out_dir,"/",out_basename1,".len",sep="")
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
  if (!file.exists(paste("../",out_basename1,".blast",sep=""))){
    cmd=paste("diamond makedb",
              "--db",paste(proteins.faa2,".dmdb",sep=""),
              "--in",proteins.faa2,
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("diamond blastp",
              "--threads",as.character(threads),
              "--db",paste(proteins.faa2,".dmdb",sep=""),
              "--query",proteins.faa1,
              "--outfmt 6",
              "--evalue 1e-5",
              "--max-target-seqs 20",
              "--out",paste(out_basename1,".blast",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    blast1=paste(out_basename1,".blast",sep="")
    system(paste("mv",blast1,"../",sep=" "),wait=TRUE)
  }
  if (!self2self){
    if (!file.exists(paste("../",out_basename2,".blast",sep=""))){
      cmd=paste("diamond makedb",
                "--db",paste(proteins.faa1,".dmdb",sep=""),
                "--in",proteins.faa1,
                sep=" ")
      print(cmd);system(cmd,wait=TRUE)
      cmd=paste("diamond blastp",
                "--threads",as.character(threads),
                "--db",paste(proteins.faa1,".dmdb",sep=""),
                "--query",proteins.faa2,
                "--outfmt 6",
                "--evalue 1e-5",
                "--max-target-seqs 20",
                "--out",paste(out_basename2,".blast",sep=""),
                sep=" ")
      print(cmd);system(cmd,wait=TRUE)
      blast2=paste(out_basename2,".blast",sep="")
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
                 p_block=0.2, # 0-1,
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

lst=readLines("/flash/BourguignonU/Cong/Luan/pal/pal.lst")
for (g in lst){
  pal.faa=paste("/flash/BourguignonU/Cong/Luan/pal/findHomo/",
                g,".faa",sep="")
  cmd=paste("seqkit grep -n -p ",g," /flash/BourguignonU/Cong/Luan/pal/pal.faa",
            " > ",pal.faa,sep="")
  print(cmd);system(cmd,wait=T)
  findHomo(pep.faa=pal.faa,
           nr.dmdb="/apps/unit/BioinfoUgrp/DB/diamondDB/ncbi/2022-07/nr.dmnd",
           threads=20,
           out_prefix=paste("/flash/BourguignonU/Cong/Luan/pal/findHomo/",
                            g,sep=""))
  mafft(in.fa=paste("/flash/BourguignonU/Cong/Luan/pal/findHomo/",
                    g,"_nr.faa",sep=""),
        align.fa=paste("/flash/BourguignonU/Cong/Luan/pal/findHomo/",
                       g,"_mafft.faa",sep=""),
        threads=20)
  trimAL(inMSA.fa=paste("/flash/BourguignonU/Cong/Luan/pal/findHomo/",
                                 g,"_mafft.faa",sep=""),
         outMSA.fa=paste("/flash/BourguignonU/Cong/Luan/pal/findHomo/",
                                  g,"_trimal.faa",sep=""))
        
  iqtree(msa.fa=paste("/flash/BourguignonU/Cong/Luan/pal/findHomo/",
                      g,"_trimal.faa",sep=""),
         type="protein", # dna/protein
         out_prefix=paste("/flash/BourguignonU/Cong/Luan/pal/findHomo/",
                          g,"_iqtree.faa",sep=""),
         threads=10)
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
              "-format fasta_cds_na","|",
              "seqkit grep -r",
              "-p",prot,
              ">>", 
              paste(out_prefix,"_cds.fna",sep=""),
              sep=" ")
    print(cmd)
    system(cmd,wait=TRUE)
    Sys.sleep(5)
    
    cmd=paste("elink",
              "-target nuccore",
              "-db protein",
              "-name protein_nuccore",
              "-id",prot,"|",
              "efetch",
              "-format fasta_cds_aa","|",
              "seqkit grep -r",
              "-p",prot,
              ">>", 
              paste(out_prefix,"_prot.faa",sep=""),
              sep=" ")
    print(cmd)
    system(cmd,wait=TRUE)
    Sys.sleep(5)
  }
  
  library(Biostrings)
  proteins.faa=readAAStringSet(paste(out_prefix,"_prot.faa",sep=""))
  proteins.IDs=names(proteins.faa)
  proteins.IDs=gsub(pattern=".*\\[protein_id=(.*?)\\].*",replacement = "\\1",proteins.IDs)
  names(proteins.faa)=proteins.IDs
  proteins.faa=proteins.faa[names(proteins.faa) %in% protID.lst]
  proteins.faa=proteins.faa[!duplicated(names(proteins.faa))]
  proteins.IDs=names(proteins.faa)
  print("NO protein sequence found for:")
  print(protID.lst[!(protID.lst %in% proteins.IDs)])
  #[1] "EDU0222588.1"   "EFG6479230.1"   "EDR4320347.1"   "MBS1204601.1"  
  #[5] "ONI56852.1"     "KZV81323.1"     "Q81835.2"       "CAQ16927.1"    
  #[9] "WP_252153986.1" "AFM74331.1"     "GEX71379.1"     "MCH83329.1"    
  #[13] "WP_252154240.1" "MBV0899642.1"   "MCL2453738.1"   "HJD67128.1"    
  #[17] "CJV99103.1"     "MBV2146500.1"   "MCM1001045.1"   "MBV0899429.1"  
  #[21] "MCM1001449.1"   "MCL2287834.1"   "MCM1002635.1"   "MCM1001130.1"   
  writeXStringSet(proteins.faa,paste(out_prefix,"_prot.faa",sep=""))
  
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
  #[1] "EDU0222588.1"   "EFG6479230.1"   "EDR4320347.1"   "MBS1204601.1"  
  #[5] "ONI56852.1"     "WP_169807269.1"        "WP_025263969.1"
  #[9] "AFM74331.1"     "GEX71379.1"     "MCH83329.1"     "WP_063630740.1"
  #[13] "MBV0899642.1"   "MCL2453738.1"   "HJD67128.1"     "MBV2146500.1"  
  #[17] "MCM1001045.1"   "MBV0899429.1"   "MCM1001449.1"   "MCL2287834.1"  
  #[21] "MCM1002635.1"   "WP_180807986.1" "MCM1001130.1"
  rownames(df)=df$cds.ID
  cds.fna=cds.fna[names(cds.fna) %in% df$cds.ID]
  cds.fna=cds.fna[!duplicated(names(cds.fna))]
  names(cds.fna)=paste(df[names(cds.fna),"protein.ID"],df[names(cds.fna),"genome.ID"],sep=" ")
  writeXStringSet(cds.fna,paste(out_prefix,"_cds.fna",sep=""))
  
  # system(paste("rm ",out_prefix,"_cds.fna",sep=""))
  # 
  # genome.ID=df$genome.ID
  # genome.ID=genome.ID[!duplicated(genome.ID)]
  # library(taxonomizr)
  # library(parallel)
  # clus=makeCluster(as.numeric(threads))
  # clusterExport(cl=clus,varlist = list("genome.ID",
  #                                      "accessionToTaxa","getRawTaxonomy","taxonomizr.sql"))
  # txid=parSapply(clus,genome.ID,function(i){return(accessionToTaxa(i,taxonomizr.sql))}) 
  # names(txid)=genome.ID
  # taxaPath=parSapply(clus,txid,function(i){i=getRawTaxonomy(i,taxonomizr.sql);return(paste(i[[1]],collapse=","))})
  # names(taxaPath)=genome.ID
  # 
  # df$txid=txid[df$genome.ID]
  # df$taxaPath=taxaPath[df$genome.ID]
  # write.table(df,paste(out_prefix,"_taxon.tsv",sep=""),
  #             sep="\t",row.names = FALSE,quote = FALSE)
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

# MAFFT in R
# Dependencies: ips (R)
mafftr=function(seqs=seqs # An object of class DNAbin or AAbin.
                ){
  res=ips::mafft(x=seqs,
             method="auto",
             thread=1)
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

# Extract 4dtv alignment
# dependencies: ape (R)
get4dtvAlign=function(in.fa=in.fa,
                      out.fa=out.fa){
  inDNA=ape::read.dna(file=in.fa,format="fasta",as.matrix=TRUE,as.character=TRUE)
  firstP=inDNA[,seq(1,ncol(inDNA),3)]
  secondP=inDNA[,seq(2,ncol(inDNA),3)]
  thirdP=inDNA[,seq(3,ncol(inDNA),3)]
  
  dtv4=sapply(seq(1,ncol(inDNA)/3,1),
              function(i){
                r=FALSE
                if (all(firstP[,i]=="c") & all(secondP[,i]=="t")){r=TRUE}
                if (all(firstP[,i]=="g") & all(secondP[,i]=="t")){r=TRUE}
                if (all(firstP[,i]=="t") & all(secondP[,i]=="c")){r=TRUE}
                if (all(firstP[,i]=="c") & all(secondP[,i]=="c")){r=TRUE}
                if (all(firstP[,i]=="a") & all(secondP[,i]=="c")){r=TRUE}
                if (all(firstP[,i]=="g") & all(secondP[,i]=="c")){r=TRUE}
                if (all(firstP[,i]=="c") & all(secondP[,i]=="g")){r=TRUE}
                if (all(firstP[,i]=="g") & all(secondP[,i]=="g")){r=TRUE}
                return(r)
              })
  dtv4.dna=thirdP[,dtv4]
  dtv4.dna=toupper(dtv4.dna)
  ape::write.dna(x=dtv4.dna,file=out.fa,format="fasta",colsep="",indent=0,blocksep=0)
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

# setwd("/flash/BourguignonU/Cong/termite_pca/geneTrees_Yiming/notung/nodeInfo/")
# HOG_tree=readLines("/bucket/BourguignonU/Cong/termite_pca/geneTrees_Yiming/HOG_treeBuilt.lst")
# nodeInfoFile=paste("/flash/BourguignonU/Cong/termite_pca/geneTrees_Yiming/notung/nodeInfo/",
#       HOG_tree,"_nodeInfo.tsv",sep="")
# geneTreeNode2TermiteSp=expand.grid(1:45,1:45)
# geneTreeNode2TermiteSp=unique(t(apply(geneTreeNode2TermiteSp, 1, sort)))
# geneTreeNode2TermiteSp=as.data.frame(geneTreeNode2TermiteSp)
# colnames(geneTreeNode2TermiteSp)=c("n1_termiteSp","n2_termiteSp")
# f=function(n){
#   i=geneTreeNode2TermiteSp[n,1]
#   j=geneTreeNode2TermiteSp[n,2]
#   cmd1=paste("awk '{if ($18==",as.character(i)," && $19==",as.character(j),") print $0}' *_nodeInfo.tsv | wc -l",sep="")
#   cmd2=paste("awk '{if ($18==",as.character(j)," && $19==",as.character(i),") print $0}' *_nodeInfo.tsv | wc -l",sep="")
#   
#   n1=as.numeric(system(cmd1,intern=TRUE))
#   n2=as.numeric(system(cmd2,intern=TRUE))
#   
#   if (i==j){res=n1}else{res=n1+n2}
#   return(res)
# }
# library(parallel)
# clus=makeCluster(10)
# clusterExport(clus,list("geneTreeNode2TermiteSp"))
# geneTreeNode2TermiteSp$n_geneTreeNode=parSapply(clus,1:nrow(geneTreeNode2TermiteSp),f)
# write.table(geneTreeNode2TermiteSp,
#             "/flash/BourguignonU/Cong/termite_pca/geneTrees_Yiming/notung/nodeInfo/termiteSpInGeneTreeNode.tsv",
#             sep="\t",row.names=FALSE,quote=FALSE)


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
# ParaAT and KaKs_calculator
# Compute pairwise codon alignments for input sequences, and pairwise KaKs
# Dependencies: MAFFT, ParaAT, Kaks_calculator
# Ks>5 means genetic saturation? 
# https://www.nature.com/articles/s42003-023-05044-1
# 
ParaAT=function(protein.faa=protein.faa,
                cds.fna=cds.fna,
                threads=threads,
                geneticCode=1,
                out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd()
  setwd(out_dir)
  
  # protein.faa="/bucket/BourguignonU/Cong/termite_pca/seqs.trimmed_HOG/proteins/HOG0000080.faa"
  # cds.fna="/bucket/BourguignonU/Cong/termite_pca/seqs.trimmed_HOG/cds/HOG0000080.fna"
  # out_dir="."
  # threads=1
  # geneticCode=1
  
  cmd=paste("grep","'>'",protein.faa,"|",
            "sed 's/>//'",
            sep=" ")
  homolog.lst=system(cmd,intern=TRUE)
  df=expand.grid(homolog.lst,homolog.lst,stringsAsFactors=F)
  df=unique(t(apply(df, 1, sort)))
  df=as.data.frame(df)
  df=df[df$V1!=df$V2,]
  write.table(df,paste(out_dir,"/",basename(protein.faa),".homolog",sep=""),
              sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  write(as.character(threads),paste(out_dir,"/",basename(protein.faa),".threads",sep=""))
  
  cmd=paste("ParaAT.pl",
            "-homolog",paste("./",basename(protein.faa),".homolog",sep=""),
            "-aminoacid",protein.faa,
            "-nuc",cds.fna,
            "-processor",paste("./",basename(protein.faa),".threads",sep=""),
            "-output",paste("./",basename(protein.faa),sep=""),
            "-code",as.character(geneticCode),
            "-format axt -kaks -m mafft",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  res=read.table(paste(out_dir,"/",basename(protein.faa),".homolog",sep=""),
                 sep="\t",header=FALSE,quote="")
  colnames(res)=c("Gene.1","Gene.2")
  res[,c("Ka","Ks","KaABOVEKs","p.null_KaEqualKS","Length","S.sites","N.sites","Substitutions","Syn.subs","Nonsyn.subs",
         "Divergence.Distance","ML.Score")]=
  t(sapply(1:nrow(res),
           function(i){
             kaks=paste(out_dir,"/",basename(protein.faa),
                        "/",res[i,"Gene.1"],"-",res[i,"Gene.2"],".cds_aln.axt.kaks",sep="")
             kaks=read.table(kaks,sep="\t",header=TRUE,quote="")
             return(
             c(kaks[1,"Ka"],kaks[1,"Ks"],kaks[1,"Ka.Ks"],
             kaks[1,"P.Value.Fisher."],kaks[1,"Length"],kaks[1,"S.Sites"],kaks[1,"N.Sites"],kaks[1,"Substitutions"],
             kaks[1,"Syn.Subs"],kaks[1,"Nonsyn.Subs"],kaks[1,"Divergence.Distance"],kaks[1,"ML.Score"])
             )
           }))
  write.table(res,paste(out_dir,"/",basename(protein.faa),".pairwise.dNdS",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
  
  system(paste("rm ",out_dir,"/",basename(protein.faa),".homolog",sep=""))
  system(paste("rm ",out_dir,"/",basename(protein.faa),".threads",sep=""))
  system(paste("rm -r ",out_dir,"/",basename(protein.faa),sep=""))
  setwd(wd)
}

# Extract a pair of sequences, align them (Protein, DNA, Codon), and compute pairwise identity, similarity, dN, dS
# Same logic with ParaAT
# Dependencies: seqkit, mafft, pal2nal, seqinr (R), parseFastaIntoAXT.pl, KaKs_calculator
# seqinr::dist.alignment: the "identity" and "similarity" are actually DISTANCE, and the "similarity" only makes senes for peps (A,T,G,C are also amino acids)
pairwise=function(protein.faa=protein.faa,
                  cds.fna=cds.fna,
                  geneticCode=geneticCode,
                  out.tsv=out.tsv){
  gene.lst=system(paste("grep '>' ",protein.faa,sep=""),intern=TRUE)
  gene.lst=sub(">","",gene.lst)
  
  df=expand.grid(gene.lst,gene.lst)
  colnames(df)=c("Gene.1","Gene.2")
  df=df[df$Gene.1!=df$Gene.2,]
  df=as.data.frame(t(apply(df[,c(1,2)], 1, sort)))
  df=df[!duplicated(df),]
  colnames(df)=c("Gene.1","Gene.2")
  
  sapply(1:nrow(df),
         function(i){
           for (gene in df[i,]){
             cmd=paste("seqkit grep -p ",gene," ",protein.faa,
                       " >> ",out.tsv,".tmp.pep.",as.character(i),sep="")
             system(cmd,wait=TRUE)
             cmd=paste("seqkit grep -p ",gene," ",cds.fna,
                       " >> ",out.tsv,".tmp.cds.",as.character(i),sep="")
             system(cmd,wait=TRUE)
           }
           cmd=paste("mafft --auto --thread 1 ",
                     out.tsv,".tmp.pep.",as.character(i)," > ",
                     out.tsv,".tmp.mafft.",as.character(i),
                     sep="")
           system(cmd,wait=TRUE)
           cmd=paste("mafft --auto --thread 1 ",
                     out.tsv,".tmp.cds.",as.character(i)," > ",
                     out.tsv,".tmp.mafft.dna.",as.character(i),
                     sep="")
           system(cmd,wait=TRUE)
           cmd=paste("pal2nal.pl ",
                     out.tsv,".tmp.mafft.",as.character(i)," ",
                     out.tsv,".tmp.cds.",as.character(i)," ",
                     "-output fasta ",
                     "-codontable ",as.character(geneticCode),
                     " > ",out.tsv,".tmp.pal2nal.",as.character(i),
                     sep="")
           system(cmd,wait=TRUE)
           cmd=paste("parseFastaIntoAXT.pl ",
                     out.tsv,".tmp.pal2nal.",as.character(i),
                     sep="")
           system(cmd,wait=TRUE)
           cmd=paste("KaKs -i ",out.tsv,".tmp.pal2nal.",as.character(i),".axt ",
                     "-o ",out.tsv,".tmp.pal2nal.",as.character(i),".axt.kaks ",
                     "-c ",as.character(geneticCode),sep="")
           system(cmd,wait=TRUE)
         })
  
  df$identity.pep=rep(NA,nrow(df));df$identity.codon=rep(NA,nrow(df))
  df$identity.dna=rep(NA,nrow(df))
  df$similarity.pep=rep(NA,nrow(df));df$similarity.codon=rep(NA,nrow(df))
  df$similarity.dna=rep(NA,nrow(df))
  df$dN=rep(NA,nrow(df));df$dS=rep(NA,nrow(df));df$omega=rep(NA,nrow(df))
  df$N.sites=rep(NA,nrow(df));df$S.sites=rep(NA,nrow(df))
  df$p.null_KaEqualKS=rep(NA,nrow(df));df$divergence.distance=rep(NA,nrow(df))
  library(seqinr)
  for (i in 1:nrow(df)){
    pep.align=read.alignment(paste(out.tsv,".tmp.mafft.",as.character(i),sep=""),
                             format="fasta")
    cds.align=read.alignment(paste(out.tsv,".tmp.pal2nal.",as.character(i),sep=""),
                             format="fasta")
    dna.align=read.alignment(paste(out.tsv,".tmp.mafft.dna.",as.character(i),sep=""),
                             format="fasta")
    KAKS=read.table(paste(out.tsv,".tmp.pal2nal.",as.character(i),".axt.kaks",sep=""),
                    sep="\t",header=TRUE,quote="")
    df[i,"identity.pep"]=as.vector(dist.alignment(x=pep.align,matrix="identity"))
    df[i,"identity.codon"]=as.vector(dist.alignment(x=cds.align,matrix="identity"))
    df[i,"identity.dna"]=as.vector(dist.alignment(x=dna.align,matrix="identity"))
    df[i,"similarity.pep"]=as.vector(dist.alignment(x=pep.align,matrix="similarity"))
    df[i,"similarity.codon"]=as.vector(dist.alignment(x=cds.align,matrix="similarity"))
    df[i,"similarity.dna"]=as.vector(dist.alignment(x=dna.align,matrix="similarity"))
    df[i,"dN"]=KAKS[1,"Ka"]
    df[i,"dS"]=KAKS[1,"Ks"]
    df[i,"omega"]=KAKS[1,"Ka.Ks"]
    df[i,"p.null_KaEqualKS"]=KAKS[1,"P.Value.Fisher."]
    df[i,"divergence.distance"]=KAKS[1,"Divergence.Distance"]
    df[i,"N.sites"]=KAKS[1,"N.Sites"]
    df[i,"S.sites"]=KAKS[1,"S.Sites"]
  }
  
  for (i in 1:nrow(df)){
    system(paste("rm ",out.tsv,".tmp.pep.",as.character(i),sep=""))
    system(paste("rm ",out.tsv,".tmp.cds.",as.character(i),sep=""))
    system(paste("rm ",out.tsv,".tmp.mafft.",as.character(i),sep=""))
    system(paste("rm ",out.tsv,".tmp.mafft.dna.",as.character(i),sep=""))
    system(paste("rm ",out.tsv,".tmp.pal2nal.",as.character(i),sep=""))
    #system(paste("rm ",out.tsv,".tmp.pal2nal.",as.character(i),".axt",sep=""))
    #system(paste("rm ",out.tsv,".tmp.pal2nal.",as.character(i),".axt.kaks",sep=""))
  }
  
  write.table(df,out.tsv,sep="\t",row.names = FALSE,quote = FALSE)
}

# Compute pairwise Ka, Ks with seqinr::kaks, using Li (1993)
# Multiple substitution fixed with K80
kaks.seqinr=function(codon.align=codon.align,
                     out.tsv=out.tsv){
  codon.align=seqinr::read.alignment(codon.align,format="fasta")
  K=seqinr::kaks(codon.align,forceUpperCase=TRUE,rmgap=FALSE)
  res=expand.grid(codon.align$nam,codon.align$nam)
  res=as.data.frame(unique(t(apply(res,1,sort))))
  colnames(res)=c("seq1","seq2")
  res=res[res$seq1!=res$seq2,]
  ks.mat=as.matrix(K$ks)
  res$Ks=sapply(1:nrow(res),
                function(i){
                  return( ks.mat[res[i,"seq1"],res[i,"seq2"]] )
                })
  ka.mat=as.matrix(K$ka)
  res$Ka=sapply(1:nrow(res),
                function(i){
                  return( ka.mat[res[i,"seq1"],res[i,"seq2"]] )
                })
  vks.mat=as.matrix(K$vks)
  res$variance.Ks=sapply(1:nrow(res),
                         function(i){
                           return( vks.mat[res[i,"seq1"],res[i,"seq2"]] )
                         })
  vka.mat=as.matrix(K$vka)
  res$variance.Ka=sapply(1:nrow(res),
                         function(i){
                           return( vka.mat[res[i,"seq1"],res[i,"seq2"]] )
                         })
  write.table(res,out.tsv,
              row.names=FALSE,quote=FALSE,sep="\t")
}

# Multiple sequence alignment -> pairwise alignment ->
# 4-fold degenerate synonymous sites of the third codons (4dtv) alignment ->
# alignment length, K80 divergence and raw divergence
# Dependencies: seqkit, ape (R)
distFrom4dtv=function(in.fa=in.fa, # codon alignment
                      out.tsv=out.tsv){
  seqID=system(paste("grep '>' ",in.fa,sep=""),intern=TRUE)
  seqID=sub("^>","",seqID)
  res=expand.grid(seqID,seqID)
  res=as.data.frame(unique(t(apply(res,1,sort))))
  colnames(res)=c("seq1","seq2")
  res=res[res$seq1!=res$seq2,]
  
  inDNA=ape::read.dna(file=in.fa,format="fasta",as.matrix=TRUE,as.character=TRUE)
  
  res[,c("length.4dtv","divergence","rawDiv")]=
    t(sapply(1:nrow(res),
                function(i){ # 1:nrow(res)
                  d=inDNA[c(res[i,"seq1"],res[i,"seq2"]),]
                  firstP=d[,seq(1,ncol(d),3)]
                  secondP=d[,seq(2,ncol(d),3)]
                  thirdP=d[,seq(3,ncol(d),3)]
                  
                  # 4d codon sites in first seq
                  dtv4.1=sapply(seq(1,ncol(d)/3,1),
                              function(i){
                                p1=firstP[1,i]
                                p2=secondP[1,i]
                                dtv=FALSE
                                if (p1%in%c("a","t") & p2=="c"){dtv=TRUE} # acN,tcN
                                if (p1=="c" & p2%in%c("t","c","g")){dtv=TRUE} # ctN,ccN,cgN
                                if (p1=="g" & p2%in%c("t","c","g")){dtv=TRUE} # gtN,gcN,ggN
                                return(dtv)
                              })
                  dtv4.1=which(dtv4.1)
                  # 4d codon sites in second seq
                  dtv4.2=sapply(seq(1,ncol(d)/3,1),
                                function(i){
                                  p1=firstP[2,i]
                                  p2=secondP[2,i]
                                  dtv=FALSE
                                  if (p1%in%c("a","t") & p2=="c"){dtv=TRUE} # acN,tcN
                                  if (p1=="c" & p2%in%c("t","c","g")){dtv=TRUE} # ctN,ccN,cgN
                                  if (p1=="g" & p2%in%c("t","c","g")){dtv=TRUE} # gtN,gcN,ggN
                                  return(dtv)
                                })
                  dtv4.2=which(dtv4.2)
                  
                  # 4d codon position shared by two seq
                  dtv4=intersect(dtv4.1,dtv4.2)
                  if (length(dtv4)<2){return(c(NA,NA,NA))}else{
                    dtv4.1=firstP[,dtv4]
                    dtv4.2=secondP[,dtv4]
                    dtv4.3=thirdP[,dtv4]
                    
                    sites=intersect(which(dtv4.1[1,]==dtv4.1[2,]),
                                    which(dtv4.2[1,]==dtv4.2[2,]))
                    if (length(sites)<20){return(c(NA,NA,NA))}else{
                      x=dtv4.3[,]
                      length.4dtv=ncol(x)
                      
                      divergence=ape::dist.dna(ape::as.DNAbin.character(x),
                                               model="K80",as.matrix=TRUE,
                                               pairwise.deletion=TRUE)
                      divergence=divergence[1,2]
                      if (is.infinite(divergence)){divergence=NA}
                      
                      rawDiv=ape::dist.dna(ape::as.DNAbin.character(x),
                                           model="raw",as.matrix=TRUE,
                                           pairwise.deletion=TRUE)
                      rawDiv=rawDiv[1,2]
                      return(c(length.4dtv,divergence,rawDiv))
                    }
                    
                  }
                }))
  write.table(res,out.tsv,sep="\t",row.names=FALSE,quote=FALSE)
}

# absrel of hyphy: estimate omega and test positive selection
# Dependencies: hyphy
absrel=function(hyphy.res="/apps/unit/BourguignonU/hyphy/res",
                codon.align=codon.align,
                geneTree=geneTree,
                out_prefix=out_prefix){
  cmd=paste("hyphy",
            paste("LIBPATH=",hyphy.res,sep=""),
            "absrel",
            "--alignment",codon.align,
            "--tree",geneTree,
            "--output",paste(out_prefix,".json",sep=""))
  print(cmd);system(cmd,wait=TRUE)
}

# Format json from absrel
# Dependencies: jsonlite (R), treeio (R), ggtree (R)
json2tsv_absrel=function(json=json,
                         tsv=tsv){
  json=jsonlite::read_json(json)
  
  tree=json$input$trees[[1]]
  tree=paste(tree,";",sep="")
  tree=ape::read.tree(text=tree)
  tree.df=ggtree::ggtree(tree)$data
  tree.df=as.data.frame(tree.df)
  
  tree.df$node.offsprings=sapply(tree.df$node,
                                 function(i){
                                   j=tree.df[tidytree::offspring(tree,i,type="tips"),"label"]
                                   if (length(j)==0){j=tree.df[i,"label"]}
                                   j=sort(j)
                                   return(paste(j,collapse=";"))
                                 })
  tree.df$parent.offsprings=sapply(tree.df$parent,
                                 function(i){
                                   j=tree.df[tidytree::offspring(tree,i,type="tips"),"label"]
                                   if (length(j)==0){j=tree.df[i,"label"]}
                                   j=sort(j)
                                   return(paste(j,collapse=";"))
                                 })
  
  branch_attributes=json$`branch attributes`[[1]]
  tree.df$Uncorrected_P_value=sapply(tree.df$label,
                                     function(label){
                                       if (label==""){Uncorrected_P_value=NA}else{
                                       Uncorrected_P_value=(branch_attributes[[label]]$`Uncorrected P-value`)}
                                       return(Uncorrected_P_value)
                                     })
  tree.df$omega=sapply(tree.df$label,
                       function(label){
                         if (label==""){Baseline_MG94xREV_omega_ratio=NA}else{
                           Baseline_MG94xREV_omega_ratio=branch_attributes[[label]]$`Baseline MG94xREV omega ratio`}
                         return(Baseline_MG94xREV_omega_ratio)
                       })
  
  # tree.df[,c("Baseline_MG94xREV","Baseline_MG94xREV_omega_ratio","Full_adaptive_model","LRT","Nucleotide_GTR",
  #         "Rate_Distributions_dN","Rate_Distributions_dS","Rate_Distributions_omega","Uncorrected_P_value")]=
  # t(sapply(tree.df$label,
  #        function(label){
  #          
  #            Baseline_MG94xREV=branch_attributes[[label]]$`Baseline MG94xREV`
  #            Baseline_MG94xREV_omega_ratio=branch_attributes[[label]]$`Baseline MG94xREV omega ratio`
  #            Full_adaptive_model=branch_attributes[[label]]$`Full adaptive model`
  #            LRT=branch_attributes[[label]]$`LRT`
  #            Nucleotide_GTR=branch_attributes[[label]]$`Nucleotide GTR`
  #            
  #            Rate_Distributions=unlist(branch_attributes[[label]]$`Rate Distributions`)
  #            Rate_Distributions_dN=Rate_Distributions[seq(1,length(Rate_Distributions),2)]
  #            Rate_Distributions_dS=Rate_Distributions[seq(2,length(Rate_Distributions),2)]
  #            Rate_Distributions_dN.above.dS=Rate_Distributions_dN/Rate_Distributions_dS
  #            
  #            Rate_Distributions_dN=paste(as.character(Rate_Distributions_dN),collapse=";")
  #            Rate_Distributions_dS=paste(as.character(Rate_Distributions_dS),collapse=";")
  #            Rate_Distributions_dN.above.dS=paste(as.character(Rate_Distributions_dN.above.dS),collapse=";")
  #            
  #            Uncorrected_P_value=branch_attributes[[label]]$`Uncorrected P-value`
  #            res=c(Baseline_MG94xREV,Baseline_MG94xREV_omega_ratio,Full_adaptive_model,LRT,Nucleotide_GTR,
  #                  Rate_Distributions_dN,Rate_Distributions_dS,Rate_Distributions_dN.above.dS,Uncorrected_P_value)
  #        }))
  # 
  write.table(tree.df,tsv,sep="\t",row.names=FALSE,quote=FALSE)
}

# hyphy LIBPATH=/apps/unit/BourguignonU/hyphy/res FitMG94.bf --alignment CD2.nex
# FitMG94 of hyphy: estimate dN and dS on gene tree branches
# Obtain branch lengths under the nucleotide GTR model and remove zero lengths
# internal nodes might be deleted
fitMG94=function(FitMG94.bf="/apps/unit/BourguignonU/hyphy-analyses/FitMG94/FitMG94.bf",
                 hyphy.res="/apps/unit/BourguignonU/hyphy/res",
                 codon.align=codon.align,
                 geneTree=geneTree,
                 out_prefix=out_prefix){
    cmd=paste("hyphy",
              paste("LIBPATH=",hyphy.res,sep=""),
              FitMG94.bf,
              "--type local",
              "--lrt Yes",
              "--alignment",codon.align,
              "--tree",geneTree,
              "--output",paste(out_prefix,".json",sep=""))
    print(cmd);system(cmd,wait=TRUE)
}

json2tsv_fitMG94=function(json=json,
                          tsv=tsv){
  json=jsonlite::read_json(json)
  
  tree=json$input$trees[[1]]
  tree=paste(tree,";",sep="")
  tree=ape::read.tree(text=tree)
  tree.df=as.data.frame(ggtree::ggtree(tree)$data)
  
  tree.df$node.offsprings=sapply(tree.df$node,
                                 function(i){
                                   j=tree.df[tidytree::offspring(tree,i,type="tips"),"label"]
                                   if (length(j)==0){j=tree.df[i,"label"]}
                                   j=sort(j)
                                   return(paste(j,collapse=";"))
                                 })
  tree.df$parent.offsprings=sapply(tree.df$parent,
                                   function(i){
                                     j=tree.df[tidytree::offspring(tree,i,type="tips"),"label"]
                                     if (length(j)==0){j=tree.df[i,"label"]}
                                     j=sort(j)
                                     return(paste(j,collapse=";"))
                                   })
  
  branch_attributes=json$`branch attributes`[[1]]
  tree.df[,c("omega","omega_lowerBound","omega_higherBound",
          "Nucleotide_GTR","Standard_MG94","dN","dS",
          "nonsynonymous_substitution_per_codon",
          "synonymous_substitution_per_codon",
          "original_name")]=
  t(sapply(tree.df$label,
         function(label){
           if (label!=""){
             omega=branch_attributes[[label]]$`Confidence Intervals`$`MLE`
             omega_lowerBound=branch_attributes[[label]]$`Confidence Intervals`$`LB`
             omega_higherBound=branch_attributes[[label]]$`Confidence Intervals`$`UB`
             
             Nucleotide_GTR=branch_attributes[[label]]$`Nucleotide GTR`
             Standard_MG94=branch_attributes[[label]]$`Standard MG94`
             dN=branch_attributes[[label]]$`dN`
             dS=branch_attributes[[label]]$`dS`
             nonsynonymous_substitution_per_codon=branch_attributes[[label]]$`nonsynonymous`
             synonymous_substitution_per_codon=branch_attributes[[label]]$`synonymous`
             original_name=branch_attributes[[label]]$`original name`
             if (is.null(original_name)){original_name=NA}
             
             res=c(omega,omega_lowerBound,omega_higherBound,
                   Nucleotide_GTR,Standard_MG94,dN,dS,
                   nonsynonymous_substitution_per_codon,
                   synonymous_substitution_per_codon,
                   original_name)
           }else{
             res=rep(NA,10)
           }
           return(res)
         }))
  
  write.table(tree.df,tsv,sep="\t",row.names=TRUE,quote=FALSE)
}

# Label tree
labelTree=function(tree.text=tree.text,
                   node=list(c("Pada","Svic")), # label all branches descend from mrca
                   #states=c("ancestral"), # ancestral: only label ancestral branch
                                          # extant: only label extant species
                                          # all: all branch descended from the ancestral node
                   label="{Foreground}"){
  # tree.text="(Bori:0.2292580098,Cmer:0.2162881506,(Mdar:0.08741383512,(((Znev:0.08171119038,Hsjo:0.03143893694):0.02383815538):0.0174824005,((Kfla:0.04170015857,(PAsim:0.2307028781,(Gfus:0.07117802103,(Ncas:0.01840461903,((Rebo:0.01124430592,Mhub:0.02559844988):0.011116978,Cbre:0.0188960985):0.0004293834907):0.02294515693):1e-08):0.0439100372):0.05654769181,(Shal:0.112909926,((Gocu:0.1775427291,Dlon:0.0642505522):0.01172992014,(PRsim:0.4236401747,((Rfla:0.08005310728,(Hten:0.0323401336,Ctes:0.0491751227):0.03418668314):0.007710339714,((Ssph:0.08480561274,(Aaca:0.04038995652,(Ofor:0.0606636153,Mnat:0.04757076576):0.003797311233):0.07700947803):1e-08,(Fval:0.06784377311,((Aunk:0.1306170067,(Eunk:0.04223800603,(Apac:0.02508881891,Aban:0.0257989308):0.02700286682):0.0244634145):0.07372669573,((Munk:0.03865183289,(Shey:0.03763464222,(Llab:0.04259219138,Cwal:0.08379718508):0.006317183079):0.004748950186):0.008857657192,(((Punk:0.03122390861,Abea:0.02396904008):0.007952933317,(Pred:0.03509906124,Iunk:0.07815701173):0.008672879662):0.006641972461,(Cpar:0.113431962,(Ntar:0.1496951552,(Lunk:0.04040163269,((Nluj:0.01490959933,Hunk:0.04141330914):0.01058250769,(Ccav:0.03193998868,Csp4:0.1151284084):0.0003519319698):0.01406213512):0.02108581795):1e-08):1e-08):1e-08):0.04553573054):1e-08):0.004260697919):0.05288211324):0.04203912462):0.01203726929):0.04712678639):0.06609269134):1e-08):0.02961437174):0.001485681023);"
  # node=list(c("Pada","Svic"))
  # label="{Foreground}"
  tree=ape::read.tree(text=tree.text)
  tree$node.label=1:tree$Nnode+length(tree$tip.label)
  tree.df=as.data.frame(ggtree::ggtree(tree)$data)
  
  node=sapply(node,function(i){return(i[i %in% tree$tip.label])})
  mrca=sapply(node,
              function(n){
                res=ape::getMRCA(tree,n)
                if (is.null(res)){return(n)}else{return(res)}
              })
  offsprings=sapply(mrca,
                    function(m){
                      res=tidytree::offspring(tree,m,tiponly=F)
                      if (length(res)==0){return(m)}else{return(tree.df[tree.df$node %in% res,"label"])}
                    })
  offsprings=unlist(offsprings)
  offsprings=offsprings[offsprings %in% tree.df$label]
  
  tree$tip.label=sapply(tree$tip.label,
                        function(i){
                          if (i %in% offsprings){
                            return(paste(i,label,sep=""))
                          }else{
                            return(i)
                          }
                        })
  tree$node.label=sapply(tree$node.label,
                         function(i){
                           if (as.character(i) %in% offsprings){
                             return(paste(i,label,sep=""))
                           }else{
                             return(i)
                           }
                         })
  return( ape::write.tree(tree) )
}


BUSTED=function(hyphy.res="/apps/unit/BourguignonU/hyphy/res",
                codon.align=codon.align,
                geneTree=geneTree,
                foreground="Foreground",
                out_prefix=out_prefix){
  cmd=paste("hyphy",
            paste("LIBPATH=",hyphy.res,sep=""),
            "BUSTED",
            "--alignment",codon.align,
            "--tree",geneTree,
            "--branches",foreground,
            "--output",paste(out_prefix,".json",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  
  json=jsonlite::read_json(paste(out_prefix,".json",sep=""))
  return(json$`test results`$`p-value`)
}

RELAX=function(hyphy.res="/apps/unit/BourguignonU/hyphy/res",
               codon.align=codon.align,
               geneTree=geneTree,
               test="test",
               reference="reference",
               out_prefix=out_prefix){
  cmd=paste("hyphy",
            paste("LIBPATH=",hyphy.res,sep=""),
            "RELAX",
            "--alignment",codon.align,
            "--tree",geneTree,
            "--test",test,
            "--reference",reference,
            "--output",paste(out_prefix,".json",sep=""))
  system(cmd,wait=TRUE)
  json=jsonlite::read_json(paste(out_prefix,".json",sep=""))
  p=json$`test results`$`p-value`
  k=json$`test results`$`relaxation or intensification parameter`
  
  return(list(p=p,k=k))
}

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
                  start_val=0, # 0/1
                  estimation="ML", # ML: maximum likelyhood
                  rateModel="BD", # BD: Birth (per gene and per million of years), Death (per gene and per million of years)
                                  # BDI: Birth, Death, and Innovation (per million of years)
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
  # tree.nwk="/flash/BourguignonU/Cong/termite_pca/geneContent/droso.12sp.tamura.nwk"
  # famSize.tsv="/flash/BourguignonU/Cong/termite_pca/geneContent/4FAMs.12sp.tsv"
  # out_dir="/flash/BourguignonU/Cong/termite_pca/geneContent/BDI_ML"
  # threads=1
  # estimation="ML"
  # rateModel="BDI"
  
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  out_dir=sub("/$","",out_dir)
  
  # split gene families
  famSize=read.table(famSize.tsv,header=TRUE,sep="\t",quote="")
  if (!file.exists(paste(out_dir,"/og.split",sep=""))){
    system(paste("mkdir ",out_dir,"/og.split",sep=""))
  }
  sapply(1:nrow(famSize), 
         function(i){
           if (!file.exists(paste(out_dir,"/og.split/",famSize[i,1],".tsv",sep=""))){
             write.table(famSize[i,],
                         paste(out_dir,"/og.split/",famSize[i,1],".tsv",sep=""),
                         sep="\t",row.names = FALSE,quote=FALSE)
           }})
  
  # branch IDs of badirate
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
  # originalTree=ggtree(read.newick(tree.nwk))$data
  # badirateTree=ggtree(read.newick(paste(out_dir,"/branchID.nwk",sep="")))$data
  # originalTree=as.data.frame(originalTree)
  # originalTree$badirateID=str_extract(badirateTree$label,"[0-9]*$")
  # originalTree$badirate.parentID=sapply(1:nrow(originalTree),
  #                                       function(i){
  #                                         parent=originalTree[i,"parent"]
  #                                         return(originalTree[originalTree$node==parent,"badirateID"])
  #                                       })
  # originalTree$badirate.branchID=paste(as.character(originalTree$badirate.parentID),
  #                                      as.character(originalTree$badirateID),
  #                                      sep="->")
  # originalTree$estimation=rep(estimation,nrow(originalTree))
  # originalTree$rateModel=rep(rateModel,nrow(originalTree))
  
  # Setup badirate
  setupTab.0=prepare_badirate(og_path=paste(out_dir,"/og.split/",sep=""),
                              tree=tree.nwk,
                              branch_models=c(gr="GR",fr="FR"),
                              rate_model=rateModel,
                              estimation=estimation,
                              out_dir=paste(out_dir,"/rawOutput.0",sep=""),
                              script_dir=paste(out_dir,"/scripts.0",sep=""),
                              replicates=1,
                              ancestral=TRUE,
                              outlier=TRUE,
                              seed=20231212,
                              start_value=start_val,
                              pbs_q="smps",
                              badirate_path=badirate.pl,
                              create_scripts="pbs")
  if (!file.exists(paste(out_dir,"/setupTab.0.tsv",sep=""))){
    write.table(setupTab.0,paste(out_dir,"/setupTab.0.tsv",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
  }
  
  # badirate cmd
  og.lst=system(paste("ls ",out_dir,"/og.split/",sep=""),intern=TRUE)
  og.lst=paste(out_dir,"/og.split/",og.lst,sep="")
  #og.lst=paste(out_dir,"/og.split/",famSize[,1],".tsv",sep="")
  if (!file.exists(paste(out_dir,"/rawOutput.0/gr/",sep=""))){
    system(paste("mkdir ",out_dir,"/rawOutput.0/gr/",sep=""))
    system(paste("mkdir ",out_dir,"/rawOutput.0/fr/",sep=""))
  }
  cmd.0.GR.1=paste(perl.5.16,badirate.pl,
                   "-seed 20231212",
                   "-start_val",as.character(start_val),
                   "-unobs",
                   "-rmodel",rateModel,
                   #"-bmodel GR",
                   "-ep", estimation,
                   "-treefile",tree.nwk,
                   "-sizefile",og.lst,
                   "-anc -outlier",
                   "-out",paste(out_dir,"/rawOutput.0/gr/",
                                basename(og.lst),".gr01.bd",sep=""),
                   sep=" ")
  cmd.0.FR.1=paste(perl.5.16,badirate.pl,
                   "-seed 20231212",
                   "-start_val",as.character(start_val),
                   "-unobs",
                   "-rmodel",rateModel,
                   "-bmodel FR",
                   "-ep",estimation,
                   "-treefile",tree.nwk,
                   "-sizefile",og.lst,
                   "-anc -outlier",
                   "-out",paste(out_dir,"/rawOutput.0/fr/",
                                basename(og.lst),".fr01.bd",sep=""),
                   sep=" ")
  cmd=c(cmd.0.GR.1,cmd.0.FR.1)
  out.file=c(paste(out_dir,"/rawOutput.0/gr/",
                   basename(og.lst),".gr01.bd",sep=""),
             paste(out_dir,"/rawOutput.0/fr/",
                   basename(og.lst),".fr01.bd",sep=""))
  err.file=paste(out.file,".stderr",sep="")
  cmd=paste(cmd,">",err.file,sep=" ")
  tsv=paste(out_dir,"/fam_results/",basename(out.file),".tsv",sep="")
  
  # run badirate
  cmd=cmd[which(!file.exists(tsv))]
  cmd=sample(cmd,length(cmd),replace = FALSE)
  if (!file.exists(paste(out_dir,"/fam_results/",sep=""))){system(paste("mkdir ",out_dir,"/fam_results/",sep=""))}
  if (length(cmd)!=0){
    print("Start badirate jobs...")
    library(parallel)
    clus=makeCluster(as.numeric(threads))
    parSapply(clus,cmd,
              function(i){
                bd=unlist(strsplit(i," "))
                bd=bd[grepl("bd$",bd)]
                tsv=paste(out_dir,"/fam_results/",basename(bd),".tsv",sep="")
                
                bd.finished=TRUE
                if (!file.exists(bd)){bd.finished=FALSE}
                if (file.exists(bd)){if (file.size(bd)==0){bd.finished=FALSE}}
                if (!bd.finished){
                  system(i,wait=TRUE)
                }
                
                if (!file.exists(tsv)){
                  library(treeio)
                  library(ggtree)
                  library(stringr)
                  #bd="/flash/BourguignonU/Cong/termite_pca/seqs.HOG/badirate/start0_BDI/rawOutput.0/fr/HOG0000012.tsv.fr01.bd"
                  bd=readLines(bd)
                  
                  originalTree.nwk=bd[grepl("\ttreefile = ",bd)]
                  originalTree.nwk=sub("\ttreefile = ","",originalTree.nwk)
                  originalTree=ggtree(read.newick(originalTree.nwk))$data
                  originalTree=as.data.frame(originalTree)
                  badirateTree=ggtree(read.tree(text=bd[2]))$data
                  originalTree$badirateID=sub("_","",str_extract(badirateTree$label,"[0-9]*$"))
                  originalTree$badirate.parentID=sapply(1:nrow(originalTree),
                                                        function(i){
                                                          parent=originalTree[i,"parent"]
                                                          return(originalTree[originalTree$node==parent,"badirateID"])
                                                        })
                  originalTree$badirate.branchID=paste(as.character(originalTree$badirate.parentID),
                                                       as.character(originalTree$badirateID),
                                                       sep="->")
                  badirate.branchID2badirate.branchCode=bd[which(bd=="\tbmodel= "):which(grepl("\tstart_val",bd))]
                  badirate.branchID2badirate.branchCode=badirate.branchID2badirate.branchCode[-c(1,length(badirate.branchID2badirate.branchCode))]
                  originalTree$badirate.branchCode=sapply(originalTree$badirate.branchID,
                                                          function(i){
                                                            i=paste("\t",i,"\t",sep="")
                                                            r=badirate.branchID2badirate.branchCode[grepl(i,badirate.branchID2badirate.branchCode)]
                                                            r=str_extract(r,"\t[0-9]*$")
                                                            r=sub("\t","",r)
                                                            if (length(r)==0){r=NA}else{r=as.numeric(r)}
                                                            return(r)
                                                          })
                  seed=bd[grepl("\tseed = ",bd)];seed=sub("\tseed = ","",seed)
                  originalTree$seed=rep(as.numeric(seed),nrow(originalTree))
                  estimation=bd[grepl("\tep = ",bd)];estimation=sub("\tep = ","",estimation)
                  originalTree$estimation=rep(estimation,nrow(originalTree))
                  rateModel=bd[grepl("\trmodel = ",bd)];rateModel=sub("\trmodel = ","",rateModel)
                  originalTree$rateModel=rep(rateModel,nrow(originalTree))
                  start_val=bd[grepl("\tstart_val = ",bd)];start_val=sub("\tstart_val = ","",start_val)
                  originalTree$start_val=rep(as.numeric(start_val),nrow(originalTree))
                  sizefile=bd[grepl("\tsizefile = ",bd)];sizefile=sub("\tsizefile = ","",sizefile)
                  originalTree$sizefile=rep(sizefile,nrow(originalTree))
                  
                  res=originalTree
                  branchGroups=res[,"badirate.branchCode"];branchGroups=branchGroups[!is.na(branchGroups)]
                  branchGroups=length(branchGroups[!duplicated(branchGroups)])
                  
                  if (rateModel=="BD"){Parameters=1+2*branchGroups}
                  if (rateModel=="BDI"){Parameters=1+3*branchGroups}
                  
                  res$Parameters=rep(Parameters,nrow(originalTree))
                  
                  res$fam=rep(NA,nrow(res))
                  res$log.likelihood=rep(NA,nrow(res))
                  res$AIC=rep(NA,nrow(res))
                  res$ancestralState=rep(NA,nrow(res))
                  res$birthRate.bd=rep(NA,nrow(res))
                  res$deathRate.bd=rep(NA,nrow(res))
                  res$innovationRate.bd=rep(NA,nrow(res))
                  res$outlier=rep(FALSE,nrow(res))
                  if (bd[length(bd)]=="END OUTPUT"){
                    log.likelihood=bd[grepl("\t\t#Likelihood: ",bd)]
                    log.likelihood=as.numeric(sub("\t\t#Likelihood: ","",log.likelihood))
                    res$log.likelihood=rep(log.likelihood,nrow(res))
                    res$AIC=2*res$Parameters-2*res$log.likelihood # the preferred model is the one with the minimum AIC value
                    
                    anc=bd[which(bd=="\t\t#Family\tAncestral Family Size Tree")+1]
                    res$fam=rep(unlist(strsplit(anc,"\t"))[3],nrow(res))
                    anc=unlist(strsplit(anc,"\t"))[4]
                    anc=as.data.frame(ggtree(read.tree(text=anc))$data)
                    res$ancestralState=as.numeric(sub("^.*_","",anc$label))
                    
                    rates=bd[which(grepl("#Branch_Group\tBirth\tDeath",bd)):which(bd=="\t##Ancestral Family Size")]
                    rates=rates[c(-1,-length(rates),-length(rates)+1)]
                    rates=sub("\t\t","",rates)
                    rates=strsplit(rates,"\t")
                    for (j in 1:length(rates)){
                      res[res$badirate.branchCode==as.numeric(rates[[j]][1]) & !is.na(res$badirate.branchCode),
                          "birthRate.bd"]=as.numeric(rates[[j]][2])
                      res[res$badirate.branchCode==as.numeric(rates[[j]][1]) & !is.na(res$badirate.branchCode),
                          "deathRate.bd"]=as.numeric(rates[[j]][3])
                      res[res$badirate.branchCode==as.numeric(rates[[j]][1]) & !is.na(res$badirate.branchCode),
                          "innovationRate.bd"]=as.numeric(rates[[j]][4])
                    }
                    
                    outlier=bd[which(bd=="\t##Outlier Families per Branch"):which(bd=="END OUTPUT")]
                    outlier=outlier[-c(1,2,length(outlier),length(outlier)-1,length(outlier)-2)]
                    outlier=sub("\t\t","",outlier)
                    outlier=stringr::str_extract(outlier,"[0-9]*->[0-9]*")
                    if (length(outlier)!=0){res[res$badirate.branchID%in%outlier,"outlier"]=TRUE}
                  }
                  write.table(res,tsv,
                              row.names=FALSE,quote=FALSE,sep="\t")
                }
              })
  }
}

# Dependencies: ggtree (R), treeio (R), stringr (R)
bd2tsv=function(bd=bd, # output of badirate
                #originalTree.nwk=originalTree.nwk,
                #badirateTree.nwk=badirateTree.nwk, # This file from `perl.5.16 badirate.pl -print_ids -treefile ${originalTree.nwk} -sizefile ${famSize.tsv} > ${badirateTree.nwk}`
                tsv=tsv){
  library(treeio)
  library(ggtree)
  library(stringr)
  #bd="/flash/BourguignonU/Cong/termite_pca/seqs.HOG/badirate/start0_BDI/rawOutput.0/fr/HOG0000012.tsv.fr01.bd"
  bd=readLines(bd)
  
  originalTree.nwk=bd[grepl("\ttreefile = ",bd)]
  originalTree.nwk=sub("\ttreefile = ","",originalTree.nwk)
  originalTree=ggtree(read.newick(originalTree.nwk))$data
  originalTree=as.data.frame(originalTree)
  badirateTree=ggtree(read.tree(text=bd[2]))$data
  originalTree$badirateID=sub("_","",str_extract(badirateTree$label,"[0-9]*$"))
  originalTree$badirate.parentID=sapply(1:nrow(originalTree),
                                        function(i){
                                          parent=originalTree[i,"parent"]
                                          return(originalTree[originalTree$node==parent,"badirateID"])
                                        })
  originalTree$badirate.branchID=paste(as.character(originalTree$badirate.parentID),
                                       as.character(originalTree$badirateID),
                                       sep="->")
  badirate.branchID2badirate.branchCode=bd[which(bd=="\tbmodel= "):which(grepl("\tstart_val",bd))]
  badirate.branchID2badirate.branchCode=badirate.branchID2badirate.branchCode[-c(1,length(badirate.branchID2badirate.branchCode))]
  originalTree$badirate.branchCode=sapply(originalTree$badirate.branchID,
                                          function(i){
                                            i=paste("\t",i,"\t",sep="")
                                            r=badirate.branchID2badirate.branchCode[grepl(i,badirate.branchID2badirate.branchCode)]
                                            r=str_extract(r,"\t[0-9]*$")
                                            r=sub("\t","",r)
                                            if (length(r)==0){r=NA}else{r=as.numeric(r)}
                                            return(r)
                                          })
  seed=bd[grepl("\tseed = ",bd)];seed=sub("\tseed = ","",seed)
  originalTree$seed=rep(as.numeric(seed),nrow(originalTree))
  estimation=bd[grepl("\tep = ",bd)];estimation=sub("\tep = ","",estimation)
  originalTree$estimation=rep(estimation,nrow(originalTree))
  rateModel=bd[grepl("\trmodel = ",bd)];rateModel=sub("\trmodel = ","",rateModel)
  originalTree$rateModel=rep(rateModel,nrow(originalTree))
  start_val=bd[grepl("\tstart_val = ",bd)];start_val=sub("\tstart_val = ","",start_val)
  originalTree$start_val=rep(as.numeric(start_val),nrow(originalTree))
  sizefile=bd[grepl("\tsizefile = ",bd)];sizefile=sub("\tsizefile = ","",sizefile)
  originalTree$sizefile=rep(sizefile,nrow(originalTree))
  
  res=originalTree
  branchGroups=res[,"badirate.branchCode"];branchGroups=branchGroups[!is.na(branchGroups)]
  branchGroups=length(branchGroups[!duplicated(branchGroups)])
  
  if (rateModel=="BD"){Parameters=1+2*branchGroups}
  if (rateModel=="BDI"){Parameters=1+3*branchGroups}
  
  res$Parameters=rep(Parameters,nrow(originalTree))
  
  res$fam=rep(NA,nrow(res))
  res$log.likelihood=rep(NA,nrow(res))
  res$AIC=rep(NA,nrow(res))
  res$ancestralState=rep(NA,nrow(res))
  res$birthRate.bd=rep(NA,nrow(res))
  res$deathRate.bd=rep(NA,nrow(res))
  res$innovationRate.bd=rep(NA,nrow(res))
  res$outlier=rep(FALSE,nrow(res))
  if (bd[length(bd)]=="END OUTPUT"){
    log.likelihood=bd[grepl("\t\t#Likelihood: ",bd)]
    log.likelihood=as.numeric(sub("\t\t#Likelihood: ","",log.likelihood))
    res$log.likelihood=rep(log.likelihood,nrow(res))
    res$AIC=2*res$Parameters-2*res$log.likelihood # the preferred model is the one with the minimum AIC value
    
    anc=bd[which(bd=="\t\t#Family\tAncestral Family Size Tree")+1]
    res$fam=rep(unlist(strsplit(anc,"\t"))[3],nrow(res))
    anc=unlist(strsplit(anc,"\t"))[4]
    anc=as.data.frame(ggtree(read.tree(text=anc))$data)
    res$ancestralState=as.numeric(sub("^.*_","",anc$label))
    
    rates=bd[which(grepl("#Branch_Group\tBirth\tDeath",bd)):which(bd=="\t##Ancestral Family Size")]
    rates=rates[c(-1,-length(rates),-length(rates)+1)]
    rates=sub("\t\t","",rates)
    rates=strsplit(rates,"\t")
    for (i in 1:length(rates)){
      res[res$badirate.branchCode==as.numeric(rates[[i]][1]) & !is.na(res$badirate.branchCode),
              c("birthRate.bd","deathRate.bd","innovationRate.bd")]=as.numeric(rates[[i]][2:4])
    }
    
    outlier=bd[which(bd=="\t##Outlier Families per Branch"):which(bd=="END OUTPUT")]
    outlier=outlier[-c(1,2,length(outlier),length(outlier)-1,length(outlier)-2)]
    outlier=sub("\t\t","",outlier)
    outlier=stringr::str_extract(outlier,"[0-9]*->[0-9]*")
    if (length(outlier)!=0){res[res$badirate.branchID%in%outlier,"outlier"]=TRUE}
  }
  write.table(res,tsv,
              row.names=FALSE,quote=FALSE,sep="\t")
} 



# DupliPHY-ML: Rates specified for families, but not for branches
# 2 threads
# Dependencies: Dupliphy-ML, parallel (R)
DupliPHYml=function(DupliPHYml.jar="~/Softwares/Dupliphy-ML/DupliphyML.jar",
                    tree.nwk=tree.nwk,
                    famSize.tsv=famSize.tsv, # Fields: FAM_ID sp1 sp2
                    out_prefix=out_prefix){
  cmd1=paste("java -Xmx509952m -jar",DupliPHYml.jar,
            famSize.tsv,"BDI+G",tree.nwk,paste(out_prefix,".BDI",sep=""),"100000",
            ">",paste(out_prefix,".BDI.stderr",sep=""),
            sep=" ")
  #print(cmd);system(cmd,wait=TRUE)
  cmd2=paste("java -Xmx509952m -jar",DupliPHYml.jar,
            famSize.tsv,"Parsimony+G",tree.nwk,paste(out_prefix,".Parsimony",sep=""),"100000",
            ">",paste(out_prefix,".Parsimony.stderr",sep=""),
            sep=" ")
  library(parallel)
  clus=makeCluster(2)
  parSapply(clus,c(cmd2,cmd1),function(i){system(i,wait=TRUE)})
  #print(cmd);system(cmd,wait=TRUE)
  
}

# Count: duplication (lambda), loss (mu) and transfer (kappa), Rates specified for families & branches
# Dependencies: Count
Count=function(Count.jar="~/Softwares/Count/Count.jar",
               tree.nwk=tree.nwk,
               famSize.tsv=famSize.tsv, # Fields: FAM_ID sp1 sp2
               out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  out_dir=sub("/$","",out_dir)
  wd=getwd()
  setwd(out_dir)
  
  # Simple model
  cmd=paste("java -Xmx204800m -jar",Count.jar,"ML",
            "-gain_k 1 -loss_k 1 -duplication_k 1", # Uniform rates between families
            "-uniform_duplication true -uniform_duplication true", # Uniform rates between branches
            tree.nwk,famSize.tsv,"> model.1.txt",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # lineage specific rates
  cmd=paste("java -Xmx204800m -jar",Count.jar,"ML",
            "-gain_k 1 -loss_k 1 -duplication_k 1", # Uniform rates between families
            tree.nwk,famSize.tsv," model.1.txt > model.2.txt",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # lineage and family specific rates
  cmd=paste("java -Xmx204800m -jar",Count.jar,"ML",
            "-gain_k 2 -loss_k 2 -duplication_k 2", # Uniform rates between families
            tree.nwk,famSize.tsv," model.2.txt > model.3.txt",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  # lineage and family specific rates (more complex)
  cmd=paste("java -Xmx204800m -jar",Count.jar,"ML",
            "-gain_k 3 -loss_k 3 -duplication_k 3", # Uniform rates between families
            tree.nwk,famSize.tsv," model.3.txt > model.4.txt",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Ancestral state
  cmd=paste("java -Xmx204800m -jar",Count.jar,"Posteriors",
            tree.nwk,famSize.tsv," model.4.txt > ancestral.tsv",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# Ancestral states via ape
# Dependencies: ape (R), ggtree (R), treeio (R), reshape (R)
ape.ace=function(tree.nwk=tree.nwk,
                 famSize.tsv=famSize.tsv, # Fields: FAM_ID sp1 sp2
                 #type="continuous", # continuous/discrete
                 out.tsv=out.tsv){
  #tree.nwk="/bucket/BourguignonU/Cong/termite_pca/timeTree_label.nwk"
  #famSize.tsv="/bucket/BourguignonU/Cong/termite_pca/AI/true_geneANDpgeneMat_copyMat.tsv"
  library(ape)
  library(ggtree)
  library(reshape2)
  
  tree=read.tree(tree.nwk)
  famSize=read.table(famSize.tsv,sep="\t",header=TRUE,quote="")
  rownames(famSize)=famSize[,1]
  df=as.data.frame(ggtree(tree)$data)
  node2lab=df[!df$isTip,c("node","label")]
  rownames(node2lab)=node2lab$label
  
  mat=matrix(rep(0,nrow(famSize)*nrow(node2lab)),
             nrow=nrow(famSize),ncol=nrow(node2lab))
  rownames(mat)=rownames(famSize);colnames(mat)=rownames(node2lab)
  res=reshape2::melt(mat)
  colnames(res)=c("Fam","nodeLab","Node")
  res[,"Node"]=node2lab[(res[,"nodeLab"]),"node"]
  res[,"ace"]=rep(NA,nrow(res))
  res[,"95CI.low"]=rep(NA,nrow(res))
  res[,"95CI.high"]=rep(NA,nrow(res))
  rownames(res)=paste(res$Fam,res$nodeLab,sep="-")
  
  # library(parallel)
  # clus=makeCluster(as.numeric(threads))
  # clusterExport(cl=clus,varlist = list("famSize","ace","tree","res"))
  for (i in 1:nrow(famSize)){
    fam=famSize[i,1]
    x=famSize[i,-1]
    x=as.numeric(x)
    names(x)=colnames(famSize)[-1]
    ancestral=ace(x=x,
                  phy=tree,
                  type="continuous",
                  method="ML",
                  CI=TRUE)
    
    state=ancestral$ace
    names(state)=paste(fam,as.character(names(state)),sep="-")
    res[names(state),"ace"]=state
    
    CI=ancestral$CI95
    rownames(CI)=names(state)
    res[names(state),"95CI.low"]=CI[,1]
    res[names(state),"95CI.high"]=CI[,2]
  }
  write.table(res,out.tsv,sep="\t",row.names = FALSE,quote=FALSE)
}

# Count gain/loss/innovation events from ape.ace
# Gene conversion
# ../../geneconv -Seqfile=test.fna -Outfile=test.frags -Logfile=test.log -Seqtype=SIL -Include_monosites  -ListPair

#####
# selective pressure
#####
RERconverge_readTree=function(tree.tsv=tree.tsv, # fields: gene, tree (nwk, no node label, same topology)
                              useSpecies=useSpecies, # comma-lst
                              out_prefix=out_prefix){
  library(RERconverge)
  geneTrees=readTrees(tree.tsv)
  saveRDS(geneTrees,paste(out_prefix,"_readTrees.rds",sep=""))
  
  RER=getAllResiduals(geneTrees,
                      useSpecies=unlist(strsplit(useSpecies,",")),
                      transform="sqrt",weighted=T,scale=T,plot=F)
  saveRDS(RER,paste(out_prefix,"_getAllResiduals.rds",sep=""))
}

#pseudergate=c("Svic","Pada","Znev","Hsjo","Kfla","PAsim","Gfus","Ncas","Rebo","Mhub","Cbre","Isch","Shal",
#              "Gocu","PRsim")
#worker=c("Mdar","Dlon","Rfla","Hten","Cges","Ctes","Ssph","Aaca","Ofor","Mnat","Fval","Aunk","Eunk","Apac",
#         "Aban","Munk","Shey","Llab","Cwal","Punk","Abea","Pred","Iunk","Cpar","Ntar","Lunk","Nluj","Hunk",
#         "Ccav","Csp4")

#click_select_foreground_branches(toyTrees$masterTree)
RERconverge_binary=function(geneTrees.RDS=geneTrees.RDS,
                            RER.RDS=RER.RDS,
                            binaryTree.txt="",
                            out_prefix=out_prefix){
  geneTrees=readRDS(geneTrees.RDS)
  RER=readRDS(RER.RDS)
  
  library(RERconverge)
  foreground.tree=read.tree(text=binaryTree.txt)
  
  geneTreesWithForeground=tree2Paths(foreground.tree,geneTrees)
  saveRDS(geneTreesWithForeground,paste(out_prefix,"_tree2Paths.rds",sep=""))
  
  correlatons=correlateWithBinaryPhenotype(RER,geneTreesWithForeground,
                                           min.sp=10, min.pos=2,
                                           weighted="auto",
                                           winsorizeRER=NULL,winsorizetrait=NULL,
                                           bootstrap=FALSE, bootn=10000)
  saveRDS(correlatons,paste(out_prefix,"_correlateWithBinaryPhenotype.rds",sep=""))
}
  
RERconverge_continuous=function(geneTrees.RDS=geneTrees.RDS,
                               RER.RDS=RER.RDS,
                               trait_val=trait_val, # numeric vector named with spp
                               out_prefix=out_prefix){
  geneTrees=readRDS(geneTrees.RDS)
  RER=readRDS(RER.RDS)
  
  library(RERconverge)
  charpaths=char2Paths(trait_val,geneTrees)
  saveRDS(charpaths,paste(out_prefix,"_char2Paths.rds",sep=""))
  
  correlatons=correlateWithContinuousPhenotype(RER,charpaths, 
                                               min.sp=10,winsorizeRER=3,winsorizetrait=3)
  saveRDS(charpaths,paste(out_prefix,"_correlateWithContinuousPhenotype.rds",sep=""))
}


