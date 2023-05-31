# Genomic analysis of protein-coding genes: 
# genome collinearity, HGT, phylogeny, condon/AA usage, exon-intron, protein domain, gene family, selection, etc.

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
# Incomplete
# AvP
# Detect horizontal gene transfer (HGT)
# Dependencies: AvP
Avp=function(query.faa=query.faa,
             blast_filename=blast_filename,
             dmbd=dmbd, # not used if blast_filename exists
             group.yaml=group.yaml,
             config.yaml=config.yaml,
             out_dir=out_dir,
             threads=threads
){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  threads=as.character(threads)
  wd=getwd()
  setwd(out_dir)
  
  if (!file.exists(blast_filename)){
    cmd=paste("diamond blastp",
              "--threads",threads,
              "-d",dmbd,
              "--max-target-seqs 500",
              "-q",query.faa,
              "--evalue 1e-5",
              "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids",
              ">",
              blast_filename)
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("calculate_ai.py",
            "-i",blast_filename,
            "-x",group.yaml,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("avp prepare",
            "-a",paste(blast_filename,"_ai.out",sep=""),
            "-o",".",
            "-f",query.faa,
            "-b",blast_filename,
            "-x",group.yaml,
            "-c",config.yaml,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("avp detect",
            "-i","./mafftgroups/",
            "-o",".",
            "-g","./groups.tsv",
            "-t","./tmp/taxonomy_nexus.txt",
            "-c",config.yaml,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
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
                 outAlign.fa=outAlign.fa){
  cmd=paste("grep '>' ",protAlign.fa," | sed 's/>//' > ",outAlign.fa,".lst",sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("seqkit faidx ",cds.fna," --infile-list ",outAlign.fa,".lst > ",outAlign.fa,"_CDS.fa",sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("pal2nal.pl",
            protAlign.fa,
            paste(outAlign.fa,"_CDS.fa",sep=""),
            "-output fasta",">",
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
            "-keepheader",
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
            # "-mset GTR", # test GTR+... models only
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



