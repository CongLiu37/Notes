# Genomic analysis based on genome annotation (protein-coding)
# Gene function annotation, genome collinearity, horizontal gene transfer

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

# OrthoFinder
# OrthoFinder creates a results directory called ‘OrthoFinder’ inside the input proteome directory and puts the results here.
# Orthofinder with "-M msa" infers phylogeny by mafft & fasttree
# Dependencies: Orthofinder, mafft, fasttree
Orthofinder=function(in_dir=in_dir, # Input proteome directory. 
                                    # One file per species with extension '.faa'
                     threads=threads){
  threads=as.character(threads)
  
  cmd=paste("orthofinder.py",
            "-f",in_dir,
            "-M msa",
            "-t",threads,
            "-a",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Re-run orthofinder with designated phylogenetic tree
Re_orthofinder=function(previous_orthofinder_result_dir=previous_orthofinder_result_dir,
                        tree=tree,
                        threads=threads){
  threads=as.character(threads)
  
  cmd=paste("orthofinder",
            "-t",threads,
            "-a",threads,
            "-ft",previous_orthofinder_result_dir,
            "-s",tree,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Best protein-protein hit by blast
best_blastp=function(query.faa=query.faa,
                     reference.faa=reference.faa,
                     out_dir=out_dir,
                     out_basename=out_basename,
                     threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",query.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  query.faa=unlist(strsplit(query.faa,"/"));query.faa=query.faa[length(query.faa)]
  
  cmd=paste("cp",reference.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  reference.faa=unlist(strsplit(reference.faa,"/"));reference.faa=reference.faa[length(reference.faa)]
  
  cmd=paste("makeblastdb",
            "-in",reference.faa,
            "-dbtype","prot",
            "-parse_seqids",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("blastp",
            "-num_threads",threads,
            "-db",reference.faa,
            "-query",query.faa,
            "-outfmt 6",
            "-evalue 1e-5",
            "-num_alignments 1",
            "-out",paste(out_dir,"/",out_basename,".blast",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",query.faa,sep=" "),wait=TRUE)
  system(paste("rm",reference.faa,sep=" "),wait=TRUE)
  system(paste("rm"," ",reference.faa,".*",sep=""),wait=TRUE)
  setwd(wd)
}

# Interproscan
interpro=function(proteins.faa=proteins.faa,
                  out_dir=out_dir,
                  out_basename=out_basename,
                  threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",proteins.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  proteins.faa=unlist(strsplit(proteins.faa,"/"));proteins.faa=proteins.faa[length(proteins.faa)]
  
  cmd=paste("sed -i 's/*/X/g'",proteins.faa,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("interproscan.sh",
            "-dp",
            "-b",out_basename,
            "-cpu",threads,
            "-f TSV",
            "-goterms",
            "-i",proteins.faa,
            "--pathways",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",proteins.faa,sep=" "),wait=TRUE)
  setwd(wd)
}

# WGDI input files
# Dependencies: blastp
wgdi_input=function(genome.fna1=genome.fna1,gff1=gff1,proteins.faa1=proteins.faa1,cds.fna1=cds.fna1,
                    genome.fna2=genome.fna2,gff2=gff2,proteins.faa2=proteins.faa2,cds.fna2=cds.fna2,
                    threads=threads,
                    out_dir=out_dir,
                    out_basename=out_basename){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  system("mkdir tmp",wait=TRUE)
  
  self2self=genome.fna1==genome.fna2
  setwd("tmp/")
  
  cmd=paste("cp"," ",proteins.faa1," ",".",sep="")
  print(cmd);system(cmd,wait=TRUE)
  proteins.faa1=unlist(strsplit(proteins.faa1,"/"));proteins.faa1=proteins.faa1[length(proteins.faa1)]
  if (!self2self){
    cmd=paste("cp"," ",proteins.faa2," ",".",sep="")
    print(cmd);system(cmd,wait=TRUE)
  }
  proteins.faa2=unlist(strsplit(proteins.faa2,"/"));proteins.faa2=proteins.faa2[length(proteins.faa2)]
  
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
  if (!self2self){
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
  }else{
    #blast2=paste(out_dir,"/",out_basename,".blast",sep="")
    blast2="false"
  }
  
  system(paste("mv",proteins.faa1,"../",sep=" "),wait=TRUE)
  if (!self2self){system(paste("mv",proteins.faa2,"../",sep=" "),wait=TRUE)}
  system(paste("mv",blast1,"../",sep=" "),wait=TRUE)
  if (!self2self){system(paste("mv",blast2,"../",sep=" "),wait=TRUE)}
  setwd("../")
  system("rm -r tmp",wait=TRUE)
  system(paste("cp"," ",cds.fna1," ",".",sep=""),wait=TRUE)
  system(paste("cp"," ",cds.fna2," ",".",sep=""),wait=TRUE)
  
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
    
    cmd=paste("awk -F '\t' -v OFS='\t' '{$2=$7; print $0}'",paste(out_prefix,".gff",sep=""),
              ">",paste(out_prefix,"_tmp.gff",sep=""),sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    
    system(paste("rm",paste(out_prefix,".gff",sep=""),sep=" "),wait=TRUE)
    system(paste("mv",
                 paste(out_prefix,"_tmp.gff",sep=""),
                 paste(out_prefix,".gff",sep=""),
                 sep=" "),wait=TRUE)
  }
  genome_prefix1=unlist(strsplit(genome.fna1,"/"));genome_prefix1=genome_prefix1[length(genome_prefix1)]
  wgdi_input(genome.fna=genome.fna1,
             gff=gff1,
             out_prefix=genome_prefix1)
  fake_gff1=paste(out_dir,"/",genome_prefix1,".gff",sep="")
  len1=paste(out_dir,"/",genome_prefix1,".len",sep="")
  if (!self2self){
    genome_prefix2=unlist(strsplit(genome.fna2,"/"));genome_prefix2=genome_prefix2[length(genome_prefix2)]
    wgdi_input(genome.fna=genome.fna2,
               gff=gff2,
               out_prefix=genome_prefix2)
    fake_gff2=paste(out_dir,"/",genome_prefix2,".gff",sep="")
    len2=paste(out_dir,"/",genome_prefix2,".len",sep="")
  }else{
    genome_prefix2=genome_prefix1
    fake_gff2=paste(out_dir,"/",genome_prefix1,".gff",sep="")
    len2=paste(out_dir,"/",genome_prefix1,".len",sep="")
  }
  
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
  conf[17]="markersize = 0.5"
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
  conf[10]="markersize = 1"  
  conf[11]="area = 0,10"
  conf[12]="block_length = 5"
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


