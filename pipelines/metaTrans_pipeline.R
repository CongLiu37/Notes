arg=commandArgs(trailingOnly = TRUE)

main.R=arg[1]
conf=arg[2]

source(main.R)
source(conf)

out_dir=sub("/$","",out_dir)
save_dir=sub("/$","",save_dir)
threads=as.character(threads)
if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
tab=read.table(tab,sep="\t",header=FALSE,quote="")

#####
# fq_noHost
#####
if (1 %in% Step){
  print("Step 1: fq_noHost")
  if (!file.exists(paste(out_dir,"/fq_noHost",sep=""))){system(paste("mkdir ",out_dir,"/fq_noHost",sep=""))}
  for (i in 1:nrow(tab)){
    sp=tab[i,1]
    fq1.lst=tab[i,2]
    fq2.lst=tab[i,3]
    hostGenome=tab[i,4]
    bam.lst=tab[i,5]
    Stamp1=paste(out_dir,"/fq_noHost/",sp,".fq_noHost.finished",sep="")
    Stamp2=paste(save_dir,"/fq_noHost/",sp,".fq_noHost.finished",sep="")
    if (file.exists(Stamp1) | file.exists(Stamp2)){
      print(paste(sp," fq_noHost FINISHED"))
    }else{
      print(paste(sp," fq_noHost START"))
      if (bam.lst!="none"){
        print(paste(sp," BAMs provided"))
        MergeBAM(BAMs=paste(unlist(strsplit(bam.lst,",")),collapse=" "),# SPACE-separated list of bam files.
                 out_prefix=paste(out_dir,"/fq_noHost/",sp,sep=""),
                 Threads=threads)
        Extract_fq(bam=paste(out_dir,"/fq_noHost/",sp,".bam",sep=""),
                   paired=TRUE, # logical. T for paired. 
                   mapped=FALSE, # logical. T for mapped
                   out_prefix=paste(out_dir,"/fq_noHost/",sp,sep=""),
                   reads_format="fq", # "fq" or "fa"
                   threads=threads)
      }else{
        fq1=unlist(strsplit(fq1.lst,","));fq2=unlist(strsplit(fq2.lst,","))
        for (j in 1:length(fq1)){
          Hisat(fq1=fq1[j],fq2=fq2[j], # Input fq files. Make fq2="None" if single-end.
                         # Comma-separated list.
               fna=hostGenome, # genome (soft masked for training AUGUSTUS)
               index=paste(out_dir,"/fq_noHost/",sp,sep=""), # Basename of Hisat2 index of reference genome.
               out_prefix=paste(out_dir,"/fq_noHost/",sp,".",as.character(j),sep=""), # Prefix of output BAM file.
               threads=threads)
        }
        MergeBAM(BAMs=paste(out_dir,"/fq_noHost/",sp,".",as.character(1:length(fq1)),".bam",sep=""),# SPACE-separated list of bam files.
                 out_prefix=paste(out_dir,"/fq_noHost/",sp,sep=""),
                 Threads=threads)
        system(paste("rm",
                     paste(out_dir,"/fq_noHost/",sp,".",as.character(1:length(fq1)),".bam",sep=""),
                     sep=" "))
        Extract_fq(bam=paste(out_dir,"/fq_noHost/",sp,".bam",sep=""),
                   paired=TRUE, # logical. T for paired. 
                   mapped=FALSE, # logical. T for mapped
                   out_prefix=paste(out_dir,"/fq_noHost/",sp,sep=""),
                   reads_format="fq", # "fq" or "fa"
                   threads=threads)
      }
      system(paste("touch",Stamp1,sep=" "))
    }
  }
}  

#####
# fq_noHost_noRrna
#####
if (2 %in% Step){
  print("Step 2: fq_noHost_noRrna")
  if (!file.exists(paste(out_dir,"/fq_noHost_noRrna",sep=""))){system(paste("mkdir ",out_dir,"/fq_noHost_noRrna",sep=""))}
  for (i in 1:nrow(tab)){
    sp=tab[i,1]
    Stamp1=paste(out_dir,"/fq_noHost_noRrna/",sp,".fq_noHost_noRrna.finished",sep="")
    Stamp2=paste(save_dir,"/fq_noHost_noRrna/",sp,".fq_noHost_noRrna.finished",sep="")
    if (file.exists(Stamp1) | file.exists(Stamp2)){
      print(paste(sp," fq_noHost_noRrna FINISHED"),sep="")
    }else{
      print(paste(sp," fq_noHost_noRrna START"),sep="")
      sortmerna(fq1=paste(save_dir,"/fq_noHost/",sp,".1.fq.gz",sep=""),
                fq2=paste(save_dir,"/fq_noHost/",sp,".2.fq.gz",sep=""),
                reference=sortmerna.db, # smr_v4.3_default_db.fasta
                out_dir=paste(out_dir,"/fq_noHost_noRrna/",sp,sep=""),
                threads=threads)
      system(paste("touch",Stamp1,sep=" "))
    }
  }
}

#####
# SPAdes
#####
if (3 %in% Step){
  print("Step 3: SPAdes")
  if (!file.exists(paste(out_dir,"/SPAdes",sep=""))){system(paste("mkdir ",out_dir,"/SPAdes",sep=""))}
  Stamp1=paste(out_dir,"/SPAdes/","SPAdes.finished",sep="")
  Stamp2=paste(save_dir,"/SPAdes/","SPAdes.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 3: SPAdes FINISHED")
  }else{
    print("Step 3: SPAdes START")
    read1=system(paste("ls ",save_dir,"/fq_noHost_noRrna/*/noRrna_fwd.fq.gz",sep=""),wait=TRUE,intern=TRUE)
    read2=system(paste("ls ",save_dir,"/fq_noHost_noRrna/*/noRrna_rev.fq.gz",sep=""),wait=TRUE,intern=TRUE)
    SPAdes(fq1=paste(read1,collapse=","),
           fq2=paste(read2,collapse=","), 
           contigs.fa="none", # Reliable contigs of the same genome.
           pacbio_clr="none",nanopore="none",sanger="none", # Long reads
           meta=FALSE, # Logical. If TRUE, run metaSPAdes.
           # metaSPAdes supports only a single short-read library which has to be paired-end.
           rna=TRUE, # Logical. TRUE for rnaSPAdes
           bio=FALSE, # Logical. TRUE for biosyntheticSPAdes
           custom_hmms="none", # directory with custom hmms for biosyntheticSPAdes
           out_dir=paste(out_dir,"/SPAdes",sep=""),
           threads=64,
           memory=500)
    system(paste("touch",Stamp1,sep=" "))
  }
}

#####
# seqNR
#####
if (4 %in% Step){
  print("Step 4: seqNR")
  if (!file.exists(paste(out_dir,"/seqNR",sep=""))){system(paste("mkdir ",out_dir,"/seqNR",sep=""))}
  Stamp1=paste(out_dir,"/seqNR/","seqNR.finished",sep="")
  Stamp2=paste(save_dir,"/seqNR/","seqNR.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 4: seqNR FINISHED")
  }else{
    print("Step 4: seqNR START")
    seqNR(in.fasta=paste(save_dir,"/SPAdes/hard_filtered_transcripts.fasta",sep=""),
          out_dir=paste(out_dir,"/seqNR",sep=""),
          Identity=0.75, # [0.0,1.0]
          cov_mode=1, # 0: alignment covers ${coverage} of target and of query
          # 1: alignment covers ${coverage} of target
          # 2: alignment covers ${coverage} of query
          # 3: target is of ${coverage} query length
          coverage=0.75, # [0.0,1.0]
          threads=threads)
    SimplifyID(fna=paste(out_dir,"/seqNR/rep.fasta",sep=""),
               common_pattern="contig", # common pattern in simplified sequence IDs
               out_dir=paste(out_dir,"/seqNR/",sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}

#####
# cdsInTranscripts
#####
if (5 %in% Step){
  print("Step 5: cdsInTranscripts")
  if (!file.exists(paste(out_dir,"/cdsInTranscripts",sep=""))){system(paste("mkdir ",out_dir,"/cdsInTranscripts",sep=""))}
  Stamp1=paste(out_dir,"/cdsInTranscripts/","cdsInTranscripts.finished",sep="")
  Stamp2=paste(save_dir,"/cdsInTranscripts/","cdsInTranscripts.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 5: cdsInTranscripts FINISHED")
  }else{
    print("Step 5: cdsInTranscripts START")
    cdsInTranscripts(transcripts.fna=paste(save_dir,"/seqNR/rep.fasta_SimpleIDs",sep=""),
                     out_dir=paste(out_dir,"/cdsInTranscripts",sep=""),
                     dmdb=uniprot.dmdb, # DIAMOND protein db, uniprot
                     pfam=pfam, # Pfam-A.hmm
                     threads=threads)
    cmd=paste("maker_map_ids",
              "--iterate 0",
              "--prefix pep",
              paste(out_dir,"/cdsInTranscripts/TransDecoder.gff3",sep=""),">",
              paste(out_dir,"/cdsInTranscripts/MAKER.name.map",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("map_gff_ids",
              paste(out_dir,"/cdsInTranscripts/MAKER.name.map",sep=""),
              paste(out_dir,"/cdsInTranscripts/TransDecoder.gff3",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("gffread",
              "-O",paste(out_dir,"/cdsInTranscripts/TransDecoder.gff3",sep=""),
              "-S",
              "-g",paste(out_dir,"/cdsInTranscripts/transcripts.fna",sep=""),
              "-x",paste(out_dir,"/cdsInTranscripts/cds.fna",sep=""),
              "-y",paste(out_dir,"/cdsInTranscripts/pep.faa",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system(paste("touch",Stamp1,sep=" "))
  }
}

#####
# diamond_p_megan
#####
if (6 %in% Step){
  print("Step 6: diamond_p_megan")
  if (!file.exists(paste(out_dir,"/diamond_p_megan",sep=""))){system(paste("mkdir ",out_dir,"/diamond_p_megan",sep=""))}
  Stamp1=paste(out_dir,"/diamond_p_megan/","diamond_p_megan.finished",sep="")
  Stamp2=paste(save_dir,"/diamond_p_megan/","diamond_p_megan.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 6: diamond_p_megan FINISHED")
  }else{
    print("Step 6: diamond_p_megan START")
    diamond_p_megan(query.faa=paste(save_dir,"/cdsInTranscripts/pep.faa",sep=""),
                    diamond.db=nr.dmdb,
                    megan.db=megan.db,
                    out_prefix=paste(out_dir,"/diamond_p_megan/mixAssembly",sep=""),
                    threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}

#####
# minimap2
#####
if (7 %in% Step){
  print("Step 7: minimap2")
  if (!file.exists(paste(out_dir,"/minimap2",sep=""))){system(paste("mkdir ",out_dir,"/minimap2",sep=""))}
  for (i in 1:nrow(tab)){
    sp=tab[i,1]
    fq1.lst=tab[i,2]
    fq2.lst=tab[i,3]
    Stamp1=paste(out_dir,"/minimap2/",sp,".minimap2.finished",sep="")
    Stamp2=paste(save_dir,"/minimap2/",sp,".minimap2.finished",sep="")
    if (file.exists(Stamp1) | file.exists(Stamp2)){
      print(paste(sp," minimap2 FINISHED",sep=""))
    }else{
      print(paste(sp," minimap2 START"),sep="")
      fq1=unlist(strsplit(fq1.lst,","))
      fq2=unlist(strsplit(fq2.lst,","))
      for (i in 1:length(fq1)){
        cmd=paste("gzip -c -d ",fq1[i]," >> ",out_dir,"/minimap2/",sp,".read1.fq",sep="")
        print(cmd);system(cmd,wait=TRUE)
        cmd=paste("gzip -c -d ",fq2[i]," >> ",out_dir,"/minimap2/",sp,".read2.fq",sep="")
        print(cmd);system(cmd,wait=TRUE)
      }
      system(paste("gzip ",out_dir,"/minimap2/",sp,".read1.fq",sep=""))
      system(paste("gzip ",out_dir,"/minimap2/",sp,".read2.fq",sep=""))
      minimap2(long_reads=paste(out_dir,"/minimap2/",sp,".read1.fq.gz"," ",
                                out_dir,"/minimap2/",sp,".read2.fq.gz",
                                sep=""), # space-separated list for PE
               lr_type="sr", # long read type. 
               # "map-pb" for PacBio
               # "map-hifi" for HiFi
               # "map-ont" for ONT reads.
               # "sr" for NGS
               # "asm5" for accurate reads diverging <5% to assembly
               assembly=paste(save_dir,"/seqNR/rep.fasta_SimpleIDs",sep=""),
               out_prefix=paste(out_dir,"/minimap2/",sp,sep=""),
               threads=threads)
      coverage(bam=paste(out_dir,"/minimap2/",sp,".bam",sep=""),
               output=paste(out_dir,"/minimap2/",sp,".tsv",sep=""))
      cmd=paste("seqkit stats ",out_dir,"/minimap2/",sp,".read*.fq.gz  -T > ",out_dir,"/minimap2/",sp,".depth.tsv",sep="")
      print(cmd);system(cmd,wait=TRUE)
      system(paste("rm ",out_dir,"/minimap2/",sp,".read1.fq.gz",sep=""))
      system(paste("rm ",out_dir,"/minimap2/",sp,".read2.fq.gz",sep=""))
      system(paste("touch",Stamp1,sep=" "))
    }
  }
}

#####
# Tables
#####
# paste(save_dir,"/diamond_p_megan/mixAssembly_taxon.tsv",sep="")
# paste(save_dir,"/minimap2/",sp,".tsv",sep="")
# paste(save_dir,"/cdsInTranscripts/TransDecoder.gff3",sep="")
# paste(save_dir,"/minimap2/",sp,".depth.tsv",sep="")


