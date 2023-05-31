# Metagenome/metatranscriptome profiling

# Map DNA reads to reference
# Dependencies: Minimap2, SAMtools
minimap2=function(long_reads=long_reads, # space-separated list for PE
                  lr_type=lr_type, # long read type. 
                  # "map-pb" for PacBio
                  # "map-hifi" for HiFi
                  # "map-ont" for ONT reads.
                  # "sr" for NGS
                  # "asm5" for accurate reads diverging <5% to assembly
                  assembly=assembly,
                  out_prefix=out_prefix,
                  threads=threads){
  cmd=paste("minimap2",
            "-ax",lr_type,
            "-t",threads,
            "--secondary=no","--MD","-L",
            assembly,long_reads,"|",
            "samtools","view","-@",threads,"-bS","|",
            "samtools","sort",
            "-@",threads,
            "-o",paste(out_prefix,".bam",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("samtools","index",paste(out_prefix,".bam",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Map short RNA reads to reference genome.
# Compress SAM to BAM, sort BAM, index BAM.
# Dependencies: Hisat2, SAMtools
Hisat = function(fq1=fq1,fq2=fq2, # Input fq files. Make fq2="None" if single-end.
                 # Comma-separated list.
                 fna=fna, # genome (soft masked for training AUGUSTUS)
                 index=index, # Basename of Hisat2 index of reference genome.
                 out_prefix=out_prefix, # Prefix of output BAM file.
                 threads=threads){
  threads = as.character(threads)
  
  if (!file.exists(paste(index,".1.ht2",sep=""))){
    cmd = paste("hisat2-build",
                fna,
                index,sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  if (fq2!="None"){ # pair-end
    cmd = paste("hisat2",
                "--dta",
                "-x",index,
                "-p",threads,
                "-1",fq1,
                "-2",fq2,
                "|",
                "samtools","view",
                "-@",threads,
                "-bS",
                "|",
                "samtools","sort",
                "-@",threads,
                "-o",paste(out_prefix,".bam",sep=""),
                sep=" ")
  }else{ # single pair
    cmd = paste("hisat2",
                "--dta",
                "-x",index,
                "-p",threads,
                "-U",fq1,
                "|",
                "samtools","view",
                "-@",threads,
                "-bS",
                "|",
                "samtools","sort",
                "-@",threads,
                "-o",paste(out_prefix,".bam",sep=""),
                sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("samtools","index",
            paste(out_prefix,".bam",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(paste(out_prefix,".bam",sep=""))
}

# Compute coverage of each scaffold
# Dependencies: SAMtools
coverage=function(bam=bam,
                  output=output){
  # tabular
  # #rname  Reference name / chromosome
  # startpos	Start position
  # endpos	End position (or sequence length)
  # numreads	Number reads aligned to the region (after filtering)
  # covbases	Number of covered bases with depth >= 1
  # coverage	Percentage of covered bases [0..100]
  # meandepth	Mean depth of coverage
  # meanbaseq	Mean baseQ in covered region
  # meanmapq	Mean mapQ of selected reads
  
  cmd=paste("samtools","coverage",
            "-o",output,
            bam,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# GC content & length of each sequence
# Dependencies: seqkit
GC_length=function(fna=fna,out=out,threads=threads){
  cmd=paste("printf","'rname\tlength\tGC.content'",">",out)
  system(cmd,wait=TRUE)
  
  cmd=paste("seqkit","fx2tab",
            "--name --only-id --gc --length",
            "--threads",as.character(threads),
            fna,">>",out,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Extract reads mapped/unmapped to an assembly via SAMtools.
# Dependencies: SAMtools
Extract_fq=function(bam=bam,
                    paired=paired, # logical. T for paired. 
                    mapped=mapped, # logical. T for mapped
                    out_prefix=out_prefix,
                    reads_format=reads_format, # "fq" or "fa"
                    threads=threads){
  threads=as.character(threads)
  if (paired){
    if (mapped){flag="-f 2 -F 256"}else{flag="-f 12 -F 256"}
    if (reads_format=="fq"){
      samtools="samtools fastq";reads_format="fq.gz"}else{samtools="samtools fasta"}
    cmd=paste(samtools,
              flag,
              "-@",threads,
              "-1",paste(out_prefix,".1.",reads_format,sep=""),
              "-2",paste(out_prefix,".2.",reads_format,sep=""),
              bam,
              sep=" ")
  }
  if (!paired){
    if (mapped){flag="-f 2"}else{flag="-G 2"}
    if (reads_format=="fq"){
      samtools="samtools fastq";reads_format="fq.gz"}else{samtools="samtools fasta"}
    cmd=paste(samtools,
              flag,
              "-@",threads,
              bam,">",
              paste(out_prefix,".",reads_format,sep=""),
              sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
}

# Remove rRNA from metatranscriptome by SortMeRNA
# Dependencies: SortMeRNA
sortmerna=function(fq1=fq1,fq2="none",
                   reference=reference, # smr_v4.3_default_db.fasta
                   out_dir=out_dir,
                   threads=threads){
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("sortmerna",
            "--ref",reference,
            "--reads",fq1,
            sep=" ")
  if (fq2!="none"){
    cmd=paste(cmd,
              "--reads",fq2,
              "--out2 True",
              "--paired_in True",# paired-end reads as Aligned when either of them is Aligned.
              sep=" ")
  }
  cmd=paste(cmd,
            "--workdir",out_dir,
            "--threads",threads,
            "--fastx True",
            "--aligned",paste(getwd(),"/rRNA",sep=""),
            "--other",paste(getwd(),"/noRrna",sep=""),
            "--otu_map True",
            sep=" ") 
  print(cmd);system(cmd,wait=TRUE)
  
  system("rm -r idx/ kvdb/ readb/")
  
  setwd(wd)
}

# Meta-genome/transcriptome assembly by IDBA-UD
# Dependencies: IDBA-UD
idba_ud=function(fq1=fq1,fq2=fq2, # comma-list, paired
                 out_dir=out_dir,
                 threads=threads){
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  system("mkdir tmp")
  system("mkdir idba_ud")
  
  fq1=unlist(strsplit(fq1,","))
  fq2=unlist(strsplit(fq2,","))
  for (i in 1:length(fq1)){
    cmd=paste("gzip -c -d",fq1[i],">","tmp/read1.fq",sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("gzip -c -d",fq2[i],">","tmp/read2.fq",sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("fq2fa --merge --filter tmp/read1.fq tmp/read2.fq tmp/tmp.fa",sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system("cat tmp/tmp.fa >> tmp/read.fa")
    system("rm tmp/read1.fq tmp/read2.fq tmp/tmp.fa")
  }
  
  cmd=paste("idba_ud",
            "-r tmp/read.fa",
            "-o",paste(getwd(),"/idba_ud",sep=""),
            "--num_threads",threads,
            "--min_count 1",
            "--mink 20",
            "--maxk 120",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("rm -r tmp/")
  setwd(wd)
}

# Assemble meta-genome/transcriptome into proteins by Plass
# Dependencies: plass
plass=function(fq1=fq1,fq2=fq2, # comma-list
               threads=threads,
               out_dir=out_dir){
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  system("mkdir tmp")
  
  fq1=unlist(strsplit(fq1,","))
  fq2=unlist(strsplit(fq2,","))
  fq=""
  for (i in 1:length(fq1)){
    fq=paste(fq,fq1[i],fq2[i],sep=" ")
  }
  
  cmd=paste("plass assemble",fq,"plass.fa tmp",
            "--threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# SPAdes: Genome/Metagenome assembly.
# long reads can be provided for scaffolding NGS assembly
# Dependencies: SPAdes
SPAdes=function(fq1=fq1,fq2=fq2, 
                # comma-list or R vector
                # Input NGS reads (fq). Set fq2="none" if single-end.
                # PacBio CCS (HiFi) reads should be treated as fq1.
                contigs.fa="none", # Reliable contigs of the same genome.
                pacbio_clr="none",nanopore="none",sanger="none", # Long reads
                meta=meta, # Logical. If TRUE, run metaSPAdes.
                           # metaSPAdes supports only a single short-read library which has to be paired-end.
                rna=rna, # Logical. TRUE for rnaSPAdes
                bio=bio, # Logical. TRUE for biosyntheticSPAdes
                custom_hmms="none", # directory with custom hmms for biosyntheticSPAdes
                out_dir=out_dir,
                threads=threads,
                memory=500){
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  if (fq2!="none"){ # pair-end
    fq1=unlist(strsplit(fq1,","))
    fq2=unlist(strsplit(fq2,","))
    for (i in 1:length(fq1)){
      cmd=paste("gzip -c -d",fq1[i],">>","read1.fq",sep=" ")
      print(cmd);system(cmd,wait=TRUE)
      cmd=paste("gzip -c -d",fq2[i],">>","read2.fq",sep=" ")
      print(cmd);system(cmd,wait=TRUE)
    }
    system("gzip read1.fq")
    system("gzip read2.fq")
    
    cmd=paste("spades.py",
              "-t",threads,
              "-m",as.character(memory),
              "-1 read1.fq.gz",
              "-2 read2.fq.gz",
              "-o",out_dir,
              sep=" ")
  }else{
    fq1=unlist(strsplit(fq1,","))
    for (i in 1:length(fq1)){
      cmd=paste("gzip -c -d",fq1[i],">>","read1.fq",sep=" ")
      print(cmd);system(cmd,wait=TRUE)
    }
    system("gzip read1.fq")
    cmd=paste("spades.py",
              "-t",threads,
              "-m",as.character(memory),
              "-s read1.fq.gz",
              "-o",out_dir,
              sep=" ")
  }
  
  if (contigs.fa!="none"){cmd=paste(cmd,"--trusted-contigs",contigs.fa,sep=" ")}
  if (pacbio_clr!="none"){cmd=paste(cmd,"--pacbio",pacbio_clr,sep=" ")}
  if (nanopore!="none"){cmd=paste(cmd,"--nanopore",nanopore,sep=" ")}
  if (sanger!="none"){cmd=paste(cmd,"--sanger",sanger,sep=" ")}
  if (meta){cmd=paste(cmd,"--meta",sep=" ")}
  if (rna){cmd=paste(cmd,"--rna",sep=" ")}
  if (bio){cmd=paste(cmd,"--bio",sep=" ")}
  if (custom_hmms!="none"){cmd=paste(cmd,"--custom-hmms",custom_hmms,sep=" ")}
  
  print(cmd);system(cmd,wait=TRUE)
  
  system("rm read1.fq.gz read2.fq.gz")
  return(out_dir)
}

# Trinity: de novo assembly of transcriptome
trinity=function(fq1=fq1,fq2=fq2, # comma-list/R vector
                 threads=threads,
                 max_memory=max_memory, # 500G
                 out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  fq1=paste(fq1,collapse=",")
  fq2=paste(fq2,collapse=",")
  
  system("mkdir trinity")
  if (!file.exists("normalization.finished")){
  cmd=paste("Trinity",
            "--seqType","fq",
            "--left",fq1,
            "--right",fq2,
            "--CPU",threads,
            "--output",paste(out_dir,"/trinity",sep=""),
            "--max_memory",max_memory,
            "--no_run_inchworm",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system("touch normalization.finished")
  }
  if (!file.exists("inchworm.finished")){
    cmd=paste("Trinity",
              "--seqType","fq",
              "--left",fq1,
              "--right",fq2,
              "--CPU",threads,
              "--output",paste(out_dir,"/trinity",sep=""),
              "--max_memory",max_memory,
              "--no_run_chrysalis",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system("touch inchworm.finished")
  }
  if (!file.exists("chrysalis.finished")){
    cmd=paste("Trinity",
              "--seqType","fq",
              "--left",fq1,
              "--right",fq2,
              "--CPU",threads,
              "--output",paste(out_dir,"/trinity",sep=""),
              "--max_memory",max_memory,
              "--no_distributed_trinity_exec",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system("touch chrysalis.finished")
  }
  if (!file.exists("assembly.finished")){
  cmd=paste("Trinity",
            "--seqType","fq",
            "--left",fq1,
            "--right",fq2,
            "--CPU",threads,
            "--output",paste(out_dir,"/trinity",sep=""),
            "--max_memory",max_memory,
            "--full_cleanup",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system("touch assembly.finished")
  }
  setwd(wd)
}

# MMseq2: redundancy filter sequences with identical length and 100% length overlap
# Dependencies: mmseq2
seqNR=function(in.fasta=in.fasta,
               out_dir=out_dir,
               Identity=Identity, # [0.0,1.0]
               cov_mode=cov_mode, # 0: alignment covers ${coverage} of target and of query
                                  # 1: alignment covers ${coverage} of target
                                  # 2: alignment covers ${coverage} of query
                                  # 3: target is of ${coverage} query length
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
  
  setwd(wd)
}

# seqNR(in.fasta="/flash/BourguignonU/Cong/termite_genome_annotation/proteins.faa",
#       out_dir="/flash/BourguignonU/Cong/termite_genome_annotation/proteins.seqNR",
#       Identity=0.95, # [0.0,1.0]
#       cov_mode=0, # 0: alignment covers ${coverage} of target and of query
#                # 1: alignment covers ${coverage} of target
#                # 2: alignment covers ${coverage} of query
#                # 3: target is of ${coverage} query length
#       coverage=0.95, # [0.0,1.0]
#       threads=64)
  
# Find Coding Regions within Transcripts
# Dependencies: TransDecoder, diamond, hmmer, AGAT
cdsInTranscripts=function(transcripts.fna=transcripts.fna,
                          out_dir=out_dir,
                          dmdb=dmdb, # DIAMOND protein db, uniprot
                          pfam=pfam, # Pfam-A.hmm
                          threads=threads){
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  system(paste("cp",transcripts.fna,"transcripts.fna",sep=" "))
  
  cmd="TransDecoder.LongOrfs -t transcripts.fna"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("diamond blastp",
            "--query","transcripts.fna.transdecoder_dir/longest_orfs.pep",
            "--db",dmdb,
            "--max-target-seqs 1",
            "--outfmt 6",
            "--evalue 1e-5",
            "--threads",threads,
            "--out","blastp.outfmt6",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("hmmsearch",
            "--cpu",threads,
            "--domtblout pfam.domtblout",
            "-o","hmmsearch.out",
            pfam,
            "transcripts.fna.transdecoder_dir/longest_orfs.pep",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("TransDecoder.Predict",
            "-t transcripts.fna",
            "--retain_pfam_hits pfam.domtblout",
            "--retain_blastp_hits blastp.outfmt6",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system("mv transcripts.fna.transdecoder.gff3 TransDecoder.gff3")
  
  system("rm blastp.outfmt6")
  system("rm hmmsearch.out")
  system("rm pfam.domtblout")
  system("rm pipeliner.*.cmds")
  #system("rm transcripts.fna")
  system("rm -r transcripts.fna.transdecoder_dir")
  system("rm -r transcripts.fna.transdecoder_dir.__checkpoints")

  setwd(wd)
  return("TransDecoder.gff3")
}

# Protein-protein search by DIAMOND and assign to taxa by MEGAN
# Dependencies: DIAMOND
diamond_p_megan=function(query.faa=query.faa,
                         diamond.db=diamond.db, # nr
                         megan.db=megan.db,
                         out_prefix=out_prefix,
                         threads=threads){
  threads=as.character(threads)
  cmd=paste("diamond blastp",
            "-p",threads,
            "-d",diamond.db,
            "-q",query.faa,
            "-f 100",
            "--out",paste(out_prefix,".blast.daa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("daa-meganizer",
            "-i",paste(out_prefix,".blast.daa",sep=""),
            "-mdb",megan.db,
            "-t",threads,
            "-ram readCount",
            "-supp 0",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("daa2info",
            "-i",paste(out_prefix,".blast.daa",sep=""),
            "-o",paste(out_prefix,"_taxon.tsv",sep=""),
            "-r2c Taxonomy",
            "-n true",
            "-p true",
            "-r true",
            "-mro","true",
            "-u false",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}




# # IDBA-MT for metatranscriptome
# idba_mt=function(fq1=fq1,fq2=fq2, # comma-list
#                  contigs=contigs, # IDBA-UD
#                  threads=threads,
#                  out_dir=out_dir){
#   threads=as.character(threads)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
#   wd=getwd();setwd(out_dir)
#   system("mkdir tmp")
#   system("mkdir idba_mt")
#   
#   fq1=unlist(strsplit(fq1,","))
#   fq2=unlist(strsplit(fq2,","))
#   for (i in 1:length(fq1)){
#     cmd=paste("gzip -c -d",fq1[i],">","tmp/read1.fq",sep=" ")
#     print(cmd);system(cmd,wait=TRUE)
#     cmd="fq2fa --filter tmp/read1.fq tmp/read1.fa" # tool from IDBA-UD
#     print(cmd);system(cmd,wait=TRUE)
#     system("rm tmp/read1.fq")
#     
#     cmd=paste("gzip -c -d",fq2[i],">","tmp/read2.fq",sep=" ")
#     print(cmd);system(cmd,wait=TRUE)
#     cmd="fq2fa --filter tmp/read2.fq tmp/read2.fa"
#     print(cmd);system(cmd,wait=TRUE)
#     system("rm tmp/read2.fq")
#     
#     ID1=system("grep '>' tmp/read1.fa | sed 's/>//'",intern=TRUE)
#     ID2=system("grep '>' tmp/read2.fa | sed 's/>//'",intern=TRUE)
#     IDs=intersect(ID1,ID2)
#     writeLines(IDs,"tmp/IDs.lst")
#     
#     cmd=paste("seqkit grep --threads ",threads," -n -w0 -f tmp/IDs.lst tmp/read1.fa >> tmp/1.fa",sep="")
#     print(cmd);system(cmd,wait=TRUE)
#     cmd=paste("seqkit grep --threads ",threads," -n -w0 -f tmp/IDs.lst tmp/read2.fa >> tmp/2.fa",sep="")
#     print(cmd);system(cmd,wait=TRUE)
#     
#     system("rm tmp/read1.fa tmp/read2.fa tmp/IDs.lst")
#   }
#   
#   cmd=paste("idba-mt",
#             "-t tmp/1.fa",
#             "-f tmp/2.fa",
#             "-O",paste(getwd(),"/idba_mt",sep=""),
#             "-c",contigs,
#             "-r 300",# read length
#             "-i 200",# insert size
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   system("rm -r tmp")
#   setwd(wd)
# }
# # # IDBA-MTP for metatranscriptome
# # Dependencies: ISBA-MTP, seqkit
# idba_mtp=function(fq1=fq1,fq2=fq2, # comma-list
#                   proteins=proteins,
#                   score_mat=score_mat,# https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#                   out_dir=out_dir){
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
#   wd=getwd();setwd(out_dir)
#   system("mkdir tmp")
#   system("mkdir idba_mtp")
#   
#   fq1=unlist(strsplit(fq1,","))
#   fq2=unlist(strsplit(fq2,","))
#   for (i in 1:length(fq1)){
#     cmd=paste("gzip -c -d",fq1[i],">","tmp/read1.fq",sep=" ")
#     print(cmd);system(cmd,wait=TRUE)
#     cmd="fq2fa tmp/read1.fq tmp/read1.fa" # tool from IDBA-UD
#     print(cmd);system(cmd,wait=TRUE)
#     system("rm tmp/read1.fq")
#     
#     cmd=paste("gzip -c -d",fq2[i],">","tmp/read2.fq",sep=" ")
#     print(cmd);system(cmd,wait=TRUE)
#     cmd="fq2fa tmp/read2.fq tmp/read2.fa"
#     print(cmd);system(cmd,wait=TRUE)
#     system("rm tmp/read2.fq")
#     
#     N1=system("seqkit grep -s -r -p 'N' tmp/read1.fa | grep '>' | sed 's/>//'",intern=TRUE)
#     N2=system("seqkit grep -s -r -p 'N' tmp/read2.fa | grep '>' | sed 's/>//'",intern=TRUE)
#     N=union(N1,N2)
#     
#     IDs=system("grep '>' tmp/read1.fa | sed 's/>//'",intern=TRUE)
#     IDs=IDs[!(IDs %in% N)]
#     writeLines(IDs,"tmp/IDs.lst")
#     
#     cmd="seqkit grep -f tmp/IDs.lst tmp/read1.fa >> tmp/1.fa"
#     print(cmd);system(cmd,wait=TRUE)
#     cmd="seqkit grep -f tmp/IDs.lst tmp/read2.fa >> tmp/2.fa"
#     print(cmd);system(cmd,wait=TRUE)
#     
#     system("rm tmp/read1.fa tmp/read2.fa tmp/IDs.lst")
#   }
#   
#   cmd=paste("idba_mtp",
#             "-t tmp/1.fa",
#             "-f tmp/2.fa",
#             "-O idba_mtp/",
#             "-p",proteins,
#             "-m",score_mat,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   system("rm -r tmp")
#   setwd(wd)
# }