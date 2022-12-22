# Genome annotation (protein-coding genes)

# Split paired fq
# Dependencies: seqkit
SplitFQ=function(fq1=fq1,fq2=fq2,
                 out_dir=out_dir,
                 part=part,
                 threads=threads){
  threads=as.character(threads)
  part=as.character(part)
  
  cmd=paste("seqkit split2",
            "--read1",fq1,
            "--read2",fq2,
            "--by-part",part,
            "--out-dir",out_dir,
            "--threads",threads)
  # out_dir/*part_001*
  print(cmd);system(cmd,wait=TRUE)
}

# Split fasta
# Dependencies: seqkit
SplitFA=function(fasta=fasta,
                 out_dir=out_dir,
                 part=part,
                 threads=threads){
  threads=as.character(threads)
  part=as.character(part)
  
  cmd=paste("seqkit split2",
            fasta,
            "--by-part",part,
            "--out-dir",out_dir,
            "--threads",threads)
  print(cmd);system(cmd,wait=TRUE)
}

# Simplify sequence ID in genome
# Dependencies: simplifyFastaHeaders.pl from AUGUSTUS
SimplifyID=function(fna=fna,
                    common_pattern=common_pattern # common pattern in simplified sequence IDs
                    ){
  cmd=paste("simplifyFastaHeaders.pl",
            fna,
            common_pattern,
            paste(fna,"_SimpleIDs",sep=""), # FASTA with simplified IDs
            paste(fna,"_IDconvert.tsv",sep="")) # old ID to new ID (tabular)
  print(cmd);system(cmd,wait=TRUE)
  return(paste(fna,"_SimpleIDs",sep=""))
}

# RepeatModeler
# Dependencies: RepeatModeler, Singularity
repeatmodeler=function(fna=fna,
                       out_dir=out_dir,
                       out_prefix=out_prefix,
                       Threads=Threads){
  pwd_begin=getwd();setwd(out_dir)
  Threads=as.character(Threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  
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

# RepeatMasker
# Dependencies: Singularity, RepeatMasker
repeatmasker=function(fna=fna,# Fasta file of genome.
                    out_dir=out_dir,
                    RepeatLib.fa=RepeatLib.fa,
                    Threads=Threads){
  pwd_begin=getwd()
  setwd(out_dir)
  
  Threads=as.character(Threads)
  out_dir=sub("/$","",out_dir)
  
  system(paste("cp",fna,out_dir,sep=" "),wait=TRUE)
  fna_name=unlist(strsplit(fna,"/"));fna_name=fna_name[length(fna_name)]
  fna=paste(out_dir,"/",fna_name,sep="")
  
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
            "-dir",out_dir,
            "-gff",
            fna,
            sep=" ")
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
  
  setwd(pwd_begin)
  
  # masked genome
  return(paste(out_dir,"/",fna,".masked",sep=""))
}

# GenomeThreader: protein-genome spliced alignments
# Dependencies: GenomeThreader, scripts from EvidenceModeler
gth=function(fna=fna, # masked genome
             faa=faa, # protein sequences
             out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  cmd=paste("cp",fna,".",sep=" ");system(cmd,wait=TRUE)
  fna=unlist(strsplit(fna,"/"));fna_file=fna[length(fna)]
  cmd=paste("cp",faa,".",sep=" ");system(cmd,wait=TRUE)
  faa=unlist(strsplit(faa,"/"));faa_file=faa[length(faa)]
  
  cmd=paste("gth",
            "-genomic",fna_file,
            "-protein",faa_file,
            "-gff3out -intermediate",
            "-o","protein_alignment.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",fna_file,sep=" "),wait=TRUE)
  system(paste("rm",faa_file,sep=" "),wait=TRUE)
  system(paste("rm"," ",fna_file,"*",sep=""),wait=TRUE)
  system(paste("rm"," ",faa_file,"*",sep=""),wait=TRUE)
  # EvidenceModeler
  cmd="genomeThreader_to_evm_gff3.pl protein_alignment.gff3 >protein_alignment_evm.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
  return(out_dir)
}

# Compute gene models (gff3) by miniprot
# Get gene models (gff3) to protein alignments (gff3)
# Dependencies: miniprot
miniprot=function(fna=fna,
                  faa=faa,
                  threads=threads,
                  out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  # Compute gene models by spliced alignments
  cmd=paste("miniprot",
            "-t",threads,
            "-d","genome.mpi",
            fna,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("miniprot",
            "-t",threads,
            "--gff",
            "genome.mpi",
            faa,
            ">","miniprot.gff3", # mRNA, CDS, stip_codon features
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  # separate genes by blank lines
  cmd="sed 's/.*##PAF.*//' miniprot.gff3 > gene_models.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  # Convert gene models to protein alignments
  system("echo '##gff-version 3' > protein_align.gff3")
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3==\"CDS\") $3=\"match\"; if ($3==\"match\") $8=\".\"; if ($3==\"match\") print $0}'",
            "gene_models.gff3",">>","protein_align.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="sed -i 's/Parent=/ID=/' protein_align.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
}

# Prothints+GeneMark-EP
# Dependencies: GeneMark
gm_ep = function(genome=genome,
                 faa=faa, # reference proteins
                 threads=threads,
                 out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  cmd=paste("gmes_petap.pl","--verbose",
            "--format GFF3",
            "--EP",
            "--dbep",faa,
            "--soft_mask","auto",
            "--cores",threads,
            "--sequence",genome,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
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
            "-@",threads,
            sep=" ")
  print(cmd);system(cmd)

  return(paste(out_prefix,".bam",sep=""))
}

# Compute coverage of each scaffold from BAM
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

# SAMtools: Merge BAM.
MergeBAM=function(BAMs=BAMs,# SPACE-separated list of bam files.
                  out_prefix=out_prefix,
                  Threads=Threads){
  Threads=as.character(Threads)
  
  cmd=paste("samtools","merge",
            "-@",Threads,
            "-",
            BAMs,"|",
            "samtools","sort",
            "-@",threads,
            "-o",paste(out_prefix,".bam",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("samtools","index",
            paste(out_prefix,".bam",sep=""),
            "-@",threads,
            sep=" ")
  print(cmd);system(cmd)
  
  return(0)
}

# Assemble transcripts from BAM.
# Dependencies: StringTie, script from EvidenceModeler
StringTie = function(input_bam=input_bam, # Input BAM.
                     output_prefix=out_prefix,
                     threads=threads){
  threads = as.character(threads)
  
  cmd = paste("stringtie",
              "-p", threads,
              "-o", paste(output_prefix,".gtf",sep=""),
              input_bam,
              sep=" ")
  print(cmd);system(cmd,wait = T)
  
  cmd=paste("cufflinks_gtf_to_alignment_gff3.pl", # evidencemodeler
            paste(output_prefix,".gtf",sep=""),">",
            paste(output_prefix,"_tmp.gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($2==\"Cufflinks\")  $2=\"StringTie\";print$0}'",
            paste(output_prefix,"_tmp.gff3",sep=""),">",
            paste(output_prefix,".gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system(paste("rm"," ",output_prefix,"_tmp.gff3",sep=""),wait=TRUE)
  return(paste(output_prefix,".gtf",sep=""))
}

# TransDecoder: Find Coding Regions within Transcripts
# Dependencies: TransDecoder, diamond, hmmer, AGAT
TransDecoder=function(gtf=gtf, # gtf from StringTie
                      fna=fna,
                      out_dir=out_dir,
                      dmdb=dmdb, # DIAMOND protein db, uniprot
                      pfam=pfam, # Pfam-A.hmm
                      threads=threads){
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("gtf_genome_to_cdna_fasta.pl",
            gtf,
            fna,
            ">",
            "transcripts.fna",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)

  cmd=paste("gtf_to_alignment_gff3.pl",
            gtf,
            ">",
            "temp.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)

  cmd=paste("TransDecoder.LongOrfs",
            "-t","transcripts.fna",
            sep=" ")
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
  
  cmd=paste("cdna_alignment_orf_to_genome_orf.pl",
            "transcripts.fna.transdecoder.gff3",
            "temp.gff3",
            "transcripts.fna",
            ">",
            "TransDecoder.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("agat_convert_sp_gff2gtf.pl",
            "--gff TransDecoder.gff3",
            "-o TransDecoder.gtf",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
  return("TransDecoder.gff3")
}

# Extract exon coordinates from transcript alignment (gff3) and randomly select 1000 genes
# Run GlimmerHMM
# Dependencies: glimmerhmm, scripts from EvidenceModeler
glimmerhmm_RNA=function(transcript_align.gff3=transcript_align.gff3, # get exon/match features
                                                                     # empty line separate genes
                        genome=genome,
                        out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  # Get exon coordinates
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($0==\"\" || $3==\"exon\" || $3==\"match\") print $0}'",
            transcript_align.gff3,"> tmp.gff3")
  print(cmd);system(cmd,wait=TRUE)
  transcript_align=readLines("tmp.gff3")
  training=sapply(transcript_align,
                  function(st){
                    if (st==""){
                      return("")
                    }else{
                      st=unlist(strsplit(st,"\t"))
                      seq=st[1];strand=st[7]
                      if (strand=="+"){start=st[4];end=st[5]}
                      if (strand=="-"){start=st[5];end=st[4]}
                      return(paste(seq,start,end,sep=" "))
                    }
                  })
  
  # How many genes provided?
  j=0
  for (i in 1:length(training)){
    if (training[i]==""){
      j=j+1
      cat(paste(as.character(i),"\n",sep=""),file="list",append=TRUE)
    }
  }
  
  if (j<=1000){
    writeLines(training,"exon_file")
  }else{ # randomly select 1000 genes for training
    index=sort(runif(1000,1,j))
    is=readLines("list")[index]
    sapply(is,function(i){
      i=as.numeric(i)
      for (k in (i-1):1){
        if (training[k]==""){break}
      }
      gene=training[(k-1):i]
      for (m in gene){
        cat(paste(m,"\n",sep=""),file="exon_file",append=TRUE)
      }
    })
  }
  
  cmd=paste("trainGlimmerHMM",
            genome,"exon_file",
            "-d training",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("glimmerhmm_linux_x86_64",
            genome,"training","-g",
            "> GlimmerHMM.gff3",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="glimmerHMM_to_GFF3.pl GlimmerHMM.gff3 > GlimmerHMM_evm.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  system("rm list",wait=TRUE)
  
  setwd(wd_begin)
}

# Get hints (introns) from RNA-mapping BAM
# Extract and filter hints (introns) from RNA-mapping BAM 
# GeneMark-ET with hints (not robust to termite genomes, not converging for Incistermes)
# Get redundant training set (AUGUSTUS)
# Dependencies: GeneMark-ET, AUGUSTUS, SAMtools, scripts from BRAKER, stringr (R)
gm_et = function(genome=genome, # soft masked
                 bam=bam,
                 PE=TRUE, # logical. TRUE for pair-end RNA reads
                 out_dir=out_dir,
                 threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  # Filter BAM
  cmd=paste("samtools sort -n",
            "-@",threads,
            bam,
            # "|","samtools sort","-@",threads,
            ">","Aligned.out.ss.bam",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  if (PE){st="--paired --pairwiseAlignment"}else{st=""}
  cmd=paste("filterBam","--uniq",st,
            "--in","Aligned.out.ss.bam",
            "--out Aligned.out.ssf.bam",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("samtools sort",
            "-@",threads,
            "Aligned.out.ssf.bam",
            ">","Aligned.out.ss.bam",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)

  # Extract and filter hints
  cmd=paste("bam2hints","--intronsonly", # AUGUSTUS
            "--in=Aligned.out.ss.bam",
            "--out=introns.gff",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("filterIntronsFindStrand.pl", # BRAKER
            genome,"introns.gff","--score",
            ">","introns.f.gff",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # GeneMark-ET
  # generating trained gmhmm.mod for GeneMarkHMM
  cmd=paste("gmes_petap.pl","--verbose",  
            "--ET","introns.f.gff",      
            "--soft_mask","auto",
            "--cores",threads,
            "--sequence",genome,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("filterGenemark.pl", # BRAKER
            "--genemark=genemark.gtf",
            "--hints=introns.f.gff",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system("mv genemark.f.good.gtf training_AUGUSTUS.gtf")
  
  cmd="computeFlankingRegion.pl training_AUGUSTUS.gtf" # AUGUSTUS
  print(cmd)
  flanking_DNA=system(cmd,wait=TRUE,intern=TRUE)[4]
  library(stringr)
  flanking_DNA=str_extract(flanking_DNA,": [0-9]*")
  flanking_DNA=sub(": ","",flanking_DNA)
  
  cmd=paste("gff2gbSmallDNA.pl genemark.gtf", # AUGUSTUS
            genome,flanking_DNA,"tmp.gb",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("filterGenesIn_mRNAname.pl", # AUGUSTUS
            "training_AUGUSTUS.gtf",
            "tmp.gb",">","training_AUGUSTUS.gb")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
}

# SNAP for gene prediction
# Dependencies: agat, SNAP, scripts from EvidenceModeler
snap=function(fna=fna,
              train.gff3=train.gff3,
              out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  cmd=paste("cp",fna,".",sep=" ");system(cmd,wait=TRUE)
  fna=unlist(strsplit(fna,"/"));fna_file=fna[length(fna)]
  cmd=paste("cp",train.gff3,".",sep=" ");system(cmd,wait=TRUE)
  train.gff3=unlist(strsplit(train.gff3,"/"));train.gff3_file=train.gff3[length(train.gff3)]
  
  cmd=paste("agat_convert_sp_gff2zff.pl",
            "--fasta",fna_file,
            "--gff",train.gff3_file,
            "-o train",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system("mv train.ann train.zff",wait=TRUE)
  
  cmd=paste("fathom -validate train.zff",fna_file,"> train.validate",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("fathom -categorize 100 train.zff",fna_file,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="fathom -export 100 -plus uni.*"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="fathom -validate export.ann export.dna"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="forge export.ann export.dna"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="hmm-assembler.pl model . > model.hmm"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("snap model.hmm",fna_file,">","SNAP.zff",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="zff2gff3.pl SNAP.zff > SNAP.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="SNAP_to_GFF3.pl SNAP.gff3 > SNAP_evm.gff3" # EvidenceModeler
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",fna_file,sep=" "),wait=TRUE)
  system(paste("rm",train.gff3_file,sep=" "),wait=TRUE)
  setwd(wd_begin)
}

# GeneMarkerHMM for gene prediction
# Dependencies: GeneMark
gmhmm=function(genome=genome,
               model_file=model_file, # gmhmm.mod from GeneMark-ET/EP
               threads=threads,
               out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  cmd=paste("gmes_petap.pl",
            "--predict_with",model_file,
            "--format","GFF3",
            "--soft_mask auto",
            "--cores",threads,
            "--sequence",genome,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
}

# Get redundant training set (AUGUSTUS) from gene models (gff3, CDS/exon features)
# Dependencies: stringr (R), scripts from AUGUSTUS, BRAKER 
augustus_trainset=function(gene_models.gff3=gene_models.gff3,
                           fna=fna,
                           out_dir=out_dir){
  library(stringr)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  cmd=paste("gth2gtf.pl", # AUGUSTUS 
            gene_models.gff3,
            "training_AUGUSTUS.gtf",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Compute flanking region (half of average mRNA size)?
  cmd="computeFlankingRegion.pl training_AUGUSTUS.gtf" # AUGUSTUS
  print(cmd)
  flanking_DNA=system(cmd,wait=TRUE,intern=TRUE)[4]
  flanking_DNA=str_extract(flanking_DNA,": [0-9]*")
  flanking_DNA=sub(": ","",flanking_DNA)
  
  # gb for training AUGUSTUS
  cmd=paste("gff2gbSmallDNA.pl training_AUGUSTUS.gtf", # AUGUSTUS
            fna,flanking_DNA,"training_AUGUSTUS.gb",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
  return(0)
}

# Remove redundant gene structures in training data of AUGUSTUS
# Dependencies: blast+, scripts from AUGUSTUS, stringr (R)
non_redundant=function(training_AUGUSTUS.gb=training_AUGUSTUS.gb, # GenBank 
                       training_AUGUSTUS.gtf=training_AUGUSTUS.gtf, # gtf 
                       genome=genome,
                       threads=threads,
                       out_dir=out_dir){
  library(stringr)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  # Depending on how training_AUGUSTUS.gb and training_AUGUSTUS.gtf are generated,
  # there might be fewer genes in training_AUGUSTUS.gb than training_AUGUSTUS.gtf
  
  # gene ID in training_AUGUSTUS.gb
  GenBank=readLines(training_AUGUSTUS.gb)
  gb=GenBank[grepl("/gene=",GenBank)]
  gb=sub("                     /gene=\"","",gb)
  gb=sub("\"","",gb)
  gb=unique(gb)
  write.table(gb,"traingenes.lst",sep="\t",col.names=FALSE,row.names=FALSE,quote=TRUE)
  
  # Gene ID to locus ID in training_AUGUSTUS.gb
  LOCUS=GenBank[grepl("LOCUS",GenBank)]
  LOCUS=sub("LOCUS       ","",LOCUS)
  LOCUS=sub("   [0-9]* bp  DNA","",LOCUS)
  write.table(data.frame(gb,LOCUS),"loci.lst",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  # Extract gtf of genes in training_AUGUSTUS.gb
  cmd=paste("grep",
            "-f","traingenes.lst",
            "-F",training_AUGUSTUS.gtf,
            ">",
            "bonafide.f.gtf",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("gtf2aa.pl", # AUGUSTUS
            genome,
            "bonafide.f.gtf",
            "prot.aa",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Simplify gene IDs
  cmd=paste("simplifyFastaHeaders.pl", # AUGUSTUS
            "prot.aa",
            "Gene",
            "prot.aa_SimpleIDs", # FASTA with simplified IDs
            "prot.aa_IDconvert.tsv", # new ID to old ID (tabular)
            sep=" ") 
  print(cmd);system(cmd,wait=TRUE)
  
  # remove genes share >80% similarity at protein level
  cmd=paste("aa2nonred.pl",
            paste("--cores=",threads,sep=""),
            "prot.aa_SimpleIDs","prot.nr.aa",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # non-redundant gene ID and locus ID
  cmd="grep '>' prot.nr.aa | perl -pe 's/>//' > nonred.lst"
  print(cmd);system(cmd,wait=TRUE)
  cmd="grep -f nonred.lst prot.aa_IDconvert.tsv | cut -f2 | perl -pe 's/>//' > original_geneID.lst"
  print(cmd);system(cmd,wait=TRUE)
  cmd="grep -f original_geneID.lst loci.lst | cut -f2 > nonred.loci.lst"
  print(cmd);system(cmd,wait=TRUE)
  
  # non-redundant training set for AUGUSTUS
  cmd=paste("filterGenesIn.pl",
            "nonred.loci.lst",
            training_AUGUSTUS.gb,
            ">",
            "training_AUGUSTUS_nr.gb",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
  return(paste(out_dir,"/training_AUGUSTUS_nr.gb",sep=""))
}

# AUGUSTUS training and gene prediction.
# Dependencies: AUGUSTUS, parallel (R), scripts from EvidenceModeler
augustus=function(fna=fna, # genome
                  species=species, # Species name for the trained model
                  training.gb=training.gb, # non-redundant training genes in GenBank format
                                           # training_AUGUSTUS_nr.gb
                  out_dir=out_dir,
                  threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("autoAug.pl",
            paste("--cpus=",threads,sep=""),
            paste("--workingdir=",out_dir,sep=""),
            paste("--species=",species,sep=""),
            paste("--genome=",fna,sep=""),
            paste("--trainingset=",training.gb,sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  parSapply(clus,
            1:(as.numeric(threads)-1),
            function(i){
              cmd=paste("sed -i 's/augustus/augustus --gff3=on/' autoAug/autoAugPred_abinitio/shells/aug",
                        as.character(i),sep="")
              print(cmd);system(cmd,wait=TRUE)
              cmd=paste("autoAug/autoAugPred_abinitio/shells/aug",
                        as.character(i),sep="")
              print(cmd);system(cmd,wait=TRUE)
            })
  
  cmd=paste("cat",
            "autoAug/autoAugPred_abinitio/shells/aug*.out",
            ">",
            "AUGUSTUS.gff3")
  print(cmd);system(cmd)
  
  cmd="cat AUGUSTUS.gff3 | join_aug_pred.pl > AUGUSTUS_nr.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  # EvidenceModeler
  cmd="augustus_GFF3_to_EVM_GFF3.pl AUGUSTUS_nr.gff3 > AUGUSTUS_evm.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# EvidenceModeler
evm=function(protein_alignments.gff3="none", # Absolute path. GFF3 for protein-genome spliced alignments
             gene_predictions.gff3="none", # Absolute path. GFF3 from ab initio gene prediction
             transcript_alignments.gff3="none", # Absolute path. GFF3 for transcript-genome spliced alignments
             evm_weights=evm_weights, # Absolute path
                                      # tabular table with fields:
                                      # class: ABINITIO_PREDICTION, PROTEIN, TRANSCRIPT, OTHER_PREDICTION
                                      # type: the 2nd column value of gff3 files provided.
                                      # weight: int. TRANSCRIPT >> PROTEIN > ABINITIO_PREDICTION
                                      # An example:
                                      # PROTEIN gth 2
                                      # ABINITIO_PREDICTION AUGUSTUS 1
                                      # ABINITIO_PREDICTION Glimmer 1
                                      # TRANSCRIPT  StringTie 10
                                      # OTHER_PREDICTIONS TransDecoder  10
                                      # gff3 of AUGUSTUS, Glimmer and TransDecoder should be merged 
                                      # into a single gff3 whose 2nd column is AUGUSTUS/Glimmer/TransDecoder
             genome=genome,
             out_dir=out_dir,
             threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  system("mkdir temp_dir")
  setwd("temp_dir")
  
  cmd_part=""
  if (protein_alignments.gff3!="none"){cmd_part=paste(cmd_part,
                                                      "--protein_alignments"," ",
                                                      protein_alignments.gff3,
                                                      sep="")}
  if (gene_predictions.gff3!="none"){cmd_part=paste(cmd_part," ",
                                                    "--gene_predictions"," ",
                                                    gene_predictions.gff3,
                                                    sep="")}
  if (transcript_alignments.gff3!="none"){cmd_part=paste(cmd_part," ",
                                                         "--transcript_alignments"," ",
                                                         transcript_alignments.gff3,
                                                         sep="")}
  
  cmd=paste("partition_EVM_inputs.pl",
            "--genome",genome,
            cmd_part,
            "--segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("write_EVM_commands.pl",
            "--genome",genome,
            "--weights",evm_weights,
            cmd_part,
            "--output_file_name evm.out  --partitions partitions_list.out >  commands.list",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("split",
            "-d",
            "-n",paste("l/",as.character(as.numeric(threads)-1),sep=""),
            "commands.list",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  parSapply(clus,
            1:(as.numeric(threads)-1),
            function(i){
              scr=system("ls x*",wait=TRUE,intern=TRUE)[i]
              cmd=paste("bash"," ",scr,sep="")
              print(cmd);system(cmd,wait=TRUE)})
  
  cmd="recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("convert_EVM_outputs_to_GFF3.pl",
            "--partitions partitions_list.out --output evm.out",
            "--genome",genome,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="find . -regex '.*evm.out.gff3' -exec cat {} \\; > EvidenceModeler.gff3"
  cat(cmd,file="merge_gff3.sh")
  system("bash merge_gff3.sh",wait=TRUE)
  
  cmd="mv evm.out ../"
  print(cmd);system(cmd,wait=TRUE)
  cmd="mv partitions_list.out ../"
  print(cmd);system(cmd,wait=TRUE)
  cmd="mv commands.list ../"
  print(cmd);system(cmd,wait=TRUE)
  cmd="mv EvidenceModeler.gff3 ../"
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(out_dir)
  system("rm -r temp_dir")
  
  setwd(wd)
}

# Update gff3 by PASA
# Dependencies: pasa,sqlite3,scripts from PASAPipeline
pasa=function(genome=genome,
              transcripts=transcripts,
              original.gff3=original.gff3, # only protein-coding gene
              pasa.alignAssembly.conf=pasa.alignAssembly.conf, # Path to PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt
              pasa.annotationCompare.conf=pasa.annotationCompare.conf, # Path to PASApipeline/pasa_conf/pasa.annotationCompare.Template.txt
              threads=threads,
              out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",genome,".",sep=" ");system(cmd,wait=TRUE)
  genome=unlist(strsplit(genome,"/"))
  genome=paste(getwd(),"/",genome[length(genome)],sep="")
  cmd=paste("cp",transcripts,".",sep=" ");system(cmd,wait=TRUE)
  transcripts=unlist(strsplit(transcripts,"/"))
  transcripts=paste(getwd(),"/",transcripts[length(transcripts)],sep="")
  cmd=paste("cp",original.gff3,".",sep=" ");system(cmd,wait=TRUE)
  original.gff3=unlist(strsplit(original.gff3,"/"))
  original.gff3=paste(getwd(),"/",original.gff3[length(original.gff3)],sep="")
  
  #cmd="sqlite3 pasa.sqlite.db"
  #print(cmd);system(cmd,wait=TRUE)
  if (!file.exists(paste(getwd(),"/pasa.sqlite.db",sep=""))){
  cmd=paste("cp",pasa.alignAssembly.conf,"./alignAssembly.config",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  a=readline("./alignAssembly.config")
  a[5]=paste("DATABASE=",getwd(),"/pasa.sqlite.db",sep="")
  writeLines(a,"./alignAssembly.config")
  
  cmd=paste("Launch_PASA_pipeline.pl",
            "-c",paste(out_dir,"/alignAssembly.config",sep=""),
            "-C -R",
            "-g",genome,
            "-t",transcripts,
            "--ALIGNERS minimap2",
            "--CPU",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("Load_Current_Gene_Annotations.dbi",
            "-c",paste(out_dir,"/alignAssembly.config",sep=""),
            "-g",genome,
            "-P",original.gff3,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("cp",pasa.annotationCompare.conf,"./annotationCompare.config",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  a=readline("./annotationCompare.config")
  a[5]=paste("DATABASE=",getwd(),"/pasa.sqlite.db",sep="")
  writeLines(a,"./annotationCompare.config")

  cmd=paste("Launch_PASA_pipeline.pl",
            "-c",paste(out_dir,"/annotationCompare.config",sep=""),
            "-A",
            "-g",genome,
            "-t",transcripts,
            "--CPU",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  re=system("ls pasa.sqlite.db.gene_structures_post_PASA_updates.[0-9]*.gff3",wait=TRUE,intern=TRUE)
  cmd=paste("mv",re,"PASA.gff3",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",genome,sep=" "),wait=TRUE)
  system(paste("rm",transcripts,sep=" "),wait=TRUE)
  system(paste("rm",original.gff3,sep=" "),wait=TRUE)
}

# Generate a comprehensive genomic database including:
# genes.gff3, transcripts_rep.fna, transcripts_iso.fna, cds_rep.fna, cds_iso.fna, proteins_rep.faa, proteins_iso.faa
# Dependencies: maker, gffread, seqkit
pasa_more=function(species=species,
                   genome=genome,
                   PASA.gff3=PASA.gff3,
                   out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",PASA.gff3,
            paste("./",species,"_genes.gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  gff3=paste("./",species,"_genes.gff3",sep="")
  
  cmd=paste("maker_map_ids",
            "--iterate 0",
            "--prefix",species,
            gff3,"> MAKER.name.map",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("map_gff_ids","MAKER.name.map",gff3,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system("rm MAKER.name.map",wait=TRUE)
  
  cmd=paste("gffread","-O",
            gff3,"-S",
            "-g",genome,
            "-w",paste(species,"_transcripts.fna",sep=""),
            "-x",paste(species,"_cds.fna",sep=""),
            "-y",paste(species,"_proteins.faa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  f=function(fasta,rep.fa,iso.fa){
    iso=system(paste("grep '>'",fasta,sep=" "),wait=TRUE,intern=TRUE)
    iso=sub(">","",iso)
    iso=sub(" .*$","",iso)
    writeLines(iso,"IDs")
    system("grep '.*-R[1-9]*$' IDs > iso_IDs",wait=TRUE)
    system("grep '.*-R0' IDs > rep_IDs",wait=TRUE)
    system("rm IDs",wait=TRUE)
    
    cmd=paste("seqkit grep",
              "-f rep_IDs",
              fasta,">",rep.fa,
              sep=" ")
    system(cmd,wait=TRUE)
    system("rm rep_IDs",wait=TRUE)
    
    cmd=paste("seqkit grep",
              "-f iso_IDs",
              fasta,">",iso.fa,
              sep=" ")
    system(cmd,wait=TRUE)
    system("rm iso_IDs",wait=TRUE)
  }
  
  f(paste(species,"_transcripts.fna",sep=""),
    paste(species,"_transcripts_rep.fna",sep=""),
    paste(species,"_transcripts_iso.fna",sep=""))
  f(paste(species,"_cds.fna",sep=""),
    paste(species,"_cds_rep.fna",sep=""),
    paste(species,"_cds_iso.fna",sep=""))
  f(paste(species,"_proteins.faa",sep=""),
    paste(species,"_proteins_rep.faa",sep=""),
    paste(species,"_proteins_iso.faa",sep=""))
  
  system(paste("rm",paste(species,"_transcripts.fna",sep=""),sep=" "),wait=TRUE)
  system(paste("rm",paste(species,"_cds.fna",sep=""),sep=" "),wait=TRUE)
  system(paste("rm",paste(species,"_proteins.faa",sep=""),sep=" "),wait=TRUE)
}

# Trinity: de novo assembly of transcriptome
trinity=function(fq1=fq1,fq2=fq2,
                 threads=threads,
                 max_memory=max_memory, # 500G
                 out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  wd=getwd()
  
  cmd=paste("Trinity",
            "--seqType","fq",
            "--left",fq1,
            "--right",fq2,
            "--CPU",threads,
            "--output",out_dir,
            "--max_memory",max_memory,
            "--full_cleanup",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}


















# MAKER
# Dependencies: MAKER, SNAP, RepeatMasker, blast, AUGUSTUS, exonerate, blast+, GAAS, stringr (R)
#               scripts from AUGUSTUS
maker=function(fna=fna, # masked
               transcript.file=transcript.file, # gff3/gtf with transcript/mRNA and exon/CDS features e.g. gtf from stringtie
               protein.file=protein.file, # gff3/gtf with transcript/mRNA and exon/CDS features
               threads=threads,
               augustus_species=augustus_species,
               out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  if (!file.exists("maker_1/MAKER.gff3")){
  if (file.exists("maker_1")){system("rm -r maker_1",wait=TRUE)}
  cmd=paste("agat_sp_alignment_output_style.pl",
            "-g",transcript.file,
            "-o","transcript_align.gff3")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("agat_sp_alignment_output_style.pl",
            "-g",protein.file,
            "-o","protein_align.gff3")
  print(cmd);system(cmd,wait=TRUE)
  
  transcript_align.gff3="transcript_align3.gff"
  protein_align.gff3="protein_align3.gff"
  
  system("mkdir maker_1",wait=TRUE);setwd("maker_1")
  cmd="maker -CTL"
  print(cmd);system(cmd,wait=TRUE)
  a=readLines("maker_opts.ctl")
  a[2]=paste("genome=",fna,sep="")
  a[18]=paste("est_gff=",out_dir,"/",transcript_align.gff3,sep="")
  a[23]=paste("protein_gff=",out_dir,"/",protein_align.gff3,sep="")
  a[41]="est2genome=1"
  a[42]="protein2genome=1"
  writeLines(a,"maker_opts.ctl")
  cmd=paste("maker",
            "-c",threads,
            "-RM_off","-quiet",
            "-base maker",
            "maker_opts.ctl maker_bopts.ctl maker_exe.ctl",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="gff3_merge -s -n -g -d maker.maker.output/maker_master_datastore_index.log > MAKER.gff3"
  print(cmd);system(cmd,wait=TRUE)
  setwd("..")
  }else{
  transcript_align.gff3="transcript_align3.gff"
  protein_align.gff3="protein_align3.gff"
  }
  
  if (!file.exists("snap_1/MAKER.hmm")){
  if (file.exists("snap_1")){system("rm -r snap_1",wait=TRUE)}
  system("mkdir snap_1",wait=TRUE);setwd("snap_1")
  cmd="maker2zff -x 0.25 -l 50 -d ../maker_1/maker.maker.output/maker_master_datastore_index.log"
  print(cmd);system(cmd,wait=TRUE)
  cmd="fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1"
  print(cmd);system(cmd,wait=TRUE)
  cmd="fathom genome.ann genome.dna -validate > validate.log 2>&1"
  print(cmd);system(cmd,wait=TRUE)
  cmd="fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1"
  print(cmd);system(cmd,wait=TRUE)
  cmd="fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1"
  print(cmd);system(cmd,wait=TRUE)
  cmd="forge export.ann export.dna > forge.log 2>&1"
  print(cmd);system(cmd,wait=TRUE)
  cmd="hmm-assembler.pl genome . > MAKER.hmm"
  print(cmd);system(cmd,wait=TRUE)
  setwd("..")
  }
  
  config=system("echo $AUGUSTUS_CONFIG_PATH",wait=TRUE,intern=TRUE)
  if (!file.exists(paste(config,"/species/",augustus_species,sep=""))){
  if (file.exists("augustus_1")){system("rm -r augustus_1",wait=TRUE)}
  system("mkdir augustus_1",wait=TRUE);setwd("augustus_1")
  library(stringr)
  cmd=paste("gth2gtf.pl", # AUGUSTUS 
            "../maker_1/MAKER.gff3",
            "training_AUGUSTUS.gtf",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="computeFlankingRegion.pl training_AUGUSTUS.gtf" # AUGUSTUS
  print(cmd)
  flanking_DNA=system(cmd,wait=TRUE,intern=TRUE)[4]
  flanking_DNA=str_extract(flanking_DNA,": [0-9]*")
  flanking_DNA=sub(": ","",flanking_DNA)
  cmd=paste("gff2gbSmallDNA.pl training_AUGUSTUS.gtf", # AUGUSTUS
            fna,flanking_DNA,"training_AUGUSTUS.gb",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  GenBank=readLines("training_AUGUSTUS.gb")
  gb=GenBank[grepl("/gene=",GenBank)]
  gb=sub("                     /gene=\"","",gb)
  gb=sub("\"","",gb)
  gb=unique(gb)
  write.table(gb,"traingenes.lst",sep="\t",col.names=FALSE,row.names=FALSE,quote=TRUE)
  LOCUS=GenBank[grepl("LOCUS",GenBank)]
  LOCUS=sub("LOCUS       ","",LOCUS)
  LOCUS=sub("   [0-9]* bp  DNA","",LOCUS)
  write.table(data.frame(gb,LOCUS),"loci.lst",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  cmd=paste("grep",
            "-f","traingenes.lst",
            "-F","training_AUGUSTUS.gtf",
            ">",
            "bonafide.f.gtf",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("gtf2aa.pl", # AUGUSTUS
            fna,
            "bonafide.f.gtf",
            "prot.aa",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("simplifyFastaHeaders.pl", # AUGUSTUS
            "prot.aa",
            "Gene",
            "prot.aa_SimpleIDs", # FASTA with simplified IDs
            "prot.aa_IDconvert.tsv", # new ID to old ID (tabular)
            sep=" ") 
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("aa2nonred.pl",
            paste("--cores=",threads,sep=""),
            "prot.aa_SimpleIDs","prot.nr.aa",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="grep '>' prot.nr.aa | perl -pe 's/>//' > nonred.lst"
  print(cmd);system(cmd,wait=TRUE)
  cmd="grep -f nonred.lst prot.aa_IDconvert.tsv | cut -f2 | perl -pe 's/>//' > original_geneID.lst"
  print(cmd);system(cmd,wait=TRUE)
  cmd="grep -f original_geneID.lst loci.lst | cut -f2 > nonred.loci.lst"
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("filterGenesIn.pl",
            "nonred.loci.lst",
            "training_AUGUSTUS.gb",
            ">",
            "training_AUGUSTUS_nr.gb",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("autoAug.pl",
            paste("--cpus=",threads,sep=""),
            paste("--workingdir=",out_dir,"/augustus_1",sep=""),
            paste("--species=",augustus_species,sep=""),
            paste("--genome=",fna,sep=""),
            paste("--trainingset=","training_AUGUSTUS_nr.gb",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  setwd("..")
  }
  
  if (!file.exists("maker_2/MAKER.gff3")){
  if (file.exists("maker_2")){system("rm -r maker_2",wait=TRUE)}
  system("mkdir maker_2",wait=TRUE);setwd("maker_2")
  cmd="maker -CTL"
  print(cmd);system(cmd,wait=TRUE)
  a=readLines("maker_opts.ctl")
  a[2]=paste("genome=",fna,sep="")
  a[18]=paste("est_gff=",out_dir,"/",transcript_align.gff3,sep="")
  a[23]=paste("protein_gff=",out_dir,"/",protein_align.gff3,sep="")
  a[34]=paste("snaphmm=",out_dir,"/snap_1/MAKER.hmm",sep="")
  a[36]=paste("augustus_species=",augustus_species,sep="")
  writeLines(a,"maker_opts.ctl")
  cmd=paste("maker",
            "-c",threads,
            "-RM_off","-quiet",
            "-base maker",
            "maker_opts.ctl maker_bopts.ctl maker_exe.ctl",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="gff3_merge -s -n -g -d maker.maker.output/maker_master_datastore_index.log > MAKER.gff3"
  print(cmd);system(cmd,wait=TRUE)
  setwd("..")
  }
  
  if (!file.exists("MAKER.gff3")){
  cmd=paste("maker_map_ids",
            "--prefix",augustus_species,
            "maker_2/MAKER.gff3 > MAKER.name.map",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="map_gff_ids MAKER.name.map maker_2/MAKER.gff3"
  print(cmd);system(cmd,wait=TRUE)
  system("mv maker_2/MAKER.gff3 ./MAKER.gff3")
  }
  
  setwd(wd)
}

















# Star: Map short RNA reads to reference genome. 
# SAMtools: Compress SAM to BAM, sort BAM, index BAM.
# ulimit -n
Star = function(fq1=fq1,fq2=fq2, # Input fq files. Make fq2="None" if single-end.
                # Comma-separated list.
                fna=fna, # genome (soft masked for training AUGUSTUS)
                out_dir=out_dir,
                out_basename=out_basename,
                threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  system("mkdir index")
  cmd=paste("STAR",
            "--runThreadN",threads,
            "--runMode genomeGenerate",
            "--genomeDir index",
            "--genomeFastaFiles",fna,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("awk '/^>/{if (l!=\"\") print(l);l=0;next}{l+=length($0)}END{print l}'",
            fna,sep=" ")
  print(cmd)
  size=system(cmd,wait=TRUE,intern=TRUE)
  size=sum(as.numeric(size))
  genomeSAindexNbases=14
  if (log2(size)/2-1<14){genomeSAindexNbases=floor(log2(size)/2-1)}
  
  cmd=paste("STAR",
            "--runMode alignReads",
            "--runThreadN",threads,
            "--genomeDir index",
            "--readFilesIn",fq1,fq2,
            "--outFileNamePrefix",out_basename,
            "--outSAMtype BAM SortedByCoordinate",
            "--genomeSAindexNbases",as.character(genomeSAindexNbases),
            "--readFilesCommand gunzip -c",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("samtools","index",
            paste(out_basename,"Aligned.sortedByCoord.out.bam",sep=""),
            "-@",threads,
            sep=" ")
  print(cmd);system(cmd)
  
  return(paste(out_dir,"/",out_basename,"Aligned.sortedByCoord.out.bam",sep=""))
}