# Genome annotation (protein-coding genes and pseudogenes)

# BBduk: trimming NGS reads
bbduk=function(fq1=fq1,fq2=fq2, # Input fq files. Set fq2="none" if single-end.
               clean_fq1=clean_fq1,clean_fq2=clean_fq2,
               adaptor.fa=adaptor.fa, #bbmap/resources/adapters.fa 
               threads=threads){
  if (fq2!="none"){
    cmd=paste("bbduk.sh",
              paste("-t=",threads,sep=""),
              paste("in1=",fq1,sep=""),
              paste("in2=",fq2,sep=""),
              paste("out1=",clean_fq1,sep=""),
              paste("out2=",clean_fq2,sep=""),
              "ktrim=r",
              paste("ref=",adaptor.fa,sep=""),
              "k=23 mink=7 hdist=1 tpe tbo maq=10",
              "qtrim=rl trimq=15 minlength=35",
              sep=" ")
  }else{
    cmd=paste("bbduk.sh",
              paste("-t=",threads,sep=""),
              paste("in=",fq1,sep=""),
              paste("out=",clean_fq1,sep=""),
              "ktrim=r",
              paste("ref=",adaptor.fa,sep=""),
              "k=23 mink=7 hdist=1",
              "tpe tbo maq=10 qtrim=rl trimq=15 minlength=35",
              sep=" ")
  }
  
  print(cmd);system(cmd,wait=TRUE)
}

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
#hisat2 --dta -x /flash/BourguignonU/Cong/termite_genome_annotation/Hisat/Coatitermes/Coatitermes -p 1 
#-1 /flash/BourguignonU/Cong/termite_genome_annotation/Hisat/Coatitermes/split/GT155_S90_R1_001.part_001.fastq.gz 
#-2 /flash/BourguignonU/Cong/termite_genome_annotation/Hisat/Coatitermes/split/GT155_S90_R2_001.part_001.fastq.gz | 
#samtools view -@ 1 -bS | samtools sort -@ 1 
#-o /flash/BourguignonU/Cong/termite_genome_annotation/Hisat/Coatitermes/1_GT155_S90_R1_001.fastq.gz.bam

#hisat2 --dta -x /flash/BourguignonU/Cong/termite_genome_annotation/Hisat/Coatitermes/Coatitermes -p 1 -1 /flash/BourguignonU/Cong/termite_genome_annotation/Hisat/Coatitermes/split/GT155_S90_R1_001.part_002.fastq.gz -2 /flash/BourguignonU/Cong/termite_genome_annotation/Hisat/Coatitermes/split/GT155_S90_R2_001.part_002.fastq.gz | samtools view -@ 1 -bS | samtools sort -@ 1 -o /flash/BourguignonU/Cong/termite_genome_annotation/Hisat/Coatitermes/2_GT155_S90_R1_001.fastq.gz.bam

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
# Make sure fasta header is simple and remove all space
# Dependencies: simplifyFastaHeaders.pl from AUGUSTUS
SimplifyID=function(fna=fna,
                    common_pattern=common_pattern, # common pattern in simplified sequence IDs
                    out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  fna_name=unlist(strsplit(fna,"/"));fna_name=fna_name[length(fna_name)]
  
  cmd=paste("simplifyFastaHeaders.pl", # AUGUSTUS
            fna,
            common_pattern,
            paste(out_dir,"/",fna_name,"_SimpleIDs",sep=""), # FASTA with simplified IDs
            paste(out_dir,"/",fna_name,"_IDconvert.tsv",sep="")) # old ID to new ID (tabular)
  print(cmd);system(cmd,wait=TRUE)
  return(paste(out_dir,"/",fna_name,"_SimpleIDs",sep=""))
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

# RepeatMasker
# Dependencies: Singularity, RepeatMasker
repeatmasker=function(fna=fna,# Fasta file of genome.
                    out_dir=out_dir,
                    RepeatLib.fa=RepeatLib.fa,
                    Threads=Threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  
  pwd_begin=getwd()
  setwd(out_dir)
  Threads=as.character(Threads)
  
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

# Compute gene models (gff3) by miniprot
# Get gene models (gff3) to protein alignments (gff3)
# Dependencies: miniprot
miniprot=function(fna=fna,
                  faa=faa, # comma-list
                  threads=threads,
                  out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd_begin=getwd();setwd(out_dir)
  
  # Compute gene models by spliced alignments
  if (!file.exists("genome.mpi")){
    cmd=paste("miniprot",
              "-t",threads,
              "-d","genome.mpi",
              fna,
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  prot=unlist(strsplit(faa,","))
  for (i in prot){
    system(paste("cat",i,">","proteins.faa",sep=" "),wait=TRUE)
    prefix=unlist(strsplit(i,"/"));prefix=prefix[length(prefix)]
    cmd=paste("miniprot",
              "-t",threads,
              "-P",prefix,
              "--gff",
              "genome.mpi",
              "proteins.faa","|",
              "sed '/##/d' ",
              ">>","gene_models.gff3", # mRNA, CDS, stop_codon features
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system("rm proteins.faa")
  }
  
  # Convert gene models to protein alignments
  system("echo '##gff-version 3' > protein_align.gff3")
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3==\"CDS\") $3=\"match\"; if ($3==\"match\") $8=\".\"; if ($3==\"match\") print $0}'",
            "gene_models.gff3",">>","protein_align.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="sed -i 's/Parent=/ID=/' protein_align.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="awk -F '\t' -v OFS='\t' '{if ($2==\"miniprot\")  $2=\"PROTEIN\";print$0}' protein_align.gff3 > protein.gff"
  print(cmd);system(cmd,wait=TRUE)
  system("rm protein_align.gff3")
  system("mv protein.gff protein_align.gff3")
  
  #system("rm proteins.faa",wait=TRUE)
  system("rm genome.mpi",wait=TRUE)
  system("rm miniprot.gff3",wait=TRUE)
  
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
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)

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
            "-@",Threads,
            "-o",paste(out_prefix,".bam",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("samtools","index",
            paste(out_prefix,".bam",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(0)
}

StringTie = function(input_bam=input_bam, # Input BAM.
                     output_prefix=out_prefix,
                     threads=threads){
  threads = as.character(threads)
  
  cmd = paste("stringtie",
              "-p", threads,
              "-o", paste(output_prefix,".gtf",sep=""),
              input_bam,
              sep=" ")
  print(cmd);system(cmd,wait = TRUE)
  
  cmd=paste("cufflinks_gtf_to_alignment_gff3.pl", # evidencemodeler
            paste(output_prefix,".gtf",sep=""),">",
            paste(output_prefix,"_tmp.gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($2==\"Cufflinks\")  $2=\"TRANSCRIPT\";print$0}'",
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
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
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
  
  cmd=paste("agat_convert_sp_gxf2gxf.pl",
            "--gff TransDecoder.gff3",
            "-o sorted.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("rm blastp.outfmt6")
  system("rm hmmsearch.out")
  system("rm pfam.domtblout")
  system("rm pipeliner.*.cmds")
  system("rm temp.gff3")
  system("rm transcripts.fna.transdecoder.bed")
  system("rm transcripts.fna.transdecoder.cds")
  system("rm -r transcripts.fna.transdecoder_dir")
  system("rm -r transcripts.fna.transdecoder_dir.__checkpoints")
  system("rm transcripts.fna.transdecoder.gff3")
  system("rm transcripts.fna.transdecoder.pep")
  system("rm TransDecoder.agat.log")
  setwd(wd)
  return("TransDecoder.gff3")
}

# Prothints+GeneMark-EP
# Dependencies: GeneMark
gm_ep = function(genome=genome,
                 faa=faa, # reference proteins
                 threads=threads,
                 out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd_begin=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  cmd=paste("gmes_petap.pl","--verbose",
            "--EP",
            "--dbep",faa,
            "--soft_mask","auto",
            "--cores",threads,
            "--sequence",genome,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="agat_convert_sp_gxf2gxf.pl --gff genemark.gtf -o sorted_tmp.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="awk -F '\t' -v OFS='\t' '{if ($2==\"GeneMark.hmm3\")  $2=\"GeneMark_EP\";print$0}' sorted_tmp.gff3 > sorted.gff3"
  print(cmd)
  system(cmd,wait=TRUE)
  
  system("mv prothint/prothint_augustus.gff .")
  system("rm -r prothint")
  system("rm sorted_tmp.gff3")
  system("rm -r data")
  system("rm genemark.agat.log")
  system("rm gmes.log")
  system("rm -r info")
  system("rm -r output")
  system("rm -r run")
  system("rm rm run.cfg")
  
  setwd(wd_begin)
}

# Get hints (introns) from RNA-mapping BAM
# Extract and filter hints (introns) from RNA-mapping BAM 
# GeneMark-ET with hints
# Dependencies: GeneMark-ET, AUGUSTUS, scripts from BRAKER
gm_et = function(genome=genome, # soft masked
                 bam=bam,
                 PE=TRUE, # logical. TRUE for pair-end RNA reads
                 out_dir=out_dir,
                 threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd_begin=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  # Filter BAM
  # cmd=paste("samtools sort -n",
  #           "-@",threads,
  #           bam,
  #           # "|","samtools sort","-@",threads,
  #           ">","Aligned.out.ss.bam",
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  # if (PE){st="--paired --pairwiseAlignment"}else{st=""}
  # cmd=paste("filterBam","--uniq",st,
  #           "--in","Aligned.out.ss.bam",
  #           "--out Aligned.out.ssf.bam",
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  # cmd=paste("samtools sort",
  #           "-@",threads,
  #           "Aligned.out.ssf.bam",
  #           ">","Aligned.out.ss.bam",
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  
  # Extract and filter hints
  cmd=paste("bam2hints","--intronsonly", # AUGUSTUS
            paste("--in=",bam,sep=""),
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
  
  cmd=paste("agat_convert_sp_gxf2gxf.pl",
            "--gff genemark.gtf",
            "-o sorted_tmp.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="awk -F '\t' -v OFS='\t' '{if ($2==\"GeneMark.hmm3\")  $2=\"GeneMark_ET\";print$0}' sorted_tmp.gff3 > sorted.gff3"
  print(cmd)
  system(cmd,wait=TRUE)
  
  system("rm sorted_tmp.gff3")
  system("rm -r data")
  system("rm genemark.agat.log")
  system("rm gmes.log")
  system("rm -r info")
  system("rm -r output")
  system("rm -r run")
  system("rm run.cfg")
  
  setwd(wd_begin)
}

# Get redundant training set (AUGUSTUS) from gene models (gff3, CDS/exon features)
# Dependencies: stringr (R), scripts from AUGUSTUS, BRAKER 
augustus_trainset=function(gene_models.gff3=gene_models.gff3,
                           fna=fna,
                           out_dir=out_dir){
  library(stringr)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
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
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
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
  
  system("rm *.lst")
  system("rm prot.*")
  system("rm bonafide.f.gtf")
  
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
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
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
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="cat AUGUSTUS.gff3 | join_aug_pred.pl > AUGUSTUS_nr.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  # EvidenceModeler
  cmd="augustus_GFF3_to_EVM_GFF3.pl AUGUSTUS_nr.gff3 > AUGUSTUS_evm.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("agat_convert_sp_gxf2gxf.pl",
            "--gff AUGUSTUS_evm.gff3",
            "-o sorted.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("mv"," ","$AUGUSTUS_CONFIG_PATH/species/",species," ",".",sep=""))
  
  system("rm -r autoAug")
  system("rm AUGUSTUS_evm.agat.log")
  system("rm AUGUSTUS.gff3")
  
  setwd(wd)
}

# GALBA: gene prediction trained by protein-genome alignments
# Dependencies: GALBA
galba=function(genome=genome,
               proteins=proteins, # comma list
               species=species,
               threads=threads,
               out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  prot=unlist(strsplit(proteins,","))
  for (i in prot){
    system(paste("cat",i,">>","proteins.faa",sep=" "),wait=TRUE)
  }
  
  cmd=paste("galba.pl",
            paste("--species=",species,sep=""),
            paste("--genome=",genome,sep=""),
            "--prot_seq=proteins.faa",
            paste("--workingdir=",out_dir,sep=""),
            paste("--threads=",threads,sep=""),
            "--skipGetAnnoFromFasta",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("agat_convert_sp_gxf2gxf.pl",
            "--gff augustus.hints.gtf",
            "-o sorted_tmp.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="awk -F '\t' -v OFS='\t' '{if ($2==\"AUGUSTUS\")  $2=\"GALBA\";print$0}' sorted_tmp.gff3 > sorted.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("mv"," ","$AUGUSTUS_CONFIG_PATH/species/",species," ",
               out_dir,sep=""))
  
  system("rm proteins.faa")
  system("rm sorted_tmp.gff3")
  system("rm -r errors")
  system("rm genome_header.map")
  system("rm pygustus_hints.out")
  system("rm pygustus_hints.py")
  system("rm what-to-cite.txt")
  setwd(wd)
}

# BRAKER
# if bam and ref_protein provided, run braker twice and run TSEBRA
# else: run braker once
braker=function(genome=genome,
                bam="none", # comma-list
                ref_proteins="none",
                species=species,
                tsebra.conf="none", # TSEBRA/config/default.cfg 
                out_dir=out_dir,
                threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  # BRAKER with RNA-seq
  if (bam!="none"){
    system(paste("mkdir"," ",out_dir,"/braker_rna/",sep=""),wait=TRUE)
    cmd=paste("braker.pl",
              paste("--cores=",threads,sep=""),
              paste("--workingdir=",out_dir,"/braker_rna/",sep=""),
              paste("--species=",species,"_braker_rna",sep=""),
              paste("--genome=",genome,sep=""),
              #paste("--hints=",hints_genemark_et.gff,sep=""),
              #paste("--geneMarkGtf=",genemark_et.gtf),
              paste("--bam=",bam,sep=""),
              "--softmasking",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system(paste("mv"," ","$AUGUSTUS_CONFIG_PATH/species/",species,"_braker_rna"," ",
                 out_dir,"/braker_rna/",sep=""))
  }
  
  # BRAKER with OrthoDB
  if (ref_proteins!="none"){
    system(paste("mkdir"," ",out_dir,"/braker_prot/",sep=""),wait=TRUE)
    cmd=paste("braker.pl",
              paste("--cores=",threads,sep=""),
              paste("--workingdir=",out_dir,"/braker_prot/",sep=""),
              paste("--species=",species,"_braker_prot",sep=""),
              paste("--genome=",genome,sep=""),
              #paste("--hints=",hints_genemark_ep.gff,sep=""),
              #paste("--geneMarkGtf=",genemark_ep.gtf),
              paste("--prot_seq=",ref_proteins,sep=""),
              "--softmasking",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    system(paste("mv"," ","$AUGUSTUS_CONFIG_PATH/species/",species,"_braker_prot"," ",
                 out_dir,"/braker_prot/",sep=""))
  }
  
  # TSEBRA
  if (tsebra.conf!="none"){
    cmd=paste("tsebra.py",
              "-g",paste(out_dir,"/braker_rna/augustus.hints.gtf",",",
                         out_dir,"/braker_prot/augustus.hints.gtf",
                         sep=""),
              "-c",tsebra.conf,
              "-e",paste(out_dir,"/braker_rna/hintsfile.gff",",",
                         out_dir,"/braker_prot/hintsfile.gff",
                         sep=""),
              "-o","braker_final.gtf",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  if (tsebra.conf=="none"){
    if (bam!="none"){
      system(paste("cp"," ",out_dir,"/braker_rna/augustus.hints.gtf"," ","braker_final.gtf",sep=""))
    }
    if (ref_proteins!="none"){
      system(paste("cp"," ",out_dir,"/braker_prot/augustus.hints.gtf"," ","braker_final.gtf",sep=""))
    }
  }
  
  cmd=paste("agat_convert_sp_gxf2gxf.pl",
            "--gff braker_final.gtf",
            "-o sorted_tmp.gff",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system("rm braker_final.agat.log")
  
  cmd="awk -F '\t' -v OFS='\t' '{if ($2==\"AUGUSTUS\")  $2=\"BRAKER\";print$0}' sorted_tmp.gff > sorted.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  system("rm sorted_tmp.gff")
  system("rm ./braker_rna/augustus.hints.aa")
  system("rm ./braker_rna/augustus.hints.codingseq")
  system("rm ./braker_rna/bam_header.map")
  system("rm -r ./braker_rna/errors")
  system("rm -r ./braker_rna/GeneMark-ET")
  system("rm ./braker_rna/genemark_hintsfile.gff")
  system("rm ./braker_rna/genome_header.map")
  system("rm ./braker_prot/what-to-cite.txt")
  
  system("rm ./braker_prot/augustus.hints.aa")
  system("rm ./braker_prot/augustus.hints.codingseq")
  system("rm -r ./braker_prot/errors")
  system("rm ./braker_prot/genemark_hintsfile.gff")
  system("rm ./braker_prot/genome_header.map")
  system("rm ./braker_prot/what-to-cite.txt")
  system("rm ./braker_prot/augustus.hints_iter1.gtf")
  system("rm ./braker_prot/evidence.gff")
  system("rm -r ./braker_prot/GeneMark-EP")
  system("rm -r ./braker_prot/GeneMark-ES")
  system("rm ./braker_prot/genemark_evidence.gff")
  system("rm ./braker_prot/prothint.gff")
  setwd(wd)
}

# EvidenceModeler
# Dependencies: EvidenceModeler, agat
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
             # OTHER_PREDICTION TransDecoder  10
             # gff3 of AUGUSTUS, Glimmer and TransDecoder should be merged 
             # into a single gff3 whose 2nd column is AUGUSTUS/Glimmer/TransDecoder
             genome=genome,
             out_dir=out_dir,
             threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  system("mkdir temp_dir",wait=TRUE)
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
            "1>partition_EVM_inputs.stdout 2>partition_EVM_inputs.stderr",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="tail -n50 partition_EVM_inputs.stdout";print(cmd);system(cmd,wait=TRUE)
  cmd="tail -n50 partition_EVM_inputs.stderr";print(cmd);system(cmd,wait=TRUE)
  system("rm partition_EVM_inputs.stdout");system("rm partition_EVM_inputs.stderr")
  
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
  stopCluster(clus)
  
  cmd="recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("convert_EVM_outputs_to_GFF3.pl",
            "--partitions partitions_list.out --output evm.out",
            "--genome",genome,
            "1> convert_EVM_outputs_to_GFF3.stdout 2> convert_EVM_outputs_to_GFF3.stderr",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="tail -n50 convert_EVM_outputs_to_GFF3.stdout";print(cmd);system(cmd,wait=TRUE)
  cmd="tail -n50 convert_EVM_outputs_to_GFF3.stderr";print(cmd);system(cmd,wait=TRUE)
  system("rm convert_EVM_outputs_to_GFF3.stdout");system("rm convert_EVM_outputs_to_GFF3.stderr")
  
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
  system("rm x*")
  system("rm merge_gff3.sh")
  
  setwd(out_dir)
  
  cmd="evm_evidence.py temp_dir/ | awk -F '\t' -v OFS='\t' '{print$1,$2,$3,$4,$5,$7}' | sort > evidence_summary.tsv"
  print(cmd);system(cmd,wait=TRUE)
  system("rm -r temp_dir",wait=TRUE)
  
  cmd=paste("agat_convert_sp_gxf2gxf.pl",
            "--gff EvidenceModeler.gff3",
            "-o sorted.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system("rm EvidenceModeler.agat.log")
  cmd="grep 'gene' sorted.gff3 | awk -F '\t' -v OFS='\t' '{print$1,$2,$3,$4,$5,$9}' | sed 's/ID=//' | sed 's/;Name=.*$//' | sort > gene_coordinate.tsv"
  print(cmd);system(cmd,wait=TRUE)
  
  a=read.table("evidence_summary.tsv",header=FALSE,sep="\t")
  b=read.table("gene_coordinate.tsv",header=FALSE,sep="\t")
  c=merge(a,b,by=c("V1","V2","V3","V4","V5"),all=TRUE)
  colnames(c)=c("chr","source","feature","start","end","software","geneID")
  c[,"evidence"]=sapply(c[,"software"],
                        function(i){
                          if (grepl("TRANSCRIPT",i)){
                            return("Transcript")
                          }else if(grepl("PROTEIN",i)){
                            return("Protein")
                          }else{
                            return("Hypothetical")
                          }
                        })
  write.table(c,"gene_evidence.tsv",sep="\t",row.names=FALSE,quote=FALSE)
  system("rm evidence_summary.tsv")
  system("rm gene_coordinate.tsv")
  
  system("mkdir add_evidence")
  setwd("add_evidence")
  cmd=paste("split",
            "-d",
            "-n",paste("l/",as.character(as.numeric(threads)-1),sep=""),
            "../gene_evidence.tsv",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  clus=makeCluster(as.numeric(threads))
  parSapply(clus,
            1:(as.numeric(threads)-1),
            function(i){
              scr=system("ls x*",wait=TRUE,intern=TRUE)[i]
              cmd=paste("awk -F '\t' -v OFS='\t' '{print $7}'",
                        scr,">",paste("gene_list.",as.character(i),sep=""),
                        sep=" ")
              system(cmd,wait=TRUE)
              cmd=paste("agat_sp_filter_feature_from_keep_list.pl",
                        "--gff ../sorted.gff3",
                        "--keep_list",paste("gene_list.",as.character(i),sep=""),
                        "--output",paste(as.character(i),".gff3",sep=""),
                        ">",paste("agat.out.",as.character(i),sep=""),
                        sep=" ")
              system(cmd,wait=TRUE)
              map=system(paste("awk -F '\t' -v OFS='\t' '{print $7,$8}'",
                               scr,sep=" "),
                         intern=TRUE)
              sapply(map,
                     function(j){
                       j=unlist(strsplit(j,"\t"))
                       gene=j[1];evidence=j[2]
                       cmd=paste("sed -i",
                                 paste("'s/","ID=",gene,";/ID=",gene,";Evidence=",evidence,";/'",sep=""),
                                 paste(as.character(i),".gff3",sep=""),
                                 sep=" ")
                       system(cmd,wait=TRUE)
                     })
              
              cmd=paste("cat",paste("./add_evidence/",as.character(i),".gff3",sep=""),
                        ">>","final.gff3",sep=" ")
              cat(paste(cmd,"\n",sep=""),file="../cmd.sh",append=TRUE)
            }
  )
  stopCluster(clus)
  setwd(out_dir)
  system("bash cmd.sh",wait=TRUE)
  system("mv sorted.gff3 ./add_evidence/")
  system("rm -r ./add_evidence")
  
  cmd=paste("agat_convert_sp_gxf2gxf.pl",
            "--gff final.gff3",
            "-o sorted.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system("rm final.agat.log")
  
  setwd(wd)
}

# Filter evm out, remove hypothetical genes supported by only 1 predictor
filterEvm=function(gene_evidence.tsv=gene_evidence.tsv, # geneID, software, evidence
                   genome.fna=genome.fna,
                   evm.gff3=evm.gff3,
                   out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  geneEvidence=read.table(gene_evidence.tsv,sep="\t",header=TRUE,quote="")
  geneID=sapply(1:nrow(geneEvidence),
                function(i){
                  d=geneEvidence[i,]
                  if (d[,"evidence"]!="Hypothetical"){
                    return(d[,"geneID"])
                  }else{
                    support=unlist(strsplit(d[,"software"],","))
                    if (length(support)>1){
                      return(d[,"geneID"])
                    }else{
                      return(NA)
                    }
                  }
                })
  geneID=geneID[!is.na(geneID)]
  writeLines(geneID,"./geneID.lst")
  
  cmd=paste("agat_sp_filter_feature_from_keep_list.pl",
            "--gff",evm.gff3,
            "--keep_list ./geneID.lst",
            "--out ./filtered.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}


# Update gff3 by PASA
# Dependencies: pasa,sqlite3,,agat,scripts from PASAPipeline
pasa=function(genome=genome,
              transcripts=transcripts,
              original.gff3=original.gff3, # only protein-coding gene with 'evidence' in 9th
              pasa.alignAssembly.conf=pasa.alignAssembly.conf, # Path to PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt
              pasa.annotationCompare.conf=pasa.annotationCompare.conf, # Path to PASApipeline/pasa_conf/pasa.annotationCompare.Template.txt
              threads=threads,
              out_dir=out_dir){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  system("mkdir ./tmp");setwd("tmp")
  
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
              "-c",paste(getwd(),"/alignAssembly.config",sep=""),
              "-C -R",
              "-g",genome,
              "-t",transcripts,
              "--ALIGNERS minimap2",
              "--CPU",threads,
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("Load_Current_Gene_Annotations.dbi",
            "-c",paste(getwd(),"/alignAssembly.config",sep=""),
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
            "-c",paste(getwd(),"/annotationCompare.config",sep=""),
            "-A",
            "-g",genome,
            "-t",transcripts,
            "--CPU",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  re=system("ls pasa.sqlite.db.gene_structures_post_PASA_updates.[0-9]*.gff3",wait=TRUE,intern=TRUE)
  cmd=paste("cp",re,"../PASA.gff3",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",genome,sep=" "),wait=TRUE)
  system(paste("rm",transcripts,sep=" "),wait=TRUE)
  system(paste("rm",original.gff3,sep=" "),wait=TRUE)
  
  setwd("../")
  cmd=paste("agat_convert_sp_gxf2gxf.pl",
            "--gff","PASA.gff3",
            "-o","sorted.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("rm -r tmp")
  system("rm PASA.agat.log")
  
  setwd(wd)
}

filterPasa=function(genome.fna=genome.fna,
                    pasa.gff3=pasa.gff3,
                    out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  system(paste("cp",genome.fna,"./genome.fna",sep=" "))
  cmd=paste("gffread","-O",pasa.gff3,
            "-g ./genome.fna -y ./pep.faa",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="seqkit grep -s -r -p \"[\\*\\.][A-Z]\" pep.faa | grep '>' | sed 's/>//' > ./inFrameStop.lst"
  print(cmd);system(cmd,wait=TRUE)
  system("rm ./genome.fna ./genome.fna.fai")
  
  cmd=paste("agat_sp_filter_feature_from_kill_list.pl",
            "--gff",pasa.gff3,
            "--kill_list ./inFrameStop.lst",
            "--out ./filtered.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}

# Generate a comprehensive genomic database including:
# genes.gff3, transcripts_rep.fna, transcripts_iso.fna, cds_rep.fna, cds_iso.fna, proteins_rep.faa, proteins_iso.faa
# Dependencies: maker, gffread, seqkit
pasa_more=function(species=species,
                   genome=genome,
                   PASA.gff3=PASA.gff3,
                   gene_evidence.tsv=gene_evidence.tsv, # from function evm
                   out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  # Change gene ID
  cmd=paste("cp",PASA.gff3,paste("./",species,"_genes.gff3",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  gff3=paste("./",species,"_genes.gff3",sep="")
  cmd=paste("maker_map_ids --iterate 0",
            "--prefix",species,gff3,"> MAKER.name.map",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("map_gff_ids","MAKER.name.map",gff3,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Evidence for genes
  cmd=paste("cp",gene_evidence.tsv,paste("./",species,"_gene_evidence.tsv",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  a=read.table(paste("./",species,"_gene_evidence.tsv",sep=""),sep="\t",header=TRUE,quote="")
  a=a[,c("geneID","evidence","software")]
  b=read.table("MAKER.name.map",sep="\t",header=FALSE,quote="")
  b=b[!grepl("-R[0-9]",b$V2),];b=b[!duplicated(b$V2),]
  c=merge(a,b,by.x="geneID",by.y="V1",all=TRUE)
  colnames(c)=c("oldID","evidence","software","geneID")
  c=c[!is.na(c$geneID),]
  for (i in 1:nrow(c)){
    if (is.na(c[i,"evidence"]) & is.na(c[i,"software"])){
      c[i,"evidence"]="Transcript"
      c[i,"software"]="PASA"
      print(c[i,])
    }
  }
  write.table(c,paste("./",species,"_gene_evidence.tsv",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
  
  # Extract transcripts, cds, proteins
  system(paste("cp",genome,paste(species,"_maskedGenome.fna",sep=""),sep=" "))
  genome=paste(species,"_maskedGenome.fna",sep="")
  cmd=paste("gffread","-O",
            gff3,"-S",
            "-g",genome,
            "-w",paste(species,"_transcripts.fna",sep=""),
            "-x",paste(species,"_cds.fna",sep=""),
            "-y",paste(species,"_proteins.faa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Separate alternative splicing
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
  
  # cmd=paste("sed -i 's/-R0//' ",species,"_transcripts_rep.fna",sep="")
  # print(cmd);system(cmd,wait=TRUE)
  # cmd=paste("sed -i 's/-R0//' ",species,"_cds_rep.fna",sep="")
  # print(cmd);system(cmd,wait=TRUE)
  # cmd=paste("sed -i 's/-R0//' ",species,"_proteins_rep.faa",sep="")
  # print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",paste(species,"_transcripts.fna",sep=""),sep=" "),wait=TRUE)
  system(paste("rm",paste(species,"_cds.fna",sep=""),sep=" "),wait=TRUE)
  system(paste("rm",paste(species,"_proteins.faa",sep=""),sep=" "),wait=TRUE)
  system("rm *.fai")
  setwd(wd)
}

# PseudoPipe for pseudogene identification
# Dependencies: seqkit,PseudoPipe (modified)
PseudoPipe=function(genome=genome, # soft masked
                    protein.fa=protein.fa, # space list
                    gff=gff,
                    out_dir=out_dir,
                    threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  system("mkdir ./input")
  
  system("mkdir ./input/dna")
  cmd=paste("cp",genome,"./input/dna/masked.fa",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("seqkit","seq","-u",
            "-j",threads,
            "./input/dna/masked.fa > ./input/dna/unmasked.fa",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("seqkit","split2","./input/dna/unmasked.fa","-s 1","-j",threads," 2> lst",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="grep 'write 1 sequences to file:' lst"
  split.fa=system(cmd,intern=TRUE)
  split.fa=sub("^.*write 1 sequences to file: ","",split.fa)
  
  split=sapply(1:length(split.fa),
               function(i){
                 seqID=system(paste("head -n1",split.fa[i],sep=" "),intern=TRUE)
                 seqID=sub(">","",seqID)
                 seqID=sub(" .*$","",seqID)
                 
                 cmd=paste("mv",split.fa[i],
                           paste("./input/dna/unmasked.",seqID,".fa",sep=""),
                           sep=" ")
                 system(cmd,wait=TRUE)
                 return(paste("./input/dna/unmasked.",seqID,".fa",sep=""))
               })
  system("rm -r ./input/dna/unmasked.fa.split/")
  system("rm ./input/dna/unmasked.fa")
  system("rm lst")
  
  system("mkdir ./input/pep")
  cmd=paste("cat",protein.fa,">","./input/pep/protein.fa",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("mkdir ./input/mysql")
  cmd=paste("grep","'exon'",gff,"|",
            "awk","-F '\t' -v OFS='\t' '{print $1,$2,$4,$5}'",
            ">","./input/mysql/gff",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  sapply(1:length(split),
         function(i){
           seqID=system(paste("head -n1",split[i],sep=" "),intern=TRUE)
           seqID=sub(">","",seqID)
           seqID=sub(" .*$","",seqID)
           cmd=paste("grep"," ","'",seqID,"'"," ","./input/mysql/gff",
                     " > ","./input/mysql/gff.",seqID,".tsv",
                     sep="")
           system(cmd,wait=TRUE)
         })
  #########################
  # python2
  system("module load sango-legacy-modules")
  system("module load python/2.7.3")
  system("module load fasta/35.4.12")
  #########################
  cmd=paste("pseudopipe.sh",
            getwd(),
            paste(getwd(),"/input/dna/masked.fa",sep=""),
            paste(getwd(),"/input/dna/unmasked.%s.fa",sep=""),
            paste(getwd(),"/input/pep/protein.fa",sep=""),
            paste(getwd(),"/input/mysql/gff.%s.tsv",sep=""),
            "0",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # system("rm -r input/")
  # system("rm -r blast/")
  # system("rm -r dna/")
  # system("rm -r pep/")
  system("mv ./pgenes/*.gz ./")
  system("mv ./pgenes/*.txt ./")
  # system("rm -r pgenes/")
  setwd(wd)
}

# PseudoPipe output to gff3
PseudoPipe2gff3=function(Pseudo.out.txt=Pseudo.out.txt,
                         species=species,
                         out.gff3=out.gff3){
  origin=readLines(Pseudo.out.txt)
  origin=origin[2:length(origin)]
  
  in.type=c("DUP","PSSD","FRAG")
  out.type=c("duplicated","processed","fragment")
  names(out.type)=in.type
  
  gff3=sapply(1:length(origin),
              function(i){
                txt=unlist(strsplit(origin[i],"\t"))
                g=paste(txt[1],"\t","PseudoPipe","\t","pseudogene","\t",
                        txt[2],"\t",txt[3],"\t",".","\t",txt[4],"\t",
                        ".","\t",
                        paste("ID=pgene_",species,as.character(i),";",
                              "Parent=",txt[5],";",
                              "Type=",unname(out.type[txt[14]]),
                              sep=""),
                        sep="")
                return(g)
              })
  writeLines(gff3,out.gff3)
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
  
  system("mkdir trinity")
  cmd=paste("Trinity",
            "--seqType","fq",
            "--left read1.fq.gz",
            "--right read2.fq.gz",
            "--CPU",threads,
            "--output",paste(out_dir,"/trinity",sep=""),
            "--max_memory",max_memory,
            "--full_cleanup",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("rm read1.fq.gz read2.fq.gz")
  system("rm -r ./trinity")
  setwd(wd)
}

















# # Sort gff3 
# # Dependencies: AGAT
# sort_gff3=function(in.gff3=in.gff3,
#                    out.gff3=out.gff3){
#   cmd=paste("agat_convert_sp_gxf2gxf.pl",
#             "--gff",in.gff3,
#             "-o",out.gff3,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
# }
# 
# # MAKER
# # Dependencies: MAKER, SNAP, RepeatMasker, blast, AUGUSTUS, exonerate, blast+, GAAS, stringr (R)
# #               scripts from AUGUSTUS
# maker=function(fna=fna, # masked
#                transcript.file=transcript.file, # gff3/gtf with transcript/mRNA and exon/CDS features e.g. gtf from stringtie
#                protein.file=protein.file, # gff3/gtf with transcript/mRNA and exon/CDS features
#                threads=threads,
#                augustus_species=augustus_species,
#                out_dir=out_dir){
#   threads=as.character(threads)
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   wd=getwd();setwd(out_dir)
#   
#   if (!file.exists("maker_1/MAKER.gff3")){
#   if (file.exists("maker_1")){system("rm -r maker_1",wait=TRUE)}
#   cmd=paste("agat_sp_alignment_output_style.pl",
#             "-g",transcript.file,
#             "-o","transcript_align.gff3")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("agat_sp_alignment_output_style.pl",
#             "-g",protein.file,
#             "-o","protein_align.gff3")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   transcript_align.gff3="transcript_align3.gff"
#   protein_align.gff3="protein_align3.gff"
#   
#   system("mkdir maker_1",wait=TRUE);setwd("maker_1")
#   cmd="maker -CTL"
#   print(cmd);system(cmd,wait=TRUE)
#   a=readLines("maker_opts.ctl")
#   a[2]=paste("genome=",fna,sep="")
#   a[18]=paste("est_gff=",out_dir,"/",transcript_align.gff3,sep="")
#   a[23]=paste("protein_gff=",out_dir,"/",protein_align.gff3,sep="")
#   a[41]="est2genome=1"
#   a[42]="protein2genome=1"
#   writeLines(a,"maker_opts.ctl")
#   cmd=paste("maker",
#             "-c",threads,
#             "-RM_off","-quiet",
#             "-base maker",
#             "maker_opts.ctl maker_bopts.ctl maker_exe.ctl",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="gff3_merge -s -n -g -d maker.maker.output/maker_master_datastore_index.log > MAKER.gff3"
#   print(cmd);system(cmd,wait=TRUE)
#   setwd("..")
#   }else{
#   transcript_align.gff3="transcript_align3.gff"
#   protein_align.gff3="protein_align3.gff"
#   }
#   
#   if (!file.exists("snap_1/MAKER.hmm")){
#   if (file.exists("snap_1")){system("rm -r snap_1",wait=TRUE)}
#   system("mkdir snap_1",wait=TRUE);setwd("snap_1")
#   cmd="maker2zff -x 0.25 -l 50 -d ../maker_1/maker.maker.output/maker_master_datastore_index.log"
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1"
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="fathom genome.ann genome.dna -validate > validate.log 2>&1"
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1"
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1"
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="forge export.ann export.dna > forge.log 2>&1"
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="hmm-assembler.pl genome . > MAKER.hmm"
#   print(cmd);system(cmd,wait=TRUE)
#   setwd("..")
#   }
#   
#   config=system("echo $AUGUSTUS_CONFIG_PATH",wait=TRUE,intern=TRUE)
#   if (!file.exists(paste(config,"/species/",augustus_species,sep=""))){
#   if (file.exists("augustus_1")){system("rm -r augustus_1",wait=TRUE)}
#   system("mkdir augustus_1",wait=TRUE);setwd("augustus_1")
#   library(stringr)
#   cmd=paste("gth2gtf.pl", # AUGUSTUS 
#             "../maker_1/MAKER.gff3",
#             "training_AUGUSTUS.gtf",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="computeFlankingRegion.pl training_AUGUSTUS.gtf" # AUGUSTUS
#   print(cmd)
#   flanking_DNA=system(cmd,wait=TRUE,intern=TRUE)[4]
#   flanking_DNA=str_extract(flanking_DNA,": [0-9]*")
#   flanking_DNA=sub(": ","",flanking_DNA)
#   cmd=paste("gff2gbSmallDNA.pl training_AUGUSTUS.gtf", # AUGUSTUS
#             fna,flanking_DNA,"training_AUGUSTUS.gb",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   GenBank=readLines("training_AUGUSTUS.gb")
#   gb=GenBank[grepl("/gene=",GenBank)]
#   gb=sub("                     /gene=\"","",gb)
#   gb=sub("\"","",gb)
#   gb=unique(gb)
#   write.table(gb,"traingenes.lst",sep="\t",col.names=FALSE,row.names=FALSE,quote=TRUE)
#   LOCUS=GenBank[grepl("LOCUS",GenBank)]
#   LOCUS=sub("LOCUS       ","",LOCUS)
#   LOCUS=sub("   [0-9]* bp  DNA","",LOCUS)
#   write.table(data.frame(gb,LOCUS),"loci.lst",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
#   cmd=paste("grep",
#             "-f","traingenes.lst",
#             "-F","training_AUGUSTUS.gtf",
#             ">",
#             "bonafide.f.gtf",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("gtf2aa.pl", # AUGUSTUS
#             fna,
#             "bonafide.f.gtf",
#             "prot.aa",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("simplifyFastaHeaders.pl", # AUGUSTUS
#             "prot.aa",
#             "Gene",
#             "prot.aa_SimpleIDs", # FASTA with simplified IDs
#             "prot.aa_IDconvert.tsv", # new ID to old ID (tabular)
#             sep=" ") 
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("aa2nonred.pl",
#             paste("--cores=",threads,sep=""),
#             "prot.aa_SimpleIDs","prot.nr.aa",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="grep '>' prot.nr.aa | perl -pe 's/>//' > nonred.lst"
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="grep -f nonred.lst prot.aa_IDconvert.tsv | cut -f2 | perl -pe 's/>//' > original_geneID.lst"
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="grep -f original_geneID.lst loci.lst | cut -f2 > nonred.loci.lst"
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("filterGenesIn.pl",
#             "nonred.loci.lst",
#             "training_AUGUSTUS.gb",
#             ">",
#             "training_AUGUSTUS_nr.gb",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("autoAug.pl",
#             paste("--cpus=",threads,sep=""),
#             paste("--workingdir=",out_dir,"/augustus_1",sep=""),
#             paste("--species=",augustus_species,sep=""),
#             paste("--genome=",fna,sep=""),
#             paste("--trainingset=","training_AUGUSTUS_nr.gb",sep=""),
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   setwd("..")
#   }
#   
#   if (!file.exists("maker_2/MAKER.gff3")){
#   if (file.exists("maker_2")){system("rm -r maker_2",wait=TRUE)}
#   system("mkdir maker_2",wait=TRUE);setwd("maker_2")
#   cmd="maker -CTL"
#   print(cmd);system(cmd,wait=TRUE)
#   a=readLines("maker_opts.ctl")
#   a[2]=paste("genome=",fna,sep="")
#   a[18]=paste("est_gff=",out_dir,"/",transcript_align.gff3,sep="")
#   a[23]=paste("protein_gff=",out_dir,"/",protein_align.gff3,sep="")
#   a[34]=paste("snaphmm=",out_dir,"/snap_1/MAKER.hmm",sep="")
#   a[36]=paste("augustus_species=",augustus_species,sep="")
#   writeLines(a,"maker_opts.ctl")
#   cmd=paste("maker",
#             "-c",threads,
#             "-RM_off","-quiet",
#             "-base maker",
#             "maker_opts.ctl maker_bopts.ctl maker_exe.ctl",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="gff3_merge -s -n -g -d maker.maker.output/maker_master_datastore_index.log > MAKER.gff3"
#   print(cmd);system(cmd,wait=TRUE)
#   setwd("..")
#   }
#   
#   if (!file.exists("MAKER.gff3")){
#   cmd=paste("maker_map_ids",
#             "--prefix",augustus_species,
#             "maker_2/MAKER.gff3 > MAKER.name.map",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="map_gff_ids MAKER.name.map maker_2/MAKER.gff3"
#   print(cmd);system(cmd,wait=TRUE)
#   system("mv maker_2/MAKER.gff3 ./MAKER.gff3")
#   }
#   
#   setwd(wd)
# }
# 
# # Star: Map short RNA reads to reference genome. 
# # SAMtools: Compress SAM to BAM, sort BAM, index BAM.
# # ulimit -n
# Star = function(fq1=fq1,fq2=fq2, # Input fq files. Make fq2="None" if single-end.
#                 # Comma-separated list.
#                 fna=fna, # genome (soft masked for training AUGUSTUS)
#                 out_dir=out_dir,
#                 out_basename=out_basename,
#                 threads=threads){
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   wd_begin=getwd();setwd(out_dir)
#   threads=as.character(threads)
#   
#   system("mkdir index",wait=TRUE)
#   cmd=paste("STAR",
#             "--runThreadN",threads,
#             "--runMode genomeGenerate",
#             "--genomeDir index",
#             "--genomeFastaFiles",fna,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("awk '/^>/{if (l!=\"\") print(l);l=0;next}{l+=length($0)}END{print l}'",
#             fna,sep=" ")
#   print(cmd)
#   size=system(cmd,wait=TRUE,intern=TRUE)
#   size=sum(as.numeric(size))
#   genomeSAindexNbases=14
#   if (log2(size)/2-1<14){genomeSAindexNbases=floor(log2(size)/2-1)}
#   
#   cmd=paste("STAR",
#             "--runMode alignReads",
#             "--runThreadN",threads,
#             "--genomeDir index",
#             "--readFilesIn",fq1,fq2,
#             "--outFileNamePrefix",out_basename,
#             "--outSAMtype BAM SortedByCoordinate",
#             "--genomeSAindexNbases",as.character(genomeSAindexNbases),
#             "--readFilesCommand gunzip -c",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("samtools","index",
#             paste(out_basename,"Aligned.sortedByCoord.out.bam",sep=""),
#             "-@",threads,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   return(paste(out_dir,"/",out_basename,"Aligned.sortedByCoord.out.bam",sep=""))
# }
# 
# # GenomeThreader: protein-genome spliced alignments
# # Dependencies: GenomeThreader, scripts from EvidenceModeler
# gth=function(fna=fna, # masked genome
#              faa=faa, # comma-list
#              out_dir=out_dir){
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   wd_begin=getwd();setwd(out_dir)
#   
#   prot=unlist(strsplit(faa,","))
#   for (i in prot){
#     system(paste("cat",i,">","proteins.faa",sep=" "),wait=TRUE)
#   }
#   
#   cmd=paste("cp",fna,".",sep=" ");system(cmd,wait=TRUE)
#   fna=unlist(strsplit(fna,"/"));fna_file=fna[length(fna)]
#   
#   cmd=paste("gth",
#             "-genomic",fna_file,
#             "-protein","proteins.faa",
#             "-gff3out -intermediate",
#             "-o","protein_alignment.gff3",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   system(paste("rm",fna_file,sep=" "),wait=TRUE)
#   system(paste("rm","proteins.faa",sep=" "),wait=TRUE)
#   system(paste("rm"," ",fna_file,"*",sep=""),wait=TRUE)
#   system(paste("rm"," ","proteins.faa","*",sep=""),wait=TRUE)
#   # EvidenceModeler
#   cmd="genomeThreader_to_evm_gff3.pl protein_alignment.gff3 >protein_alignment_evm.gff3"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   setwd(wd_begin)
#   return(out_dir)
# }
# 
# # Spaln for protein-genome alignments
# # Dependencies: Spaln, gffread, scripts from AUGUSTUS
# spaln=function(genome=genome,
#                faa=faa, # comma-list
#                threads=threads,
#                out_dir=out_dir){
#   threads=as.character(threads)
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   wd_begin=getwd();setwd(out_dir)
#   
#   cmd=paste("cp",genome,out_dir,sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   genome=unlist(strsplit(genome,"/"));genome=genome[length(genome)]
#   
#   prot=unlist(strsplit(faa,","))
#   for (i in prot){
#     system(paste("cat",i,">","proteins.faa",sep=" "),wait=TRUE)
#   }
#   
#   # Format genome
#   cmd=paste("spaln",
#             "-W -KP",
#             "-t",threads,
#             genome,sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   # Gene models
#   cmd=paste("spaln",
#             "-Q 7",
#             "-O 0",
#             "-t",threads,
#             "-d",genome,"proteins.faa",
#             ">","Spaln_gene.gff3",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   # Format to protein-genome alignments
#   
#   system(paste("rm",genome,sep=" "),wait=TRUE)
#   system("rm proteins.faa",wait=TRUE)
#   setwd(wd_begin)
# }
# 
# # GeMoMa: Predict genes from transcript/protein-genome alignments
# # Dependencies: GeMoMa, scripts from evidencemodeler
# gemoma=function(genome=genome,
#                 ref_gff=ref_gff,# comma-list
#                 ref_genome=ref_genome, # comma-list
#                 RNA_lib="FR_UNSTRANDED", # FR_UNSTRANDED/FR_FIRST_STRAND/FR_SECOND_STRAND
#                 bam=bam,
#                 threads=threads,
#                 out_dir=out_dir){
#   #######################################
#   #Absolute path to gemoma
#   #######################################
#   path="/home/c/c-liu/miniconda3/pkgs/gemoma-1.9-hdfd78af_0/share/gemoma-1.9-0/GeMoMa-1.9.jar"
#   system("conda activate gemoma",wait=TRUE)
#   # setwd("/home/c/c-liu/miniconda3/pkgs/gemoma-1.9-hdfd78af_0/share/gemoma-1.9-0/")
#   # system("cd /home/c/c-liu/miniconda3/pkgs/gemoma-1.9-hdfd78af_0/share/gemoma-1.9-0/")
#   #######################################
#   
#   threads=as.character(threads)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   
#   ref_genome=unlist(strsplit(ref_genome,","))
#   ref_gff=unlist(strsplit(ref_gff,","))
#   f=function(i){
#     s1=paste("s=own",
#              paste("a=",ref_gff[i],sep=""),
#              paste("g=",ref_genome[i],sep=""),
#              sep=" ")
#     return(s1)
#   }
#   s=sapply(1:length(ref_genome),f)
#   s=paste(s,collapse=" ")
#   
#   cmd=paste("java -jar",
#             path,"CLI GeMoMaPipeline",
#             paste("threads=",threads,sep=""),
#             paste("t=",genome,sep=""),
#             s,
#             "tblastn=true",
#             paste("outdir=",out_dir,sep=""),
#             "r=MAPPED",
#             paste("ERE.s=",RNA_lib,sep=""),
#             paste("ERE.m=",bam,sep=""),
#             "ERE.c=true AnnotationFinalizer.r=NO",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("GeMoMa_gff_to_gff3.pl"," ",
#             out_dir,"/final_annotation.gff"," ",
#             ">"," ",
#             out_dir,"/GeMoMa_evm.gff")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("agat_convert_sp_gxf2gxf.pl"," ",
#             "--gff"," ",out_dir,"/GeMoMa_evm.gff"," ",
#             "-o"," ",out_dir,"/sorted.gff",
#             sep="")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   wd=getwd();setwd(out_dir)
#   system("rm GeMoMa_evm.agat.log") 
#   system("rm predicted_proteins.fasta")
#   system("rm protocol_GeMoMaPipeline.txt") 
#   system("rm reference_gene_table.tabular")
#   setwd(wd)
# }
# 
# # GeneMarkerHMM for gene prediction
# # Dependencies: GeneMark
# gmhmm=function(genome=genome,
#                model_file=model_file, # gmhmm.mod from GeneMark-ET/EP
#                threads=threads,
#                out_dir=out_dir){
#   threads=as.character(threads)
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   wd_begin=getwd();setwd(out_dir)
#   
#   cmd=paste("gmes_petap.pl",
#             "--predict_with",model_file,
#             "--format","GFF3",
#             "--soft_mask auto",
#             "--cores",threads,
#             "--sequence",genome,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("agat_convert_sp_gxf2gxf.pl",
#             "--gff genemark.gff3",
#             "-o sorted.gff3",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   system("rm -r data")
#   system("rm genemark.agat.log")
#   system("rm rm gmes.log")
#   system("rm -r info")
#   system("rm -r output")
#   system("rm -r run")
#   system("rm run.cfg")
#   
#   setwd(wd_begin)
# }
# 
# # Extract exon coordinates from transcript alignment (gff3) and randomly select 1000 genes
# # Run GlimmerHMM
# # Dependencies: glimmerhmm, scripts from EvidenceModeler
# glimmerhmm=function(transcript_align.gff3=transcript_align.gff3, # get exon/match features
#                     # empty line separate genes
#                     genome=genome,
#                     out_dir=out_dir){
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   wd_begin=getwd();setwd(out_dir)
#   
#   # Get exon coordinates
#   cmd=paste("awk -F '\t' -v OFS='\t' '{if ($0==\"\" || $3==\"exon\" || $3==\"match\") print $0}'",
#             transcript_align.gff3,"> tmp.gff3")
#   print(cmd);system(cmd,wait=TRUE)
#   transcript_align=readLines("tmp.gff3")
#   training=sapply(transcript_align,
#                   function(st){
#                     if (st==""){
#                       return("")
#                     }else{
#                       st=unlist(strsplit(st,"\t"))
#                       seq=st[1];strand=st[7]
#                       if (strand=="+"){start=st[4];end=st[5]}
#                       if (strand=="-"){start=st[5];end=st[4]}
#                       return(paste(seq,start,end,sep=" "))
#                     }
#                   })
#   training=unname(training)
#   
#   # How many genes provided?
#   j=0
#   for (i in 1:length(training)){
#     if (training[i]==""){
#       j=j+1
#       cat(paste(as.character(i),"\n",sep=""),file="list",append=TRUE)
#     }
#   }
#   
#   c=5000 # how many genes for training?
#   if (j<=c){
#     writeLines(training,"exon_file")
#   }else{ # pick 1000 longest genes
#     exon_length=sapply(1:length(training),
#                        function(line){
#                          line=training[line]
#                          if (line==""){
#                            return(NA)
#                          }else{
#                            line=unlist(strsplit(line," "))
#                            start=as.numeric(line[2])
#                            end=as.numeric(line[3])
#                            return(end-start+1)
#                          }
#                        })
#     
#     # compute gene length
#     gene_length=rep(NA,j)
#     gene_index=1
#     gene_len=0
#     for (i in exon_length){
#       if (!is.na(i)){
#         gene_len=gene_len+i
#       }else{
#         gene_length[gene_index]=gene_len
#         gene_len=0
#         gene_index=gene_index+1
#       }
#     }
#     gene_length=abs(gene_length)
#     names(gene_length)=1:j
#     gene_length=sort(gene_length,decreasing=TRUE)
#     
#     # Pick long genes
#     index=names(head(gene_length,n=c))
#     index=as.numeric(index)
#     is=readLines("list")[index]
#     sapply(is,function(i){
#       i=as.numeric(i)
#       for (k in (i-1):1){
#         if (training[k]==""){break}
#       }
#       gene=training[(k-1):i]
#       for (m in gene){
#         cat(paste(m,"\n",sep=""),file="exon_file",append=TRUE)
#       }
#     })
#   }
#   
#   cmd=paste("trainGlimmerHMM",
#             genome,"exon_file",
#             "-d training",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("glimmerhmm_linux_x86_64",
#             genome,"training","-g",
#             "> GlimmerHMM.gff3",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="glimmerHMM_to_GFF3.pl GlimmerHMM.gff3 > GlimmerHMM_evm.gff3"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("agat_convert_sp_gxf2gxf.pl",
#             "--gff GlimmerHMM_evm.gff3",
#             "-o sorted.gff3",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   system("rm list",wait=TRUE)
#   system("rm GlimmerHMM_evm.agat.log")
#   system("rm tmp.gff3")
#   system("rm -r training")
#   system("rm training.log")
#   
#   setwd(wd_begin)
# }
# 
# # SNAP for gene prediction
# # Dependencies: agat, SNAP, scripts from EvidenceModeler
# snap=function(fna=fna,
#               train.gff3=train.gff3, # Gene models
#               out_dir=out_dir){
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   wd_begin=getwd();setwd(out_dir)
#   
#   cmd=paste("agat_convert_sp_gff2zff.pl",
#             "--fasta",fna,
#             "--gff",train.gff3,
#             "-o ./train",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   system("mv train.ann train.zff",wait=TRUE)
#   
#   cmd="fathom -validate train.zff train.dna > train.validate"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="fathom -categorize 100 train.zff train.dna"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="fathom -export 100 -plus uni.*"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="fathom -validate export.ann export.dna"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="forge export.ann export.dna"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="hmm-assembler.pl model . > model.hmm"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="snap model.hmm train.dna > SNAP.zff"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="zff2gff3.pl SNAP.zff > SNAP.gff3"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="SNAP_to_GFF3.pl SNAP.gff3 > SNAP_evm.gff3" # EvidenceModeler
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("agat_convert_sp_gxf2gxf.pl",
#             "--gff SNAP_evm.gff3",
#             "-o sorted.gff3",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   system("rm Acceptor*")
#   system("rm alt*")
#   system("rm Coding*")
#   system("rm Donor*")
#   system("rm E*")
#   system("rm e*")
#   system("rm I*")
#   system("rm P*")
#   system("rm Start*")
#   system("rm Stop*")
#   system("rm UTR*")
#   system("rm *.ann")
#   system("rm t*")
#   system("rm *.dna")
#   system("rm phaseprefs model.hmm")
#   system("rm SNAP_evm.agat.log")
#   setwd(wd_begin)
# }
# 
# # Extract features (gene+exon+CDS) by list of gene IDs
# SplitGFF=function(in.gff=in.gff,
#                   gene.list=gene.list,
#                   out.gff=out.gff){
#   cmd=paste("agat_sp_filter_feature_from_keep_list.pl",
#             "--gff",in.gff,
#             "--keep_list",gene.list,
#             "--out",out.gff)
#   print(cmd);system(cmd,wait=TRUE)
# }