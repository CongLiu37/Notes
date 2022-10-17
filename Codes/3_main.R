# Prepare genome assembly

# Check quality of fastq file.
# Dependencies: FastQC
QualityCheck=function(fq1=fq1, # Input fq file.
                      fq2=fq2, # Input fq file. Set "none" if single-end.
                      out_dir=out_dir, # directory in which quality check reports are saved.
                      threads=threads){
  threads = as.character(threads)
  if (fq2!="none"){ # pair-end
    cmd = paste("fastqc",
                "-o",out_dir,
                "-t",threads,
                fq1,fq2,
                sep=" ")
  }else{ # single-end
    cmd = paste("fastqc",
                "-o",out_dir,
                "-t",threads,
                fq1,
                sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  return(0)
}

# Filter NGS reads.
# Dependencies: Trimmomatic
QualityFilter = function(fq1=fq1,fq2=fq2, # Input fq files. Set fq2="none" if single-end.
                         clean_fq1=clean_fq1,clean_fq2=clean_fq2,# Names of clean fq files. clean_fq2 not used if single-end.
                         unpaired_fq1=unpaired_fq1,unpaired_fq2=unpaired_fq2, # Names of fq files for unpaired clean reads. unpaired_fq2 not used if single-end.
                         QualityFilter=QualityFilter, # String of step options of Trimmomatic.
                                                      # General steps included in quality filter:
                                                      # 1. adaptor; 2. low quality tail; 3. read length; 4. average quality score
                         threads=threads){
  threads = as.character(threads)
  
  if (fq2!="none"){ # pair-end
    cmd = paste("trimmomatic",
                "PE",
                "-threads",threads,
                "-phred33",
                fq1,fq2,
                clean_fq1,unpaired_fq1,
                clean_fq2,unpaired_fq2,
                QualityFilter,
                sep=" ")
  }else{ # single-end
    cmd = paste("trimmomatic",
                "SE",
                "-threads",threads,
                "-phred33",
                fq1,
                clean_fq1,
                QualityFilter,
                sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  
  return(0)
}

# FMLRC2: long read error correction with NGS short reads
# It corrects long reads based on short read assembly
# Dependencies: repobwt2, FMLRC2
fmlrc=function(long_reads=long_reads,
               short_reads=short_reads, # fq.gz, space-separated list
               out_dir=out_dir,
               threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("gunzip","-c",short_reads,"|",
            "awk 'NR % 4 == 2' |",
            "sort |",
            "tr NT TN |",
            "ropebwt2 -LR |",
            "tr NT TN |",
            "fmlrc2-convert comp_msbwt.npy",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("fmlrc2",
            "-t",threads,
            "comp_msbwt.npy",
            long_reads,
            "corrected_long_reads.fna",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# proovframe: frame-shift correction for long reads
#             can also be an additional polishment for assemblies
# Correct frame-shift based on alignments computed by DIAMOND
# Dependencies: proovframe, DIAMOND
proovframe=function(long_reads=long_reads, # can be assemblies
                    proteins.faa=proteins.faa, # Uniprot recommended
                    output_prefix=output_prefix,
                    threads=threads){
  threads=as.character(threads)
  
  cmd=paste("proovframe","map",
            "-t",threads,
            "-a",proteins.faa,
            "-o",paste(output_prefix,"_proovframe.tsv",sep=""),
            long_reads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("proovframe","fix",
            "-o",paste(output_prefix,"_corrected.fa",sep=""),
            long_reads,
            paste(output_prefix,"_proovframe.tsv",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Build kmc database
# Dependencies: KMC
kmc=function(fq1=fq1, # Input fq1
             fq2=fq2, # Input fq2. "none" if single-end.
             input_format=input_format, # "-fa" for FASTA, "-fq" for FASTQ
             # "-fm" for multiple FASTA
             # "-fbam" for BAM
             kmer=31,
             out_basename=out_basename, # kmc database
             out_dir=out_dir,
             threads=threads){
  threads=as.character(threads)
  out_dir=gsub("/$","",out_dir)
  kmer=as.character(kmer)
  wd=getwd()
  setwd(out_dir)
  
  cmd=paste("echo",fq1,">",
            paste(out_dir,"/",out_basename,".FILE",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  if (fq2!="none"){
    cmd=paste("echo",fq2,">>",
              paste(out_dir,"/",out_basename,".FILE",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("kmc",
            input_format,
            paste("-k",kmer,sep=""),
            paste("-t",threads,sep=""),
            "-m64","-ci1 -cs10000",
            paste("@",out_dir,"/",out_basename,".FILE",sep=""),
            out_basename, # kmc database
            ".",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",paste(out_dir,"/",out_basename,".FILE",sep=""),sep=" "))
  
  setwd(wd)
}

# Merge kmc database
merge_kmc=function(kmcdb1=kmcdb1,kmcdb2=kmcdb2,
                   out_prefix=out_prefix){
  cmd=paste("kmc_tools simple",
            kmcdb1,kmcdb2,
            "union",out_prefix,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# kmc database to kmer hist
# Dependencies: kmc
kmc_hist=function(kmcdb=kmcdb,
                  out_prefix=out_prefix){
  cmd=paste("kmc_tools transform",
            kmcdb,
            "histogram",
            paste(out_prefix,".hist",sep=""),
            "-cx10000",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Genome survey via Genomescope.
# Dependencies: GenomeScope, stringr (R)
Genomescope=function(hist=hist, # kmer histogram from kmc.
                     ploidy="none", # 1-6
                     out_dir=out_dir,
                     out_basename=out_basename,
                     kmer=31){
  kmer=as.character(kmer)
  out_dir=sub("/$","",out_dir)
  
  if (!file.exists(paste(out_dir,"/",out_basename,"_summary.txt",sep=""))){
    if (ploidy!="none"){
      out_basename=paste(out_basename,"_p",as.character(ploidy),sep="")
    }
    cmd=paste("genomescope.R",
              "-i",hist,
              "-o",out_dir,
              "-n",out_basename,
              "-k",kmer,
              sep=" ")
    if (ploidy!="none"){
      cmd=paste(cmd,"-p",as.character(ploidy))
    }
    print(cmd);system(cmd,wait=TRUE)
  }
  
  library(stringr)
  res=readLines(paste(out_dir,"/",out_basename,"_summary.txt",sep=""))
  Homozygous=unlist(str_extract_all(res[grepl("Homozygous",res)],"[0-9.]*%"))
  Heterozygous=unlist(str_extract_all(res[grepl("Heterozygous",res)],"[0-9.]*%"))
  HaploidLength=unlist(str_extract_all(res[grepl("Genome Haploid Length",res)],"[0-9,]* bp"))
  
  o=list(Homozygous=Homozygous,Heterozygous=Heterozygous,HaploidLength=HaploidLength)
  return(o)
}

# Genome survey via Smudgeplot.
# Dependencies: Smudgeplot, KMC
Smudgeplot=function(hist=hist, # kmer histogram from kmc.
                    kmcdb=kmcdb, # kmc database.
                    L="none",U="none", # integers estimated from kmer profile (Genomescope)
                    # L should be chosen so it cuts out majority of the error kmers
                    # U should be so it does contain all the k-mers that are at least in four genomic copies (but you might want to increase this number for genomes that might have higher ploidy levels).      
                    out_prefix=out_prefix){
  L=as.character(L);U=as.character(U)
  
  if (L=="none"){
    L=system(paste("smudgeplot.py cutoff ",hist," L",sep=""),wait=TRUE,intern=TRUE)
  }
  if (U=="none"){
    U=system(paste("smudgeplot.py cutoff ",hist," U",sep=""),wait=TRUE,intern=TRUE)
  }
  
  cmd=paste("kmc_tools transform",
            kmcdb,
            paste("-ci",L," ","-cx",U," ","dump",sep=""),
            "-s",paste(kmcdb,"_L",L,"_U",U,".dump",sep=""),
            sep=" ")
  print(cmd);system(cmd)
  
  cmd=paste("smudgeplot.py hetkmers -o",paste(out_prefix,"_L",L,"_U",U,sep=""),
            paste(kmcdb,"_L",L,"_U",U,".dump",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("smudgeplot.py plot","-o",paste(out_prefix,"_L",L,"_U",U,sep=""),
            paste(out_prefix,"_L",L,"_U",U,"_coverages.tsv",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# kmergenie: analyze NGS reads
# Dependencies: kmergenie
kmergenie=function(reads=reads, # vector of paths
                   out_prefix=out_prefix,
                   threads=threads){
  threads=as.character(threads)
  
  if (file.exists(paste(out_prefix,".FILE",sep=""))){
    system(paste("rm"," ",out_prefix,".FILE",sep=""),wait=TRUE)
  }
  
  for (i in reads){
    cmd=paste("echo",i,">>",paste(out_prefix,".FILE",sep=""),sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("kmergenie",
            paste(out_prefix,".FILE",sep=""),
            "-t",threads,
            "-o",out_prefix)
  print(cmd);system(cmd,wait=TRUE)
  
  return(out_prefix)
}

# Assemble organelle genome from WGS NGS reads
# Dependencies: GetOrganelle
getorganelle=function(fq1=fq1,
                      fq2=fq2,
                      target_organelle_type=target_organelle_type, # embplant_pt, embplant_mt, embplant_nr, fungus_mt, fungus_nr, animal_mt, and/or other_pt
                      seed.fna="none",
                      out_dir=out_dir,
                      threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  # GetOrganelle will make the out_dir
  
  cmd="get_organelle_from_reads.py"
  if (fq2!="none"){
    cmd=paste(cmd,"-1",fq1,"-2",fq2,sep=" ")
  }else{
    cmd=paste(cmd,"-u",fq1,sep=" ")
  }
  
  cmd=paste(cmd,
            "-t",threads,
            "-o",out_dir,
            "-F",target_organelle_type,
            sep=" ")
  
  if (seed.fna!="none"){
    cmd=paste(cmd,"-s",seed.fna,sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("summary_get_organelle_output.py",
            out_dir,
            "-o",paste(out_dir,"/summary.tsv",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Map reads to reference
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
    if (mapped){flag="-f 2"}else{flag="-G 2"}
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

# SPAdes: Genome/Metagenome assembly.
# long reads can be provided for scaffolding NGS assembly
# Dependencies: SPAdes
SPAdes=function(fq1=fq1,fq2=fq2, # Input NGS reads (fq). Set fq2="none" if single-end.
                                 # PacBio CCS (HiFi) reads should be treated as fq1.
                contigs.fa="none", # Reliable contigs of the same genome.
                pacbio_clr="none",nanopore="none",sanger="none", # Long reads
                meta=meta, # Logical. If TRUE, run metaSPAdes.
                           # metaSPAdes supports only a single short-read library which has to be paired-end.
                out_dir=out_dir,
                threads=threads
){
  threads=as.character(threads)
  
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  
  if (fq2!="none"){ # pair-end
    cmd=paste("spades.py",
              "-t",threads,
              "-1",fq1,
              "-2",fq2,
              "-o",out_dir,
              sep=" ")
  }else{
    cmd=paste("spades.py",
              "-t",threads,
              "-s",fq1,
              "-o",out_dir,
              sep=" ")
  }
  
  if (contigs.fa!="none"){cmd=paste(cmd,"--trusted-contigs",contigs.fa,sep=" ")}
  if (pacbio_clr!="none"){cmd=paste(cmd,"--pacbio",pacbio_clr,sep=" ")}
  if (nanopore!="none"){cmd=paste(cmd,"--nanopore",nanopore,sep=" ")}
  if (sanger!="none"){cmd=paste(cmd,"--sanger",sanger,sep=" ")}
  if (meta){cmd=paste(cmd,"--meta",sep=" ")}
  
  print(cmd);system(cmd,wait=TRUE)
  return(out_dir)
}

# Flye: assemble long reads
# Dependencies: flye
flye=function(long_reads=long_reads, # space-separated list
              read_type=read_type,
              meta=FALSE, # logical. TRUE for metagenomic long reads
              out_dir=out_dir,
              threads=threads){
  # read_type        Description
  # --pacbio-raw     PacBio regular CLR reads (<20% error)
  # --pacbio-corr    PacBio reads that were corrected with other methods (<3% error)
  # --pacbio-hifi    PacBio HiFi reads (<1% error)
  # --nano-raw       ONT regular reads, pre-Guppy5 (<20% error)
  # --nano-corr      ONT reads that were corrected with other methods (<3% error)
  # --nano-hq        ONT high-quality reads: Guppy5+ SUP or Q20 (<5% error)
  
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  
  cmd=paste("flye",
            read_type,long_reads,
            "--out-dir",out_dir,
            "-t",threads,
            sep=" ")
  if (meta){cmd=paste(cmd,"--meta",sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
}

# Pilon: polish long read assembly using short reads
# Dependencies: pilon
pilon=function(assembly=assembly,
               bam=bam, # Map NGS to assembly
               out_basename=out_basename,
               threads=threads,
               out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  cmd=paste("pilon",
            "--changes",
            "--genome",assembly,
            "--bam",bam,
            "--output",out_basename,
            "--outdir",".",
            "--threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# Run pilon iteratively
# Dependencies: Minimap2, SAMtools, pilon
run_pilon=function(assembly=assembly,
                   NGS=NGS, # space-separated list
                   out_dir=out_dir,
                   out_basename=out_basename,
                   iteration=4, 
                   threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  for (i in 1:iteration){
    system(paste("mkdir"," ","Pilon_",as.character(i),sep=""))
    minimap2(long_reads=NGS, # space-separated list for PE
             lr_type="sr", # long read type. 
             assembly=assembly,
             out_prefix=paste("Pilon_",as.character(i),"/",out_basename,sep=""),
             threads=threads)
    pilon(assembly=assembly,
          bam=paste("Pilon_",as.character(i),"/",out_basename,".bam",sep=""), # Map NGS to assembly
          out_basename=out_basename,
          threads=threads,
          out_dir=paste("Pilon_",as.character(i),sep=""))
    assembly=paste("Pilon_",as.character(i),"/",out_basename,".fasta",sep="")
    changes=paste("Pilon_",as.character(i),"/",out_basename,".changes",sep="")
    c=system(paste("wc -l"," ",
                   "Pilon_",as.character(i),"/",out_basename,".changes",sep=""),
             wait=TRUE,intern=TRUE)
    if (c=="0"){break}
  }
  c=system(paste("wc -l"," ",
                 "Pilon_",as.character(iteration),"/",out_basename,".changes",sep=""),
           wait=TRUE,intern=TRUE)
  if (c=="0"){
    print(out_basename)
    print("Polishment finished!")
    system(paste("mv"," ",
                 "Pilon_",as.character(iteration),"/",out_basename,"fasta",
                 " ",
                 "./",out_basename,"_pilon.fna",sep=""),wait=TRUE)
  }else{
    print(out_basename)
    print("More iteration needed")
  }
  
  setwd(wd)
}

# Hypo: polish long read assembly
# NGS paired reads or PacBio HiFi reads are required
# Long reads are optional
# Dependencies: minimap2, SAMtools, hypo
hypo=function(assembly=assembly, # draft assembly of long reads
              fq1="none",fq2="none",NGS.bam="none", # pair-end NGS reads.
              hifi="none",hifi.bam="none", # PacBio HiFi reads. "none" if not provided
              long_reads.bam="none", # BAM for long reads
              genome_size=genome_size, # taking k/m/g as unit, e.g. 3.2m
              coverage=coverage, # int. mean coverage of the short/HiFi reads 
              out_dir=out_dir,
              polished_assembly=polished_assembly, # file name
              threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd_tmp=paste("hypo",
                "-d",assembly,
                "-s",as.character(genome_size),
                "-c",as.character(coverage),
                "-t",threads,
                "-o",polished_assembly,
                sep=" ")
  if (fq1!="none"){ # NGS reads provided
    cmd=paste("echo"," ",fq1,"\n",fq2," > ","fq.file",sep="")
    print(cmd);system(cmd,wait=TRUE)
    cmd_tmp=paste(cmd_tmp,
                  "-r","fq.file",
                  "-b",NGS.bam,
                  sep=" ")
  }
  if (hifi!="none"){
    cmd_tmp=paste(cmd_tmp,
                  "-r",hifi,
                  "-b",hifi.bam,
                  sep=" ")
  }
  if (long_reads!="none"){
    cmd_tmp=paste(cmd_tmp,
                  "-B",long_reads.bam,
                  sep=" ")
  }
  
  print(cmd_tmp);system(cmd_tmp,wait=TRUE)
  setwd(wd)
}

# SALSA: scaffolding assembly with paired short reads
# Dependencies: BEDtools, SAMtools, SALSA
salsa=function(assembly=assembly,
               bam=bam, # map PE short reads to assembly
               scaffolds.fna=scaffolds.fna, # file name for scaffolds
               out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("bamToBed",
            "-i",bam,
            ">","alignment.bed",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="sort -k 4 alignment.bed > tmp && mv tmp alignment.bed"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("samtools","faidx",
            "--fai-idx","assembly.fai",
            assembly,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("run_pipeline.py",
            "-a",assembly,
            "-l","assembly.fai",
            "-b","alignment.bed",
            "-o",scaffolds.fna,
            "-m yes",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}







# Bowtie2: Map reads to genome. 
# SAMtools: Compress SAM to BAM, sort BAM, index BAM.
# Dependencies: Bowtie2, SAMtools
Bowtie2 = function(fq1=fq1,fq2=fq2, # Input fq files. Make fq2="none" if single-end.
                   # Can be comma-separated list of files if multiple libraries used
                   fna=fna, # FASTA of genome
                   index=index, # Basename of Bowtie2 index of reference genome.
                   out_prefix=out_prefix, # Prefix of output BAM file.
                   threads=threads){
  threads = as.character(threads)
  
  cmd = paste("bowtie2-build",
              fna,
              index,
              sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  bam_filename=paste(out_prefix,".bam",sep="")
  if (fq2!="none"){ # pair-end
    cmd = paste("bowtie2",
                "-x",index,
                "-p",threads,
                "-1",fq1,
                "-2",fq2,
                "|",
                "samtools","view",
                "-@",threads,"-bS",
                "|",
                "samtools","sort",
                "-@",threads,
                "-o",bam_filename,
                sep=" ")
  }else{ # single pair
    cmd = paste("bowtie2",
                "-x",index,
                "-p",threads,
                "-U",fq1,
                "|",
                "samtools","view",
                "-@",threads,"-bS",
                "|",
                "samtools","sort",
                "-@",threads,
                "-o",bam_filename,
                sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("samtools","index",bam_filename,sep=" ")
  print(cmd);system(cmd)
  
  return(bam_filename)
}