# Genome/Metagenome assembly, binning and assembly quality assessment
# Dependencies:
#    Softwares: PhyloFlash, FastQC, Trimmomatic, SPAdes, BUSCO, CheckM, Bowtie2, SAMtools
#    R packages: ape, Biostrings

# PhyloFlash: Taxon profiling of reads/assembly via reconstructing the SSU rRNAs.
PhyloFlash=function(fq1=fq1,fq2=fq2, # Input fq files. Set fq2="none" if single-end.
                    contigs.fa=contigs.fa, # If not "none", SSU rRNA will be called from contigs to screen vs reads
                    out_prefix=out_prefix,
                    threads=threads){
  threads=as.character(threads)
  if (contigs.fa=="none"){
    if (fq2!="none"){ # pair-end
      cmd=paste("phyloFlash.pl","-almosteverything","-lib",out_prefix,"-read1",fq1,"-read2",fq2,"-CPUs",threads,sep=" ")
    }else{ # single-end
      cmd=paste("phyloFlash.pl","-almosteverything","-lib",out_prefix,"-read1",fq1,"-CPUs",threads,sep=" ")
    }
  }else{
    if (fq2!="none"){ # pair-end
      cmd=paste("phyloFlash.pl","-almosteverything","-lib",out_prefix,"-read1",fq1,"-read2",fq2,"-trusted",contigs.fa,"-CPUs",threads,sep=" ")
    }else{ # single-end
      cmd=paste("phyloFlash.pl","-almosteverything","-lib",out_prefix,"-read1",fq1,"-trusted",contigs.fa,"-CPUs",threads,sep=" ")
    }
  }
  
  print(cmd);system(cmd,wait=TRUE)
  return(0)
}

# FastQC: Check quality of sequencing data.
QualityCheck=function(fq1=fq1, # Input fq file.
                      fq2=fq2, # Input fq file. Set "None" if single-end.
                      out_dir=out_dir, # directory in which quality check reports are saved.
                      threads=threads){
  threads = as.character(threads)
  if (fq2!="None"){ # pair-end
    cmd = paste("fastqc","-o",out_dir,"-t",threads,fq1,fq2,sep=" ")
  }else{ # single-end
    cmd = paste("fastqc","-o",out_dir,"-t",threads,fq1,sep=" ")
  }
  print(cmd)
  system(cmd,wait=TRUE)
  return(0)
}

# Trimmomatic: Quality control of sequencing data.
QualityFilter = function(fq1=fq1,fq2=fq2, # Input fq files. Set fq2="None" if single-end.
                         clean_fq1=clean_fq1,clean_fq2=clean_fq2,# Names of clean fq files. Set clean_fq2="None" if single-end.
                         unpaired_fq1=unpaired_fq1,unpaired_fq2=unpaired_fq2, # Names of fq files for unpaired clean reads. Set unpaired_fq2="None" if single-end.
                         QualityFilter=QualityFilter, # String of step options of Trimmomatic.
                         # General steps included in quality filter:
                         # 1. adaptor; 2. low quality tail; 3. read length; 4. average quality score
                         threads=threads){
  threads = as.character(threads)
  
  if (fq2!="None"){ # pair-end
    cmd = paste("trimmomatic","PE","-threads",threads,"-phred33",fq1,fq2,clean_fq1,unpaired_fq1,clean_fq2,unpaired_fq2,QualityFilter,sep=" ")
  }else{ # single-end
    cmd = paste("trimmomatic","SE","-threads",threads,"-phred33",fq1,clean_fq1,QualityFilter,sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  
  return(0)
}

# SPAdes: Genome/Metagenome assembly.
SPAdes=function(fq1=fq1,fq2=fq2, # Input fq files. Set fq2="none" if single-end.
                                 # PacBio CCS reads should be treated as fq1.
                contigs.fa=contigs.fa, # Reliable contigs of the same genome. Set "none" if unavailable.
                pacbio_clr=pacbio_clr,nanopore=nanopore,sanger=sanger, # File of long reads. set "none" if unavailable.
                meta=meta, # Logical. If TRUE, run metaSPAdes.
                           # metaSPAdes supports only a single short-read library which has to be paired-end.
                out_dir=out_dir,
                threads=threads
){
  threads=as.character(threads)
  
  if (fq2!="none"){ # pair-end
    cmd=paste("spades.py","-t",threads,"-1",fq1,"-2",fq2,"-o",out_dir,sep=" ")
  }else{
    cmd=paste("spades.py","-t",threads,"-s",fq1,"-o",out_dir,sep=" ")
  }
  
  if (contigs.fa!="none"){cmd=paste(cmd,"--trusted-contigs",contigs.fa,sep=" ")}
  if (pacbio_clr!="none"){cmd=paste(cmd,"--pacbio",pacbio_clr,sep=" ")}
  if (nanopore!="none"){cmd=paste(cmd,"--nanopore",nanopore,sep=" ")}
  if (sanger!="none"){cmd=paste(cmd,"--sanger",sanger,sep=" ")}
  if (meta){cmd=pastes(cmd,"--meta")}
  
  print(cmd);system(cmd,wait=TRUE)
  return(out_dir)
}

# Basic statistics of assembly.
# N50, min/max/median contig size, gap percent, GC content
BSG=function(fna=fna,
             GenomeSize=GenomeSize, # Genome size estimated by cytological technique
                                    # Set "none" if unavailable.
             Threshold=Threshold # Contigs below this threshold are defined as small contigs. 
                                 # Small contigs are not considered, and their statistics are reported separatedly.
){
  Threshold=as.numeric(Threshold)
  
  library(ape)
  data=read.dna(fna,format="fasta",as.character=TRUE)
  ContigSizes=sapply(data,length)
  
  # Statistics of contigs smaller than Threshold
  SmallContigSizes=ContigSizes[which(ContigSizes<=Threshold)]
  SmallContigCount=length(SmallContigSizes)
  TotalSmallContigSizes=sum(SmallContigSizes)
  
  data=data[which(ContigSizes>Threshold)] # Remove contigs smaller than Threshold
  ContigSizes=sapply(data,length)
  
  f1=function(sequ){return(length(sequ[which(sequ%in%c("N","n"))]))}
  GapSizes=sapply(data,f1)
  
  f2=function(sequ){return(length(sequ[which(sequ%in%c("G","C","g","c"))]))}
  GcSizes=sapply(data,f2)
  
  library(Biostrings)
  AssemblySize=sum(ContigSizes)
  ContigCount=length(ContigSizes)
  N50=N50(ContigSizes)
  Min=min(ContigSizes);Max=max(ContigSizes);Median=median(ContigSizes);Mean=mean(ContigSizes)
  GapPercent=sum(GapSizes)/sum(ContigSizes)
  if (GenomeSize!="none"){LateralCoverage=sum(ContigSizes)/GenomeSize}else{LateralCoverage="none"}
  GcContent=sum(GcSizes)/sum(ContigSizes)
  
  o=list(fna=fna,ThresholdSize=Threshold,
         SmallContigCount=SmallContigCount,TotalSmallContigSizes=TotalSmallContigSizes,
         AssemblySize=AssemblySize,
         ContigCount=ContigCount,
         N50=N50,MinContig=Min,MaxContig=Max,MedianContig=Median,MeanContig=Mean,
         GapPercent=GapPercent,
         LateralCoverage=LateralCoverage,
         GC=GcContent)
  return(o)
}

# BUSCO: Assess genome completeness via searching universal single-copy orthologue genes by BUSCO.
BUSCO=function(fna=fna, # Fasta file of nucleotide or protein.
                        # Be consistent with Mode.
               Mode=Mode, # genome/proteins/transcriptome
               Lineage=Lineage, # Lineage dataset, e.g. insecta_odb10
                                # Available datasets: https://busco-data.ezlab.org/v5/data/lineages/
                                # BUSCO will download lineage dataset automatically.
               Out_prefix=Out_prefix, # Give the analysis run a recognisable short name.
                                      # Output folders and files will be labelled with this name.
               Threads=Threads
){
  Threads=as.character(Threads)
  
  cmd=paste("busco","--in",fna,"--lineage_dataset",Lineage,"--out",Out_prefix,"--mode",Mode,"--cpu",Threads,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Remove temporaries
  # cmd=paste("rm","-r","augustus_config",sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  
  # Print result summary
  cmd=paste("cat",paste(Out_prefix,"/short_summary.specific.",Lineage,".",Out_prefix,".txt"),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # BUSCO completeness
  re=readLines(paste(Out_prefix,"/short_summary.specific.",Lineage,".",Out_prefix,".txt"))[9]
  re=gsub("\t*C:","",re)
  re=strsplit(re,"%")[[1]][1]
  re=as.numeric(re)/100
  print("BUSCO completeness:")
  print(re)
  
  return(list(BuscoCompleteness=re,Database=Lineage))
}

# CheckM: Assess completeness and contamination of genomic bins via using collocated sets of genes that are ubiquitous and single-copy within a phylogenetic lineage by CheckM.
checkm=function(bin_dir=bin_dir, # the directory in which bins are located.
                bin_basename=bin_basename, # the base name of bins, e.g. fna
                out_dir=out_dir,
                threads=threads){
  threads=as.character(threads)
  
  cmd=paste("checkm","lineage_wf","-x",bin_basename,"-t",threads,bin_dir,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(out_dir)
}

# Bowtie2; Build bowtie2 index of genome.
Bowtie2Build = function(fna=fna, # FASTA of genome
                        index_prefix=index_prefix){
  cmd = paste("bowtie2-build",fna,index_prefix,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(index_prefix)
}

# Bowtie2: Map reads to genome. 
# SAMtools: Compress SAM to BAM and sort BAM.
Bowtie2 = function(fq1=fq1,fq2=fq2, # Input fq files. Make fq2="none" if single-end.
                   index=index, # Basename of Hisat2 index of reference genome.
                   out_prefix=out_prefix, # Prefix of output BAM file.
                   threads=threads){
  threads = as.character(threads)
  
  bam_filename=paste(out_prefix,".bam",sep="")
  if (fq2!="none"){ # pair-end
    cmd = paste("bowtie2","-x",index,"-p",threads,"-1",fq1,"-2",fq2,"|",
                "samtools","view","-@",threads,"-bS","|",
                "samtools","sort","-@",threads,"-o",bam_filename,sep=" ")
  }else{ # single pair
    cmd = paste("bowtie2","-x",index,"-p",threads,"-U",fq1,"|",
                "samtools","view","-@",threads,"-bS","|",
                "samtools","sort","-@",threads,"-o",bam_filename,sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  return(bam_filename)
}

# SprayNPray: Binning based on tetranucleotide frequency, GC-content, codon usage bias, and read coverage.
Binning=function(fna=fna,
                 bam=bam,
                 out_basename=out_basename,
                 threads=threads){
  threads=as.character(threads)
  
  cmd=paste("spray-and-pray.py","-g",fna,"-bam",bam,"-out",out_basename,"-t",threads,"--bin",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(out_basename)
}
