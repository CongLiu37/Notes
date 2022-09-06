# Prepare genome assembly

# Profile contaminations via reconstructing the SSU rRNAs.
# Dependencies: PhyloFlash
PhyloFlash=function(fq1=fq1,fq2=fq2, # Input fq files. Set fq2="none" if single-end.
                    contigs.fa="none", # If provided, SSU rRNA will be called from contigs to screen vs reads
                    out_prefix=out_prefix,
                    threads=threads){
  threads=as.character(threads)
  
  cmd=paste("phyloFlash.pl",
            "-almosteverything",
            "-lib",out_prefix,
            "-CPUs",threads,
            "-read1",fq1,
            sep=" ")
  if (fq2!="none"){ # pair-end
    cmd=paste(cmd,
              "-read2",fq2,
              sep=" ")
  }
  if (contigs.fa!="none"){ # contigs provided
    cmd=paste(cmd,
              "-trusted",contigs.fa,
              sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  return(out_prefix)
}

# Check quality of short reads (fq).
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

# SPAdes: Genome/Metagenome assembly.
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

# FMLRC: long read error correction with NGS short reads
# Dependencies: repobwt2, FMLRC v.2
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

# proovframe: frame-shift correction for long read
# Dependencies: proovframe, DIAMOND
proovframe=function(long_reads=long_reads,
                    proteins.faa=proteins.faa, # Uniprot recommended
                    output_prefix=output_prefix,
                    threads=threads){
  threads=as.character(threads)
  
  cmd=paste("proovframe","map",
            "-t",threads,
            "-a",proteins.faa,
            "-o",paste(output_prefix,".tsv",sep=""),
            long_reads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("proovframe","fix",
            "-o",paste(output_prefix,"_corrected.fa"),
            long_reads,
            paste(output_prefix,".tsv",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Flye: assemble long reads
# Dependencies: flye
flye=function(long_reads=long_reads, # space-separated list
              read_type=read_type,
              meta=meta, # logical. TRUE for metagenomic long reads
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
# Dependencies: minimap2, pilon
pilon=function(assembly=assembly,
               fq1=fq1,fq2=fq2,
               polished_assembly=polished_assembly, # file name
               threads=threads,
               out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  cmd=paste("minimap2",
            "--secondary=no","--MD",
            "-ax","sr",
            "-t",threads,
            assembly,fq1,fq2,"|",
            "samtools","view","-@",threads,"-bS","|",
            "samtools","sort",
            "-@",threads,
            "-o","NGS.bam",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("pilon",
            "--genome",assembly,
            "--bam","NGS.bam",
            "--output",polished_assembly,
            "--outdir",".",
            "--threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}

# Hypo: polish long read assembly
# NGS pair reads or PacBio HiFi reads are required
# Long reads are optional
# Dependencies: minimap2, SAMtools, hypo
hypo=function(assembly=assembly, # draft assembly of long reads
              fq1="none",fq2="none", # pair-end NGS reads.
              hifi="none", # PacBio HiFi reads. "none" if not provided
              long_reads="none", # long reads for polishing
              lr_type="none", # long read type. "map-pb" for PacBio and "map-ont" for ONT reads.
              genome_size=genome_size, # int taking k/m/g as unit
              coverage=coverage, # str. mean coverage of the short/HiFi reads 
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
                "-c",coverage,
                "-t",threads,
                "-o",polished_assembly,
                sep=" ")
  
  if (fq1!="none"){ # NGS reads provided
    cmd=paste("minimap2",
              "--secondary=no","--MD",
              "-ax","sr",
              "-t",threads,
              assembly,fq1,fq2,"|",
              "samtools","view","-@",threads,"-bS","|",
              "samtools","sort",
              "-@",threads,
              "-o","NGS.bam",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("samtools","index","NGS.bam","-@",threads,sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("echo"," ",fq1,"\n",fq2," > ","fq.file",sep="")
    print(cmd);system(cmd,wait=TRUE)
    cmd_tmp=paste(cmd_tmp,
                  "-r","fq.file",
                  "-b","NGS.bam",
                  sep=" ")
  }
  if (hifi!="none"){
    cmd=paste("minimap2","--secondary=no --MD",
              "-ax","asm20",
              "-t",threads,
              assembly,hifi,"|",
              "samtools","view","-@",threads,"-bS","|",
              "samtools","sort",
              "-@",threads,
              "-o","HIFI.bam",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("samtools","index","HIFI.bam","-@",threads,sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd_tmp=paste(cmd_tmp,
                  "-r",hifi,
                  "-b","HIFI.bam",
                  sep=" ")
  }
  if (long_reads!="none"){
    cmd=paste("minimap2","--secondary=no","--MD",
              "-ax",lr_type,
              "-t",threads,
              assembly,long_reads,"|",
              "samtools","view","-@",threads,"-bS","|",
              "samtools","sort",
              "-@",threads,
              "-o","LR.bam",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("samtools","index","LR.bam","-@",threads,sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd_tmp=paste(cmd_tmp,
                  "-B","LR.bam",
                  sep=" ")
  }
  
  print(cmd_tmp);system(cmd_tmp,wait=TRUE)
  setwd(wd)
}

# SALSA: scaffolding assembly with paired short reads
# Dependencies: BEDtools, SALSA
salsa=function(assembly=assembly,
               bam=bam, # map PE short reads to assembly
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
            "-o","scaffolds.fna",
            "-m yes",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}

# Quast: Quality of assembly.
# Dependencies: QUAST
Quast=function(fna=fna, # FASTA of genome assembly.
               out_dir=out_dir,
               threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  
  if (!file.exists(paste(out_dir,"/report.tsv",sep=""))){
    if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))
  }
  
    cmd=paste("quast --min-contig 0",
              "-o",out_dir,
              "-t",threads,
              fna,
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  res=read.table(paste(out_dir,"/report.tsv",sep=""),
                 header=FALSE,row.names=1,sep="\t",quote="",comment.char="")
  o=list(Assembly=fna,
         ContigCount=res["# contigs",],
         MaxContig=res["Largest contig",],
         TotalSize=res["Total length",],
         GCpercent=res["GC (%)",],
         N50=res["N50",],
         N90=res["N90",],
         auN=res["auN",],
         L50=res["L50",],
         L90=res["L90",],
         GapPer100kb=res["# N's per 100 kbp",])
  return(o)
}

# BUSCO: Assess genome completeness via searching universal single-copy orthologue genes by BUSCO.
# Dependencies: BUSCO
BUSCO=function(fna=fna, # Fasta file of nucleotide or protein.
                        # Be consistent with Mode.
               Mode=Mode, # genome/proteins/transcriptome
               Lineage=Lineage, # Lineage dataset, e.g. insecta_odb10
                                # Available datasets: https://busco-data.ezlab.org/v5/data/lineages/
                                # BUSCO will download lineage dataset automatically.
               Out_prefix=Out_prefix, # Give the analysis run a recognisable short name.
                                      # Output folders and files will be labelled with this name.
                                      # Cannot be path
               out_dir=out_dir,
               Threads=Threads){
  Threads=as.character(Threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  
  if (!file.exists(paste(out_dir,"/",Out_prefix,"/short_summary.specific.",Lineage,".",Out_prefix,".txt",sep=""))){
    wd_begin=getwd();setwd(out_dir)
    cmd=paste("busco",
              "--in",fna,
              "--lineage_dataset",Lineage,
              "--out",Out_prefix,
              "--mode",Mode,
              "--cpu",Threads,
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    setwd(wd_begin)
  }
  
  # BUSCO results
  re=readLines(paste(out_dir,"/",Out_prefix,"/short_summary.specific.",Lineage,".",Out_prefix,".txt",sep=""))[9]
  re=sub("\t","",re);re=sub("\t   ","",re)
  
  return(re)
}

# CheckM: Assess completeness and contamination of genomic bins via using collocated sets of genes that are ubiquitous and single-copy within a phylogenetic lineage by CheckM.
# CheckM is dependent on python3, HMMER (>=3.1b1), prodigal (2.60 or >=2.6.1), pplacer (>=1.1). 
checkm=function(bin_dir=bin_dir, # the directory in which bins are located.
                bin_basename=bin_basename, # the base name of bins, e.g. fna
                out_dir=out_dir, # directory for output.
                threads=threads){
  threads=as.character(threads)
  
  cmd=paste("checkm",
            "lineage_wf",
            "-x",bin_basename,
            "-t",threads,
            bin_dir,
            out_dir,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(out_dir)
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
  
  cmd=paste("samtools","index",
            bam_filename,
            "-@",threads,
            sep=" ")
  print(cmd);system(cmd)
  
  return(bam_filename)
}

# Metabat2: Unsupervised binning by Metabat2.
# Dependencies: Metabat2
Metabat=function(fna=fna, # fna. Input DNA sequences.
                 bam=bam, # BAM. For coverage computation.
                 MinContig=MinContig, # Integer. Input DNA sequences shorter than MinContig are removed.
                 out_prefix=out_prefix,
                 threads=threads){
  threads=as.character(threads)
  
  cmd=paste("jgi_summarize_bam_contig_depths",
            "--outputDepth",paste(out_prefix,".depth.txt",sep=""),
            bam,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("metabat2",
            "-i",fna, 
            "-a",paste(out_prefix,".depth.txt",sep=""),
            "-o",out_prefix,
            "-m",as.character(MinContig),
            "-t",threads,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(out_prefix)
}

# DNA-DNA search by blastn.
# Compute contig length, coverage, GC-content, taxonomy at major ranks by blobtools.
# Dependencies: blast, blobtools
blastn_blobtools=function(fna=fna, # fna. Input DNA sequences.
                          bam=bam, # BAM. For coverage computation.
                          out_basename=out_basename,
                          blast_dir=blast_dir, # Directory for blastn output.
                          blob_dir=blob_dir, # Directory for blobtools output.
                          blastn_db=blastn_db, # blast database.
                          threads=threads){
  threads=as.character(threads)
  blast_dir=sub("/$","",blast_dir)
  blob_dir=sub("/$","",blob_dir)
  
  cmd=paste("blastn",
            "-query",fna,
            "-db",blastn_db,
            "-out",paste(blast_dir,"/",out_basename,".blast",sep=""),
            "-outfmt","'6 qseqid staxids bitscore std'",
            "-max_target_seqs 10",
            "-max_hsps 1",
            "-evalue 1e-25",
            "-num_threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("blobtools create",
            "-i",fna,
            "-b",bam,
            "-t",paste(blast_dir,"/",out_basename,".blast",sep=""),
            "-o",paste(blob_dir,"/",out_basename,sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("blobtools view",
            "-i", paste(blob_dir,"/",out_basename,".blobDB.json",sep=""),
            "-o",paste(blob_dir,"/",out_basename,sep=""),
            "-r","all",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# DNA-protein search by diamond.
# Compute taxonomy at major ranks by Megan.
# Diamond and Megan in long-read mode.
# Dependencies: DIAMOND, MEGAN
Diamond_Megan=function(fna, # fna. Input DNA sequences.
                       out_basename=out_basename,
                       blast_dir=blast_dir, # Directory for diamond output.
                       rma_dir=rma_dir, # Directory for rma output of megan.
                       assignment_dir=assignment_dir, # Directory for taxonomy table.
                       ref_diamond=ref_diamond, # Diamond database.
                       ref_megan=ref_megan, # Megan database.
                       threads=threads){
  threads=as.character(threads)
  blast_dir=sub("/$","",blast_dir)
  rma_dir=sub("/$","",rma_dir)
  assignment_dir=sub("/$","",assignment_dir)
  
  cmd=paste("diamond blastx",
            "-p",threads,
            "-d",ref_diamond,
            "-q",fna,
            "--long-reads",
            "--out",paste(blast_dir,"/",out_basename,".blast",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("blast2rma",
            "-i",paste(blast_dir,"/",out_basename,".blast",sep=""),
            "-o",paste(rma_dir,"/",out_basename,".rma",sep=""),
            "-f","BlastTab",
            "-bm","BlastX",
            "--paired","false",
            "-lg","true",
            "-mdb",ref_megan,
            "-t",threads,
            "-ram","readCount",
            "-supp","0",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("rma2info",
            "-i",paste(rma_dir,"/",out_basename,".rma",sep=""),
            "-o",paste(assignment_dir,"/",out_basename,".tsv",sep=""),
            "-r2c Taxonomy",
            "-n true",
            "-p true",
            "-r true",
            "-mro","true",
            "-u false",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(paste(assignment_dir,"/",out_basename,".tsv",sep=""))
}

# SprayNPray: Predict genes by Prodigal, and compute contig length, coding density, coverage, GC-content, average AAI, domain of closest DIAMOND hits.
# SprayNPray is dependent on DIAMOND, Prodigal, Metabat, Python3, Biopython3, Joblib.
SprayNPray=function(fna=fna, # fna. Input DNA sequences.
                    bam=bam, # BAM. For coverage computation.
                    ref=ref, # Diamond database.
                    blast="none",
                    out_basename=out_basename, # cannot be path.
                                               # output files are stored a directory named as out_basename
                    out_dir=out_dir, # the directory where out_basename are moved to.
                    threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  wd=getwd()
  setwd(out_dir)
  
  cmd=paste("spray-and-pray.py",
            "-g",fna,
            "-bam",bam,
            "-out",out_basename,
            "-lvl","Domain",
            "-t",threads,
            "-ref",ref,
            sep=" ")
  if (blast!="none"){cmd=paste(cmd,"-blast",blast,sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
  return(paste(out_dir,"/",out_basename,sep=""))
}

# TAGC plot
# Taxa supported less than 0.1% of total scaffolds are removed as false positives.
# Dependencies: ggplot2 (R), stringr (R)
TAGC=function(df, # Each row represents a scaffold.
                  # Columns include cov, GC.content (e.g. 38.1), taxon (e.g. [D] Bacteria; [P] Proteobacteria)
              rank=rank, # Major taxonomy rank to be included in plot.
              title="", # title of plot.
              out.pdf=out.pdf){
  library(ggplot2);library(stringr)
  
  d=df;d[is.na(d)]="none";d[d==""]="none"
  r=substr(rank,1,1);r=toupper(r)
  
  pattern=paste("\\[",r,"\\] [A-Za-z]*;",sep="")
  plot_taxa=str_extract(d[,"taxon"],pattern)
  plot_taxa=plot_taxa[!is.na(plot_taxa)]
  t=table(plot_taxa)
  plot_taxa=names(t[which(t>0.001*sum(t))])
  data=d[str_extract(d[,"taxon"],pattern)%in%plot_taxa,]
  
  p=ggplot()+
    geom_point(data=data,size=0.2,
               aes(x=log10(cov),y=GC.content,
                   color=str_extract(data[,"taxon"],pattern)))+
    theme_classic()+labs(title=title,color="taxonomy",x="log10(coverage)",y="GC%")
  pdf(out.pdf);print(p);dev.off()
  return(p)
}

# Evaluate contamination level of draft genome by GC-coverage plot, coverage distribution and GC distribution.
# Dependencies: ggplot2 (R), ggExtra (R)
ContaminationPlot=function(df, # Each row represents a scaffold.
                               # Columns include cov, GC.content (e.g. 38.1)
                          title="", # plot title
                          out.pdf=out.pdf){
  library(ggplot2);library(ggExtra)
  p1=ggplot(df)+geom_point(size=0.2,aes(x=cov,y=GC.content))+
    theme_classic()+
    labs(title=title,x="coverage",y="GC%")+
    scale_x_continuous(limits=c(0,100),expand=c(0,0))+
    scale_y_continuous(limits=c(0,100),expand=c(0,0))
  p=ggMarginal(p1,type="histogram",size=10,
             xparams=list(bins=100),
             yparams=list(bins=100))
  pdf(out.pdf);print(p);dev.off()
  return(p)
}

# Extract short reads (fq) mapped to an assembly via SAMtools.
# Dependencies: Bowtie2, SAMtools
Extract_fq=function(fq1=fq1,fq2=fq2, # Original fq. fq2="none" id single-end.
                    fna=fna, # target assembly
                    out_prefix=out_prefix,
                    threads=threads
){
  threads=as.character(threads)
  
  cmd = paste("bowtie2-build",
              fna,
              paste(out_prefix,"_index",sep=""),
              sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  if (fq2!="none"){ # pair-end
    cmd = paste("bowtie2",
                "-x",paste(out_prefix,"_index",sep=""),
                "-p",threads,
                "-1",fq1,
                "-2",fq2,
                "|",
                "samtools fastq",
                "-@",threads,
                "-G","12",
                "-1",paste(out_prefix,".1.fq.gz",sep=""),
                "-2",paste(out_prefix,".2.fq.gz",sep=""),
                sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }else{ # single pair
    cmd = paste("bowtie2",
                "-x",paste(out_prefix,"_index",sep=""),
                "-p",threads,
                "-U",fq1,
                "|",
                "samtools fastq",
                "-@",threads,
                "-G","4",
                "-0",paste(out_prefix,".fq.gz",sep=""),
                sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  return(0)
}

# k-mer spectrum via KMC.
# Dependencies: KMC
kmc=function(fq1=fq1, # Input fq1
             fq2=fq2, # Input fq2. "none" if single-end.
             kmer=21,
             out_basename=out_basename,
             out_dir=out_dir,
             threads=threads){
  threads=as.character(threads)
  out_dir=gsub("/$","",out_dir)
  kmer=as.character(kmer)
  wd=getwd()
  setwd(out_dir)
  
  cmd=paste("echo",
            fq1,
            ">",
            paste(out_dir,"/",out_basename,".FILE",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  if (fq2!="none"){
    cmd=paste("echo",
              fq2,
              ">>",
              paste(out_dir,"/",out_basename,".FILE",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("kmc",
            paste("-k",kmer,sep=""),
            paste("-t",threads,sep=""),
            "-m64","-ci1 -cs10000",
            paste("@",out_dir,"/",out_basename,".FILE",sep=""),
            out_basename,
            ".",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("kmc_tools transform",
            paste(out_dir,"/",out_basename,sep=""),
            "histogram",paste(out_dir,"/",out_basename,".hist",sep=""),
            "-cx10000",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",paste(out_dir,"/",out_basename,".FILE",sep=""),sep=" "))
  
  setwd(wd)
}

# Genome survey via Genomescope.
# Dependencies: GenomeScope, stringr (R)
Genomescope=function(hist=hist, # kmer histogram from kmc.
                     out_dir=out_dir,
                     out_basename=out_basename,
                     kmer=kmer){
  kmer=as.character(kmer)
  out_dir=sub("/$","",out_dir)
  
  if (!file.exists(paste(out_dir,"/",out_basename,"_summary.txt",sep=""))){
    cmd=paste("genomescope.R",
              "-i",hist,
              "-o",out_dir,
              "-n",out_basename,
              "-k",kmer,
              sep=" ")
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

# Genome survey via Sumdgeplot.
# Dependencies: Smudgeplot, KMC
Sumdgeplot=function(hist=hist, # kmer histogram from kmc.
                    kmcdb=kmcdb, # kmc database.
                    out_prefix=out_prefix
){
  L=system(paste("smudgeplot.py cutoff ",hist," L",sep=""),wait=TRUE,intern=TRUE)
  U=system(paste("smudgeplot.py cutoff ",hist," U",sep=""),wait=TRUE,intern=TRUE)
  
  cmd=paste("kmc_tools transform",
            kmcdb,
            paste("-ci",L," ","-cx",U," dump",sep=""),
            "-s",paste(kmcdb,"_L",L,"_U",U,".dump",sep=""),
            sep=" ")
  print(cmd);system(cmd)
  
  cmd=paste("smudgeplot.py hetkmers -o",out_prefix,
            paste(kmcdb,"_L",L,"_U",U,".dump",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("smudgeplot.py plot","-o",out_prefix,
            paste(out_prefix,"_coverages.tsv",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# kmergenie
# Dependencies: kmergenie
kmergenie=function(fq1=fq1, # Input fq1
                   fq2=fq2, # Input fq2. "none" if single-end.
                   out_prefix=out_prefix,
                   threads=threads){
  threads=as.character(threads)
  
  cmd=paste("echo",fq1,">",paste(out_prefix,".FILE",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  if (fq2!="none"){
    cmd=paste("echo",fq2,">>",paste(out_prefix,".FILE",sep=""),sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("kmergenie",
            paste(out_prefix,".FILE",sep=""),
            "-t",threads,
            "-o",out_prefix)
  print(cmd);system(cmd,wait=TRUE)
  
  return(out_prefix)
}
