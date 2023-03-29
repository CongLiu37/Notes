# Dependencies: 
# R, Python
# Softwares: minimap2, SAMtools, SprayNPray, diamond, megan (blast2rma & rma2info) and its database, seqkit, QUAST, BUSCO
# Python modules: numpy, pandas, sklearn
# R packages: reticulate, stringr, ggplot2, ggExtra

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
  
  system(paste("cp",fna,".",sep=" "),wait=TRUE)
  fna=unlist(strsplit(fna,"/"))[length(unlist(strsplit(fna,"/")))]
  fna=paste(out_dir,"/",fna,sep="")
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
  system(paste("rm",fna,sep=" "),wait=TRUE)
  setwd(wd)
  return(paste(out_dir,"/",out_basename,sep=""))
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

# Quast: Quality of assembly.
# Dependencies: QUAST
Quast=function(fna=fna, # FASTA of genome assembly.
               out_dir=out_dir,
               threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  
  if (!file.exists(paste(out_dir,"/report.tsv",sep=""))){
    if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
    cmd=paste("quast.py --min-contig 0",
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
  
  f1=paste(out_dir,"/",Out_prefix,"/short_summary.specific.",Lineage,".",Out_prefix,".txt",sep="")
  f2=paste(out_dir,"/",Out_prefix,"/short_summary.specific..",Out_prefix,".txt",sep="")
  if (!file.exists(f1) & !file.exists(f2)){
    wd_begin=getwd();setwd(out_dir)
    if (file.exists(Out_prefix)){system(paste("rm"," -r ",Out_prefix,sep=""))}
    cmd=paste("busco","--force",
              "--in",fna,
              "--lineage_dataset",Lineage,
              "--out",Out_prefix,
              "--mode",Mode,
              "--cpu",Threads,
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    setwd(wd_begin)
  }
  
  if (file.exists(f1)){f=f1}
  if (file.exists(f2)){f=f2}
  re=readLines(f)[8]
  re=sub("\t","",re);re=sub("\t   ","",re)
  system(paste("cat",f,sep=" "))
  return(re)
}
