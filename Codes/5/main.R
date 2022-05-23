# Genome annotation.
# Dependencies:
#    Softwares: RepeatModeler & RepeatMasker (Singularity container), 
#               GenomeThreader, 
#               Augustus, 
#               Hisat2, Stringtie, SAMtools,
#               gff3toolkit (python 3+),Cufflinks
#    R packages:

# RepeatModeler & RepeatMasker: Genome mask
GenomeMask=function(fna=fna,# Fasta file of genome.
                    out_prefix=out_prefix,
                    Threads=Threads
){
  Threads=as.character(Threads)
  
  # RepeatModeler & RepeatMasker installed in Singularity container
  path="singularity run /home/c/c-liu/Softwares/dfam-tetools-latest.sif"
  
  # RepeatModeler: de novo identification of repeats.
  cmd=paste(path,"BuildDatabase","-name",
            paste(out_prefix,"_RepeatModeler.db",sep=""),
            "-engine","ncbi",fna,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste(path,"RepeatModeler","-database",
            paste(out_prefix,"_RepeatModeler.db",sep=""),
            "-engine","ncbi","-pa",Threads,"-LTRStruct",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # RepeatMasker: Mask repeats find by RepeatModeler
  cmd=paste(path,"RepeatMasker","-lib",
            paste(out_prefix,"_RepeatModeler.db-families.fa",sep=""),
            "-pa",Threads,"-gff",fna)
  print(cmd);system(cmd,wait=TRUE)
  
  # RepeatMasker: Mask repeats by homology search
  cmd=paste(path,"RepeatMasker","-pa",Threads,"-gff",
            paste(fna,".masked",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Remove temporaries
  cmd=paste("rm",paste(out_prefix,"_RepeatModeler.db*",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("rm","-r","RM_*",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("rm",paste(fna,".masked.cat.gz",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("rm",paste(fna,".cat.gz",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  print("Output:")
  print(paste(fna,".masked.masked",sep=""))
  print(paste(fna,".out",sep=""))
  print(paste(fna,".out.gff",sep=""))
  print(paste(fna,".masked.out",sep=""))
  print(paste(fna,".masked.out.gff",sep=""))
  
  return(paste(fna,".masked.masked"))
}

# GenomeThreader: Homology-based gene prediction
GenomeThreader=function(fna=fna, # Fasta genome, masked
                        faas=faas, # A vector of paths to fasta files of homology proteins
                        Threads=Threads
){
  Threads=as.character(Threads)
  
  # Merge protein sequences for homology search
  for (faa in faas){
    cmd=paste("cat",faa,">> seed.faa",sep=" ")
    system(cmd,wait=TRUE)
  }
  
  # Gene prediction
  # gff3 output
  cmd=paste("gth","-genomic",fna,"-protein","seed.faa","-gff3out","-intermediate","-o",
            paste(fna,"_GenomeThreader.gff3"),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Remove temporaries
  cmd=paste("rm",paste(fna,".dna.*",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("rm",paste(fna,".protein.*",sep=""),sep=" ")
  cmd=paste("rm","seed.faa")
  print(cmd);system(cmd,wait=TRUE)
  
  print("Output:")
  print(paste(fna,"_GenomeThreader.gff3"))
  
  return(paste(fna,"_GenomeThreader.gff3"))
}

# Hisat2; Build hisat2 index of genome.
Hisat2Build = function(fna=fna, # FASTA of reference genome
                       index_prefix=index_prefix){
  cmd = paste("hisat2-build",fna,index_prefix,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(0)
}

# Hisat2: Map reads to reference genome. 
# SAMtools: Compress SAM to BAM and sort BAM.
Hisat = function(fq1=fq1,fq2=fq2, # Input fq files. Make fq2="None" if single-end.
                 index=index, # Basename of Hisat2 index of reference genome.
                 out_prefix=out_prefix, # Prefix of output BAM file.
                 threads=threads){
  threads = as.character(threads)
  
  bam_filename=paste(out_prefix,".bam",sep="")
  
  if (fq2!="None"){ # pair-end
    cmd = paste("hisat2","--dta","-x",index,"-p",threads,"-1",fq1,"-2",fq2,"|",
                "samtools","view","-@",threads,"-bS","|",
                "samtools","sort","-@",threads,"-o",bam_filename,sep=" ")
  }else{ # single pair
    cmd = paste("hisat2","--dta","-x",index,"-p",threads,"-U",fq1,"|",
                "samtools","view","-@",threads,"-bS","|",
                "samtools","sort","-@",threads,"-o",bam_filename,sep=" ")
  }
  
  print(cmd);system(cmd,wait=TRUE)
  
  return(bam_filename)
}

# SAMtools: Merge BAM.
MergeBAM=function(BAM_dir=BAM_dir,# directory containing BAM files
                  Threads=Threads){
  Threads=as.character(Threads)
  
  cmd=paste("samtools","merge","-@",Threads,"Merged.bam","'ls *bam'",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(0)
}

# StringTie: Assemble transcripts.
StringTie = function(input_bam=input_bam, # Input BAM.
                     threads=threads){
  threads = as.character(threads)
  
  output_prefix = gsub(pattern="\\.bam",replacement="",x=input_bam)
  gtf_filename = paste(output_prefix,".gtf",sep="")
  
  cmd = paste("stringtie","-p", threads,"-o", gtf_filename,input_bam,sep=" ")
  print(cmd);system(cmd,wait = T)
  
  return(gtf_filename)
}

# TransDecoder: Gene prediction from transcripts
TransDecoder=function(gtf=gtf, # gtf from StringTie
                      fna=fna){
  # From https://github.com/TransDecoder/TransDecoder
  path="/home/c/c-liu/Softwares/TransDecoder/util/"
  
  cmd=paste(path,"gtf_genome_to_cdna_fasta.pl",gtf,fna,">","transcripts.fna",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste(path,"gtf_to_alignment_gff3.pl",gtf,">","temp.gff3",sep=" ")
  print(cmd);system(cmd=TRUE)
  
  cmd=paste("TransDecoder.LongOrfs","-t",
            paste(fna,"_transcripts.fna",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste(path,"cdna_alignment_orf_to_genome_orf.pl","transcripts.fna.transdecoder.gff3","temp.gff3",
            paste(fna,"_transcripts.fna",sep=""),
            ">",paste(fna,"_TransDecoder.gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Remove temporary files
  cmd=paste("rm","transcripts.fna.transdecoder.gff3",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("rm","temp.gff3",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  print("Output:")
  print(paste(fna,"_TransDecoder.gff3"))
  paste(fna,"_transcripts.fna",sep="")
  
  return(paste(fna,"_TransDecoder.gff3"))
}

# gff3toolkit: Merge gff3 files
MergeGff3=function(gff3_1=gff3_1,gff3_2=gff3_2,fna=fna,out=out){
  cmd=paste("gff3_merge","-g1",gff3_1,"-g2",gff3_2,"-f",fna,"-og",out,"-r merged_report.txt",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(out)
}

# Cufflinks: Get transcripts from gff, and translate to protein. 
gtf2faa = function(gff=gff,fna=fna,faa=faa){
  cmd = paste("gffread",gff,"-g",fna,"-y",faa,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(0)
}

# Augustus: de novo gene prediction
Augustus=function(fna=fna, # Fasta genome, masked
                  species=species,
                  faa=faa, # Protein sequences used for model trainning.
                          # Set "none" to skip model trainning and use pre-trained model provided by Augustus
                          # augustus --species=help\
                  Threads=Threads
){
  Threads=as.character(Threads)
  # Path to tool scripts
  path="/home/c/c-liu/Softwares/Augustus/scripts/"
  
  if (faa!="none"){
    cmd=paste("perl",
              paste(path,"autoAug.pl",sep=""),
              paste("--genome=",fna,sep=""),
              paste("--species=",species,sep=""),
              paste("--trainingset=",faa,sep=""),
              paste("--cpus=",Threads,sep=""),
              "-v",sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }else{
    cmd=paste("augustus",
              paste("--species=",species,sep=""),
              fna,
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
}