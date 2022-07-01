# Genome annotation.

# RepeatModeler & RepeatMasker: Genome mask
# Masked genome: <file name of fna>.masked.masked
GenomeMask=function(fna=fna,# Fasta file of genome.
                    out_dir=out_dir,
                    out_prefix=out_prefix,
                    Threads=Threads){
  pwd_begin=getwd()
  setwd(out_dir)
  
  Threads=as.character(Threads)
  out_dir=sub("/$","",out_dir)
  
  system(paste("cp",fna,out_dir,sep=" "),wait=TRUE)
  fna_name=unlist(strsplit(fna,"/"));fna_name=fna_name[length(fna_name)]
  fna=paste(out_dir,"/",fna_name,sep="")
  
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
  cmd=paste(path,"RepeatMasker",
            "-xsmall", # soft masking
            "-lib",paste(out_prefix,"_RepeatModeler.db-families.fa",sep=""),
            "-pa",Threads,
            "-gff",
            fna,
            "-dir",out_dir)
  print(cmd);system(cmd,wait=TRUE)
  # RepeatMasker: Mask repeats by homology search
  cmd=paste(path,"RepeatMasker",
            "-xsmall", # soft masking
            "-pa",Threads,
            "-gff",
            paste(fna,".masked",sep=""),
            "-dir",out_dir,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Remove temporaries
  cmd=paste("rm",fna,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(pwd_begin)
  
  return(paste(out_dir,"/",fna,".masked.masked",sep=""))
}

# Protein based gene prediction with Braker, GenomeThreader
GenePrediction_protein=function(fna=fna, # masked genome
                                faa=faa,
                                out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)

  cmd=paste("startAlign.pl",
            paste("--genome=",fna,sep=""),
            paste("--prot=",faa,sep=""),
            "--prg=gth",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("gth2gtf.pl",
            "align_gth/gth.concat.aln",
            "bonafide.gtf",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # cmd="computeFlankingRegion.pl bonafide.gtf"
  # print(cmd);system(cmd,wait=TRUE)
  # cmd="gff2gbSmallDNA.pl bonafide.gtf genome.fa 10000 bonafide.gb"
  # print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
}

# GenomeThreader: Homology-based exon-intron prediction
GenomeThreader=function(fna=fna, # Fasta genome, masked
                        faa=faa, # proteins
                        out_prefix=out_prefix,
                        Threads=Threads){
  Threads=as.character(Threads)

  # Gene prediction
  # gff3 output
  cmd=paste("gth",
            "-genomic",fna,
            "-protein",faa,
            "-gff3out","-skipalignmentout","-paralogs",
            "-gcmincoverage 80 -prseedlength 20 -prminmatchlen 20 -prhdist 2",
            "-o",paste(out_prefix,".gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Remove temporaries
  cmd=paste("rm",paste(out_prefix,"_seed.faa",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  print("Output:")
  print(paste(out_prefix,".gff3"))
  
  return(paste(out_prefix,".gff3"))
}

# ProtHint: Homology-based prediction of hints in the form of introns, start and stop codons.
prothint=function(fna=fna, # Fasta genome, masked
                  faa=faa, # proteins
                  out_dir=out_dir,
                  threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  
  cmd=paste("prothint.py",
            "--workdir",out_dir,
            "--threads",threads,
            fna,
            faa,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Augustus
augustus=function(fna=fna,
                  species=species,
                  training_gff=training_gff,
                  hints_gff=hints_gff,
                  out_dir=out_dir,
                  threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  
  cmd=paste("autoAug.pl",
            paste("--workingdir=",out_dir,sep=""),
            paste("--cpus=",threads,sep=""),
            paste("--species=",species,sep=""),
            paste("--genome=",fna,sep=""),
            paste("--trainingset=",training_gff,sep=""),
            paste("--hints=",hints_gff,sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Hisat2; Build hisat2 index of genome.
Hisat2Build = function(fna=fna, # FASTA of reference genome
                       index_prefix=index_prefix){
  cmd = paste("hisat2-build",fna,index_prefix,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(0)
}

# Hisat2: Map reads to reference genome. 
# SAMtools: Compress SAM to BAM, sort BAM, index BAM.
Hisat = function(fq1=fq1,fq2=fq2, # Input fq files. Make fq2="None" if single-end.
                                  # Comma-separated list.
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
  
  cmd=paste("samtools","index",bam_filename,"-@",threads,sep=" ")
  print(cmd);system(cmd)
  
  return(bam_filename)
}

# SAMtools: Merge BAM.
MergeBAM=function(BAMs=BAMs,# space-separated list of bam files.
                  out_prefix=out_prefix,
                  Threads=Threads){
  Threads=as.character(Threads)
  
  cmd=paste("samtools","merge",
            "-@",Threads,
            paste(out_prefix,".bam"),
            BAMs,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(0)
}

# StringTie: Assemble transcripts.
StringTie = function(input_bam=input_bam, # Input BAM.
                     output_prefix=out_prefix,
                     threads=threads){
  threads = as.character(threads)
  gtf_filename = paste(output_prefix,".gtf",sep="")
  
  cmd = paste("stringtie","-p", threads,"-o", gtf_filename,input_bam,sep=" ")
  print(cmd);system(cmd,wait = T)
  
  return(gtf_filename)
}

# TransDecoder: Find Coding Regions within Transcripts
TransDecoder=function(gtf=gtf, # gtf from StringTie
                      out_prefix=out_prefix,
                      fna=fna){
  
  cmd=paste("gtf_genome_to_cdna_fasta.pl",
            gtf,
            fna,
            ">",
            paste(out_prefix,"_","transcripts.fna",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("gtf_to_alignment_gff3.pl",
            gtf,
            ">",
            paste(out_prefix,"_","temp.gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd=TRUE)
  
  cmd=paste("TransDecoder.LongOrfs",
            "-t",paste(out_prefix,"_","transcripts.fna",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("cdna_alignment_orf_to_genome_orf.pl",
            paste(out_prefix,"_","transcripts.fna.transdecoder.gff3",sep=""),
            paste(out_prefix,"_","temp.gff3",sep=""),
            paste(out_prefix,"_","transcripts.fna",sep=""),
            ">",
            paste(out_prefix,"_TransDecoder.gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(paste(out_prefix,"_TransDecoder.gff3",sep=""))
}

# gff3toolkit: Merge gff3 files
MergeGff3=function(gff3_1=gff3_1,gff3_2=gff3_2,fna=fna,out=out){
  cmd=paste("gff3_merge","-g1",gff3_1,"-g2",gff3_2,"-f",fna,"-og",out,"-r merged_report.txt",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(out)
}

# gffread: Extract sequences from gff.
gffread=function(gff=gff,
                 fna=fna,
                 exons="none",
                 cds="none",
                 pep="none"){
  cmd=paste("gffread",
            gff,
            "-g",fna,
            sep=" ")
  if (exons!="none"){cmd=paste(cmd,"-w",exons,sep=" ")}
  if (cds!="none"){cmd=paste(cmd,"-x",cds,sep=" ")}
  if (pep!="none"){cmd=paste(cmd,"-y",pep,sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
}

# OrthoFinder