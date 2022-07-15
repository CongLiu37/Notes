# Genome annotation.

# Simplify sequence ID in genome
SimplifyID=function(fna=fna,
                    common_pattern=common_pattern){
  f1=function(i){
    new_ID=paste(common_pattern,as.character(i),sep="")
    return(new_ID)
  }
  library(ape)
  data=read.dna(fna,format="fasta",as.character=TRUE)
  old_IDs=names(data)
  new_IDs=sapply(1:length(old_IDs),f1)
  
  names(data)=new_IDs
  write.dna(data,paste(fna,"_SimpleIDs",sep=""),"fasta",colsep="")
  
  ID_convert=data.frame(old_IDs,new_IDs)
  write.table(ID_convert,paste(fna,"_IDconvert.tsv",sep=""),
              sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
}

# RepeatModeler & RepeatMasker: Genome mask
# Masked genome: <file name of fna>.masked.masked
# Dependencies: Singularity, RepeatMasker, RepeatModeler
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
  
  #####################################################################
  # RepeatModeler & RepeatMasker installed in Singularity container
  path="singularity run /home/c/c-liu/Softwares/dfam-tetools-latest.sif"
  #####################################################################
  
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

# Protein based gene prediction with GenomeThreader (gff3)
# gtf and GenBank format for AUGUSTUS training
# Dependencies: GenomeThreader, scripts from AUGUSTUS
GenePrediction_protein=function(fna=fna, # masked genome
                                faa=faa,
                                out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  cmd=paste("gth",
            "-genomic",fna,
            "-protein",faa,
            "-gff3out -skipalignmentout -paralogs",
            "-o","bonafide.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("gth2gtf.pl", # AUGUSTUS Remove alternative splicement
            "bonafide.gff3",
            "training_AUGUSTUS.gtf",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="computeFlankingRegion.pl training_AUGUSTUS.gtf" # AUGUSTUS
  print(cmd);library(stringr)
  flanking_DNA=system(cmd,wait=TRUE,intern=TRUE)[4]
  flanking_DNA=str_extract(flanking_DNA,": [0-9]*")
  flanking_DNA=sub(": ","",flanking_DNA)
  cmd=paste("gff2gbSmallDNA.pl training_AUGUSTUS.gtf",
            fna,flanking_DNA,"training_AUGUSTUS.gb",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
  return(0)
}

# Remove redundant gene structures in training data of AUGUSTUS
# Dependencies: blast+
non_redundant=function(training_AUGUSTUS.gb=training_AUGUSTUS.gb, # GenBank
                       training_AUGUSTUS.gtf=training_AUGUSTUS.gtf, # bonafide
                       genome=genome,
                       out_dir=out_dir){
  library(stringr)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  # Extract gene ID
  GenBank=readLines(training_AUGUSTUS.gb)
  gb=GenBank[grepl("/gene=",GenBank)]
  gb=sub("                     /gene=\"","",gb)
  gb=sub("\"","",gb)
  gb=unique(gb)
  write.table(gb,"traingenes.lst",sep="\t",col.names=FALSE,row.names=FALSE,quote=TRUE)
  
  # Gene ID to locus ID
  LOCUS=GenBank[grepl("LOCUS",GenBank)]
  LOCUS=sub("LOCUS       ","",LOCUS)
  LOCUS=sub("   [0-9]* bp  DNA","",LOCUS)
  write.table(data.frame(gb,LOCUS),"loci.lst",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  # remove genes share >80% similarity at protein level
  cmd=paste("grep",
            "-f","traingenes.lst",
            "-F",training_AUGUSTUS.gtf,
            ">",
            "bonafide.f.gtf",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("gtf2aa.pl", # AUGUSTUS
            genome,
            "bonafide.f.gtf",
            "prot.aa",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("aa2nonred.pl","prot.aa","prot.nr.aa",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # non-redundant gene ID and locus ID
  cmd="grep '>' prot.nr.aa | perl -pe 's/>//' > nonred.lst"
  print(cmd);system(cmd,wait=TRUE)
  cmd="grep -f nonred.lst loci.lst | cut -f2 > nonred.loci.lst"
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

augustus=function(fna=fna,
                  species=species,
                  training.gb=training.gb,
                  out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("autoAug.pl",
            paste("--workingdir=",out_dir,sep=""),
            paste("--species=",species,sep=""),
            paste("--genome=",fna,sep=""),
            paste("--trainingset=",training.gb,sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

# # Training AUGUSTUS
# AUGUSTUS_training=function(training.gb=training.gb,
#                            species=species,
#                            out_dir=out_dir){
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
#   wd=getwd();setwd(out_dir)
# 
#   cmd=paste("new_species.pl"," ",
#             "--species=",species,
#             sep="")
#   print(cmd);system(cmd,wait=TRUE)
# 
#   cmd=paste("etraining"," ",
#             "--species=",species," ",
#             training.gb,
#             " &> ",
#             "training.out",
#             sep="")
#   print(cmd);system(cmd,wait=TRUE)
# 
#   error=system("grep -c 'Variable stopCodonExcludedFromCDS set right' training.out",intern=TRUE)
#   total=system(paste("grep -c 'LOCUS'",training.gb,sep=" "),intern=TRUE)
#   percent=as.numeric(error)/as.numeric(total)
#   print(paste("stopCodonExcludedFromCDS error percent:",as.character(percent),sep=" "))
#   if (percent>0.5){
#     cmd=paste("sed",
#               "-i",
#               "'s/stopCodonExcludedFromCDS false/stopCodonExcludedFromCDS true/g'",
#               paste(species,"_parameters.cfg",sep=""),
#               sep=" ")
#     print(cmd);system(cmd,wait=TRUE)
#   }
#   
#   cmd=paste("etraining --species=",species," ",training.gb,
#             " 2>&1 | grep 'in sequence' | perl -pe 's/.*n sequence (\S+):.*/$1/' | sort -u > bad.lst",
#             sep="")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("filterGenes.pl",
#             "bad.lst",
#             training.gb,
#             ">",
#             "training.f.gb",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   setwd(wd)
# }

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
  cmd=paste("gffread","-O",
            gff,
            "-g",fna,
            sep=" ")
  if (exons!="none"){cmd=paste(cmd,"-w",exons,sep=" ")}
  if (cds!="none"){cmd=paste(cmd,"-x",cds,sep=" ")}
  if (pep!="none"){cmd=paste(cmd,"-y",pep,sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
}

# OrthoFinder
# By default OrthoFinder creates a results directory called ‘OrthoFinder’ inside the 
# input proteome directory and puts the results here.
Orthofinder=function(in_dir=in_dir, # Input proteome directory. 
                                    # One file per species with extension '.faa'
                     threads=threads){
  threads=as.character(threads)
  
  cmd=paste("orthofinder",
            "-f",in_dir,
            "-t",threads,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}
