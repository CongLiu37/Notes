# Genome annotation

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
  
  # RepeatMasker: Mask repeats in de novo repeat library from RepeatModeler
  cmd=paste(path,
            "RepeatMasker",
            "-xsmall", # soft masking
            "-lib",paste(out_prefix,"_RepeatModeler.db-families.fa",sep=""),
            "-pa",Threads,
            "-gff",
            fna,
            "-dir",out_dir)
  print(cmd);system(cmd,wait=TRUE)
  
  # RepeatMasker: Mask repeats in RepBase
  cmd=paste(path,
            "RepeatMasker",
            "-xsmall", # soft masking
            "-pa",Threads,
            "-gff",
            paste(fna,".masked",sep=""),
            "-dir",out_dir,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # Remove temporaries
  system(paste("rm",fna,sep=" "),wait=TRUE)
  
  setwd(pwd_begin)
  
  return(paste(out_dir,"/",fna,".masked.masked",sep=""))
}

# GenomeThreader: protein-genome spliced alignments
# Dependencies: GenomeThreader
gth=function(fna=fna, # masked genome
             faa=faa, # protein sequences
             out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  cmd=paste("gth",
            "-genomic",fna,
            "-protein",faa,
            "-gff3out -intermediate",
            "-o","protein_alignment.gff3",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd_begin)
  return(out_dir)
}

# # Prepare AUGUSTUS training set with proteins
# # Dependencies: GenomeThreader, scripts from AUGUSTUS, stringr (R)
# trainAUGUSTUS_protein=function(fna=fna, # masked genome
#                                faa=faa, # reference proteins
#                                out_dir=out_dir){
#   library(stringr)
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
#   wd_begin=getwd();setwd(out_dir)
#   
#   cmd=paste("gth",
#             "-genomic",fna,
#             "-protein",faa,
#             "-gff3out -skipalignmentout -paralogs",
#             "-o","GenomeThreader.gff3",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("gth2gtf.pl", # AUGUSTUS Remove alternative splicement
#             "GenomeThreader.gff3",
#             "training_AUGUSTUS.gtf",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd="computeFlankingRegion.pl training_AUGUSTUS.gtf" # AUGUSTUS
#   print(cmd)
#   flanking_DNA=system(cmd,wait=TRUE,intern=TRUE)[4]
#   flanking_DNA=str_extract(flanking_DNA,": [0-9]*")
#   flanking_DNA=sub(": ","",flanking_DNA)
#   cmd=paste("gff2gbSmallDNA.pl training_AUGUSTUS.gtf",
#             fna,flanking_DNA,"training_AUGUSTUS.gb",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   setwd(wd_begin)
#   return(0)
# }

# Hisat2: Map short RNA reads to reference genome. 
# SAMtools: Compress SAM to BAM, sort BAM, index BAM.
Hisat = function(fq1=fq1,fq2=fq2, # Input fq files. Make fq2="None" if single-end.
                                  # Comma-separated list.
                 fna=fna, # genome (soft masked for training AUGUSTUS)
                 index=index, # Basename of Hisat2 index of reference genome.
                 out_prefix=out_prefix, # Prefix of output BAM file.
                 threads=threads){
  threads = as.character(threads)
  
  cmd = paste("hisat2-build",
              fna,
              index,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
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
  
  cmd=paste("awk '{if ($2==\"Cufflinks\")  $2=\"StringTie\";print$0}'",
            paste(output_prefix,"_tmp.gff3",sep=""),">",
            paste(output_prefix,".gff3",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system(paste("rm"," ",output_prefix,"_tmp.gff3",sep=""),wait=TRUE)
  return(paste(output_prefix,".gtf",sep=""))
}

# TransDecoder: Find Coding Regions within Transcripts
TransDecoder=function(gtf=gtf, # gtf from StringTie
                      out_dir=out_dir,
                      fna=fna){
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
  print(cmd);system(cmd=TRUE)
  
  cmd=paste("TransDecoder.LongOrfs",
            "-t","transcripts.fna",
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
  
  setwd(wd)
  return("TransDecoder.gff3")
}

# Prepare AUGUSTUS training set with BAM (short RNA read to genome)
# Dependencies: GeneMark-ET, scripts from BRAKER, AUGUSTUS
trainAUGUSTUS_shortRNA=function(genome=genome, # soft masked
                                bam=bam,
                                out_dir=out_dir,
                                threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  threads=as.character(threads)
  
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
  
  cmd=paste("gmes_petap.pl","--verbose", # GeneMark-ET
            "--sequence",genome,         # gmhmm.mod
            "--ET","introns.f.gff",
            "--cores",threads,
            "--soft_mask 1000",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("filterGenemark.pl", # BRAKER
            "--genemark=genemark.gtf",
            "--introns=introns.f.gff",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("mv genemark.f.good.gtf training_AUGUSTUS.gtf")
  
  cmd="computeFlankingRegion.pl training_AUGUSTUS.gtf" # AUGUSTUS
  print(cmd)
  flanking_DNA=system(cmd,wait=TRUE,intern=TRUE)[4]
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

# GeneMarker-hmm for gene prediction
gmhmm=function(genome=genome,
               model_file=model_file,
               out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  cmd=paste("gmhmme3",
            "-m",model_file,
            "-f","gff3",
            genome,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
}

# Train glimmerhmm with gff3 from stringtie+EvidenceModeler and gene prediction
# Dependencies: glimmerhmm, scripts from EvidenceModeler
glimmerhmm=function(transcript_align.gff3=transcript_align.gff3,
                    genome=genome,
                    out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd_begin=getwd();setwd(out_dir)
  
  transcript_align=readLines(transcript_align.gff3)
  training=sapply(transcript_align,
                  function(st){
                    st=unlist(strsplit(st," "))
                    if (length(st)==2){
                      return("")
                    }else{
                      seq=st[1];strand=st[7]
                      if (strand=="+"){start=st[4];end=st[5]}
                      if (strand=="-"){start=st[5];end=st[4]}
                      return(paste(seq,start,end,sep=" "))
                    }
                  })
  writeLines(training,"exon_file")
  
  cmd=paste("trainGlimmerHMM",genome,"exon_file",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("glimmerhmm",genome,".","-g",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  #cmd=paste("glimmerHMM_to_GFF3.pl","glimmerHMM.output",sep=" ") # EvidenceModeler
  #print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd_begin)
}

# Remove redundant gene structures in training data of AUGUSTUS
# Dependencies: blast+, scripts from AUGUSTUS, stringr (R)
non_redundant=function(training_AUGUSTUS.gb=training_AUGUSTUS.gb, # GenBank
                       training_AUGUSTUS.gtf=training_AUGUSTUS.gtf, # gtf
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
            "bonafide.f.gtf",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("gtf2aa.pl", # AUGUSTUS
            genome,
            "bonafide.f.gtf",
            "prot.aa",
            sep=" ")
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

# AUGUSTUS training and gene prediction.
# Dependencies: AUGUSTUS, parallel (R)
augustus=function(fna=fna, # genome
                  species=species, # Species name for the trained model
                  training.gb=training.gb, # non-redundant training genes in GenBank format
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
  
  setwd(wd)
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

# EvidenceModeler
evm=function(protein_alignments.gff3="none", # GFF3 for protein-genome spliced alignments
             gene_predictions.gff3="none", # GFF3 from ab initio gene prediction
             transcript_alignments.gff3="none", # GFF3 for transcript-genome spliced alignments
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
