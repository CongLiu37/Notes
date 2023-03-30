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

# Get redundant training set (AUGUSTUS) from gene models (gff3, CDS/exon features)
# Dependencies: stringr (R), scripts from AUGUSTUS 
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
  #system("rm -r errors")
  system("rm genome_header.map")
  #system("rm pygustus_hints.out")
  #system("rm pygustus_hints.py")
  system("rm what-to-cite.txt")
  
  system("rm augustus.hints.aa")
  system("rm augustus.hints.agat.log")
  system("rm augustus.hints.codingseq")
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
  #system("rm -r ./braker_rna/errors")
  system("rm -r ./braker_rna/GeneMark-ET")
  system("rm ./braker_rna/genemark_hintsfile.gff")
  system("rm ./braker_rna/genome_header.map")
  system("rm ./braker_prot/what-to-cite.txt")
  
  system("rm ./braker_prot/augustus.hints.aa")
  system("rm ./braker_prot/augustus.hints.codingseq")
  #system("rm -r ./braker_prot/errors")
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

gffread=function(gff=gff,
                 fna=fna, # genome
                 exons="none",
                 cds="none",
                 pep="none",
                 tmp_dir=tmp_dir){
  tmp_dir=sub("/$","",tmp_dir)
  if (!file.exists(tmp_dir)){system(paste("mkdir",tmp_dir,sep=" "))}
  system(paste("cp",fna,paste(tmp_dir,"/genome.fna",sep=""),sep=" "))
  cmd=paste("gffread","-O",
            gff,
            "-g",paste(tmp_dir,"/genome.fna",sep=""),
            sep=" ")
  if (exons!="none"){cmd=paste(cmd,"-w",exons,sep=" ")}
  if (cds!="none"){cmd=paste(cmd,"-x",cds,sep=" ")}
  if (pep!="none"){cmd=paste(cmd,"-y",pep,sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
  system(paste("rm -r",tmp_dir,sep=" "))
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

# gene model statistics
gene_model_stat=function(gff3=gff3,
                         out=out,
                         tmp_dir=tmp_dir){
  if (!file.exists(tmp_dir)){system(paste("mkdir",tmp_dir,sep=" "))}
  wd=getwd();setwd(tmp_dir)
  
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3==\"gene\") print $4,$5}'",
            gff3,"> gene_coordinate.tsv",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  gene=read.table("gene_coordinate.tsv",sep="\t",header=FALSE,quote="")
  gene$V3=gene$V2-gene$V1+1
  
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3==\"exon\") print $4,$5}'",
            gff3,"> exon_coordinate.tsv",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  exon=read.table("exon_coordinate.tsv",sep="\t",header=FALSE,quote="")
  exon$V3=exon$V2-exon$V1+1
  
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3==\"mRNA\") print $4,$5}'",
            gff3,"> mRNA_coordinate.tsv",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  mRNA=read.table("mRNA_coordinate.tsv",sep="\t",header=FALSE,quote="")
  mRNA$V3=mRNA$V2-mRNA$V1+1
  
  d=data.frame(gff3=gff3,
               gene_number=nrow(gene),
               min_gene_length=min(gene$V3),
               max_gene_length=max(gene$V3),
               average_gene_length=mean(gene$V3),
               exon_number=nrow(exon),
               exon_per_gene=nrow(exon)/nrow(gene),
               average_exon_length=mean(exon$V3),
               transcript_number=nrow(mRNA),
               transcript_per_gene=nrow(mRNA)/nrow(gene),
               gene_less_200bp=nrow(subset(gene,gene$V3<200)))
  write.table(d,out,sep="\t",row.names=FALSE,quote=FALSE)
  setwd(wd)
  system(paste("rm -r",tmp_dir,sep=" "))
}

# Update gff3 by PASA
# Dependencies: pasa,sqlite3,,agat,scripts from PASAPipeline
pasa=function(genome=genome,
              transcripts=transcripts,
              original.gff3=original.gff3, # only protein-coding gene
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
  system(paste("cp",genome,paste(species,"_genome.fna",sep=""),sep=" "))
  genome=paste(species,"_genome.fna",sep="")
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
  
  system(paste("rm",paste(species,"_transcripts.fna",sep=""),sep=" "),wait=TRUE)
  system(paste("rm",paste(species,"_cds.fna",sep=""),sep=" "),wait=TRUE)
  system(paste("rm",paste(species,"_proteins.faa",sep=""),sep=" "),wait=TRUE)
  system("rm *.fai")
  setwd(wd)
}
