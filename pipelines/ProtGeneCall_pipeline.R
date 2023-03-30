arg=commandArgs(trailingOnly = TRUE)

main=arg[1]
conf=arg[2]

source(main)
source(conf)

out_dir=sub("/$","",out_dir)
save_dir=sub("/$","",save_dir)
threads=as.character(threads)
if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}

genome_file=unlist(strsplit(genome,"/"))
genome_file=genome_file[length(genome_file)]
masked_genome=paste(save_dir,"/masking/",label,"/",genome_file,".masked",sep="")

############################################
# Repeat masking
if (1 %in% step){
  Stamp1=paste(out_dir,"/masking/",label,"/",label,"_masking.finished",sep="")
  Stamp2=paste(save_dir,"/masking/",label,"/",label,"_masking.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 1: Repeat masking FINISHED")
  }else{
    print("Step 1: Repeat masking START")
    if (!file.exists(paste(out_dir,"/masking/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/masking/",sep=""),wait=TRUE)
    }
    repeatmodeler(fna=genome,
                  out_dir=paste(out_dir,"/masking/",label,sep=""),
                  out_prefix=label,
                  Threads=threads)
    repeatmasker(fna=genome,# Fasta file of genome.
                 out_dir=paste(out_dir,"/masking/",label,sep=""),
                 RepeatLib.fa=paste(out_dir,"/masking/",label,"/",label,"_RepeatModeler.db-families.fa",sep=""),
                 Threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# miniprot
if (2 %in% step){
  Stamp1=paste(out_dir,"/miniprot/",label,"/",label,"_miniprot.finished",sep="")
  Stamp2=paste(save_dir,"/miniprot/",label,"/",label,"_miniprot.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 2: miniprot FINISHED")
  }else{
    print("Step 2: miniprot START")
  
    if (!file.exists(paste(out_dir,"/miniprot/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/miniprot/",sep=""),wait=TRUE)
    }
    miniprot(fna=masked_genome,
             faa=ref_protein, # comma-list
             threads=threads,
             out_dir=paste(out_dir,"/miniprot/",label,sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# Hisat
if (3 %in% step){
  Stamp1=paste(out_dir,"/Hisat/",label,"/",label,"_Hisat.finished",sep="")
  Stamp2=paste(save_dir,"/Hisat/",label,"/",label,"_Hisat.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 3: Hisat FINISHED")
  }else{
    print("Step 3: Hisat START")
    if (!file.exists(paste(out_dir,"/Hisat/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/Hisat/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/Hisat/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/Hisat/",label,sep=""))
    }
    RNA_f=unlist(strsplit(RNA1,","))
    RNA_r=unlist(strsplit(RNA2,","))
    for (i in 1:length(RNA_f)){
      n=unlist(strsplit(RNA_f[i],"/"));n=n[length(n)]
      if (!file.exists(paste(out_dir,"/Hisat/",label,"/",label,"_",n,".bam",sep="")) 
          & !file.exists(paste(save_dir,"/Hisat/",label,"/",label,"_",n,".bam",sep=""))){
        Hisat(fq1=RNA_f[i],fq2=RNA_r[i], # Input fq files. Make fq2="None" if single-end.
            # Comma-separated list.
              fna=masked_genome, # genome (soft masked for training AUGUSTUS)
              index=paste(out_dir,"/Hisat/",label,"/",label,sep=""), # Basename of Hisat2 index of reference genome.
              out_prefix=paste(out_dir,"/Hisat/",label,"/",label,"_",n,sep=""), # Prefix of output BAM file.
              threads=threads)
      }
    }
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# StringTie
if (4 %in% step){
  Stamp1=paste(out_dir,"/StringTie/",label,"/",label,"_StringTie.finished",sep="")
  Stamp2=paste(save_dir,"/StringTie/",label,"/",label,"_StringTie.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 4: StringTie FINISHED")
  }else{
    print("Step 4: StringTie START")
    if (!file.exists(paste(out_dir,"/MergeBAM/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/MergeBAM/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/MergeBAM/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/MergeBAM/",label,sep=""))
    }
    BAM_lst=system(paste("ls"," ",save_dir,"/Hisat/",label,"/",label,"_*.bam",sep=""),intern = TRUE)
    MergeBAM(BAMs=paste(BAM_lst,collapse=" "),# SPACE-separated list of bam files.
             out_prefix=paste(out_dir,"/MergeBAM/",label,"/",label,"_total",sep=""),
             Threads=threads)
    if (!file.exists(paste(out_dir,"/StringTie/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/StringTie/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/StringTie/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/StringTie/",label,sep=""))
    }
    StringTie(input_bam=paste(out_dir,"/MergeBAM/",label,"/",label,"_total.bam",sep=""), # Input BAM.
              output_prefix=paste(out_dir,"/StringTie/",label,"/",label,sep=""),
              threads=threads)
    system(paste("rm"," ",out_dir,"/MergeBAM/",label,"/",label,"_total.bam",sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# TransDecoder
if (5 %in% step){
  Stamp1=paste(out_dir,"/TransDecoder/",label,"/",label,"_TransDecoder.finished",sep="")
  Stamp2=paste(save_dir,"/TransDecoder/",label,"/",label,"_TransDecoder.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 5: TransDecoder FINISHED")
  }else{
    print("Step 5: TransDecoder START")
    if (!file.exists(paste(out_dir,"/TransDecoder/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/TransDecoder/",sep=""))
    }
    TransDecoder(gtf=paste(save_dir,"/StringTie/",label,"/",label,".gtf",sep=""), # gtf from StringTie
                 fna=masked_genome,
                 out_dir=paste(out_dir,"/TransDecoder/",label,sep=""),
                 dmdb=UniRef_dmnd, # DIAMOND protein db, uniprot
                 pfam=PfamA, # Pfam-A.hmm
                 threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# augustus
if (6 %in% step){
  Stamp1=paste(out_dir,"/augustus/",label,"/",label,"_augustus.finished",sep="")
  Stamp2=paste(save_dir,"/augustus/",label,"/",label,"_augustus.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 6: augustus FINISHED")
  }else{
    print("Step 6: augustus START")
    if (!file.exists(paste(out_dir,"/augustus/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/augustus/",sep=""))
    }
    augustus_trainset(gene_models.gff3=paste(save_dir,"/TransDecoder/",label,"/TransDecoder.gff3",sep=""),
                      fna=masked_genome,
                      out_dir=paste(out_dir,"/augustus/",label,sep=""))
    non_redundant(training_AUGUSTUS.gb=paste(out_dir,"/augustus/",label,"/training_AUGUSTUS.gb",sep=""), # GenBank 
                  training_AUGUSTUS.gtf=paste(out_dir,"/augustus/",label,"/training_AUGUSTUS.gtf",sep=""), # gtf 
                  genome=masked_genome,
                  threads=threads,
                  out_dir=paste(out_dir,"/augustus/",label,sep=""))
    augustus(fna=masked_genome, # genome
             species=paste(label,"_TransDecoder",sep=""), # Species name for the trained model
             training.gb=paste(out_dir,"/augustus/",label,"/training_AUGUSTUS_nr.gb",sep=""), # non-redundant training genes in GenBank format
             # training_AUGUSTUS_nr.gb
             out_dir=paste(out_dir,"/augustus/",label,sep=""),
             threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# braker
if (7 %in% step){
  Stamp1=paste(out_dir,"/braker/",label,"/",label,"_braker.finished",sep="")
  Stamp2=paste(save_dir,"/braker/",label,"/",label,"_braker.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 7: braker FINISHED")
  }else{
    print("Step 7: braker START")
    if (!file.exists(paste(out_dir,"/braker/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/braker/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/MergeBAM/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/MergeBAM/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/MergeBAM/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/MergeBAM/",label,sep=""))
    }
    BAM_lst=system(paste("ls"," ",save_dir,"/Hisat/",label,"/",label,"_*.bam",sep=""),intern = TRUE)
    MergeBAM(BAMs=paste(BAM_lst,collapse=" "),# SPACE-separated list of bam files.
             out_prefix=paste(out_dir,"/MergeBAM/",label,"/",label,"_total",sep=""),
             Threads=threads)
    braker(genome=masked_genome,
           bam=paste(out_dir,"/MergeBAM/",label,"/",label,"_total.bam",sep=""), # comma-list
           ref_proteins="none",
           species=label,
           tsebra.conf="none", # TSEBRA/config/default.cfg 
           out_dir=paste(out_dir,"/braker/",label,sep=""),
           threads=threads)
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# galba
if (8 %in% step){
  Stamp1=paste(out_dir,"/galba/",label,"/",label,"_galba.finished",sep="")
  Stamp2=paste(save_dir,"/galba/",label,"/",label,"_galba.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 8: galba FINISHED")
  }else{
    print("Step 8: galba START")
    if (!file.exists(paste(out_dir,"/galba/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/galba/",sep=""))
    }
    galba(genome=masked_genome,
          proteins=galba_protein,
          species=paste(label,"_galba",sep=""),
          threads=threads,
          out_dir=paste(out_dir,"/galba/",label,sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# EvidenceModeler
if (9 %in% step){
  Stamp1=paste(out_dir,"/evm/",label,"/",label,"_evm.finished",sep="")
  Stamp2=paste(save_dir,"/evm/",label,"/",label,"_evm.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 9: evm FINISHED")
  }else{
    print("Step 9: evm START")
    if (!file.exists(paste(out_dir,"/evm/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/evm/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/evm/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/evm/",label,sep=""))
    }
  
    weight=paste(out_dir,"/evm/",label,"/evm_weights.txt",sep="")
    cat("TRANSCRIPT\tTRANSCRIPT\t30\n",file=weight,append=TRUE)
    TRANSCRIPT=paste(save_dir,"/StringTie/",label,"/",label,".gff3",sep="")
  
    cat("OTHER_PREDICTION\ttransdecoder\t15\n",file=weight,append=TRUE)
    transdecoder=paste(save_dir,"/TransDecoder/",label,"/sorted.gff3",sep="")
  
    cat("PROTEIN\tPROTEIN\t10\n",file=weight,append=TRUE)
    PROTEIN=paste(save_dir,"/miniprot/",label,"/protein_align.gff3",sep="")
  
    cat("ABINITIO_PREDICTION\tBRAKER\t8\n",file=weight,append=TRUE)
    braker=paste(save_dir,"/braker/",label,"/sorted.gff3",sep="")
  
    cat("ABINITIO_PREDICTION\tGALBA\t6\n",file=weight,append=TRUE)
    galba=paste(save_dir,"/galba/",label,"/sorted.gff3",sep="")
  
    cat("ABINITIO_PREDICTION\tAugustus\t1\n",file=weight,append=TRUE)
    Augustus=paste(save_dir,"/augustus/",label,"/sorted.gff3",sep="")
  
    cmd=paste("cat",
              transdecoder,Augustus,braker,galba,
              ">",paste(out_dir,"/evm/",label,"/ab_initio.gff3",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  
    evm(protein_alignments.gff3=PROTEIN, # Absolute path. GFF3 for protein-genome spliced alignments
        gene_predictions.gff3=paste(out_dir,"/evm/",label,"/ab_initio.gff3",sep=""), # Absolute path. GFF3 from ab initio gene prediction
        transcript_alignments.gff3=TRANSCRIPT, # Absolute path. GFF3 for transcript-genome spliced alignments
        evm_weights=paste(out_dir,"/evm/",label,"/evm_weights.txt",sep=""), # Absolute path
        genome=masked_genome,
        out_dir=paste(out_dir,"/evm/",label,sep=""),
        threads=threads)
    system(paste("rm"," ",out_dir,"/evm/",label,"/ab_initio.gff3",sep=""))
    system(paste("touch",Stamp1))
  }
}
############################################
# Quality of evm, filter, and quality of filtered
if (10 %in% step){
  Stamp1=paste(out_dir,"/quality_filter/",label,"/",label,"_quality_filter.finished",sep="")
  Stamp2=paste(save_dir,"/quality_filter/",label,"/",label,"_quality_filter.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 10: quality_filter FINISHED")
  }else{
    print("Step 10: quality_filter START")
    
    if (!file.exists(paste(out_dir,"/quality_evm/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/quality_evm/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/quality_evm/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/quality_evm/",label,sep=""))
    }
    gffread(gff=paste(save_dir,"/evm/",label,"/sorted.gff3",sep=""),
            fna=masked_genome, # genome
            exons="none",cds="none",
            pep=paste(out_dir,"/quality_evm/",label,"/pep.faa",sep=""),
            tmp_dir=paste(out_dir,"/quality_evm/",label,"/gffread_tmp",sep=""))
    BUSCO(fna=paste(out_dir,"/quality_evm/",label,"/pep.faa",sep=""), # Fasta file of nucleotide or protein.
                 # Be consistent with Mode.
          Mode="proteins", # genome/proteins/transcriptome
          Lineage=BUSCO_db, # Lineage dataset, e.g. insecta_odb10
                 # Available datasets: https://busco-data.ezlab.org/v5/data/lineages/
                 # BUSCO will download lineage dataset automatically.
          Out_prefix=paste(label,"_busco",sep=""), # Give the analysis run a recognisable short name.
                 # Output folders and files will be labelled with this name.
                 # Cannot be path
          out_dir=paste(out_dir,"/quality_evm/",label,sep=""),
          Threads=threads)
    system(paste("rm ",out_dir,"/quality_evm/",label,"/pep.faa",sep=""))
    gene_model_stat(gff3=paste(save_dir,"/evm/",label,"/sorted.gff3",sep=""),
                    out=paste(out_dir,"/quality_evm/",label,"/gene_model_stat.tsv",sep=""),
                    tmp_dir=paste(out_dir,"/quality_evm/",label,"/tmp",sep=""))

    ## Filter hypothetical gene supported by only 1 predictor
    if (!file.exists(paste(out_dir,"/filter_evm/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/filter_evm/",sep=""))
    }
    filterEvm(gene_evidence.tsv=paste(save_dir,"/evm/",label,"/gene_evidence.tsv",sep=""),
              genome.fna=masked_genome,
              evm.gff3=paste(save_dir,"/evm/",label,"/sorted.gff3",sep=""),
              out_dir=paste(out_dir,"/filter_evm/",label,sep=""))
  
    if (!file.exists(paste(out_dir,"/quality_filter/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/quality_filter/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/quality_filter/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/quality_filter/",label,sep=""))
    }
    gffread(gff=paste(out_dir,"/filter_evm/",label,"/filtered.gff3",sep=""),
            fna=masked_genome, # genome
            exons="none",cds="none",
            pep=paste(out_dir,"/quality_filter/",label,"/pep.faa",sep=""),
            tmp_dir=paste(out_dir,"/quality_filter/",label,"/gffread_tmp",sep=""))
    BUSCO(fna=paste(out_dir,"/quality_filter/",label,"/pep.faa",sep=""), # Fasta file of nucleotide or protein.
          # Be consistent with Mode.
          Mode="proteins", # genome/proteins/transcriptome
          Lineage=BUSCO_db, # Lineage dataset, e.g. insecta_odb10
          # Available datasets: https://busco-data.ezlab.org/v5/data/lineages/
          # BUSCO will download lineage dataset automatically.
          Out_prefix=paste(label,"_busco",sep=""), # Give the analysis run a recognisable short name.
          # Output folders and files will be labelled with this name.
          # Cannot be path
          out_dir=paste(out_dir,"/quality_filter/",label,sep=""),
          Threads=threads)
    system(paste("rm ",out_dir,"/quality_filter/",label,"/pep.faa",sep=""))
    gene_model_stat(gff3=paste(out_dir,"/filter_evm/",label,"/filtered.gff3",sep=""),
                    out=paste(out_dir,"/quality_filter/",label,"/gene_model_stat.tsv",sep=""),
                    tmp_dir=paste(out_dir,"/quality_filter/",label,"/tmp",sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# pasa
if (11 %in% step){
  Stamp1=paste(out_dir,"/pasa/",label,"/",label,"_pasa.finished",sep="")
  Stamp2=paste(save_dir,"/pasa/",label,"/",label,"_pasa.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 11: pasa FINISHED")
  }else{
    print("Step 11: pasa START")
    if (!file.exists(paste(out_dir,"/pasa/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/pasa/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/pasa/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/pasa/",label,sep=""))
    }
    pasa(genome=masked_genome,
         transcripts=paste(save_dir,"/TransDecoder/",label,"/transcripts.fna",sep=""),
         original.gff3=paste(save_dir,"/filter_evm/",label,"/filtered.gff3",sep=""), # only protein-coding gene
        pasa.alignAssembly.conf=pasa.alignAssembly.conf, # Path to PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt
        pasa.annotationCompare.conf=pasa.annotationCompare.conf, # Path to PASApipeline/pasa_conf/pasa.annotationCompare.Template.txt
        threads=threads,
        out_dir=paste(out_dir,"/pasa/",label,"/1",sep=""))
    pasa(genome=masked_genome,
         transcripts=paste(save_dir,"/TransDecoder/",label,"/transcripts.fna",sep=""),
         original.gff3=paste(out_dir,"/pasa/",label,"/1/sorted.gff3",sep=""), # only protein-coding gene
        pasa.alignAssembly.conf=pasa.alignAssembly.conf, # Path to PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt
        pasa.annotationCompare.conf=pasa.annotationCompare.conf, # Path to PASApipeline/pasa_conf/pasa.annotationCompare.Template.txt
        threads=threads,
        out_dir=paste(out_dir,"/pasa/",label,"/2",sep=""))
    system(paste("mv"," ",out_dir,"/pasa/",label,"/2/sorted.gff3",
                      " ",out_dir,"/pasa/",label,"/",sep=""))
    system(paste("touch",Stamp1,sep=" "))
  }
}
############################################
# filter PASA
if (12 %in% step){
  Stamp1=paste(out_dir,"/quality_pasa/",label,"/",label,"_quality_pasa.finished",sep="")
  Stamp2=paste(save_dir,"/quality_pasa/",label,"/",label,"_quality_pasa.finished",sep="")
  if (file.exists(Stamp1) | file.exists(Stamp2)){
    print("Step 12: quality_pasa FINISHED")
  }else{
    print("Step 12: pasa START")
    if (!file.exists(paste(out_dir,"/filter_pasa/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/filter_pasa/",sep=""))
    }
    filterPasa(genome.fna=masked_genome,
               pasa.gff3=paste(save_dir,"/pasa/",label,"/sorted.gff3",sep=""),
               out_dir=paste(out_dir,"/filter_pasa/",label,sep=""))
    if (!file.exists(paste(out_dir,"/protein_1/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/protein_1/",sep=""))
    }
    pasa_more(species=label,
              genome=masked_genome,
              gene_evidence.tsv=paste(save_dir,"/evm/",label,"/gene_evidence.tsv",sep=""),
              PASA.gff3=paste(out_dir,"/filter_pasa/",label,"/filtered.gff3",sep=""),
              out_dir=paste(out_dir,"/protein_1/",label,sep=""))
    if (!file.exists(paste(out_dir,"/quality_pasa/",sep=""))){
      system(paste("mkdir"," ",out_dir,"/quality_pasa/",sep=""))
    }
    if (!file.exists(paste(out_dir,"/quality_pasa/",label,sep=""))){
      system(paste("mkdir"," ",out_dir,"/quality_pasa/",label,sep=""))
    }
    BUSCO(fna=paste(out_dir,"/protein_1/",label,"/",label,"_proteins_rep.faa",sep=""), # Fasta file of nucleotide or protein.
          # Be consistent with Mode.
          Mode="proteins", # genome/proteins/transcriptome
          Lineage=BUSCO_db, # Lineage dataset, e.g. insecta_odb10
          # Available datasets: https://busco-data.ezlab.org/v5/data/lineages/
          # BUSCO will download lineage dataset automatically.
          Out_prefix=paste(label,"_busco",sep=""), # Give the analysis run a recognisable short name.
          # Output folders and files will be labelled with this name.
          # Cannot be path
          out_dir=paste(out_dir,"/quality_pasa/",label,sep=""),
          Threads=threads)
    gene_model_stat(gff3=paste(out_dir,"/protein_1/",label,"/",label,"_genes.gff3",sep=""),
                    out=paste(out_dir,"/quality_pasa/",label,"/gene_model_stat.tsv",sep=""),
                    tmp_dir=paste(out_dir,"/quality_pasa/",label,"/tmp",sep=""))
  }
}
