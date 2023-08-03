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
  
  f=paste(out_dir,"/",Out_prefix,"/short_summary*.txt",sep="")
  if (system(paste("if [ -e ",f," ]; then echo TRUE; fi",sep=""),intern=TRUE)!="TRUE"){
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
  f=system(paste("ls",f,sep=" "),intern=TRUE)
  re=readLines(f)
  re=re[grepl("C:.*$",re)]
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
               exon_per_transcript=nrow(exon)/nrow(mRNA),
               average_exon_length=mean(exon$V3),
               transcript_number=nrow(mRNA),
               transcript_per_gene=nrow(mRNA)/nrow(gene),
               gene_less_200bp=nrow(subset(gene,gene$V3<200)))
  write.table(d,out,sep="\t",row.names=FALSE,quote=FALSE)
  setwd(wd)
  system(paste("rm -r",tmp_dir,sep=" "))
}
