# Curation and evaluation of genome and genome annotation (kmer-based genome survey, decontamination, quality of genome/genome annotation)

#######################################################################
# kmer-based genome survey
#######################################################################
# Build kmc database
# Dependencies: KMC
kmc=function(fq1=fq1, # Input fq1
             fq2=fq2, # Input fq2. "none" if single-end.
             input_format=input_format, # "-fa" for FASTA, "-fq" for FASTQ
             # "-fm" for multiple FASTA
             # "-fbam" for BAM
             kmer=31,
             out_basename=out_basename, # kmc database
             out_dir=out_dir,
             threads=threads){
  threads=as.character(threads)
  out_dir=gsub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  kmer=as.character(kmer)
  wd=getwd()
  setwd(out_dir)
  
  cmd=paste("echo",fq1,">",
            paste(out_dir,"/",out_basename,".FILE",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  if (fq2!="none"){
    cmd=paste("echo",fq2,">>",
              paste(out_dir,"/",out_basename,".FILE",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("kmc",
            input_format,
            paste("-k",kmer,sep=""),
            paste("-t",threads,sep=""),
            "-m64","-ci1 -cs10000",
            paste("@",out_dir,"/",out_basename,".FILE",sep=""),
            out_basename, # kmc database
            ".",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",paste(out_dir,"/",out_basename,".FILE",sep=""),sep=" "))
  
  setwd(wd)
}

# Merge kmc database
merge_kmc=function(kmcdb1=kmcdb1,kmcdb2=kmcdb2,
                   out_prefix=out_prefix){
  cmd=paste("kmc_tools simple",
            kmcdb1,kmcdb2,
            "union",out_prefix,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# kmc database to kmer hist
# Dependencies: kmc
kmc_hist=function(kmcdb=kmcdb,
                  out_prefix=out_prefix){
  cmd=paste("kmc_tools transform",
            kmcdb,
            "histogram",
            paste(out_prefix,".hist",sep=""),
            "-cx10000",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Genome survey via Genomescope.
# Dependencies: GenomeScope, stringr (R)
Genomescope=function(hist=hist, # kmer histogram from kmc.
                     ploidy="none", # 1-6
                     out_dir=out_dir,
                     out_basename=out_basename,
                     kmer=31){
  kmer=as.character(kmer)
  out_dir=sub("/$","",out_dir)
  
  if (!file.exists(paste(out_dir,"/",out_basename,"_summary.txt",sep=""))){
    if (ploidy!="none"){
      out_basename=paste(out_basename,"_p",as.character(ploidy),sep="")
    }
    cmd=paste("genomescope.R",
              "-i",hist,
              "-o",out_dir,
              "-n",out_basename,
              "-k",kmer,
              sep=" ")
    if (ploidy!="none"){
      cmd=paste(cmd,"-p",as.character(ploidy))
    }
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

# Genome survey via Smudgeplot.
# Dependencies: Smudgeplot, KMC
Smudgeplot=function(hist=hist, # kmer histogram from kmc.
                    kmcdb=kmcdb, # kmc database.
                    L="none",U="none", # integers estimated from kmer profile (Genomescope)
                    # L should be chosen so it cuts out majority of the error kmers
                    # U should be so it does contain all the k-mers that are at least in four genomic copies (but you might want to increase this number for genomes that might have higher ploidy levels).      
                    out_prefix=out_prefix){
  L=as.character(L);U=as.character(U)
  
  if (L=="none"){
    L=system(paste("smudgeplot.py cutoff ",hist," L",sep=""),wait=TRUE,intern=TRUE)
  }
  if (U=="none"){
    U=system(paste("smudgeplot.py cutoff ",hist," U",sep=""),wait=TRUE,intern=TRUE)
  }
  
  cmd=paste("kmc_tools transform",
            kmcdb,
            paste("-ci",L," ","-cx",U," ","dump",sep=""),
            "-s",paste(kmcdb,"_L",L,"_U",U,".dump",sep=""),
            sep=" ")
  print(cmd);system(cmd)
  
  cmd=paste("smudgeplot.py hetkmers -o",paste(out_prefix,"_L",L,"_U",U,sep=""),
            paste(kmcdb,"_L",L,"_U",U,".dump",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("smudgeplot.py plot","-o",paste(out_prefix,"_L",L,"_U",U,sep=""),
            paste(out_prefix,"_L",L,"_U",U,"_coverages.tsv",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# kmergenie: analyze NGS reads
# Dependencies: kmergenie
kmergenie=function(reads=reads, # vector of paths
                   out_prefix=out_prefix,
                   threads=threads){
  threads=as.character(threads)
  
  if (file.exists(paste(out_prefix,".FILE",sep=""))){
    system(paste("rm"," ",out_prefix,".FILE",sep=""),wait=TRUE)
  }
  
  for (i in reads){
    cmd=paste("echo",i,">>",paste(out_prefix,".FILE",sep=""),sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("kmergenie",
            paste(out_prefix,".FILE",sep=""),
            "-t",threads,
            "-o",out_prefix)
  print(cmd);system(cmd,wait=TRUE)
  
  return(out_prefix)
}

# KAT: evaluate contamination in assembly by GC & kmer coverage
kat_sect=function(assembly=assembly,
                  kmer=27,
                  threads=threads,
                  out_prefix=out_prefix){
  threads=as.character(threads)
  kmer=as.character(kmer)
  
  cmd=paste("kat sect",
            "-t",threads,
            "-m",kmer,
            "-o",out_prefix,
            assembly,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# KAT: kmer based assembly evaluation
kat_comp=function(assembly=assembly,
                  reads=reads, # space separated list
                  kmer=27,
                  threads=threads,
                  out_prefix=out_prefix){
  threads=as.character(threads)
  kmer=as.character(kmer)
  cmd=paste("kat comp",
            "-t",threads,
            "-m",kmer,
            "-o",out_prefix,
            reads,
            assembly,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Merqury: kmer based assembly evaluation
# Dependencies: merqury, meryl
merqury=function(reads=reads, # comma list
                 assembly=assembly,
                 kmer=21, # best_k.sh from merqury
                 out_dir=out_dir,
                 out_basename=out_basename){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  reads=unlist(strsplit(reads,","))
  for (i in 1:length(reads)){
    cmd=paste("meryl","count",
              paste("k=",as.character(kmer)),
              reads[i],
              "output",
              paste(out_basename,as.character(i),".meryl",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("meryl","union-sum",
            "output",paste(out_basename,".meryl",sep=""),
            paste(out_basename,"*.meryl",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("merqury.sh",
            paste(out_basename,".meryl",sep=""),
            assembly,out_basename,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)

  setwd(wd)
}

#######################################################################
# Decontamination 
#######################################################################
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

# Map reads to reference
# Dependencies: Minimap2, SAMtools
minimap2=function(long_reads=long_reads, # space-separated list for PE
                  lr_type=lr_type, # long read type. 
                  # "map-pb" for PacBio
                  # "map-hifi" for HiFi
                  # "map-ont" for ONT reads.
                  # "sr" for NGS
                  # "asm5" for accurate reads diverging <5% to assembly
                  chunk_mem="5G",
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
            "-m", chunk_mem,
            "-@",threads,
            "-o",paste(out_prefix,".bam",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("samtools","index",paste(out_prefix,".bam",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
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

# GC content & length of each sequence
# Dependencies: seqkit
GC_length=function(fna=fna,out=out,threads=threads){
  cmd=paste("printf","'rname\tlength\tGC.content'",">",out)
  system(cmd,wait=TRUE)
  
  cmd=paste("seqkit","fx2tab",
            "--name --only-id --gc --length",
            "--threads",as.character(threads),
            fna,">>",out,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# DepthSizer estimating genome size by BUSCO single copy genes
# Dependencies: python2, DepthSizer, SAMtools
depthsizer=function(genome=genome,
                    bam=bam,
                    BUSCO.tsv=BUSCO.tsv,
                    path2depthsizer.py=path2depthsizer.py,
                    out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("python2",path2depthsizer.py,
            paste("seqin=",genome,sep=""),
            paste("bam=",bam,sep=""),
            paste("busco=",BUSCO.tsvsep=""),
            sep=" ")
  
  setwd(wd)
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
#View(read.table("~/Desktop/PhD/Results/termite_pca/AI/possible_wolbachia_scaffold_barrnap.fna2nt.blast",
#           header=FALSE,sep="\t",quote=""))
# Scaffolds-protein search by diamond.
# Compute taxonomy at major ranks by Megan.
# Diamond and Megan in long-read mode.
# Dependencies: DIAMOND, MEGAN
## fna 2 daa 2 tsv
Diamond_Megan=function(fna, # fna. DNA assembly.
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
            "--outfmt 100",
            "-p",threads,
            "-d",ref_diamond,
            "-q",fna,
            "--long-reads",
            "-f 100",
            "--out",paste(blast_dir,"/",out_basename,".blast.daa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("daa-meganizer",
            "-i",paste(blast_dir,"/",out_basename,".blast.daa",sep=""),
            "-lg",
            "-mdb",ref_megan,
            "-t",threads,
            "-ram readCount",
            "-supp 0",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("daa2info",
            "-i",paste(blast_dir,"/",out_basename,".blast.daa",sep=""),
            "-o",paste(assignment_dir,"/",out_basename,".tsv",sep=""),
            "-r2c Taxonomy",
            "-n true",
            "-p true",
            "-r true",
            "-mro true",
            "-u false",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  return(paste(assignment_dir,"/",out_basename,".tsv",sep=""))
}

#######################################################################
# Quality of genome/genome annotation
#######################################################################
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
  
  f=paste(out_dir,"/",Out_prefix,"/short_summary*.txt",sep="")
  if (system(paste("if [ -e ",f," ]; then echo TRUE; else echo FALSE; fi",sep=""),intern=TRUE)!="TRUE"){
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

# Convert genbank to gff
# Dependencies: biocode
genbank2gff=function(genbank=genbank,
                     gff=gff){
  cmd=paste("convert_genbank_to_gff3.py",
            "-i",genbank,
            "-o",gff,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(gff)
}

# Extract features (gene+exon+CDS) by list of gene IDs
SplitGFF=function(in.gff=in.gff,
                  gene.list=gene.list,
                  out.gff=out.gff){
  cmd=paste("agat_sp_filter_feature_from_keep_list.pl",
            "--gff",in.gff,
            "--keep_list",gene.list,
            "--out",out.gff)
  print(cmd);system(cmd,wait=TRUE)
}

# # gFACs: filter gene models (gff3)
# gfacs=function(evm.gff3=evm.gff3, # EvidenceModeler.gff3
#                genome.fna=genome.fna,
#                out_dir=out_dir){
#   out_dir=sub("/$","",out_dir)
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
#   wd=getwd();setwd(out_dir)
#   system(paste("cp",genome.fna,"./genome.fa",sep=" "),wait=TRUE)
#   
#   cmd=paste("gFACs.pl",
#             "-f EVM_1.1.1_gff3 -p gFACs --create-gtf --statistics --fasta ./genome.fa",
#             "--statistics-at-every-step --splice-table --rem-all-incompletes",
#             "--allowed-inframe-stop-codons 0",
#             "--min-exon-size 20 --min-intron-size 20 --min-CDS-size 150",
#             #"--rem-genes-without-start-and-stop-codon",
#             #"--rem-genes-without-start-codon --rem-genes-without-stop-codon",
#             "-O ./",evm.gff3,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("awk -F '\t' -v OFS='\t' '{if ($3==\"gene\") print $9}' gFACs_out.gtf | sed 's/ID=//' > geneID.lst",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("agat_sp_filter_feature_from_keep_list.pl",
#             "--gff",evm.gff3,
#             "--keep_list geneID.lst --out gfacs.gff3",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd="agat_convert_sp_gxf2gxf.pl --gff gfacs.gff3 -o sorted.gff3"
#   print(cmd);system(cmd,wait=TRUE)
#   
#   system("rm geneID.lst")
#   system("rm gfacs.agat.log")
#   system("rm gFACs_gene_table.txt")
#   system("rm gfacs.gff3")
#   system("rm gFACs_out.gtf")
#   system("rm gfacs_report.txt")
#   system("rm gFACs_statistics.txt")
#   system("rm tempy.txt")
#   system("rm ./genome.fa")
#   system("rm ./genome.fa.idx")
#   setwd(wd)
# }

# # Blast protein against itself and extract proteins IDs hit too many times
# # Dependencies: blastp
# proteinWithManyHits=function(proteins.faa=proteins.faa,
#                              out_prefix=out_prefix,
#                              threads=threads){
#   threads=as.character(threads)
#   system(paste("cp"," ",proteins.faa," ",out_prefix,".faa",sep=""))
#   proteins=paste(out_prefix,".faa",sep="")
#   
#   cmd=paste("makeblastdb",
#             "-in",proteins.faa1,
#             "-dbtype prot",
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   cmd=paste("blastp",
#             "-num_threads",threads,
#             "-db",proteins.faa1,
#             "-query",proteins.faa2,
#             "-outfmt 6",
#             "-evalue 1e-5",
#             "-num_alignments 20",
#             "-out",paste(out_basename,".blast",sep=""),
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
# }











# Compare two gff3
# Dependencies: aegean
parseval=function(reference.gff3=reference.gff3,
                  prediction.gff3=prediction.gff3,
                  out_prefix=out_prefix){
  cmd=paste("parseval",
            "--outfile",paste(out_prefix,".txt",sep=""),
            reference.gff3,
            prediction.gff3,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("grep","-A10","'CDS structure comparison'",
            paste(out_prefix,".txt",sep=""),
            ">",
            paste(out_prefix,"_CDS.txt",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("sed -i 's/     |    Annotation edit distance:      //'",
            paste(out_prefix,"_CDS.txt",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("grep","-A10","'Exon structure comparison'",
            paste(out_prefix,".txt",sep=""),
            ">",
            paste(out_prefix,"_Exon.txt",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("sed -i 's/     |    Annotation edit distance:      //'",
            paste(out_prefix,"_Exon.txt",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("grep","-A10","'UTR structure comparison'",
            paste(out_prefix,".txt",sep=""),
            ">",
            paste(out_prefix,"_UTR.txt",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("sed -i 's/     |    Annotation edit distance:      //'",
            paste(out_prefix,"_UTR.txt",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("grep","-A6","'Nucleotide-level comparison'",
            paste(out_prefix,".txt",sep=""),
            ">",
            paste(out_prefix,"_Nucleotide.txt",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("sed -i 's/     |    Annotation edit distance:      //'",
            paste(out_prefix,"_Nucleotide.txt",sep=""),sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Evaluate predicted gene structures by comparing with reference
# Dependencies: GenomeTools
gene_structure_evaluation=function(reference.gff3=reference.gff3,
                                   prediction.gff3=prediction.gff3,
                                   out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("gt gff3 -sort -tidy",reference.gff3,"> sorted_reference.gff3",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("gt gff3 -sort -tidy",prediction.gff3,"> sorted_prediction.gff3",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("gt eval",
            "sorted_reference.gff3","sorted_prediction.gff3",
            ">","evaluation.txt",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  system(paste("echo",reference.gff3,">> evaluation.txt",sep=" "))
  system(paste("echo",prediction.gff3,">> evaluation.txt",sep=" "))
  
  setwd(wd)
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
  plot_taxa=names(t[which(t>0.0001*sum(t))])
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
  #df=data.frame(cov=rnorm(1000,50,10),GC.content=rnorm(1000,38,5))
  #title=""
  p1=ggplot(df,aes(x=cov,y=GC.content))+
    geom_point(size=0.01)+
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


# #
# gaeval=function(protein_alignment.gff3="none",
#                 transcript_alignment.gff3="none",
#                 gene_structure.gff3=gene_structure.gff3,
#                 out_prefix=out_prefix){
#   if (transcript_alignment.gff3!="none"){
#     cmd=paste("cat",transcript_alignment.gff3,"|",
#               "awk -v OFS='\t' '{$2=\"nucleotide_match\";print $0}' >>",
#               paste(out_prefix,"_alignment.gff3",sep=""),
#               sep=" ")
#     print(cmd);system(cmd,wait=TRUE)
#   }
#   if (protein_alignment.gff3!="none"){
#     cmd=paste("cat",protein_alignment.gff3,"|",
#               "awk -v OFS='\t' '{$2=\"nucleotide_match\";print $0}' >>",
#               paste(out_prefix,"_alignment.gff3",sep=""),
#               sep=" ")
#     print(cmd);system(cmd,wait=TRUE)
#   }
#   
#   cmd=paste("gaeval",
#             "--alpha 0.5",
#             "--beta 0.5",
#             "--gamma 0",
#             "--epsilon 0",
#             paste(out_prefix,"_alignment.gff3",sep=""),
#             gene_structure.gff3,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
# }
# # Compare gff
# # Dependencies: biocode
# compare_gff=function(ref.gff=ref.gff,
#                      predicted.gff=predicted.gff,
#                      feature=feature, # Exon/CDS
#                      out_dir=out_dir){
#   if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
#   
#   cmd=paste("compare_gene_structures.py",
#             "-a1",ref.gff,
#             "-a2",predicted.gff,
#             "-f",feature,
#             "-o",out_dir,
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
# }
##fna 2 daa 2 rma
# Diamond_Megan=function(fna, # fna. Input DNA sequences.
#                        out_basename=out_basename,
#                        blast_dir=blast_dir, # Directory for diamond output.
#                        rma_dir=rma_dir, # Directory for rma output of megan.
#                        assignment_dir=assignment_dir, # Directory for taxonomy table.
#                        ref_diamond=ref_diamond, # Diamond database.
#                        ref_megan=ref_megan, # Megan database.
#                        threads=threads){
#   threads=as.character(threads)
#   blast_dir=sub("/$","",blast_dir)
#   rma_dir=sub("/$","",rma_dir)
#   assignment_dir=sub("/$","",assignment_dir)
#   
#   cmd=paste("diamond blastx",
#             "-p",threads,
#             "-d",ref_diamond,
#             "-q",fna,
#             "--long-reads",
#             "-f 100",
#             "--out",paste(blast_dir,"/",out_basename,".blast.daa",sep=""),
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("daa2rma",
#             "-i",paste(blast_dir,"/",out_basename,".blast.daa",sep=""),
#             "-o",paste(rma_dir,"/",out_basename,".rma",sep=""),
#             "--paired","false",
#             "-lg","true",
#             "-mdb",ref_megan,
#             "-t",threads,
#             "-ram","readCount",
#             "-supp","0",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("rma2info",
#             "-i",paste(rma_dir,"/",out_basename,".rma",sep=""),
#             "-o",paste(assignment_dir,"/",out_basename,".tsv",sep=""),
#             "-r2c Taxonomy",
#             "-n true",
#             "-p true",
#             "-r true",
#             "-mro","true",
#             "-u false",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   return(paste(assignment_dir,"/",out_basename,".tsv",sep=""))
# }
##fna 2 blast 2 rma
# Diamond_Megan=function(fna, # fna. Input DNA sequences.
#                        out_basename=out_basename,
#                        blast_dir=blast_dir, # Directory for diamond output.
#                        rma_dir=rma_dir, # Directory for rma output of megan.
#                        assignment_dir=assignment_dir, # Directory for taxonomy table.
#                        ref_diamond=ref_diamond, # Diamond database.
#                        ref_megan=ref_megan, # Megan database.
#                        threads=threads){
#   threads=as.character(threads)
#   blast_dir=sub("/$","",blast_dir)
#   rma_dir=sub("/$","",rma_dir)
#   assignment_dir=sub("/$","",assignment_dir)
#   
#   cmd=paste("diamond blastx",
#             "-p",threads,
#             "-d",ref_diamond,
#             "-q",fna,
#             "--long-reads",
#             "--out",paste(blast_dir,"/",out_basename,".blast",sep=""),
#             sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("blast2rma",
#             "-i",paste(blast_dir,"/",out_basename,".blast",sep=""),
#             "-o",paste(rma_dir,"/",out_basename,".rma",sep=""),
#             "-f","BlastTab",
#             "-bm","BlastX",
#             "--paired","false",
#             "-lg","true",
#             "-mdb",ref_megan,
#             "-t",threads,
#             "-ram","readCount",
#             "-supp","0",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   cmd=paste("rma2info",
#             "-i",paste(rma_dir,"/",out_basename,".rma",sep=""),
#             "-o",paste(assignment_dir,"/",out_basename,".tsv",sep=""),
#             "-r2c Taxonomy",
#             "-n true",
#             "-p true",
#             "-r true",
#             "-mro","true",
#             "-u false",sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   return(paste(assignment_dir,"/",out_basename,".tsv",sep=""))
# }