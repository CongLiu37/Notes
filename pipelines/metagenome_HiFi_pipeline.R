main=commandArgs(trailingOnly = TRUE)[1]
conf=commandArgs(trailingOnly = TRUE)[2]
source(main)
source(conf)

if (1 %in% Step){
  print("Step 1: quality_check")
  dire=paste(out_dir,"/quality_check/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            fastq,
            ">",
            paste(dire,"/",out_prefix,".stats",sep=""))
  print(cmd);system(cmd,wait=TRUE)
}

if (2 %in% Step){
  print("Step 2: rmHost")
  dire=paste(out_dir,"/rmHost/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  minimap2(long_reads=fastq, # space-separated list for PE
           lr_type="map-hifi", # long read type.
           # "map-pb" for PacBio
           # "map-hifi" for HiFi
           # "map-ont" for ONT reads.
           # "sr" for NGS
           # "asm5" for accurate reads diverging <5% to assembly
           assembly=host_genome.fna,
           out_prefix=paste(dire,"/",out_prefix,sep=""),
           threads=threads)
  
  Extract_fq(bam=paste(dire,"/",out_prefix,".bam",sep=""),
             paired=F, # logical. T for paired.
             mapped=F, # logical. T for mapped
             out_prefix=paste(dire,"/",out_prefix,"_rmHost",sep=""),
             reads_format="fq", # "fq" or "fa"
             threads=threads)
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            paste(dire,"/",out_prefix,"_rmHost.fq.gz",sep=""),
            ">",
            paste(dire,"/",out_prefix,"_rmHost.stats",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  file.remove(paste(dire,"/",out_prefix,".bam",sep=""))
  file.remove(paste(dire,"/",out_prefix,".bam.bai",sep=""))
}

if (3 %in% Step){
  print("Step 3: metaMDBG")
  dire=paste(out_dir,"/metaMDBG/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  metaMDBG(hifi.fq=paste(out_dir,"/rmHost/",out_prefix,"_rmHost.fq.gz",sep=""),
           nanopore.fq=NA,
           out_dir=paste(dire,"/",out_prefix,sep=""),
           threads=threads)
  cmd=paste("gzip -d ",dire,"/",out_prefix,"/contigs.fasta.gz",sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            paste(dire,"/",out_prefix,"/contigs.fasta",sep=""),
            ">",
            paste(dire,"/",out_prefix,"/contigs.fasta.stats",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  minimap2(long_reads=paste(out_dir,"/rmHost/",out_prefix,"_rmHost.fq.gz",sep=""), # space-separated list for PE
           lr_type="map-hifi", # long read type.
           # "map-pb" for PacBio
           # "map-hifi" for HiFi
           # "map-ont" for ONT reads.
           # "sr" for NGS
           # "asm5" for accurate reads diverging <5% to assembly
           assembly=paste(dire,"/",out_prefix,"/contigs.fasta",sep=""),
           out_prefix=paste(dire,"/",out_prefix,"/minimap",sep=""),
           threads=as.character(threads))
  
  Extract_fq(bam=paste(dire,"/",out_prefix,"/minimap.bam",sep=""),
             paired=F, # logical. T for paired.
             mapped=F, # logical. T for mapped
             out_prefix=paste(dire,"/",out_prefix,"/minimap_unmap",sep=""),
             reads_format="fq", # "fq" or "fa"
             threads=as.character(threads))
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            paste(dire,"/",out_prefix,"/minimap_unmap.fq.gz",sep=""),
            ">",
            paste(dire,"/",out_prefix,"/minimap_unmap.stats",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  file.remove(paste(dire,"/",out_prefix,"/minimap.bam",sep=""))
  file.remove(paste(dire,"/",out_prefix,"/minimap.bai",sep=""))
  file.remove(paste(dire,"/",out_prefix,"/minimap_unmap.fq.gz",sep=""))
}

if (4 %in% Step){
  print("Step 4: hifiasm_meta")
  dire=paste(out_dir,"/hifiasm_meta/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  if (!dir.exists(paste(dire,"/",out_prefix,sep=""))){dir.create(paste(dire,"/",out_prefix,sep=""))}
  hifiasm_meta(hifi.fq=paste(out_dir,"/rmHost/",out_prefix,"_rmHost.fq.gz",sep=""),
              out_prefix=paste(dire,"/",out_prefix,"/",out_prefix,sep=""),
              threads=threads)
  
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            paste(dire,"/",out_prefix,"/",out_prefix,".p_ctg.fna",sep=""),
            ">",
            paste(dire,"/",out_prefix,"/",out_prefix,".p_ctg.stats",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)

  minimap2(long_reads=paste(out_dir,"/rmHost/",out_prefix,"_rmHost.fq.gz",sep=""), # space-separated list for PE
           lr_type="map-hifi", # long read type.
           # "map-pb" for PacBio
           # "map-hifi" for HiFi
           # "map-ont" for ONT reads.
           # "sr" for NGS
           # "asm5" for accurate reads diverging <5% to assembly
           assembly=paste(dire,"/",out_prefix,"/",out_prefix,".p_ctg.fna",sep=""),
           out_prefix=paste(dire,"/",out_prefix,"/minimap",sep=""),
           threads=as.character(threads))

  Extract_fq(bam=paste(dire,"/",out_prefix,"/minimap.bam",sep=""),
             paired=F, # logical. T for paired.
             mapped=F, # logical. T for mapped
             out_prefix=paste(dire,"/",out_prefix,"/minimap_unmap",sep=""),
             reads_format="fq", # "fq" or "fa"
             threads=as.character(threads))
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            paste(dire,"/",out_prefix,"/minimap_unmap.fq.gz",sep=""),
            ">",
            paste(dire,"/",out_prefix,"/minimap_unmap.stats",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  file.remove(paste(dire,"/",out_prefix,"/minimap.bam",sep=""))
  file.remove(paste(dire,"/",out_prefix,"/minimap.bai",sep=""))
  file.remove(paste(dire,"/",out_prefix,"/minimap_unmap.fq.gz",sep=""))
}

if (5 %in% Step){
  print("Step 5: flye_meta")
  dire=paste(out_dir,"/flye_meta/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  flye(hifi=paste(out_dir,"/rmHost/",out_prefix,"_rmHost.fq.gz",sep=""),
       out_dir=paste(dire,"/",out_prefix,"/",sep=""),
       threads=threads,
       meta=T)
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            paste(dire,"/",out_prefix,"/assembly.fasta",sep=""),
            ">",
            paste(dire,"/",out_prefix,"/assembly.stats",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  minimap2(long_reads=paste(out_dir,"/rmHost/",out_prefix,"_rmHost.fq.gz",sep=""), # space-separated list for PE
           lr_type="map-hifi", # long read type.
           # "map-pb" for PacBio
           # "map-hifi" for HiFi
           # "map-ont" for ONT reads.
           # "sr" for NGS
           # "asm5" for accurate reads diverging <5% to assembly
           assembly=paste(dire,"/",out_prefix,"/assembly.fasta",sep=""),
           out_prefix=paste(dire,"/",out_prefix,"/minimap",sep=""),
           threads=as.character(threads))
  
  Extract_fq(bam=paste(dire,"/",out_prefix,"/minimap.bam",sep=""),
             paired=F, # logical. T for paired.
             mapped=F, # logical. T for mapped
             out_prefix=paste(dire,"/",out_prefix,"/minimap_unmap",sep=""),
             reads_format="fq", # "fq" or "fa"
             threads=as.character(threads))
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            paste(dire,"/",out_prefix,"/minimap_unmap.fq.gz",sep=""),
            ">",
            paste(dire,"/",out_prefix,"/minimap_unmap.stats",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  file.remove(paste(dire,"/",out_prefix,"/minimap.bam",sep=""))
  file.remove(paste(dire,"/",out_prefix,"/minimap.bai",sep=""))
  file.remove(paste(dire,"/",out_prefix,"/minimap_unmap.fq.gz",sep=""))
}

if (6 %in% Step){
  print("Step 6: mmseqs")
  dire=paste(out_dir,"/mmseqs/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  if (!dir.exists(paste(dire,"/",out_prefix,sep=""))){dir.create(paste(dire,"/",out_prefix,sep=""))}
  
  metaMDBG_asm=paste(out_dir,"/metaMDBG/",out_prefix,"/contigs.fasta",sep="")
  cmd=paste("grep '>' ",metaMDBG_asm," | sed 's/>//' | sed 's/ .*$//'",sep="")
  metaMDBG.id=system(cmd,intern=TRUE)
  metaMDBG.new=paste(out_prefix,"_metaMDBG",sprintf("%08d",1:length(metaMDBG.id)),sep="" )
  write.table(data.frame(metaMDBG.id,metaMDBG.new),
              paste(dire,"/",out_prefix,"/metaMDBG_IDmapper",sep=""),
              sep="\t",row.names = F,col.names = F,quote=F)
  cmd=paste("seqkit replace",
            "-p '^(\\S+)'",
            "-r '{kv}$2'",
            "-k",paste(dire,"/",out_prefix,"/metaMDBG_IDmapper",sep=""),
            metaMDBG_asm,">",
            paste(dire,"/",out_prefix,"/redundant.fasta",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=T)
  
  hifiasm_meta=paste(out_dir,"/hifiasm_meta/",out_prefix,"/",out_prefix,".p_ctg.fna",sep="")
  cmd=paste("grep '>' ",hifiasm_meta," | sed 's/>//' | sed 's/ .*$//'",sep="")
  hifiasm_meta.id=system(cmd,intern=TRUE)
  hifiasm_meta.new=paste(out_prefix,"_hifiasm_meta",sprintf("%08d",1:length(hifiasm_meta.id)),sep="" )
  write.table(data.frame(hifiasm_meta.id,hifiasm_meta.new),
              paste(dire,"/",out_prefix,"/hifiasm_meta_IDmapper",sep=""),
              sep="\t",row.names = F,col.names = F,quote=F)
  cmd=paste("seqkit replace",
            "-p '^(\\S+)'",
            "-r '{kv}$2'",
            "-k",paste(dire,"/",out_prefix,"/hifiasm_meta_IDmapper",sep=""),
            hifiasm_meta,">>",
            paste(dire,"/",out_prefix,"/redundant.fasta",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=T)
  
  flye_meta=paste(out_dir,"/flye_meta/",out_prefix,"/assembly.fasta",sep="")
  cmd=paste("grep '>' ",flye_meta," | sed 's/>//' | sed 's/ .*$//'",sep="")
  flye_meta.id=system(cmd,intern=TRUE)
  flye_meta.new=paste(out_prefix,"_flye_meta",sprintf("%08d",1:length(flye_meta.id)),sep="" )
  write.table(data.frame(flye_meta.id,flye_meta.new),
              paste(dire,"/",out_prefix,"/flye_meta_IDmapper",sep=""),
              sep="\t",row.names = F,col.names = F,quote=F)
  cmd=paste("seqkit replace",
            "-p '^(\\S+)'",
            "-r '{kv}$2'",
            "-k",paste(dire,"/",out_prefix,"/flye_meta_IDmapper",sep=""),
            flye_meta,">>",
            paste(dire,"/",out_prefix,"/redundant.fasta",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=T)
  
  seqNR(in.fasta=paste(dire,"/",out_prefix,"/redundant.fasta",sep=""),
        out_dir=paste(dire,"/",out_prefix,sep=""),
        Identity=0.95, # [0.0,1.0]
        cov_mode=1, # 0: alignment covers ${coverage} of target and of query
                 # 1: alignment covers ${coverage} of target
                 # 2: alignment covers ${coverage} of query
                 # 3: target is of ${coverage} query length
                 # query as representative
        coverage=0.95, # [0.0,1.0]
        threads=threads)
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            paste(dire,"/",out_prefix,"/rep.fasta",sep=""),
            ">",
            paste(dire,"/",out_prefix,"/rep.stats",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  minimap2(long_reads=paste(out_dir,"/rmHost/",out_prefix,"_rmHost.fq.gz",sep=""), # space-separated list for PE
           lr_type="map-hifi", # long read type.
           # "map-pb" for PacBio
           # "map-hifi" for HiFi
           # "map-ont" for ONT reads.
           # "sr" for NGS
           # "asm5" for accurate reads diverging <5% to assembly
           assembly=paste(dire,"/",out_prefix,"/rep.fasta",sep=""),
           out_prefix=paste(dire,"/",out_prefix,"/minimap",sep=""),
           threads=as.character(threads))
  cmd=paste("samtools","coverage",
            "-o",paste(dire,"/",out_prefix,"/rep_cov.tsv",sep=""),
            paste(dire,"/",out_prefix,"/minimap.bam",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("seqkit","fx2tab",
            "--name --only-id --gc --length",
            paste(dire,"/",out_prefix,"/rep.fasta",sep=""),">>",
            paste(dire,"/",out_prefix,"/rep.gclen.tsv",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  Extract_fq(bam=paste(dire,"/",out_prefix,"/minimap.bam",sep=""),
             paired=F, # logical. T for paired.
             mapped=F, # logical. T for mapped
             out_prefix=paste(dire,"/",out_prefix,"/minimap_unmap",sep=""),
             reads_format="fq", # "fq" or "fa"
             threads=as.character(threads))
  cmd=paste("seqkit stats --all -T",
            "-j",as.character(threads),
            paste(dire,"/",out_prefix,"/minimap_unmap.fq.gz",sep=""),
            ">",
            paste(dire,"/",out_prefix,"/minimap_unmap.stats",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  file.remove(paste(dire,"/",out_prefix,"/minimap_unmap.fq.gz",sep=""))
}

if (7 %in% Step){
  print("Step 7: viralVerify")
  dire=paste(out_dir,"/viralVerify/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  viralVerify(in.fna=paste(out_dir,"/mmseqs/",out_prefix,"/rep.fasta",sep=""),
              viralVerify.db=viralVerify.db,
              out_dir=paste(dire,"/",out_prefix,sep=""),
              threads=threads)
}

if (8 %in% Step){
  print("Step 8: metabat")
  dire=paste(out_dir,"/metabat/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  metabat(in.fna=paste(out_dir,"/mmseqs/",out_prefix,"/rep.fasta",sep=""),
          bam.lst=paste(out_dir,"/mmseqs/",out_prefix,"/minimap.bam",sep=""), # space-lst
          threads=threads,
          out_dir=paste(dire,"/",out_prefix,sep=""))
}

if (9 %in% Step){
  print("Step 9: barrnap")
  dire=paste(out_dir,"/barrnap/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  barrnap(in.fna=paste(out_dir,"/mmseqs/",out_prefix,"/rep.fasta",sep=""),
          kingdom="bac", # Kingdom: bac euk arc mito (default 'bac')
          threads=threads, 
          out.gff3=paste(dire,"/",out_prefix,"_bac_rRNA.gff3",sep=""),
          out.fna=paste(dire,"/",out_prefix,"_bac_rRNA.fna",sep=""))
  barrnap(in.fna=paste(out_dir,"/mmseqs/",out_prefix,"/rep.fasta",sep=""),
          kingdom="euk", # Kingdom: bac euk arc mito (default 'bac')
          threads=threads, 
          out.gff3=paste(dire,"/",out_prefix,"_euk_rRNA.gff3",sep=""),
          out.fna=paste(dire,"/",out_prefix,"_euk_rRNA.fna",sep=""))
  barrnap(in.fna=paste(out_dir,"/mmseqs/",out_prefix,"/rep.fasta",sep=""),
          kingdom="arc", # Kingdom: bac euk arc mito (default 'bac')
          threads=threads, 
          out.gff3=paste(dire,"/",out_prefix,"_arc_rRNA.gff3",sep=""),
          out.fna=paste(dire,"/",out_prefix,"_arc_rRNA.fna",sep=""))
}

if (10 %in% Step){
  print("Step 10: metaProdigal")
  dire=paste(out_dir,"/metaProdigal/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  prodigal.meta(fna=paste(out_dir,"/mmseqs/",out_prefix,"/rep.fasta",sep=""),
                out_prefix=paste(dire,"/",out_prefix,sep=""),
                threads=threads)
  gff3=read.table(paste(dire,"/",out_prefix,".prodigal.gff3",sep=""),
                  sep="\t",header=F,quote="")
  gff3$V9=sub("^ID=[0-9]*_","",gff3$V9)
  gff3$V9=paste("ID=",gff3$V1,"_",gff3$V9,sep="")
  write.table(gff3,
              paste(dire,"/",out_prefix,".prodigal_IDformat.gff3",sep=""),
              sep="\t",row.names = F,quote=F,col.names = F)
    
}

if (11 %in% Step){
  print("Step 11: diamond_p_megan")
  dire=paste(out_dir,"/diamond_p_megan/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  diamond_p_megan(query.faa=paste(out_dir,"/metaProdigal/",out_prefix,".prodigal.faa",sep=""),
                  diamond.db=nr, # nr
                  megan.db=megan.db,
                  out_prefix=paste(dire,"/",out_prefix,sep=""),
                  threads=threads)
}

blastn -query $ -db $ -out $ -evalue 1e-5 -outfmt 6 -num_threads $
# qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
sort -k1,1 -k2,2n in.bed > in.sorted.bed
bedtools merge $

if (12 %in% Step){
  print("Step 12: summmary")
  dire=paste(out_dir,"/summmary/",sep="")
  if (!dir.exists(dire)){dir.create(dire)}
  rep.fasta=paste(out_dir,"/mmseqs/",out_prefix,"/rep.fasta",sep="")
  file.copy(from = rep.fasta,
            to = paste(out_dir,"/summmary/",out_prefix,".fna",sep=""))
  
  tbl=read.table(paste(out_dir,"/mmseqs/",out_prefix,"/rep.gclen.tsv",sep=""),
                 sep="\t",header=FALSE,quote="")
  colnames(tbl)=c("contig","length","gc")
  
  rep_cov=read.table(paste(out_dir,"/mmseqs/",out_prefix,"/rep_cov.tsv",sep=""),
                     sep="\t",header=FALSE,quote="")
  rownames(rep_cov)=rep_cov$V1
  tbl$meandepth=rep_cov[tbl$contig,7]
  tbl$coverage=rep_cov[tbl$contig,6]
  
  metaMDBG_IDmapper=read.table(paste(out_dir,"/mmseqs/",out_prefix,"/metaMDBG_IDmapper",sep=""),
                               sep="\t",header=FALSE,quote="")
  rownames(metaMDBG_IDmapper)=metaMDBG_IDmapper$V1
  metaMDBG.c=system(paste("grep '>' ",out_dir,"/metaMDBG/",out_prefix,"/contigs.fasta",sep=""),
                    intern=TRUE)
  metaMDBG.ID=stringr::str_extract(metaMDBG.c,">ctg[0-9]*");metaMDBG.ID=sub(">","",metaMDBG.ID)
  metaMDBG.cir=stringr::str_extract(metaMDBG.c,"circular=.*");metaMDBG.cir= (metaMDBG.cir=="circular=yes")
  df1=data.frame(oldID=metaMDBG.ID,circular=metaMDBG.cir)
  df1$newID=metaMDBG_IDmapper[df1$oldID,2]
  
  flye_meta_IDmapper=read.table(paste(out_dir,"/mmseqs/",out_prefix,"/flye_meta_IDmapper",sep=""),
                                sep="\t",header=FALSE,quote="")
  colnames(flye_meta_IDmapper)=c("oldID","newID")
  assembly_info=read.table(paste(out_dir,"/flye_meta/",out_prefix,"/assembly_info.txt",sep=""),
                           sep="\t",header=FALSE,quote="")
  rownames(assembly_info)=assembly_info$V1
  df2=flye_meta_IDmapper
  df2$circular=assembly_info[df2$oldID,4];df2$circular=df2$circular=="Y"
  
  hifiasm_meta_IDmapper=read.table(paste(out_dir,"/mmseqs/",out_prefix,"/hifiasm_meta_IDmapper",sep=""),
                                   sep="\t",header=FALSE,quote="")
  colnames(hifiasm_meta_IDmapper)=c("oldID","newID")
  hifiasm_meta_IDmapper$circular=grepl("c$",hifiasm_meta_IDmapper$oldID)
  
  df=rbind(df1,df2,hifiasm_meta_IDmapper)
  rownames(df)=df$newID
  tbl$circular=df[tbl$contig,"circular"]
  
  pv=read.csv(paste(out_dir,"/viralVerify/",out_prefix,"/rep_result_table.csv",sep=""),
                sep=",",header=TRUE)
  rownames(pv)=pv$Contig.name 
  tbl$viral=pv[tbl$contig,"Prediction"]
  
  metabat=system(paste("ls ",out_dir,"/metabat/",out_prefix,"/rep.fasta.metabat-bins*/*",sep=""),intern = T)
  metabat=read.table(metabat,sep="\t",header=F,quote="")
  rownames(metabat)=metabat$V1;metabat$V2=paste("bin",as.character(metabat$V2),sep="")
  tbl$metabat=metabat[tbl$contig,2]
  
  for (i in c("bac","arc","euk")){
    barrnap.bac=read.table(paste(out_dir,"/barrnap/filtered_",out_prefix,"_",i,"_rRNA.gff3",sep=""),
                           sep="\t",header=FALSE,quote="")
    barrnap.bac$name=stringr::str_extract(barrnap.bac$V9,"Name=.*;p")
    barrnap.bac$name=sub("Name=","",barrnap.bac$name);barrnap.bac$name=sub("_rRNA;p","",barrnap.bac$name)
    tbl[,paste("barrnap",i,sep=".")]=sapply(tbl$contig,
                           function(c){
                             l=barrnap.bac[barrnap.bac$V1==c,"name"]
                             l=l[!duplicated(l)]
                             return(paste(l,collapse = ";"))
                           })
    tbl[tbl[,paste("barrnap",i,sep=".")]=="",paste("barrnap",i,sep=".")]=NA
  }
  
  # if (file.exists(paste(out_dir,"/diamond_d_megan/",out_prefix,"_taxon.tsv",sep=""))){
  #   taxon=read.table(paste(out_dir,"/diamond_d_megan/",out_prefix,"_taxon.tsv",sep=""),
  #                    sep="\t",header=FALSE,quote="",fill = T)
  #   rownames(taxon)=taxon$V1
  #   tbl$taxon=taxon[tbl$contig,3]
  # }
  
  prodigal.gff3=paste(out_dir,"/metaProdigal/",out_prefix,".prodigal_IDformat.gff3",sep="")
  prodigal.gff3=read.table(prodigal.gff3,sep="\t",header=F,quote="")
  prodigal.gff3$V10=prodigal.gff3$V5-prodigal.gff3$V4+1
  tbl$CDS.n=table(prodigal.gff3$V1)[tbl$contig]
  tbl$CDS.len=tapply(prodigal.gff3$V10,prodigal.gff3$V1,sum)[tbl$contig]
  
  taxon=read.table(paste(out_dir,"/diamond_p_megan/",out_prefix,"_taxon.tsv",sep=""),
                   sep="\t",header=FALSE,quote="",fill = T)
  taxon=taxon[taxon$V3!="",]
  taxon$V4=sub("_[0-9]*$","",taxon$V1)
  tbl$CDS.taxa=tapply(taxon$V3,taxon$V4,
                     function(v){
                       j=sort(-table(v))[1]
                       return(names(j))
                     })[tbl$contig]
  tbl$CDSinTaxa=tapply(taxon$V3,taxon$V4,
                      function(v){
                        j=sort(-table(v))[1]
                        return(abs(unname(j)))
                      })[tbl$contig]
  write.table(tbl,paste(dire,"/",out_prefix,".tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
  
}



