
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

# Extract reads mapped/unmapped to an assembly via SAMtools.
# Dependencies: SAMtools
Extract_fq=function(bam=bam,
                    paired=paired, # logical. T for paired. 
                    mapped=mapped, # logical. T for mapped
                    out_prefix=out_prefix,
                    reads_format=reads_format, # "fq" or "fa"
                    threads=threads){
  threads=as.character(threads)
  if (paired){
    if (mapped){
      flag="-f 2" #only include reads with PROPER_PAIR (each segment properly aligned according to the aligner)
    }else{
      flag="-f 12" #only include reads with UNMAP & MUNMAP (segment unmapped & next segment in the template unmapped)
    }
    if (reads_format=="fq"){
      samtools="samtools fastq"}else{samtools="samtools fasta"}
    cmd=paste(samtools,
              flag,
              "-@",threads,
              "-1",paste(out_prefix,".1.",reads_format,sep=""),
              "-2",paste(out_prefix,".2.",reads_format,sep=""),
              bam,
              sep=" ")
  }
  if (!paired){
    if (mapped){
      flag="-G 4" # only EXCLUDE reads with UNMAP (segment unmapped)
    }else{
      flag="-f 4" # only include reads with UNMAP (segment unmapped)
    }
    if (reads_format=="fq"){
      samtools="samtools fastq"}else{samtools="samtools fasta"}
    cmd=paste(samtools,
              flag,
              "-@",threads,
              bam,">",
              paste(out_prefix,".",reads_format,sep=""),
              sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  if (paired){
    cmd=paste("gzip ",out_prefix,".1.",reads_format,sep="")
    print(cmd);system(cmd,wait=T)
    cmd=paste("gzip ",out_prefix,".2.",reads_format,sep="")
    print(cmd);system(cmd,wait=T)
  }else{
    cmd=paste("gzip ",out_prefix,".",reads_format,sep="")
    print(cmd);system(cmd,wait=T)
  }
}

# metaMDBG: Nanopore/Hifi
metaMDBG=function(hifi.fq=hifi.fq,
                  nanopore.fq=NA,
                  out_dir=out_dir,
                  threads=threads){
  cmd=paste("metaMDBG asm",
            "--out-dir",out_dir,
            "--threads",as.character(threads),
            sep=" ")
  if (!is.na(hifi.fq)){cmd=paste(cmd,"--in-hifi",hifi.fq,sep=" ")}
  if (!is.na(nanopore.fq)){cmd=paste(cmd,"--in-ont",nanopore.fq,sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
}

hifiasm_meta=function(hifi.fq=hifi.fq,
                      out_prefix=out_prefix,
                      threads=threads){
  cmd=paste("hifiasm_meta",
            "-o",out_prefix,
            "-t",as.character(threads),
            hifi.fq,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("awk '/^S/{print \">\"$2;print $3}'",
            paste(out_prefix,".p_ctg.gfa",sep=""),
            ">",
            paste(out_prefix,".p_ctg.fna",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=T)
}

flye=function(hifi=hifi,
              out_dir=out_dir,
              threads=threads,
              meta=T){
  if (!dir.exists(out_dir)){dir.create(out_dir)}
  cmd=paste("flye",
            "--pacbio-hifi",hifi,
            "--out-dir",out_dir,
            "--threads",as.character(threads),
            sep=" ")
  if (meta){cmd=paste(cmd,"--meta",sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
}

seqNR=function(in.fasta=in.fasta,
               out_dir=out_dir,
               Identity=Identity, # [0.0,1.0]
               cov_mode=cov_mode, # 0: alignment covers ${coverage} of target and of query
               # 1: alignment covers ${coverage} of target
               # 2: alignment covers ${coverage} of query
               # 3: target is of ${coverage} query length
               # quert as representative
               coverage=coverage, # [0.0,1.0]
               threads=threads){
  threads=as.character(threads)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("mmseqs createdb",in.fasta,"sequenceDB",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cov_mode=as.character(cov_mode)
  coverage=as.character(coverage)
  Identity=as.character(Identity)
  cmd=paste("mmseqs cluster sequenceDB clusterDB", out_dir,
            "--cov-mode",cov_mode,
            "-c",coverage,
            "--min-seq-id",Identity,
            "--threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="mmseqs createsubdb clusterDB sequenceDB clusterDB_rep"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="mmseqs convert2fasta clusterDB_rep rep.fasta"
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="mmseqs createtsv sequenceDB sequenceDB clusterDB clusters.tsv"
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}

# # best blastn of contigs (lowest evalue, longest alignment)
# blastn2nt=function(query.fna=query.fna,
#                    nt="/bucket/BourguignonU/Cong/public_db/ncbi.blastdb_20250610/nt",
#                    out_prefix=out_prefix,
#                    threads=threads){
#   cmd=paste("blastn",
#             "-num_threads",as.character(threads),
#             "-query",query.fna,
#             "-db",nt,
#             "-out",paste(out_prefix,".blastn",sep=""),
#             "-outfmt '6 qaccver saccver qcovs pident length evalue bitscore staxid'")
#   print(cmd);system(cmd,wait=TRUE)
# 
#   blastn=read.table(paste(out_prefix,".blastn",sep=""),sep="\t",header=FALSE,quote="")
#   colnames(blastn)=c("qaccver","saccver","qcovs","pident","length","evalue","bitscore","staxid")
# 
#   res=data.frame(query=blastn$qaccver,subject=NA)
#   res=res[!duplicated(res$query),]
#   library(parallel)
#   clus=makeCluster(threads)
#   clusterExport(clus,"blastn")
#   res[,c("subject","qcovs","pident","staxid")]=
#     t(parSapply(clus,res$query,
#                 function(q){
#                   df=blastn[blastn$qaccver==q,]
#                   df=df[df$evalue==min(df$evalue),]
#                   df=df[df$length==max(df$length),]
#                   return( df[1,c("saccver","qcovs","pident","staxid")] )
#                 }))
# 
#   taxid.lst=res$staxid;taxid.lst=taxid.lst[!duplicated(taxid.lst)]
#   writeLines(taxid.lst,paste(out_prefix,"_taxID.lst",sep=""))
#   cmd=paste("taxonkit lineage -R",
#             paste(out_prefix,"_taxID.lst",sep=""),
#             ">",paste(out_prefix,"_taxID.tsv",sep=""))
#   print(cmd);system(cmd,wait=TRUE)
#   taxID=read.table(paste(out_prefix,"_taxID.tsv",sep=""),header=FALSE,quote="",sep="\t")
#   rownames(taxID)=taxID$V1
#   res[,c("staxid.lineage","staxid.rank")]=taxID[res$staxid,c(2,3)]
# 
#   write.table(res,paste(out_prefix,".blastn2nt",sep=""),
#               sep="\t",row.names=FALSE,quote=FALSE)
# }

# # best blastn of contigs (lowest evalue, longest alignment)
# blastn2nt=function(query.fna=query.fna,
#                    nt="/bucket/BourguignonU/Cong/public_db/ncbi.blastdb_20250610/nt",
#                    out_prefix=out_prefix,
#                    threads=threads){
#   if (!dir.exists(paste(out_prefix,"_tmp",sep=""))){dir.create(paste(out_prefix,"_tmp",sep=""))}
#   cmd=paste("seqkit fx2tab --no-qual --only-id ",query.fna,
#             " > ",out_prefix,"_tmp/sequences.tsv",sep="")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   seqs=read.table(paste(out_prefix,"_tmp/sequences.tsv",sep=""),
#                   sep="\t",header=FALSE,quote="")
#   seqs.all=seqs
#   
#   seqs$stamp=paste(out_prefix,"_tmp/",seqs$V1,".finished",sep="")
#   seqs=seqs[!file.exists(seqs$stamp),]
#   seqs=seqs[sample(1:nrow(seqs),nrow(seqs),replace=F),]
#   
#   library(parallel)
#   clus=makeCluster(threads)
#   clusterExport(clus,"seqs",envir=environment())
#   clusterExport(clus,"out_prefix",envir=environment())
#   clusterExport(clus,"nt",envir=environment())
#   parSapply(clus,1:nrow(seqs),
#             function(i){
#               id=seqs[i,1]
#               s=seqs[i,2]
#               write(paste(">",id,sep=""),
#                     paste(out_prefix,"_tmp/",id,".fna",sep=""))
#               write(s,
#                     paste(out_prefix,"_tmp/",id,".fna",sep=""),
#                     append = T)
#               
#               stamp=paste(out_prefix,"_tmp/",id,".finished",sep="")
#               if (!file.exists(stamp)){
#                 cmd=paste("blastn",
#                           "-query",paste(out_prefix,"_tmp/",id,".fna",sep=""),
#                           "-db",nt,
#                           "-out",paste(out_prefix,"_tmp/",id,".blastn",sep=""),
#                           "-outfmt '6 qaccver saccver qcovs pident length evalue bitscore staxid'")
#                 print(cmd);system(cmd,wait=TRUE)
#                 if (file.size(paste(out_prefix,"_tmp/",id,".blastn",sep=""))!=0){
#                   blastn=read.table(paste(out_prefix,".blastn",sep=""),sep="\t",header=FALSE,quote="")
#                   colnames(blastn)=c("qaccver","saccver","qcovs","pident","length","evalue","bitscore","staxid")
#                   df=blastn[blastn$evalue==min(blastn$evalue),]
#                   df=df[df$length==max(df$length),]
#                   df=df[1,]
#                   write.table(df,paste(out_prefix,"_tmp/",id,".bestHit",sep=""),
#                               sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
#                 }
#                 
#                 file.create(stamp)
#               }
#             })
#   blast.lst=paste(out_prefix,"_tmp/",seqs.all$V1,".bestHit",sep="")
#   blast.lst=blast.lst[file.exists(blast.lst)]
#   for (b in blast.lst){
#     cmd=paste("cat",b,">>",
#               paste(out_prefix,".blastn2nt",sep=""))
#   }
#   blastn=read.table(paste(out_prefix,".blastn2nt",sep=""),
#                     sep="\t",header=FALSE,quote="")
#   rownames(blastn)=c("qaccver","saccver","qcovs","pident","length","evalue","bitscore","staxid")
#   
#   taxid.lst=blastn$staxid;taxid.lst=taxid.lst[!duplicated(taxid.lst)]
#   writeLines(taxid.lst,paste(out_prefix,"_taxID.lst",sep=""))
#   cmd=paste("taxonkit lineage -R",
#             paste(out_prefix,"_taxID.lst",sep=""),
#             ">",paste(out_prefix,"_taxID.tsv",sep=""))
#   print(cmd);system(cmd,wait=TRUE)
#   taxID=read.table(paste(out_prefix,"_taxID.tsv",sep=""),header=FALSE,quote="",sep="\t")
#   rownames(taxID)=taxID$V1
#   blastn[,c("staxid.lineage","staxid.rank")]=taxID[res$staxid,c(2,3)]
#   
#   write.table(blastn,paste(out_prefix,".blastn2nt",sep=""),
#               sep="\t",row.names=FALSE,quote=FALSE)
#   cmd=paste("rm -r ",out_prefix,"_tmp",sep="")
#   print(cmd);system(cmd,wait = T)
# }

viralVerify=function(in.fna=in.fna,
                     viralVerify.db="/bucket/BourguignonU/Cong/public_db/virusVerify/nbc_hmms.hmm",
                     out_dir=out_dir,
                     threads=threads){
  if (!dir.exists(out_dir)){dir.create(out_dir)}
  cmd=paste("viralverify",
            "-f",in.fna,
            "-o",out_dir,
            "--hmm",viralVerify.db,
            "-t",threads,
            "-p",
            sep=" ")
  print(cmd);system(cmd,wait=T)
}

metabat=function(in.fna=in.fna,
                 bam.lst=bam.lst, # space-lst
                 threads=threads,
                 out_dir=out_dir){
  wd=getwd()
  
  if (!dir.exists(out_dir)){dir.create(out_dir)}
  setwd(out_dir)
  cmd=paste("runMetaBat.sh",
            "-t",as.character(threads),
            "--saveCls --noBinOut",
            in.fna,bam.lst,
            sep=" ")
  print(cmd);system(cmd,wait=T)
  setwd(wd)
}

# Barrnap: rRNA
barrnap=function(in.fna=in.fna,
                 kingdom="bac", # Kingdom: bac euk arc mito (default 'bac')
                 threads=threads, 
                 out.gff3=out.gff3,
                 out.fna=out.fna){
  cmd=paste("barrnap",
            "--kingdom",kingdom,
            "--threads",as.character(threads),
            "--outseq",out.fna,
            in.fna,">",
            out.gff3,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("grep -v partial ",out.gff3," > ",dirname(out.gff3),"/filtered_",basename(out.gff3),sep="")
  print(cmd);system(cmd,wait=TRUE)
}


# prodigal for metagenome
# Dependencies: seqkit, prodigal, parallel (R)
prodigal.meta=function(fna=fna,
                       out_prefix=out_prefix,
                       threads=threads){
  cmd=paste("prodigal",
            "-i",fna,
            "-o",paste(out_prefix,".prodigal.gff3",sep=""),
            "-a",paste(out_prefix,".prodigal.faa",sep=""),
            "-d",paste(out_prefix,".prodigal.fna",sep=""),
            "-p meta -f gff",sep=" ")
  system(cmd,wait=TRUE)
  # threads=as.character(threads)
  # 
  # cmd=paste("seqkit split2",
  #           "-s 1",
  #           "-j",threads,
  #           "-O",paste(out_prefix,"_tmp",sep=""),
  #           fna,
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  # 
  # fna.lst=system(paste("ls ",out_prefix,"_tmp/*",sep=""),intern=TRUE)
  # library(parallel)
  # clus=makeCluster(threads)
  # parSapply(clus,
  #           fna.lst,
  #           function(fna.i){
  #             cmd=paste("prodigal",
  #                       "-i",fna.i,
  #                       "-o",paste(fna.i,".prodigal",sep=""),
  #                       "-a",paste(fna.i,".prodigal.faa",sep=""),
  #                       "-d",paste(fna.i,".prodigal.fna",sep=""),
  #                       "-p meta -f gff",sep=" ")
  #             system(cmd,wait=TRUE)
  #           })
  # prodigal=system(paste("ls ",out_prefix,"_tmp/*.prodigal",sep=""),intern=TRUE)
  # for (i in prodigal){
  #   system(paste("cat ",i," >> ",out_prefix,".prodigal.gff3",sep=""))
  #   system(paste("cat ",i,".faa"," >> ",out_prefix,".prodigal.faa",sep=""))
  #   system(paste("cat ",i,".fna"," >> ",out_prefix,".prodigal.fna",sep=""))
  # }
  # #system(paste("rm -r ",out_prefix,"_tmp",sep=""))
}

# Protein-protein search by DIAMOND and assign to taxa by MEGAN
# Dependencies: DIAMOND
diamond_p_megan=function(query.faa=query.faa,
                         diamond.db=diamond.db, # nr
                         megan.db=megan.db,
                         out_prefix=out_prefix,
                         threads=threads){
  threads=as.character(threads)
  cmd=paste("diamond blastp",
            "-p",threads,
            "-d",diamond.db,
            "-q",query.faa,
            "--outfmt 100",
            "--out",paste(out_prefix,".blast.daa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("daa-meganizer",
            "-i",paste(out_prefix,".blast.daa",sep=""),
            "-mdb",megan.db,
            "-t",threads,
            "-ram readCount",
            "-supp 0",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("daa2info",
            "-i",paste(out_prefix,".blast.daa",sep=""),
            "-o",paste(out_prefix,"_taxon.tsv",sep=""),
            "-r2c Taxonomy",
            "-n true",
            "-p true",
            "-r true",
            "-mro","true",
            "-u false",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}


# DNA-protein search by DIAMOND and assign to taxa by MEGAN
# Dependencies: DIAMOND
diamond_d_megan=function(query.fna=query.fna,
                         diamond.db="/bucket/BourguignonU/Cong/public_db/ncbi.blastdb_20250610/nr.dmdb.dmnd", # nr
                         megan.db="/bucket/BourguignonU/Cong/public_db/megan.db/megan-nr-r2.mdb",
                         geneticCode=1, # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
                         out_prefix=out_prefix,
                         threads=threads){
  threads=as.character(threads)
  cmd=paste("diamond blastx",
            "-p",threads,
            "-d",diamond.db,
            "-q",query.fna,
            "--outfmt 100",
            "--range-culling -k 25 -F 15",
            "--query-gencode",as.character(geneticCode),
            "--out",paste(out_prefix,".blast.daa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("daa-meganizer",
            "-i",paste(out_prefix,".blast.daa",sep=""),
            "-mdb",megan.db,
            "-t",threads,
            "-ram readCount",
            "-supp 0",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("daa2info",
            "-i",paste(out_prefix,".blast.daa",sep=""),
            "-o",paste(out_prefix,"_taxon.tsv",sep=""),
            "-r2c Taxonomy",
            "-n true",
            "-p true",
            "-r true",
            "-mro","true",
            "-u false",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}