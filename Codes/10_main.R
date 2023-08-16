# Ortholog groups, gene functions & metabolic network

#####
# Ortholog groups
#####
# OrthoFinder
# OrthoFinder creates a results directory called ‘OrthoFinder’ inside the input proteome directory and puts the results here.
# Orthofinder with "-M msa" infers phylogeny by mafft & fasttree (too slow)
# Dependencies: Orthofinder, mafft, fasttree
Orthofinder=function(in_dir=in_dir, # Input proteome directory. 
                     # One file per species with extension '.faa'
                     threads=threads){
  threads=as.character(threads)
  
  cmd=paste("orthofinder.py",
            "-f",in_dir,
            "-t",threads,
            "-a",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Re-run orthofinder with designated phylogenetic tree
Re_orthofinder=function(previous_orthofinder_result_dir=previous_orthofinder_result_dir,
                        tree=tree,
                        threads=threads){
  threads=as.character(threads)
  
  cmd=paste("orthofinder.py",
            "-t",threads,
            "-a",threads,
            "-ft",previous_orthofinder_result_dir,
            "-s",tree,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# OMA
# 不好用。。。

# Best protein-protein hit by blast
best_blastp=function(query.faa=query.faa,
                     reference.faa=reference.faa,
                     out_dir=out_dir,
                     out_basename=out_basename,
                     threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",query.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  query.faa=unlist(strsplit(query.faa,"/"));query.faa=query.faa[length(query.faa)]
  
  cmd=paste("cp",reference.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  reference.faa=unlist(strsplit(reference.faa,"/"));reference.faa=reference.faa[length(reference.faa)]
  
  cmd=paste("makeblastdb",
            "-in",reference.faa,
            "-dbtype","prot",
            "-parse_seqids",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("blastp",
            "-num_threads",threads,
            "-db",reference.faa,
            "-query",query.faa,
            "-outfmt 6",
            "-evalue 1e-5",
            "-num_alignments 1",
            "-out",paste(out_dir,"/",out_basename,".blast",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",query.faa,sep=" "),wait=TRUE)
  system(paste("rm",reference.faa,sep=" "),wait=TRUE)
  system(paste("rm"," ",reference.faa,".*",sep=""),wait=TRUE)
  setwd(wd)
}

#####
# Gene functions
#####
# Interproscan
# Dependencies: Interproscan, seqkit, parallel (R)
# AT LEAST 16 THREADS
interpro=function(proteins.faa=proteins.faa,
                  out_dir=out_dir,
                  out_basename=out_basename,
                  threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",proteins.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  proteins.faa=unlist(strsplit(proteins.faa,"/"));proteins.faa=proteins.faa[length(proteins.faa)]
  cmd=paste("sed -i 's/*/X/g'",proteins.faa,sep=" ")
  print(cmd);system(cmd,wait=TRUE)

  cmd=paste("seqkit split2",proteins.faa,"-s 1000",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  prot.lst=system(paste("ls ",proteins.faa,".split/*",sep=""),intern=TRUE)
  chunks=as.numeric(threads)-1
  
  for (prot in prot.lst){
    cmd=paste("interproscan.sh",
              "-dp",
              "-b",basename(prot),
              "-cpu",1,
              "-f TSV",
              "-goterms",
              "-i",prot,
              "--pathways",
              sep=" ")
    system(paste("echo '",cmd,"' >> jobs",sep=""))
  }
  cmd=paste("split",
            "-d",
            "-n",paste("l/",chunks,sep=""),
            "jobs",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads)-1)
  parSapply(clus,
            1:chunks,
            function(i){
              scr=system("ls x*",wait=TRUE,intern=TRUE)[i]
              
              system(paste("sed -i '1i module load Other\\/interproscan\\/5.60-92.0'",scr,sep=" "))
              system(paste("sed -i '1i module load bioinfo-ugrp-modules'",scr,sep=" "))
              
              cmd=paste("bash"," ",scr,sep="")
              print(cmd);system(cmd,wait=TRUE)})
  stopCluster(clus)
  
  
  for (prot in prot.lst){
    system(paste("cat ",basename(prot),".tsv >> ",out_basename,".interpro.tsv",sep=""))
    system(paste("rm ",basename(prot),".tsv",sep=""))
  }
  
  system("rm x*")
  system(paste("rm",proteins.faa,sep=" "),wait=TRUE)
  system(paste("rm -r ",proteins.faa,".split/",sep=""))
  system("rm -r temp")
  setwd(wd)
}

# eggNOG-mapper
# Dependencies: eggNOG-mapper, seqkit, parallel (R)
eggNOGmapper=function(proteins.faa=proteins.faa,
                      out_dir=out_dir,
                      out_basename=out_basename,
                      threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",proteins.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  proteins.faa=basename(proteins.faa)
  
  cmd=paste("seqkit split2",proteins.faa,"-s 1000",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  prot.lst=system(paste("ls ",proteins.faa,".split/*",sep=""),intern=TRUE)
  chunks=floor( as.numeric(threads)-1 )
  
  for (prot in prot.lst){
    cmd=paste("emapper.py",
              "-i",prot,
              "-o",basename(prot),
              "--cpu",1,
              sep=" ")
    system(paste("echo",cmd,">> jobs",sep=" "))
  }
  cmd=paste("split",
            "-d",
            "-n",paste("l/",chunks,sep=""),
            "jobs",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads)-1)
  parSapply(clus,
            1:chunks,
            function(i){
              scr=system("ls x*",wait=TRUE,intern=TRUE)[i]
              cmd=paste("bash"," ",scr,sep="")
              print(cmd);system(cmd,wait=TRUE)})
  stopCluster(clus)
  
  system("rm x*")
  system(paste("rm",proteins.faa,sep=" "),wait=TRUE)
  system(paste("rm -r ",proteins.faa,".split/",sep=""))
  
  for (prot in prot.lst){
    system(paste("cat ",basename(prot),".emapper.annotations >> ",out_basename,".emapper.annotations.tsv",sep=""))
    system(paste("rm ",basename(prot),".emapper.annotations",sep=""))
    system(paste("cat ",basename(prot),".emapper.hits >> ",out_basename,".emapper.hits.tsv",sep=""))
    system(paste("rm ",basename(prot),".emapper.hits",sep=""))
    system(paste("cat ",basename(prot),".emapper.seed_orthologs >> ",out_basename,".emapper.seed_orthologs.tsv",sep=""))
    system(paste("rm ",basename(prot),".emapper.seed_orthologs",sep=""))
  }
  setwd(wd)
}
  
# Gene function by best blast hits
# Dependencies: DIAMOND+nr, parallel (R)
geneFun_bbh=function(pep.faa=pep.faa,
                     nr.dmdb="/apps/unit/BioinfoUgrp/DB/diamondDB/ncbi/2022-07/nr.dmnd",
                     out_prefix=out_prefix,
                     threads=threads){
  threads=as.character(threads)
  
  # blast search
  blast=paste(out_prefix,".blast",sep="")
  if (!file.exists(blast)){
    cmd=paste("diamond blastp",
              "--threads",threads,
              "--db",nr.dmdb,
              "--query",pep.faa,
              "--out",blast,
              "--min-score 50",
              "--query-cover 75",
              "--outfmt 6 qseqid sseqid length evalue bitscore stitle",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  blast=read.table(blast,sep="\t",header=FALSE,quote="")
  colnames(blast)=c("qseqid","sseqid","length","evalue","bitscore","stitle")
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  clusterExport(clus,list("blast"),envir=environment())
  blast[,"retainHit"]=parSapply(clus,
                                1:nrow(blast),
                                function(i){
                                  query=blast[i,"qseqid"]
                                  evalue=blast[i,"evalue"]
                                  
                                  d=blast[blast[,"qseqid"]==query,]
                                  bestEval=min(d[,"evalue"])
                                  if (evalue!=bestEval){return(FALSE)}else{return(TRUE)}
                                })
  blast=blast[blast[,"retainHit"],]
  clusterExport(clus,list("blast"),envir=environment())
  blast[,"retainHit"]=parSapply(clus,
                                1:nrow(blast),
                                function(i){
                                  query=blast[i,"qseqid"]
                                  length=blast[i,"length"]
                                  
                                  d=blast[blast[,"qseqid"]==query,]
                                  longest=max(d[,"length"])
                                  if (length!=longest){return(FALSE)}else{return(TRUE)}
                                })
  stopCluster(clus)
  
  blast=blast[blast[,"retainHit"],]
  blast=blast[blast[,"evalue"]<1e-5,]
  write.table(blast,paste(out_prefix,"_bbh.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
}

# KOfamScan
KOfamScan=function(pep.faa=pep.faa,
                   KOfamScan.ko_list="/bucket/BourguignonU/Cong/public_db/kofamscan/ko_list",
                   KOfamScan.profiles="/bucket/BourguignonU/Cong/public_db/kofamscan/profiles",
                   out_prefix=out_prefix,
                   threads=threads){
  threads=as.character(threads)
  
  cmd=paste("exec_annotation",
            "-o",paste(out_prefix,".tsv",sep=""),
            pep.faa,
            "-p",KOfamScan.profiles,
            "-k",KOfamScan.ko_list,
            paste("--cpu=",threads,sep=""),
            "-f mapper")
  print(cmd);system(cmd,wait=TRUE)
}
  
  
  
  
