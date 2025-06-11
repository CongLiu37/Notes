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

# Dependencies: MAFFT, PAL2NAL,treeio (R), ape (R), ggtree (R), tidytree (R), parallel (R)
formatHOGs_orthofinder=function(N0.tsv=N0.tsv,
                                spTree_4orthofinder.nwk=spTree_4orthofinder.nwk, # only tree topology is used
                                pep.faa.lst=pep.faa.lst, # comma-lst
                                cds.fna.lst=cds.fna.lst, # comma-lst
                                out_dir=out_dir,
                                out_prefix=out_prefix,
                                threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  
  N0=read.table(N0.tsv,sep="\t",header=TRUE,quote="")
  N0=N0[,c(1,4:ncol(N0))]
  N0$HOG=sub("^N0\\.","",N0$HOG)
  HOG=N0$HOG
  N0=cbind(HOG,
           data.frame(apply(N0[,-1],c(1,2),
                            function(i){return(gsub(" ","",i))})))
  
  HOGxSpecies4copy=cbind(HOG,
                         data.frame(apply(N0[,-1],c(1,2),
                                          function(i){return(length(unlist(strsplit(i,","))))})))
  write.table(HOGxSpecies4copy,
              paste(out_dir,"/",out_prefix,"_HOGxSpecies4copy.tsv",sep=""),
              row.names=FALSE,quote=FALSE,sep="\t")

  library(parallel)
  clus=makeCluster(threads)
  
  clusterExport(cl=clus,varlist=list("HOGxSpecies4copy"),envir=environment())
  single_copy=parSapply(clus,HOGxSpecies4copy$HOG,
                        function(hog){
                          n=unlist(HOGxSpecies4copy[HOGxSpecies4copy$HOG==hog,-1])
                          return(all(n==1))
                        })
  single_copy_HOG=HOGxSpecies4copy[which(single_copy),"HOG"]
  writeLines(single_copy_HOG,
             paste(out_dir,"/",out_prefix,"_single_copy_HOG.lst",sep=""))
  
  if (!file.exists(paste(out_dir,"/",out_prefix,"_HOG_stat.tsv",sep=""))){
    tree=treeio::read.newick(spTree_4orthofinder.nwk)
    tree.data=ggtree::ggtree(tree)$data
    tree.data=as.data.frame(tree.data)
    node2label=tree.data$label
    names(node2label)=as.character(tree.data$node)
    HOG_stat=data.frame(HOG=HOG,
                        num_genes=apply(HOGxSpecies4copy[,-1],1,sum))
    clusterExport(cl=clus,varlist=list("HOGxSpecies4copy","tree","node2label"),envir=environment())
    HOG_stat[,c("sp.lst","sp.lst.len","mrca.node","mrca.offsprings",
                "mrca.offsprings.len","lost.sp")]=t(parSapply(clus,HOG_stat$HOG,
                                                           function(hog){
                                                             copyN=HOGxSpecies4copy[HOGxSpecies4copy$HOG==hog,-1]
                                                             copyN=unlist(copyN)
                                                             sp.lst=names(which(copyN!=0))
                                                             
                                                             sp.lst.len=length(sp.lst)
                                                             sp.lst=paste(sp.lst,collapse=",")
                                                             
                                                             if (sp.lst.len==1){mrca.node=sp.lst}else{mrca.node=ape::getMRCA(phy=tree,tip=unlist(strsplit(sp.lst,",")))}
                                                             if (sp.lst.len==1){mrca.offsprings=sp.lst}else{
                                                               lst=tidytree::offspring(tree,mrca.node,tiponly=TRUE)
                                                               mrca.offsprings=paste(node2label[as.character(lst)],collapse=",")
                                                             }
                                                             mrca.offsprings.len=length(unlist(strsplit(mrca.offsprings,",")))
                                                             lost.sp=setdiff(unlist(strsplit(mrca.offsprings,",")),
                                                                             unlist(strsplit(sp.lst,",")))
                                                             lost.sp=paste(lost.sp,collapse=",")
                                                             return(c(sp.lst,sp.lst.len,mrca.node,mrca.offsprings,
                                                                      mrca.offsprings.len,lost.sp))
                                                           }))
    write.table(HOG_stat,
                paste(out_dir,"/",out_prefix,"_HOG_stat.tsv",sep=""),
                row.names=FALSE,quote=FALSE,sep="\t")
  }
  
  lst.dir=paste(out_dir,"/",out_prefix,"_lst",sep="")
  system(paste("mkdir",lst.dir,sep=" "))
  sapply(N0$HOG,
         function(hog){
           i=unlist(N0[N0$HOG==hog,-1])
           i=i[i!=""]
           i=paste(i,collapse=",")
           i=unlist(strsplit(i,","))
           writeLines(i,paste(lst.dir,"/",hog,".lst",sep=""))
         })
  
  pep.faa.lst=unlist(strsplit(pep.faa.lst,","))
  all.pep=paste(out_dir,"/",out_prefix,"_all.pep.faa",sep="")
  cds.fna.lst=unlist(strsplit(cds.fna.lst,","))
  all.cds=paste(out_dir,"/",out_prefix,"_all.cds.fna",sep="")
  for (i in 1:length(pep.faa.lst)){
    cmd=paste("cat",pep.faa.lst[i],">>",all.pep,sep=" ");system(cmd,wait=TRUE)
    cmd=paste("cat",cds.fna.lst[i],">>",all.cds,sep=" ");system(cmd,wait=TRUE)
  }
  
  pep.dir=paste(out_dir,"/",out_prefix,"_pep",sep="")
  system(paste("mkdir",pep.dir,sep=" "))
  cds.dir=paste(out_dir,"/",out_prefix,"_cds",sep="")
  system(paste("mkdir",cds.dir,sep=" "))
  mafft.dir=paste(out_dir,"/",out_prefix,"_mafft",sep="")
  system(paste("mkdir",mafft.dir,sep=" "))
  pal2nal.dir=paste(out_dir,"/",out_prefix,"_pal2nal",sep="")
  system(paste("mkdir",pal2nal.dir,sep=" "))
  clusterExport(cl=clus,varlist=list("lst.dir","pep.dir","cds.dir","mafft.dir","pal2nal.dir","all.pep","all.cds"),envir=environment())
  parSapply(clus,N0$HOG,
            function(hog){
              lst=paste(lst.dir,"/",hog,".lst",sep="")
              pep=paste(pep.dir,"/",hog,".faa",sep="")
              cds=paste(cds.dir,"/",hog,".fna",sep="")
              mafft=paste(mafft.dir,"/",hog,"_mafft.faa",sep="")
              pal2nal=paste(pal2nal.dir,"/",hog,"_pal2nal.fna",sep="")
              
              cmd=paste("seqkit grep",
                        "-f",lst,
                        all.pep,">",pep,
                        sep=" ")
              system(cmd,wait=TRUE)
              cmd=paste("seqkit grep",
                        "-f",lst,
                        all.cds,">",cds,
                        sep=" ")
              system(cmd,wait=TRUE)
              
              cmd=paste("mafft --auto",
                        pep,">",mafft,
                        sep=" ")
              system(cmd,wait=TRUE)
              cmd=paste("pal2nal.pl",
                        mafft,
                        cds,
                        "-output fasta",
                        "-codontable 1",
                        ">",pal2nal,
                        sep=" ")
              system(cmd,wait=TRUE)
            })
  
  system(paste("rm",all.pep,sep=" "))
  system(paste("rm",all.cds,sep=" "))

}
  
# Orthofinder N0.tsv statistics
N0_Stat=function(N0.tsv=N0.tsv,
                 orthofinder.in_dir=orthofinder.in_dir, # all .faa for orthofinder
                 out_prefix=out_prefix){
  d=read.table(N0.tsv,sep="\t",header=TRUE,quote="")
  row.names(d)=d[,"HOG"]
  d=d[,4:ncol(d)]
  mat=as.matrix(d)
  mat=apply(mat,c(1,2),
            function(i){ return(length(unlist(strsplit(i,", ")))) })
  copyNumber=as.data.frame(mat)
  copyNumber=cbind(data.frame(HOG=rownames(copyNumber)),copyNumber)
  write.table(copyNumber,paste(out_prefix,"_copyNumber.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
  
  species=colnames(d)
  res=data.frame(species=species,
                 gene=rep(NA,length(species)), # #genes
                 geneInHog=rep(NA,length(species)), # #genes assigned to HOG
                 pctGeneInHog=rep(NA,length(species)), # %genes assigned to HOG
                 hogSpeciesContained=rep(NA,length(species)), # #HOG containing species 
                 pctHogSpecies=rep(NA,length(species)) # %HOG containing species
                 )
  res[,2:6]=
    t(sapply(species,
           function(sp){
             gene=system(paste("grep '>' ",orthofinder.in_dir,"/",sp,".faa | wc -l",sep=""),
                         intern=TRUE)
             gene=as.numeric(gene)
             geneInHog=sum(copyNumber[,sp])
             pctGeneInHog=geneInHog/gene
             hogSpeciesContained=length(copyNumber[,sp][copyNumber[,sp]!=0])
             pctHogSpecies=hogSpeciesContained/nrow(d)
             return(c(gene,geneInHog,pctGeneInHog,hogSpeciesContained,pctHogSpecies))
           }))
  write.table(res,paste(out_prefix,"_spStat.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
}
#N0.tsv="/bucket/BourguignonU/Cong/termite_pca/orthofinder/OrthoFinder/Results_Jul24_1/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
#orthofinder.in_dir="/bucket/BourguignonU/Cong/termite_pca/orthofinder"
#out_prefix="/bucket/BourguignonU/Cong/termite_pca/orthofinder/OrthoFinder/Results_Jul24_1/Phylogenetic_Hierarchical_Orthogroups/N0"
# OMA
# 不好用。。。

# Best protein-protein hit (lowest evalue, longest alignment) by blast 
# Dependencies: diamond, parallel (R)
best_blastp=function(query.faa=query.faa,
                     reference.faa=reference.faa,
                     out_dir=out_dir,
                     out_basename=out_basename,
                     threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  tmp=paste(out_basename,"_tmp_",sep="")
  system(paste("mkdir",tmp,sep=" "))
  
  cmd=paste("cp",query.faa,tmp,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  query.faa=paste(out_dir,"/",tmp,"/",basename(query.faa),sep="")
  
  cmd=paste("cp",reference.faa,tmp,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  reference.faa=paste(out_dir,"/",tmp,"/",basename(reference.faa),sep="")
  
  cmd=paste("diamond makedb",
            "--in",reference.faa,
            "--db",paste(reference.faa,".dmdb",sep=""),
            "--threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("diamond blastp",
            "--threads",threads,
            "--db",paste(reference.faa,".dmdb",sep=""),
            "--query",query.faa,
            "--out",paste(out_dir,"/",out_basename,".blast",sep=""),
            "--min-score 50",
            "--query-cover 75",
            "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm -r",tmp,sep=" "),wait=TRUE)
  
  blast=read.table(paste(out_dir,"/",out_basename,".blast",sep=""),
                   sep="\t",header=FALSE,quote="")
  colnames(blast)=c("qseqid","sseqid","pident","length","mismatch","gapopen",
                    "qstart","qend","sstart","send","evalue","bitscore",
                    "qcovhsp","scovhsp")
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
  write.table(blast,paste(out_dir,"/",out_basename,".best.blast.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
  setwd(wd)
}

bbh=function(blast1=blast1,blast2=blast2,
             # BLAST tab with header, only need qseqid & sseqid
             out_prefix=out_prefix){
  b1=read.table(blast1,header=TRUE,sep="\t",quote="")
  b1.pairs=paste(b1$qseqid,b1$sseqid,sep=" ")
  
  b2=read.table(blast2,header=TRUE,sep="\t",quote="")
  b2.pairs=paste(b2$sseqid,b2$qseqid,sep=" ")
  
  pairs=intersect(b1.pairs,b2.pairs)
  
  gene1=sapply(pairs,function(i){return(unlist(strsplit(i," "))[1])})
  gene2=sapply(pairs,function(i){return(unlist(strsplit(i," "))[2])})
  write.table(data.frame(gene1=gene1,gene2=gene2),
              paste(out_prefix,"_best_bidirectional_blast.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
}


# /flash/BourguignonU/Cong/HGT_Anna/blastBetweenLineages/
# Chitiniovibrionales.faa  Fibromonas.faa  Leadbettera.faa
#####
# Gene functions
#####
# Interproscan
# Dependencies: Interproscan, seqkit, parallel (R)
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
  
# Gene function by best blast hits (lowest evalue, longest alignment)
# Dependencies: DIAMOND+nr, parallel (R)
geneFun_bbh=function(pep.faa=pep.faa,
                     nr.dmdb="/apps/unit/BioinfoUgrp/DB/diamondDB/ncbi/2022-07/nr.dmnd",
                     out_prefix=out_prefix,
                     threads=threads){
  threads=as.character(threads)
  
  # blast search
  blast=paste(out_prefix,".blast",sep="")
  #if (!file.exists(blast)){ # the task here is to find one best hits
    cmd=paste("diamond blastp",
              "--threads",threads,
              "--db",nr.dmdb,
              "--query",pep.faa,
              "--out",blast,
              "--min-score 50",
              "--evalue 1e-5",
              "--query-cover 75",
              "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp stitle skingdoms sphylums sscinames staxids",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  #}
  
  blast=read.table(blast,sep="\t",header=FALSE,quote="",comment.char="")
  colnames(blast)=c("qseqid","sseqid","pident","length","mismatch","gapopen",
                    "qstart","qend","sstart","send","evalue","bitscore",
                    "qcovhsp","scovhsp","stitle","skingdoms","sphylums","sscinames","staxids")
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
  write.table(blast,paste(out_prefix,"_bbh.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
}

# Gene function by best blast hits (lowest evalue, longest alignment)
# Dependencies: blast+nr+NCBI Gene, parallel (R)
geneFun_blast2nr=function(pep.faa=pep.faa,
                          nr="/bucket/BourguignonU/Cong/public_db/ncbi.blastdb_20240708/nr",
                          gene2accession="/bucket/BourguignonU/Cong/public_db/ncbi_gene.20250610/gene2accession",
                          out_prefix=out_prefix,
                          threads=threads){
  threads=as.character(threads)
  
  # blast search
  blast=paste(out_prefix,".blast",sep="")
  cmd=paste("blastp",
            "-query",pep.faa,
            "-db",nr,
            "-out",blast,
            "-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid'",
            "-evalue 1e-5",
            "-num_alignments 20",
            "-num_threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  blast=read.table(blast,sep="\t",header=FALSE,quote="",comment.char="")
  colnames(blast)=c("qaccver","saccver","pident","length","mismatch",
                    "gapopen","qstart","qend","sstart","send","evalue","bitscore",
                    "staxid")
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
  blast=blast[blast[,"retainHit"],]
  
  blast[,"Symbol"]=parSapply(clus,
                             blast$saccver,
                             function(i){
                               cmd=paste("awk -F '\t' -v OFS='\t' -v i=1 ",
                                         "'{if ($6==\"",i,"\") print $16}'")
                               return(system(cmd,intern=TRUE))
                             })
  stopCluster(clus)
  
  
  write.table(blast,paste(out_prefix,"_bbh.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
}







# bitacora: gene family
# Dependencies: BITACORA, BLAST, HMMER
bitacora=function(mode_="full", # full/genome/protein
                  query_dir=query_dir, # Folder containing the query database or multiple databases, named as Example1_db.fasta and Example1_db.hmm (Mandatory)
                  genome.fna=genome.fna,
                  gff3=gff3,
                  prot.faa=prot.faa,
                  species="Out",
                  bitacora_scripts_folder="/bucket/BourguignonU/Cong/Softwares/bitacora/Scripts/",
                  GeMoMa.jar="/bucket/BourguignonU/Cong/Softwares/bitacora/GeMoMa-1.7.1/GeMoMa-1.7.1.jar",
                  eval=1e-3,
                  threads=threads,
                  out_dir=out_dir){
  if (!file.exists(out_dir)){dir.create(out_dir)}
  wd=getwd()
  setwd(out_dir)
  
  cmd=paste("runBITACORA_command_line.sh",
            "-m",mode_,
            "-q",query_dir,
            "-g",genome.fna,
            "-f",gff3,
            "-p",prot.faa,
            "-n",species,
            "-sp",bitacora_scripts_folder,
            "-gp",GeMoMa.jar,
            "-e",as.character(eval),
            "-t",as.character(threads),
            "-b T")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
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
  
#####
# Metabolic network based on KEGG
#####
# Assemble KOs into draft metabolic network
# Dependencies: KEGGREST (R), stringr(r)
KO2Network=function(gene2ko=gene2ko, # KOfamScan output, mapper format
                    kegg.compound="/bucket/BourguignonU/Cong/public_db/kofamscan/compound", # https://rest.kegg.jp/list/compound
                    kegg.reaction="/bucket/BourguignonU/Cong/public_db/kofamscan/reaction", # https://rest.kegg.jp/list/reaction
                    out_prefix=out_prefix){
  gene2ko=readLines(gene2ko)
  d=data.frame(geneID=rep(NA,length(gene2ko)),KO=rep(NA,length(gene2ko)))
  d[,c(1,2)]=t(sapply(gene2ko,
                      function(i){
                        i=unlist(strsplit(i,"\t"))
                        return(c(i[1],i[2]))
                      }))
  d=d[!is.na(d[,"KO"]),] 
  gene2ko=d # geneID,KO

  library(KEGGREST)
  ko.lst=gene2ko[,"KO"][!duplicated(gene2ko[,"KO"])]
  reactions=c()
  for (i in seq(0,length(ko.lst),100)){
    KOs=ko.lst[(i+1):(i+100)]
    KOs=KOs[!is.na(KOs)]
    if (length(KOs)!=0){
      reactions=c(reactions,keggLink("reaction",KOs))
    }
  }
  ko2reaction=data.frame(KO=names(reactions),reactionID=unname(reactions))
  ko2reaction[,"KO"]=sub("^ko:","",ko2reaction[,"KO"])
  ko2reaction[,"reactionID"]=sub("^rn:","",ko2reaction[,"reactionID"])
  
  gene.lst=gene2ko[,"geneID"][!duplicated(gene2ko[,"geneID"])]
  ko.lst=unname(sapply(gene.lst,
                       function(gene){return(paste(gene2ko[gene2ko[,"geneID"]==gene,"KO"],collapse=","))}))
  reaction.lst=unname(sapply(ko.lst,
                             function(ko){
                              ko=unlist(strsplit(ko,","))
                              return(paste(ko2reaction[ko2reaction[,"KO"] %in% ko,"reactionID"],collapse=","))
                      }))
  gene2reaction=data.frame(geneID=gene.lst,KO=ko.lst,reactionID=reaction.lst)
  gene2reaction=gene2reaction[gene2reaction[,"reactionID"]!="",]
  write.table(gene2reaction,
              paste(out_prefix,"_gene2reaction.tsv",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
  
  reaction.lst=ko2reaction[,"reactionID"][!duplicated(ko2reaction[,"reactionID"])]
  equations=rep(NA,length(reaction.lst))
  for (i in seq(0,length(reaction.lst),10)){
    reac=reaction.lst[(i+1):(i+10)]
    reac=reac[!is.na(reac)]
    if (length(reac)!=0){
      eq=keggGet(reac)
      for (j in 1:length(eq)){
        equations[j+i]=eq[[j]]$EQUATION
      }
    }
  }
  reactants=unname(sapply(equations,function(equ){return( unlist(strsplit(equ,"<=>"))[1] )}))
  products=unname(sapply(equations,function(equ){return( unlist(strsplit(equ,"<=>"))[2] )}))
  network=data.frame(reactantID=c(),reactantName=c(),
                     productID=c(),productName=c(),
                     reactionID=c(),reactionName=c(),
                     equation=c())
  for (i in 1:length(reaction.lst)){
    r=reaction.lst[i] # R13178
    rname=system(paste("awk -F '\t' -v OFS='\t' '{if ($1==\"",r,"\") print $2}' ",kegg.reaction,sep=""),intern=TRUE)
    if (length(rname)!=0){
      rname=unlist(strsplit(rname,";"))[1]
      e=equations[i]
      # print(c(r,rname,e))
      library(stringr)
      from=unlist(strsplit(reactants[i]," "));from=str_extract(from,"[CG][0-9]{5}");from=from[!is.na(from)]
      to=unlist(strsplit(products[i]," "));to=str_extract(to,"[CG][0-9]{5}");to=to[!is.na(to)]
      f=function(ids){
        id_c=ids[grepl("C[0-9]{5}",ids)]
        id_g=ids[!grepl("C[0-9]{5}",ids)]
        id_g2c=c()
        if (length(id_g)!=0){
          id_g2c=keggLink("compound",id_g);id_g2c=unname(sub("^cpd:","",id_g2c))
        }
        return(c(id_c,id_g2c))
      }
      from=f(from);to=f(to)
      d=merge(from,to);colnames(d)=c("reactantID","productID")
      d[,"reactantName"]=sapply(d[,"reactantID"],
                                function(id){
                                  return( unlist(strsplit(system(paste("awk -F '\t' -v OFS='\t' '{if ($1==\"",id,"\") print $2}' ",kegg.compound,sep=""),
                                                 intern=TRUE),";"))[1] )
                                })
      d[,"productName"]=sapply(d[,"productID"],
                               function(id){
                                 return( unlist(strsplit(system(paste("awk -F '\t' -v OFS='\t' '{if ($1==\"",id,"\") print $2}' ",kegg.compound,sep=""),
                                                                intern=TRUE),";"))[1] )
                               })
      d[,"reactionID"]=rep(r,nrow(d))
      d[,"reactionName"]=rep(rname,nrow(d))
      d[,"equation"]=rep(e,nrow(d))
      network=rbind(network,d)
    }
  }
  network[,"KO"]=sapply(network[,"reactionID"],
                        function(r){
                          return(paste(ko2reaction[ko2reaction[,"reactionID"]==r,"KO"],collapse=","))
                        })
  network[,"geneID"]=sapply(network[,"KO"],
                            function(ko){
                              ko=unlist(strsplit(ko,","))
                              return(paste(gene2ko[gene2ko[,"KO"] %in% ko,"geneID"],collapse = ","))
                            })
  write.table(network,paste(out_prefix,"_network.tsv",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
}

# Curate metabolic network
# Dependencies: KEGGREST (R)
curateNetwork=function(network.tsv=network.tsv, # from KO2Network
                       metaboliteConcentrations="/bucket/BourguignonU/Cong/public_db/Park2016/metConcentrationRanges.tsv",
                       # tsv with header
                       # Fields: compound (Compound ID of KEGG)
                       #         low (lower bound of absolute concentration in mol/L)
                       #         high (higher bound of absolute concentration in mol/L)
                       standardGibbs="/bucket/BourguignonU/Cong/public_db/eQuilibrator/standardGibbs.tsv",
                       # tsv with header
                       # Fields: reaction (Reaction ID of KEGG)
                       #         Gibbs (Standard Gibbs free energy of reaction in J/mol, temperature=298.15K, all reactants and products are of 1 mol/m^3)
                       out_prefix=out_prefix){
  # network.tsv="/Users/congliu/Desktop/PhD/Results/termite_pca/KO2Network/Nluj_network.tsv"
  # metaboliteConcentrations="~/Desktop/PhD/Results/public_db/Park2016/metConcentrationRanges.tsv"
  # standardGibbs="~/Desktop/PhD/Results/public_db/eQuilibrator/standardGibbs.tsv"
  # out_prefix="/Users/congliu/Desktop/PhD/Results/termite_pca/KO2Network/Nluj"
  # 
  network=read.table(network.tsv,sep="\t",header=TRUE,quote="")
  metaboliteConcentrations=read.table(metaboliteConcentrations,header=TRUE,sep="\t",quote ="")
  rownames(metaboliteConcentrations)=metaboliteConcentrations[,"compound"]
  metaboliteConcentrations=metaboliteConcentrations[metaboliteConcentrations$low!=0 & metaboliteConcentrations$high!=0,]
  standardGibbs=read.table(standardGibbs,header=TRUE,sep="\t",quote="")
  rownames(standardGibbs)=standardGibbs[,"reaction"]
  
  library(KEGGREST)
  compounds=c(network$reactantID,network$productID);compounds=compounds[!duplicated(compounds)]
  formula=rep(NA,length(compounds))
  for (i in seq(0,length(compounds),10)){
    compound=compounds[(i+1):(i+10)]
    compound=compound[!is.na(compound)]
    if (length(compound)!=0){
      f=keggGet(compound)
      for (j in 1:length(f)){
        if (!is.null(f[[j]]$FORMULA)){
          formula[j+i]=f[[j]]$FORMULA
        }
      }
    }
  }
  compound2formula=data.frame(compound=compounds,formula=formula)
  rownames(compound2formula)=compound2formula$compound
  network$reactantFormula=sapply(1:nrow(network),function(i){ return(compound2formula[network[i,"reactantID"],"formula"]) })
  network$productFormula=sapply(1:nrow(network),function(i){ return(compound2formula[network[i,"productID"],"formula"]) })
  
  # I. filter DNA/RNA/Protein reactions
  network=network[!(network[,"reactionID"] %in% 
                     network[network[,"reactantName"]==network[,"productName"],"reactionID"]),]
  network=network[!(network[,"reactionID"] %in% 
                      network[grepl("[DR]NA",network[,"reactantName"]) | 
                                grepl("[DR]NA",network[,"productName"]),"reactionID"]),]
  network=network[!(network[,"reactionID"] %in% 
                     network[grepl("[Pp]rotein",network[,"reactantName"]) | 
                               grepl("[Pp]rotein",network[,"productName"]),"reactionID"]),]
  ambigousCompounds=compound2formula[is.na(compound2formula$formula),"compound"]
  network=network[!(network[,"reactionID"] %in% 
                      network[network[,"reactantID"] %in% ambigousCompounds | 
                                network[,"productID"] %in% ambigousCompounds,"reactionID"]),]
  # II. filter compounds:
  # ambigously defined compounds
  RCompounds=compound2formula[grepl("R",compound2formula$formula),"compound"]
  network=network[!(network[,"reactantID"] %in% RCompounds) &
                    !(network[,"productID"] %in% RCompounds),]
  # water, hydrogen peroxide, oxygen, hydrogen
  network=network[network[,"reactantName"]!="H2O" & network[,"productName"]!="H2O",]
  network=network[network[,"reactantName"]!="Hydrogen peroxide" & network[,"productName"]!="Hydrogen peroxide",]
  network=network[network[,"reactantName"]!="Oxygen" & network[,"productName"]!="Oxygen",]
  network=network[network[,"reactantName"]!="H+" & network[,"productName"]!="H+",]
  # Phosphate
  network=network[network[,"reactantName"]!="Diphosphate" & network[,"productName"]!="Diphosphate",]
  network=network[network[,"reactantName"]!="Orthophosphate" & network[,"productName"]!="Orthophosphate",]
  network=network[network[,"reactantName"]!="Polyphosphate" & network[,"productName"]!="Polyphosphate",]
  # NTP, NMP, dNTP, dNMP
  NTP=c("ATP","TTP","GTP","CTP","UTP","Nucleoside triphosphate");dNTP=paste("d",NTP,sep="")
  NDP=c("ADP","TDP","GDP","CDP","UDP","NDP");dNDP=paste("d",NDP,sep="")
  # NMP=c("AMP","TMP","GMP","CMP","UMP");dNMP=paste("d",NMP,sep="")
  network=network[!(network[,"reactantName"] %in% c(NTP,dNTP,NDP,dNDP)) 
                  & !(network[,"productName"] %in% c(NTP,dNTP,NDP,dNDP)),]
  # H+/e- donors
  network=network[network[,"reactantName"]!="NAD+" & network[,"productName"]!="NAD+",]
  network=network[network[,"reactantName"]!="NADP+" & network[,"productName"]!="NADP+",]
  network=network[network[,"reactantName"]!="NADH" & network[,"productName"]!="NADH",]
  network=network[network[,"reactantName"]!="NADPH" & network[,"productName"]!="NADPH",]
  network=network[network[,"reactantName"]!="FAD" & network[,"productName"]!="FAD",]
  network=network[network[,"reactantName"]!="FADH2" & network[,"productName"]!="FADH2",]
  network=network[network[,"reactantName"]!="FMN" & network[,"productName"]!="FMN",]
  network=network[network[,"reactantName"]!="Reduced FMN" & network[,"productName"]!="Reduced FMN",]
  # CoA
  network=network[network[,"reactantName"]!="CoA" & network[,"productName"]!="CoA",]
  # III. direction of reactions
  reaction2equation=network[,c("reactionID","equation")]
  reaction2equation=reaction2equation[!duplicated(reaction2equation[,"reactionID"]),]
  R=8.3144598 # universal gas constatnt
  Temperature=298.15 # temperature
  reaction2equation[,"direction"]=
    sapply(1:nrow(reaction2equation),
           function(i){
             direction=0
             
             reac=reaction2equation[i,"reactionID"]
             stan.Gibbs=standardGibbs[reac,"Gibbs"]
             
             equa=reaction2equation[i,"equation"]
             equa=unlist(strsplit(equa,"<=>"))
             reactant=equa[1];product=equa[2]
             reactant=unlist(strsplit(reactant,"\\+"));product=unlist(strsplit(product,"\\+"))
             reactant=sub(" $","",sub("^ ","",reactant));product=sub(" $","",sub("^ ","",product))
             d.reactant=data.frame(reactant=rep(NA,length(reactant)),sto=rep(NA,length(reactant)),
                                   low=rep(NA,length(reactant)),high=rep(NA,length(reactant)))
             for (i in 1:nrow(d.reactant)){
               a=unlist(strsplit(reactant[i]," "))
               if (a[length(a)]!="C00001"){ # water
                 d.reactant[i,"reactant"]=a[length(a)]
                 if (a[1]!=a[length(a)]){d.reactant[i,"sto"]=as.numeric(a[1])}else{d.reactant[i,"sto"]=1}
                 l=metaboliteConcentrations[a[length(a)],"low"];h=metaboliteConcentrations[a[length(a)],"high"]
                 if (a[length(a)]=="C00080"){l=1e-7;h=1e-7} # H+: 1e-7 mol/L, pH=7
                 if (is.na(l)){l=1e-5;h=1e+5}else{l=l*1e-3/100;h=h*1e-3*100}
                 d.reactant[i,c("low","high")]=c(l,h)
               }
             }
             d.reactant=d.reactant[!is.na(d.reactant[,"reactant"]),]
             d.product=data.frame(product=rep(NA,length(product)),sto=rep(NA,length(product)),
                                  low=rep(NA,length(product)),high=rep(NA,length(product)))
             for (i in 1:nrow(d.product)){
               a=unlist(strsplit(product[i]," "))
               if (a[length(a)]!="C00001"){ # water
                 d.product[i,"product"]=a[length(a)]
                 if (a[1]!=a[length(a)]){d.product[i,"sto"]=as.numeric(a[1])}else{d.product[i,"sto"]=1}
                 l=metaboliteConcentrations[a[length(a)],"low"];h=metaboliteConcentrations[a[length(a)],"high"]
                 if (a[length(a)]=="C00080"){l=1e-7;h=1e-7} # H+
                 if (is.na(l)){l=1e-5;h=1e+5}else{l=l*1e-3/100;h=h*1e-3*100}
                 d.product[i,c("low","high")]=c(l,h)
               }
             }
             d.product=d.product[!is.na(d.product[,"product"]),]
             
             if (!is.na(stan.Gibbs)){
               Q.min=prod(d.product[,"low"]^d.product[,"sto"])/prod(d.reactant[,"high"]^d.reactant[,"sto"])
               Q.max=prod(d.product[,"high"]^d.product[,"sto"])/prod(d.reactant[,"low"]^d.reactant[,"sto"])
               deltaG.min=stan.Gibbs+R*Temperature*log(Q.min)
               deltaG.max=stan.Gibbs+R*Temperature*log(Q.max)
               if (deltaG.min>0){direction=-1}
               if (deltaG.max<0){direction=1}
               if (deltaG.min<0 & deltaG.max>0){direction=0}
             }
             # CO2 can only be synthesized
             if ("C00011" %in% d.product[,"product"]){direction=1}
             if ("C00011" %in% d.reactant[,"reactant"]){direction=-1}
             # O2 can only be consumed
             if ("C00007" %in% d.product[,"product"]){direction=-1}
             if ("C00007" %in% d.reactant[,"reactant"]){direction=1}
             # ATP can only be consumed
             # C00002, C00131, C00459, C00044, C00286, C00063, C00458, C00201
             if ("C00002" %in% d.product[,"product"]){direction=-1}
             if ("C00002" %in% d.reactant[,"reactant"]){direction=1}
             return(direction)
           })
  # direction:
  #   0: reversible
  #   1: left to right
  #   -1: right to left
  
  # Correct directions
  wrongDir=reaction2equation[reaction2equation[,"direction"]==-1,"reactionID"]
  wrongNetwork=network[network$reactionID %in% wrongDir,]
  wrongNetwork=wrongNetwork[,c("productID","reactantID","productName","reactantName",
                               "reactionID","reactionName","equation","KO","geneID",
                               "productFormula","reactantFormula")]
  colnames(wrongNetwork)=colnames(network)
  correctNetwork=network[!(network$reactionID %in% wrongDir),]
  network=rbind(wrongNetwork,correctNetwork)
  network[,"reversible"]=sapply(network[,"reactionID"],
                                function(ID){
                                  direction=reaction2equation[reaction2equation$reactionID==ID,"direction"]
                                  if (direction==0){return(TRUE)}else{return(FALSE)}
                                })
  write.table(network,paste(out_prefix,"_curateNetwork.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
}

#######
# Enrichment
#######
# topGO
topGO=function(Gene2GO.tsv=Gene2GO.tsv, # no header
               target.gene.lst=target.gene.lst, # one gene per line
               go_category="BP", # BP/MF/CC,
               out.tsv=out.tsv){
  mapfile=Gene2GO.tsv
  input=target.gene.lst
  
  library(topGO)
  gene_id=readMappings(file=mapfile,sep="\t",IDsep=",")
  gene_names=names(gene_id)
  my_genes=readLines(input)
  
  gene_list=rep(1,length(gene_id))
  names(gene_list)=names(gene_id)
  
  gene_list[match(my_genes,names(gene_list))]=0
  
  topGo_data=new("topGOdata",
                 nodeSize=5,
                 ontology=go_category,
                 allGenes=gene_list,
                 annot=annFUN.gene2GO,
                 gene2GO=gene_id,
                 geneSel=function(allScore){return(allScore==0)}) 
  result_KS.elim=runTest(topGo_data,
                         algorithm="elim",
                         statistic="ks")
  
  allres=GenTable(topGo_data,
                  KS=result_KS.elim,
                  ranksOf="classic",
                  topNodes=attributes(result_KS.elim)$geneData[4])
  library(GO.db)
  allres$term.full=Term(allres$GO.ID)
  
  write.table(allres,out.tsv,sep="\t",row.names=FALSE,quote=FALSE)
}

# genome=read.table("~/quality_genomePeptide.tsv",sep="\t",header=TRUE,quote="")
# for (sp in genome$Label){
#   df=paste("/bucket/BourguignonU/Cong/termite_pca/geneFunction/sum_geneFun/",
#            sp,"_geneFunction.tsv",sep="")
#   df=read.table(df,header=TRUE,sep="\t",quote="")
#   GO_Universe=df[,c("Gene","GO_eggNOG")]
#   write.table(GO_Universe,
#               paste("/flash/BourguignonU/Cong/termite_pca/paralogs/topGO/",sp,"_GO_Universe.tsv",sep=""),
#               sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
#   writeLines(df[df$HOG.category=="single-copy","Gene"],
#              paste("/flash/BourguignonU/Cong/termite_pca/paralogs/topGO/",sp,"_singlecopy.lst",sep=""))
#   writeLines(df[df$HOG.category=="paralogous","Gene"],
#              paste("/flash/BourguignonU/Cong/termite_pca/paralogs/topGO/",sp,"_paralogous.lst",sep=""))
#   for (go_category in c("BP","MF","CC")){
#     topGO(Gene2GO.tsv=paste("/flash/BourguignonU/Cong/termite_pca/paralogs/topGO/",sp,"_GO_Universe.tsv",sep=""), # no header
#           target.gene.lst=paste("/flash/BourguignonU/Cong/termite_pca/paralogs/topGO/",sp,"_singlecopy.lst",sep=""), # one gene per line
#           go_category=go_category, # BP/MF/CC,
#           out.tsv=paste("/flash/BourguignonU/Cong/termite_pca/paralogs/topGO/",sp,"_singlecopy_",go_category,".tsv",sep=""))
#     
#     topGO(Gene2GO.tsv=paste("/flash/BourguignonU/Cong/termite_pca/paralogs/topGO/",sp,"_GO_Universe.tsv",sep=""), # no header
#           target.gene.lst=paste("/flash/BourguignonU/Cong/termite_pca/paralogs/topGO/",sp,"_paralogous.lst",sep=""), # one gene per line
#           go_category=go_category, # BP/MF/CC,
#           out.tsv=paste("/flash/BourguignonU/Cong/termite_pca/paralogs/topGO/",sp,"_paralogous_",go_category,".tsv",sep=""))
#     
#   }
#   
# }










# Gene function enrichment
# Dependencies: clusterProfiler (R)
enrichment=function(genes=genes, # a vector of gene id
                    gene2term.tsv=gene2term.tsv, 
                    # with header and fields:
                    # GeneID: gene id
                    # Terms: comma-list or ""
                    out.tsv=out.tsv){
  hgtGene=read.table("~/Desktop/PhD/Results/termite_pca/AI/hgtGene.tsv",
                     sep="\t",header=TRUE,quote="")
  hgtGene=hgtGene[hgtGene$Donor!="0",]
  genes=hgtGene$HOG
  genes=genes[!duplicated(genes)]
  gene2term.tsv="~/Desktop/PhD/Results/termite_pca/hogFunction/hog2go.tsv"
  
  gene2term=read.table(gene2term.tsv,sep="\t",header=TRUE,quote="");colnames(gene2term)=c("gene","term")
  universe=gene2term$gene
  gene2term=gene2term[gene2term$term!="",]
  gene2term=do.call(rbind.data.frame,
                    lapply(1:nrow(gene2term),
                           function(idx){
                             data.frame(
                             gene = gene2term[idx, "gene"],
                             term = strsplit(gene2term[idx, "term"], ",")[[1]],
                             stringsAsFactors = FALSE)}))
  term2gene=gene2term[,c("term","gene")]
  
  library(clusterProfiler)
  a=enricher(gene=genes,
             universe=universe,
             TERM2GENE=term2gene)
  write.table(a@result,out.tsv,sep="\t",row.names=FALSE,quote=FALSE)
}

# Remove parent GO IDs if its child appears
# Dependencies: annotate (R)
goTrim=function(GO.IDs=GO.IDs # vector of GO terms
                ){
  #GO.IDs="GO:0000003,GO:0003674,GO:0003824,GO:0004653,GO:0005575,GO:0005622,GO:0005623,GO:0005737,GO:0005794,GO:0005795,GO:0006464,GO:0006486,GO:0006493,GO:0006807,GO:0008150,GO:0008152,GO:0008194,GO:0008376,GO:0009058,GO:0009059,GO:0009100,GO:0009101,GO:0009987,GO:0012505,GO:0016020,GO:0016021,GO:0016740,GO:0016757,GO:0016758,GO:0018193,GO:0018210,GO:0018243,GO:0019538,GO:0031224,GO:0031984,GO:0031985,GO:0032501,GO:0032504,GO:0034645,GO:0036211,GO:0043170,GO:0043226,GO:0043227,GO:0043229,GO:0043231,GO:0043412,GO:0043413,GO:0044237,GO:0044238,GO:0044249,GO:0044260,GO:0044267,GO:0044422,GO:0044424,GO:0044425,GO:0044431,GO:0044444,GO:0044446,GO:0044464,GO:0070085,GO:0071704,GO:0098791,GO:0140096,GO:1901135,GO:1901137,GO:1901564,GO:1901566,GO:1901576"
  #GO.IDs=unlist(strsplit(GO.IDs,","))
  library(annotate)
  parentsGO = sapply(GO.IDs,
                     function(ID){
                       parents=try(getGOParents(ID),silent=TRUE)[[1]][2][[1]]
                       return(unname(parents))
                     })
  parentsGO = unlist(parentsGO)
  parentsGO = parentsGO[!duplicated(parentsGO)]
  trimmedGO.IDs=setdiff(GO.IDs,parentsGO)
  return(trimmedGO.IDs)
}

# Term2Desc=getGOTerm(o$Term)
# BP=Term2Desc$BP;BP_id=names(BP);BP_desc=unname(BP)
# category=rep("Biological Processes",length(BP))
# BP = data.frame(Term=BP_id,Description=BP_desc,Category=category)
# MF=Term2Desc$MF;MF_id=names(MF);MF_desc=unname(MF)
# category=rep("Molecular Functions",length(MF))
# MF = data.frame(Term=MF_id,Description=MF_desc,Category=category)
# CC=Term2Desc$CC;CC_id=names(CC);CC_desc=unname(CC)
# category=rep("Cellular Components",length(CC))
# CC = data.frame(Term=CC_id,Description=CC_desc,Category=category)
# Term2Desc = rbind(BP,MF,CC)
# p=merge(o,Term2Desc,by.x="Term",by.y="Term")
# p=subset(p,!is.na(p$Description))

#network=network[!grepl("G[0-9]{5}",network[,"equation"]),]
#network=network[!grepl("C[0-9]{5}\\(n\\)",network[,"equation"]),]
#Alkyl=network[grepl("[Aa]lkyl",network[,"reactantName"]) | grepl("[Aa]lkyl",network[,"productName"]),"reactionID"]
# # 金属
# network=network[network[,"reactantName"]!="Fe2+" & network[,"productName"]!="Fe2+",] # Fe2+
# network=network[network[,"reactantName"]!="Fe3+" & network[,"productName"]!="Fe3+",] # Fe3+
# network=network[network[,"reactantName"]!="Zinc cation" & network[,"productName"]!="Zinc cation",] # Zn2+
# network=network[network[,"reactantName"]!="Calcium cation" & network[,"productName"]!="Calcium cation",] # Ca2+
# network=network[network[,"reactantName"]!="Cobalt ion" & network[,"productName"]!="Cobalt ion",] # Co2+
# network=network[network[,"reactantName"]!="Potassium cation" & network[,"productName"]!="Potassium cation",] # K+
# network=network[network[,"reactantName"]!="Magnesium cation" & network[,"productName"]!="Magnesium cation",] # Mg2+
# network=network[network[,"reactantName"]!="Mercury(2+)" & network[,"productName"]!="Mercury(2+)",] # Hg+
# network=network[network[,"reactantName"]!="Sodium cation" & network[,"productName"]!="Sodium cation",] # Na+
# network=network[network[,"reactantName"]!="Barium cation" & network[,"productName"]!="Barium cation",] # Ba2+
# network=network[network[,"reactantName"]!="Cu2+" & network[,"productName"]!="Cu2+",] # Cu2+
# network=network[network[,"reactantName"]!="Cu+" & network[,"productName"]!="Cu+",] # Cu+
# network=network[network[,"reactantName"]!="Manganese(3+)" & network[,"productName"]!="Manganese(3+)",] # Mn3+
# network=network[network[,"reactantName"]!="Manganese(2+)" & network[,"productName"]!="Manganese(2+)",] # Mn2+
# network=network[network[,"reactantName"]!="Nickel(2+)" & network[,"productName"]!="Nickel(2+)",] # Ni2+
# network=network[network[,"reactantName"]!="Rubidium cation" & network[,"productName"]!="Rubidium cation",] # Rb+
# network=network[network[,"reactantName"]!="Strontium cation" & network[,"productName"]!="Strontium cation",] # Sr2+