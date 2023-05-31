arg=commandArgs(trailingOnly = TRUE)

main=arg[1]
conf=arg[2]

source(main)
source(conf)

threads=as.character(threads)
out_dir=sub("/$","",out_dir)
if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}

#####
# getBuscoSeq
#####
if (!file.exists(paste(out_dir,"/getBuscoSeq",sep=""))){system(paste("mkdir ",out_dir,"/getBuscoSeq",sep=""))}
if (file.exists(paste(out_dir,"/getBuscoSeq/getBuscoSeq.finished",sep=""))){
  print("Step 1: getBuscoSeq FINISHED")
}else{
  print("Step 1: getBuscoSeq START")
  getBuscoSeq(tab=tab, # tsv. 
                       # First column is species label and second column is run_/busco_sequences/single_copy_busco_sequences/
              out_dir=paste(out_dir,"/getBuscoSeq/seqs",sep=""))
  system(paste("touch",paste(out_dir,"/getBuscoSeq/getBuscoSeq.finished",sep=""),sep=" "))
}

#####
# mafft
#####
if (!file.exists(paste(out_dir,"/mafft",sep=""))){system(paste("mkdir ",out_dir,"/mafft",sep=""))}
if (!file.exists(paste(out_dir,"/mafft/seqs",sep=""))){system(paste("mkdir ",out_dir,"/mafft/seqs",sep=""))}
if (file.exists(paste(out_dir,"/mafft/mafft.finished",sep=""))){
  print("Step 2: mafft FINISHED")
}else{
  print("Step 2: mafft START")
  genes=system(paste("ls ",out_dir,"/getBuscoSeq/seqs/",sep=""),intern=TRUE)
  genes=paste(out_dir,"/getBuscoSeq/seqs/",genes,sep="")
  copyNum=sapply(1:length(genes),
                 function(i){
                   cmd=paste("grep '>' ",genes[i]," | wc -l",sep="")
                   return(system(cmd,intern=TRUE))
                 })
  copyNum=as.numeric(copyNum)
  genes=genes[which(copyNum>4)]
  for (gene in genes){
    ID=unlist(strsplit(gene,"/"));ID=ID[length(ID)]
    mafft(in.fa=gene,
          align.fa=paste(out_dir,"/mafft/seqs/",ID,sep=""),
          threads=threads)
  }
  system(paste("touch",paste(out_dir,"/mafft/mafft.finished",sep=""),sep=" "))
}

#####
# trimAL
#####
if (!file.exists(paste(out_dir,"/trimAL",sep=""))){system(paste("mkdir ",out_dir,"/trimAL",sep=""))}
if (!file.exists(paste(out_dir,"/trimAL/seqs",sep=""))){system(paste("mkdir ",out_dir,"/trimAL/seqs",sep=""))}
if (file.exists(paste(out_dir,"/trimAL/trimAL.finished",sep=""))){
  print("Step 3: trimAL FINISHED")
}else{
  print("Step 3: trimAL START")
  genes=system(paste("ls ",out_dir,"/getBuscoSeq/seqs/",sep=""),intern=TRUE)
  for (gene in genes){
    trimAL(inMSA.fa=paste(out_dir,"/mafft/seqs/",gene,sep=""),
           outMSA.fa=paste(out_dir,"/trimAL/seqs/",gene,sep=""))
  }
  system(paste("touch",paste(out_dir,"/trimAL/trimAL.finished",sep=""),sep=" "))
}

#####
# iqtree
#####
if (!file.exists(paste(out_dir,"/iqtree",sep=""))){system(paste("mkdir ",out_dir,"/iqtree",sep=""))}
if (!file.exists(paste(out_dir,"/iqtree/trees",sep=""))){system(paste("mkdir ",out_dir,"/iqtree/trees",sep=""))}
if (file.exists(paste(out_dir,"/iqtree/iqtree.finished",sep=""))){
  print("Step 4: iqtree FINISHED")
}else{
  print("Step 4: iqtree START")
  genes=system(paste("ls ",out_dir,"/trimAL/seqs/",sep=""),intern=TRUE)
  genes=paste(out_dir,"/trimAL/seqs/",genes,sep="")
  library(parallel)
  index=splitIndices(length(genes),as.numeric(threads)-1)
  clus=makeCluster(as.numeric(threads))
  clusterExport(clus, c("genes","index","iqtree","out_dir"), envir = .GlobalEnv)
  parSapply(clus,
            1:(as.numeric(threads)-1),
            function(i){
              for (gene in genes[index[[i]]]){
                ID=unlist(strsplit(gene,"/"));ID=ID[length(ID)]
                iqtree(msa.fa=gene, 
                       type="protein",
                       out_prefix=paste(out_dir,"/iqtree/trees/",ID,sep=""),
                       threads=1)
              }
            })
  stopCluster(clus)
  system(paste("touch",paste(out_dir,"/iqtree/iqtree.finished",sep=""),sep=" "))
}

#####
# ASTRAL
#####
if (!file.exists(paste(out_dir,"/astral",sep=""))){system(paste("mkdir ",out_dir,"/astral",sep=""))}
if (file.exists(paste(out_dir,"/astral/astral.finished",sep=""))){
  print("Step 5: astral FINISHED")
}else{
  print("Step 5: astral START")
  astral(trees=paste(system(paste("ls ",out_dir,"/iqtree/trees/*.treefile",sep=""),intern=TRUE),collapse=","), # comma-list, nwk trees
         path2astral="/home/c/c-liu/Softwares/ASTRAL/astral.5.7.8.jar",
         out_prefix=paste(out_dir,"/astral/Phylo",sep=""))
  system(paste("touch",paste(out_dir,"/astral/astral.finished",sep=""),sep=" "))
}