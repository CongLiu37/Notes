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
  genes=genes[which(copyNum>4)] # only retain busco genes present in at least 5 species
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
  astral(trees=paste(system(paste("ls ",out_dir,"/iqtree/trees/*.contree",sep=""),intern=TRUE),collapse=","), # comma-list, nwk trees
         path2astral="/home/c/c-liu/Softwares/ASTRAL/astral.5.7.8.jar",
         out_prefix=paste(out_dir,"/astral/Phylo",sep=""))
  system(paste("touch",paste(out_dir,"/astral/astral.finished",sep=""),sep=" "))
}

#####
# supermatrix
#####
if (!file.exists(paste(out_dir,"/supermatrix",sep=""))){system(paste("mkdir ",out_dir,"/supermatrix",sep=""))}
genes=system(paste("ls ",out_dir,"/trimAL/seqs/*",sep=""),intern=TRUE)
copyNum=sapply(1:length(genes),
               function(i){
                 cmd=paste("grep '>' ",genes[i]," | wc -l",sep="")
                 return(system(cmd,intern=TRUE))
               })
copyNum=as.numeric(copyNum)
spNum=system(paste("wc -l",tab,sep=" "),intern=TRUE)
spNum=as.numeric(unlist(strsplit(spNum," "))[1])
sp.lst=system(paste("awk '{print $1}' ",tab,sep=""),intern=TRUE)
for (percent in seq(0.85,1,0.025)){ #
  if (!file.exists(paste(out_dir,"/supermatrix/",as.numeric(percent),sep=""))){
    system(paste("mkdir ",out_dir,"/supermatrix/",as.numeric(percent),sep=""))}
  stamp=paste("mkdir ",out_dir,"/supermatrix/",as.numeric(percent),".finished",sep="")
  if (!file.exists(stamp)){
    gene.lst=genes[which(copyNum>=floor(percent*spNum))]
    if (length(gene.lst)>=75){
      sapply(gene.lst,
             function(gene){
               system(paste("seqkit seq -w0 ",gene," > ",out_dir,"/supermatrix/",as.numeric(percent),"/",basename(gene),sep=""))
               length=as.numeric(system(paste("seqkit fx2tab ",gene," -l -n -i -H | awk '{print $2}'",sep=""),intern=TRUE)[2])
               for (sp in sp.lst){
                 if (paste(">",sp,sep="") %in% 
                     system(paste("grep '>' ",out_dir,"/supermatrix/",as.numeric(percent),"/",basename(gene),sep=""),intern=TRUE)){
                }else{
                 system(paste("echo '>",sp,"' >> ",out_dir,"/supermatrix/",as.numeric(percent),"/",basename(gene),sep=""))
                 system(paste("echo '",paste(rep("-",length),collapse=""),"' >> ",out_dir,"/supermatrix/",as.numeric(percent),"/",basename(gene),sep=""))
                }
               }
             })
      catMSA(align.fa=paste(system(paste("ls ",out_dir,"/supermatrix/",as.numeric(percent),"/*",sep=""),
                                   intern=TRUE),collapse=","), # comma list
             cat.fa=paste(out_dir,"/supermatrix/",as.numeric(percent),"/catMSA_raw.faa",sep=""))
      cmd=paste("seqkit grep -r -p [A-Z] -s ",
                out_dir,"/supermatrix/",as.numeric(percent),"/catMSA_raw.faa > ",
                out_dir,"/supermatrix/",as.numeric(percent),"/catMSA.faa",sep="")
      system(cmd,wait=TRUE)
      iqtree(msa.fa=paste(out_dir,"/supermatrix/",as.numeric(percent),"/catMSA.faa",sep=""), 
             type="protein", # dna/protein
             out_prefix=paste(out_dir,"/supermatrix/",as.numeric(percent),"/",as.numeric(percent),"_iqtree",sep=""),
             threads=threads)
    }
  }
  system(paste("touch",stamp,sep=" "))
}

