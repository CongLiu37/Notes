# Collect single-copy busco sequences from busco runs
# Dependencies:
getBuscoSeq=function(tab=tab, # tsv. 
                     # First column is species label and second column is run_/busco_sequences/single_copy_busco_sequences/
                     out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  
  tab=read.table(tab,sep="\t",header=FALSE,quote="")
  for (i in 1:nrow(tab)){
    sp=tab[i,1]
    seqs=system(paste("ls",tab[i,2],sep=" "),intern=TRUE)
    seqs=paste(tab[i,2],"/",seqs,sep="")
    for (j in seqs){
      out=unlist(strsplit(j,"/"));out=paste(out_dir,"/",out[length(out)],sep="")
      sed=paste("'s/>.*/>",sp,"/'",sep="")
      cmd=paste("sed",sed,j,">>",out,sep=" ")
      system(cmd,wait=TRUE)
    }
  }
}

# MAFFT: multiple sequence alignment
# Dependencies: MAFFT
mafft=function(in.fa=in.faa,
               align.fa=align.fa,
               threads=threads){
  threads=as.character(threads)
  cmd=paste("mafft",
            "--auto",
            "--thread",threads,
            in.fa,">",
            align.fa,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# trimAL: curate multiple sequence alignments for phylogenetic analysis
# Dependencies: trimAL, seqkit
trimAL=function(inMSA.fa=inMSA.fa,
                outMSA.fa=outMSA.fa){
  cmd=paste("trimal",
            "-in",inMSA.fa,
            "-automated1",
            "-keepheader",
            "| seqkit seq -u",
            ">",outMSA.fa,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Phylogenetic tree
# Dependencies: iqtree
iqtree=function(msa.fa=msa.fa, 
                type="dna", # dna/protein
                out_prefix=out_prefix,
                threads=threads){
  threads=as.character(threads)
  cmd=paste("iqtree",
            "-s",msa.fa,
            "--prefix",out_prefix,
            "-T",threads,
            # "-m TEST", # Standard model selection followed by tree inference
            # "-m MFP", # Extended model selection followed by tree inference
            # "-mset GTR", # test GTR+... models only
            # "--msub nuclear", # Amino-acid model source (nuclear, mitochondrial, chloroplast or viral)
            "-B 1000",
            sep=" ")
  if (type=="dna"){cmd=paste(cmd,"-mset GTR",sep=" ")}
  if (type=="protein"){cmd=paste(cmd,"--msub nuclear",sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
}

# ASTRAL: supertree
# Dependencies: ASTRAL
astral=function(trees=trees, # comma-list, nwk trees
                path2astral="/home/c/c-liu/Softwares/ASTRAL/astral.5.7.8.jar",
                out_prefix=out_prefix){
  trees=unlist(strsplit(trees,","))
  for (tree in trees){
    cmd=paste("cat ",tree," >> ",out_prefix,"_inTrees.nwk",sep="")
    system(cmd,wait=TRUE)
  }
  cmd=paste("java -jar",path2astral,
            "-i",paste(out_prefix,"_inTrees.nwk",sep=""),
            "-o",paste(out_prefix,"_astralTree.nwk",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}