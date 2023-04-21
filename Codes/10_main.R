# Gene functions & metabolic network

#####
# Gene function annotation
#####
# OrthoFinder
# OrthoFinder creates a results directory called ‘OrthoFinder’ inside the input proteome directory and puts the results here.
# Orthofinder with "-M msa" infers phylogeny by mafft & fasttree
# Dependencies: Orthofinder, mafft, fasttree
Orthofinder=function(in_dir=in_dir, # Input proteome directory. 
                     # One file per species with extension '.faa'
                     threads=threads){
  threads=as.character(threads)
  
  cmd=paste("orthofinder.py",
            "-f",in_dir,
            "-M msa",
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
  
  cmd=paste("orthofinder",
            "-t",threads,
            "-a",threads,
            "-ft",previous_orthofinder_result_dir,
            "-s",tree,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

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

# Interproscan
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
  
  cmd=paste("interproscan.sh",
            "-dp",
            "-b",out_basename,
            "-cpu",threads,
            "-f TSV",
            "-goterms",
            "-i",proteins.faa,
            "--pathways",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",proteins.faa,sep=" "),wait=TRUE)
  setwd(wd)
}
