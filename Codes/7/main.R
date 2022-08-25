# Downstream analysis in genomics

# Retain the longest protein isoform of each gene by header lines of fasta.
# Proteins are considered as isoforms only when the header line declares isoform.
# e.g. these are the same gene as they declare isoform and share same name "pyruvate decarboxylase 2-like"
# >XP_024356342.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
# >XP_024356343.1 pyruvate decarboxylase 2-like isoform X1 [Physcomitrella patens]
# these are different genes:
# >XP_024390255.1 40S ribosomal protein S12-like [Physcomitrella patens]
# >XP_024399722.1 40S ribosomal protein S12-like [Physcomitrella patens]
# Dependencies: seqkit
Isoform_filter=function(faa_in=faa_in,
                        faa_out=faa_out,
                        threads=threads){
  threads=as.character(threads)
  tmp=paste(faa_out,"_tmp",sep="")
  system(paste("mkdir",tmp,sep=" "))
  
  # Sequence header (everything except >) to length (tabular)
  cmd=paste("seqkit","fx2tab",
            "--threads",threads,
            "--length --name --header-line",
            faa_in,
            ">",
            paste(tmp,"/Header2Length",sep=""))
  system(cmd,wait=TRUE)
  
  Header2Length=read.table(paste(tmp,"/Header2Length",sep=""),sep="\t",header=FALSE,quote="",
                           col.names=c("Header","Length"))
  Header2Length[,"ID"]=sapply(Header2Length[,"Header"],
                              function(Header){
                                Header=unlist(strsplit(Header," "))
                                return(Header[1])
                              })

  # Select longest isoform
  Isoforms=Header2Length[grepl("[Ii]soform",Header2Length[,"Header"]),] # Protein records that declare isoform
  Isoforms[,"Header"]=sub("[Ii]soform [0-9, A-Z]*","",Isoforms[,"Header"])
  if (nrow(Isoforms)!=0){
  Isoforms[,"Header"]=sapply(1:nrow(Isoforms),
                             function(i){ # convert header to name
                               return(sub(Isoforms[i,"ID"],"",Isoforms[i,"Header"]))
                             })
  IDs_LongestIso=sapply(names(table(Isoforms[,"Header"])),
                        function(Header){
                          d=Isoforms[Isoforms[,"Header"]==Header,]
                          d=d[d[,"Length"]==max(d[,"Length"]),]
                          ID=sample(d[,"ID"],1)
                          return(ID)
                        })
  IDs_LongestIso=unname(IDs_LongestIso)
  }else{IDs_LongestIso=c()}
  IDs_noIso=Header2Length[!grepl("[Ii]soform",Header2Length[,"Header"]),][,"ID"] # Protein IDs that do not declare isoform
  
  target_IDs=c(IDs_noIso,IDs_LongestIso)
  target_IDs=data.frame(target_IDs)
  write.table(target_IDs,paste(tmp,"/target_IDs",sep=""),
              sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  # Extract fasta by IDs
  cmd=paste("seqkit","grep",
            "--threads",threads,
            "-f",paste(tmp,"/target_IDs",sep=""),
            faa_in,
            ">",
            faa_out,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm","-r",tmp,sep=" "),wait=TRUE)
}

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

# Incomplete
# AvP
# Detect horizontal gene transfer (HGT)
# Dependencies: AvP
Avp=function(query.faa=query.faa,
             blast_filename=blast_filename,
             dmbd=dmbd,
             group.yaml=group.yaml,
             config.yaml=config.yaml,
             out_dir=out_dir,
             threads=threads
){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  threads=as.character(threads)
  wd=getwd()
  setwd(out_dir)
  
  if (!file.exists(blast_filename)){
    cmd=paste("diamond blastp",
              "--threads",threads,
              "-d",dmbd,
              "--max-target-seqs 500",
              "-q",query.faa,
              "--evalue 1e-5",
              "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids",
              ">",
              blast_filename)
    print(cmd);system(cmd,wait=TRUE)
  }
  cmd=paste("calculate_ai.py",
            "-i",blast_filename,
            "-x",group.yaml,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("avp prepare",
            "-a",paste(blast_filename,"_ai.out",sep=""),
            "-o",".",
            "-f",query.faa,
            "-b",blast_filename,
            "-x",group.yaml,
            "-c",config.yaml,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("avp detect",
            "-i","./mafftgroups/",
            "-o",".",
            "-g","./groups.tsv",
            "-t","./tmp/taxonomy_nexus.txt",
            "-c",config.yaml,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}









