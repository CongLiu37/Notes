# Downstream analysis of gene prediction

# Retain the longest protein isoform
# Dependencies: seqkit, parrell (R)
Isoform_filter=function(faa_in=faa_in,
                        faa_out=faa_out,
                        threads=threads){
  threads=as.character(threads)
  tmp=paste(faa_out,"_tmp",sep="")
  system(paste("mkdir",tmp,sep=" "))
  
  # remove line breakers in sequences
  cmd=paste("seqkit","seq","--threads",threads,"-w 0",faa_in,
            ">",
            paste(tmp,"/no_wrap.faa",sep=""),
            sep=" ")
  system(cmd,wait=TRUE)
  
  # Sequence header (everything except >) to length (tabular)
  cmd=paste("seqkit","fx2tab","--threads",threads,"--length --name --header-line",faa_in,
            ">",
            paste(tmp,"/Header2Length",sep=""))
  system(cmd,wait=TRUE)
  
  Header2Length=read.table(paste(tmp,"/Header2Length",sep=""),
                           sep="\t",header=FALSE,quote="",
                           col.names=c("Header","Length"))
  Header2Length[,"ID"]=sapply(Header2Length[,"Header"],
                              function(Header){
                                Header=unlist(strsplit(Header," "))
                                return(Header[1])
                              })
  
  # Select longest isoform
  Isoforms=Header2Length[grepl("[Ii]soform",Header2Length[,"Header"]),]
  Isoforms[,"Header"]=sub("[Ii]soform [0-9, A-Z]*","",Isoforms[,"Header"])
  if (nrow(Isoforms)!=0){
  Isoforms[,"Header"]=sapply(1:nrow(Isoforms),
                             function(i){
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
  IDs_noIso=Header2Length[!grepl("[Ii]soform",Header2Length[,"Header"]),][,"ID"]
  target_IDs=c(IDs_noIso,IDs_LongestIso)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads)-1)
  clusterExport(cl=clus,c("target_IDs","tmp"),envir=environment())
  parSapply(cl=clus,target_IDs,
         function(ID){
           cmd=paste("grep","-A1",ID,paste(tmp,"/no_wrap.faa",sep=""),
                     ">>",
                     paste(faa_out,";",sep=""),sep=" ")
           system(paste("echo",paste("'",cmd,"'",sep=""),
                        ">>",
                        paste(tmp,"/commands.sh",sep=""),sep=" "))
         })
  system(paste("bash",paste(tmp,"/commands.sh",sep=""),sep=" "),wait=TRUE)
  #system(paste("rm","-r",tmp,sep=" "),wait=TRUE)
}

# species=readLines("/flash/HusnikU/Cong/Orthofinder_noIsoforms/SpeciesList.txt")
# for (line in species){
#   name=unlist(strsplit(line,"/"))
#   name=name[length(name)]
#   faa_in=line
#   faa_out=paste("/flash/HusnikU/Cong/Orthofinder_noIsoforms/",name,sep="")
#   print(faa_out)
#   Isoform_filter(faa_in=faa_in,faa_out=faa_out,threads=64)
# }

# OrthoFinder
# By default OrthoFinder creates a results directory called ‘OrthoFinder’ inside the 
# input proteome directory and puts the results here.
# Dependencies: Orthofinder, mafft, fasttree
Orthofinder=function(in_dir=in_dir, # Input proteome directory. 
                     # One file per species with extension '.faa'
                     threads=threads){
  threads=as.character(threads)
  
  cmd=paste("orthofinder.py",
            "-f",in_dir,
            "-M msa",
            "-t",threads,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# AvP
# Detect horizontal gene transfer (HGT)
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











