# PseudoPipe for pseudogene identification
# Dependencies: seqkit,PseudoPipe (modified)
PseudoPipe=function(genome=genome, # soft masked
                    protein.fa=protein.fa, # space list
                    gff=gff,
                    out_dir=out_dir,
                    threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  system("mkdir ./input")
  
  system("mkdir ./input/dna")
  cmd=paste("cp",genome,"./input/dna/masked.fa",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  # cmd=paste("seqkit","seq","-u",
  #           "-j",threads,
  #           "./input/dna/masked.fa > ./input/dna/unmasked.fa",
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  cmd=paste("seqkit","split2","./input/dna/masked.fa","-s 1","-j",threads," 2> lst",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd="grep 'write 1 sequences to file:' lst"
  split.fa=system(cmd,intern=TRUE)
  split.fa=sub("^.*write 1 sequences to file: ","",split.fa)
  
  split=sapply(1:length(split.fa),
               function(i){
                 seqID=system(paste("head -n1",split.fa[i],sep=" "),intern=TRUE)
                 seqID=sub(">","",seqID)
                 seqID=sub(" .*$","",seqID)
                 
                 cmd=paste("mv",split.fa[i],
                           paste("./input/dna/masked.",seqID,".fa",sep=""),
                           sep=" ")
                 system(cmd,wait=TRUE)
                 return(paste("./input/dna/masked.",seqID,".fa",sep=""))
               })
  system("rm -r ./input/dna/masked.fa.split/")
  system("rm lst")
  
  system("mkdir ./input/pep")
  cmd=paste("cat",protein.fa,">","./input/pep/protein.fa",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("mkdir ./input/mysql")
  cmd=paste("grep","'exon'",gff,"|",
            "awk","-F '\t' -v OFS='\t' '{print $1,$2,$4,$5}'",
            ">","./input/mysql/gff",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  sapply(1:length(split),
         function(i){
           seqID=system(paste("head -n1",split[i],sep=" "),intern=TRUE)
           seqID=sub(">","",seqID)
           seqID=sub(" .*$","",seqID)
           cmd=paste("grep"," ","'",seqID,"'"," ","./input/mysql/gff",
                     " > ","./input/mysql/gff.",seqID,".tsv",
                     sep="")
           system(cmd,wait=TRUE)
         })
  #########################
  # python2
  system("module load sango-legacy-modules")
  system("module load python/2.7.3")
  system("module load fasta/35.4.12")
  #########################
  cmd=paste("pseudopipe.sh",
            getwd(),
            paste(getwd(),"/input/dna/masked.fa",sep=""),
            paste(getwd(),"/input/dna/masked.%s.fa",sep=""),
            paste(getwd(),"/input/pep/protein.fa",sep=""),
            paste(getwd(),"/input/mysql/gff.%s.tsv",sep=""),
            "0",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # system("rm -r input/")
  # system("rm -r blast/")
  # system("rm -r dna/")
  # system("rm -r pep/")
  system("mv ./pgenes/*.gz ./")
  system("mv ./pgenes/*.txt ./")
  # system("rm -r pgenes/")
  setwd(wd)
}

# PseudoPipe output to gff3
PseudoPipe2gff3=function(Pseudo.out.txt=Pseudo.out.txt,
                         species=species,
                         out.gff3=out.gff3){
  origin=readLines(Pseudo.out.txt)
  origin=origin[2:length(origin)]
  
  in.type=c("DUP","PSSD","FRAG")
  out.type=c("duplicated","processed","fragment")
  names(out.type)=in.type
  
  gff3=sapply(1:length(origin),
              function(i){
                txt=unlist(strsplit(origin[i],"\t"))
                g=paste(txt[1],"\t","PseudoPipe","\t","pseudogene","\t",
                        txt[2],"\t",txt[3],"\t",".","\t",txt[4],"\t",
                        ".","\t",
                        paste("ID=pgene_",species,as.character(i),";",
                              "Parent=",txt[5],";",
                              "Type=",unname(out.type[txt[14]]),
                              sep=""),
                        sep="")
                return(g)
              })
  writeLines(gff3,out.gff3)
}