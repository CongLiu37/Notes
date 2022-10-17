# Search for specific pattern in promoters of a genome
# Promoter: 1000 bp upstream + 500 bp downstream of the transcription start site

# Extract bed of promoters from gff file.
gff2bed_promoters=function(gff=gff,bed=bed){
  gff=read.table(gff,sep="\t",header=FALSE,quote="")
  gene=subset(gff,gff[,3]=="gene")
  bed_promoters=data.frame(gene[,1],
                           gene[,4]-1501,
                           gene[,4]+500,
                           gene[,9],gene[,6],gene[,7])
  bed_promoters[,2][bed_promoters[,2]<0]=0
  
  # Wrangle gene IDs
  func=function(str){
    str=unlist(strsplit(str,";"))[1]
    str=gsub(pattern="ID=gene-",replacement="",x=str)
  }
  
  bed_promoters[,4]=sapply(bed_promoters[,4],func)
  write.table(bed_promoters,bed,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  return(0)
}

# Extract fasta by bed.
# Dependencies: BEDtools
bed2fasta_promoter=function(bed=bed,
                            fna=fna,
                            promoters.fna=promoters.fna){
  cmd=paste("bedtools getfasta -s -name -fi",fna,"-bed",bed,">",promoters.fna,sep=" ")
  system(cmd,wait=TRUE)
  return(0)
}

# Search patterns.
# Dependencies: Biostrings (R)
PatternInPromoters=function(promoters.fna=promoters.fna,
                            pattern=pattern, # IUPAC ambiguity code.
                                             # e.g. hypoxia response element
                                             # RCGTGVBB
                                             # means [AG][C][G][T][G][ACG][AG][AG]
                            out=out){
  library(Biostrings)
  promoters=readBStringSet(promoters.fna, format="fasta",
                           nrec=-1L, skip=0L, seek.first.rec=FALSE,
                           use.names=TRUE, with.qualities=FALSE)
  
  PatternCount=vcountPattern(pattern=pattern,subject=DNAStringSet(promoters),
                             max.mismatch=0,min.mismatch=0,with.indels=FALSE,
                             fixed=FALSE,algorithm="auto")
  f1=function(n){if (n==0){return(FALSE)}else{return(TRUE)}}
  PatternPresent=sapply(PatternCount,f1)
  
  promoters=as.data.frame(promoters)
  
  TranscriptID=rownames(promoters)
  
  # Wrangle gene IDs
  func=function(string){string=unlist(strsplit(string,":"))[1];return(string)}
  TranscriptID=sapply(TranscriptID,func);TranscriptID=unname(TranscriptID)
  o=data.frame(Transcript=TranscriptID,Pattern=PatternPresent)
  write.table(o,out,sep="\t",row.names=FALSE,quote=FALSE)
  return(0)
}