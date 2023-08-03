# Genomic analysis of protein-coding genes:
#  condon/AA usage, exon-intron, protein domain, gene family, etc.

#####
# condon/AA usage
#####
# Amino acid 2 nitrogen
AA2N=data.frame(AA=c("R","H","K","D","E",
                     "S","T","N","Q",
                     "C","G","P",
                     "A","V","I","L","M","F","Y","W"),
                N=c(4,3,2,1,1,
                    1,1,2,2,
                    1,1,1,
                    1,1,1,1,1,1,1,2))

# codon, nitrogen and amino acid
codonDNA=c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA",
           "CTG","CCT","CCC","CCA","CCG","CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT","ACC",
           "ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG","GAT",
           "GAC","GAA","GAG","GGT","GGC","GGA","GGG")
N_DNA=sapply(codonDNA,
             function(dnaCodon){
               dnaCodon=unlist(strsplit(dnaCodon,""))
               n=rep(NA,3)
               for (i in 1:3){
                 if (dnaCodon[i]=="A"){n[i]=5}
                 if (dnaCodon[i]=="T"){n[i]=2}
                 if (dnaCodon[i]=="C"){n[i]=3}
                 if (dnaCodon[i]=="G"){n[i]=5}
               }
               return(sum(n))
             })
codonRNA=sapply(codonDNA,
                function(dnaCodon){
                  dnaCodon=unlist(strsplit(dnaCodon,""))
                  rnaCodon=rep(NA,3)
                  for (i in 1:3){
                    if (dnaCodon[i]=="A"){rnaCodon[i]="U"}
                    if (dnaCodon[i]=="T"){rnaCodon[i]="A"}
                    if (dnaCodon[i]=="C"){rnaCodon[i]="G"}
                    if (dnaCodon[i]=="G"){rnaCodon[i]="C"}
                  }
                  return(paste(rnaCodon,collapse=""))
                })
N_RNA=sapply(codonRNA,
             function(rnaCodon){
               rnaCodon=unlist(strsplit(rnaCodon,""))
               n=rep(NA,3)
               for (i in 1:3){
                 if (rnaCodon[i]=="A"){n[i]=5}
                 if (rnaCodon[i]=="U"){n[i]=2}
                 if (rnaCodon[i]=="C"){n[i]=3}
                 if (rnaCodon[i]=="G"){n[i]=5}
               }
               return(sum(n))
             })
codon2N=data.frame(codonDNA=codonDNA,
                   N_DNA=N_DNA,
                   codonRNA=codonRNA,
                   N_RNA=N_RNA,
                   AA=c("F","F","L","L","S","S","S","S","Y","Y","*","*","C","C","*","W","L","L","L","L","P","P","P","P","H","H","Q","Q","R",
                        "R","R","R","I","I","I","M","T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A","D","D",
                        "E","E","G","G","G","G"))
codon2N=codon2N[order(codon2N$AA),]

# synonymous codon usage
# Dependencies: seqkit, parallel (R)
codonUsage=function(cds.fna=cds.fna,
                    threads=threads,
                    output.tsv=output.tsv){
  
  d=codon2N
  d[,"unifCodonFreq"]=sapply(d[,"AA"],
                             function(aa){return( 1/nrow(d[d[,"AA"]==aa,]) )})
  d=d[order(d[,"codonDNA"]),]
  cds.lst=system(paste("grep '>' ",cds.fna," | sed 's/>//'",sep=""),intern=TRUE)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads)-1)
  clusterExport(clus,list("d","cds.fna"),envir=environment())
  res=parSapply(clus,
                cds.lst,
                function(cds){
                  cmd=paste("seqkit grep -n -p",cds,cds.fna,"| seqkit seq -u -w 0",sep=" ")
                  seq=system(cmd,intern=TRUE)[2]
                  codonSeq=substring(seq,seq(1,nchar(seq),3),seq(3,nchar(seq),3))
                  codonUsage=as.data.frame(table(codonSeq));colnames(codonUsage)=c("Codon","codonUsage")
                  codonUsage=merge(codonUsage,d[,c("codonDNA","AA")],by.x="Codon",by.y="codonDNA",all=TRUE)
                  codonUsage[is.na(codonUsage)]=0
                  codonUsage=codonUsage[order(codonUsage[,"Codon"]),]
                  codonUsage[,"codonFreq"]=rep(NA,nrow(codonUsage))
                  for (i in 1:nrow(codonUsage)){
                    codonUsage[i,"codonFreq"]=codonUsage[i,"codonUsage"]/sum(codonUsage[codonUsage[,"AA"]==codonUsage[i,"AA"],"codonUsage"])
                  }
                  codonUsage[is.na(codonUsage)]=0
                  if (nrow(codonUsage)==64){ # Some genes covers gaps in the genome, leading to Ns in cds
                    return(codonUsage[,"codonFreq"])}else{
                      return(rep(NA,64))
                    }
                })
  stopCluster(clus)
  res=as.data.frame(res)
  res=res[,colSums(is.na(res))==0]
  d=cbind(d,res)
  d=d[order(d[,"AA"]),]
  write.table(d,output.tsv,sep="\t",row.names=FALSE,quote=FALSE)
}

# aa usage
# Dependencies: seqkit, parallel (R)
aaUsage=function(pep.faa=pep.faa,
                 threads=threads,
                 output.tsv=output.tsv){
  d=AA2N
  d=d[order(d[,"AA"]),]
  pep.lst=system(paste("grep '>' ",pep.faa," | sed 's/>//'",sep=""),intern=TRUE)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads)-1)
  clusterExport(clus,list("d","pep.faa"),envir=environment())
  res=parSapply(clus,
                pep.lst,
                function(pep){
                  cmd=paste("seqkit grep -n -p",pep,pep.faa,"| seqkit seq -u -w 0",sep=" ")
                  seq=system(cmd,intern=TRUE)[2]
                  aaSeq=unlist(strsplit(seq,""))
                  aaUsage=as.data.frame(table(aaSeq));colnames(aaUsage)=c("AA","aaUsage")
                  aaUsage=merge(d[,c("AA","N")],aaUsage,by="AA",all=TRUE)
                  aaUsage[is.na(aaUsage)]=0
                  aaUsage=aaUsage[order(aaUsage[,"AA"]),]
                  if (nrow(aaUsage)==20){
                    return(aaUsage[,"aaUsage"]/sum(aaUsage[,"aaUsage"]))}else{
                      return(rep(NA,20))
                    }
                })
  stopCluster(clus)
  res=as.data.frame(res)
  res=res[,colSums(is.na(res))==0]
  d=cbind(d,res)
  write.table(d,output.tsv,sep="\t",row.names=FALSE,quote=FALSE)
}

# # nitrogenPerAa
# # Dependencies: seqkit,parallel (R)
# nitrogenPerAa=function(proteins.faa=proteins.faa,
#                        threads=threads,
#                        output.tsv=output.tsv){
#   cmd=paste("seqkit fx2tab",proteins.faa,"-l -n -i -H >",output.tsv,sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   d=read.table(output.tsv,sep="\t",header=FALSE,quote="")
#   colnames(d)=c("ID","protLength")
#   
#   library(parallel)
#   clus=makeCluster(as.numeric(threads)-1)
#   clusterExport(clus,list("d","proteins.faa","AA2N"),envir=environment())
#   d$nitrogenInProt=
#     parSapply(clus,
#              1:nrow(d),
#              function(i){
#                ID=d[i,"ID"]
#                cmd=paste("seqkit grep -n -p",ID,proteins.faa,"| seqkit seq -u -w 0",sep=" ")
#                AAseq=system(cmd,intern=TRUE)[2]
#                AAseq=unlist(strsplit(AAseq,""))
#                nitrogen=sapply(AAseq,
#                                function(aa){
#                                  rownames(AA2N)=AA2N[,"AA"]
#                                  return( AA2N[aa,"N"] ) 
#                                })
#                return(sum(nitrogen))
#              })
#   stopCluster(clus)
#   
#   d$nitrogenPerAA=d$nitrogenInProt/d$protLength
#   
#   write.table(d,output.tsv,sep="\t",row.names=FALSE,quote=FALSE)
#   return(d)
# }
# 
# # nitrogenPerCodon
# # Dependencies: seqkit,parallel (R)
# nitrogenPerCodon=function(cds.fna=cds.fna,
#                           threads=threads,
#                           output.tsv=output.tsv){
#   cmd=paste("seqkit fx2tab",cds.fna,"-l -n -i -H >",output.tsv,sep=" ")
#   print(cmd);system(cmd,wait=TRUE)
#   
#   d=read.table(output.tsv,sep="\t",header=FALSE,quote="")
#   colnames(d)=c("ID","cdsLength")
#   d$codonNum=d$cdsLength/3
#   
#   library(parallel)
#   clus=makeCluster(as.numeric(threads)-1)
#   clusterExport(clus,list("d","cds.fna","codon2N"),envir=environment())
#   d$totalNitrogenDNA=
#     parSapply(clus,
#               1:nrow(d),
#               function(i){
#                 ID=d[i,"ID"]
#                 cmd=paste("seqkit grep -n -p",ID,cds.fna,"| seqkit seq -u -w 0",sep=" ")
#                 CDSseq=system(cmd,intern=TRUE)[2]
#                 CDSseq=substring(CDSseq,seq(1,nchar(CDSseq),3),seq(3,nchar(CDSseq),3))
#                 nitrogen=sapply(CDSseq,
#                                 function(codon){
#                                   rownames(codon2N)=codon2N[,"codonDNA"]
#                                   return( codon2N[codon,"N_DNA"] ) 
#                                 })
#                 return(sum(nitrogen))
#               })
#   d$totalNitrogenRNA=
#     parSapply(clus,
#               1:nrow(d),
#               function(i){
#                 ID=d[i,"ID"]
#                 cmd=paste("seqkit grep -n -p",ID,cds.fna,"| seqkit seq -u -w 0",sep=" ")
#                 CDSseq=system(cmd,intern=TRUE)[2]
#                 CDSseq=substring(CDSseq,seq(1,nchar(CDSseq),3),seq(3,nchar(CDSseq),3))
#                 nitrogen=sapply(CDSseq,
#                                 function(codon){
#                                   rownames(codon2N)=codon2N[,"codonDNA"]
#                                   return( codon2N[codon,"N_RNA"] ) 
#                                 })
#                 return(sum(nitrogen))
#               })
#   stopCluster(clus)
#   
#   d$nitrogenPerDnaCodon=d$totalNitrogenDNA/d$codonNum
#   d$nitrogenPerRnaCodon=d$totalNitrogenRNA/d$codonNum
#   
#   write.table(d,output.tsv,sep="\t",row.names=FALSE,quote=FALSE)
#   return(d)
# }



