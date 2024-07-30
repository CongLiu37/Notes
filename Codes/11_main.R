# Genomic analysis of protein-coding genes:
#  condon/AA usage, exon-intron, protein domain, gene family, etc.

#####
# condon/AA usage (standard genetic code)
#####
# Codon usage indices
# Dependencies: seqinr (R), Biostrings (R), parallel (R)
uco.seqinr=function(cds.fna,
                    threads=threads,
                    out.tsv=out.tsv){
  # cds.fna="/bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/Aaca/Aaca_cds_rep.fna"
  
  library(Biostrings)
  library(seqinr)
  cds=readDNAStringSet(cds.fna)
  cds=Biostrings::as.data.frame(cds)
  colnames(cds)="seq"
  cds$ID=rownames(cds)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  clusterExport(cl=clus,varlist=list("cds","uco","s2c"),envir=environment())
  res=parSapply(clus,1:nrow(cds),
                function(i){
                  index=uco(seq=s2c(cds[i,1]),frame=0,as.data.frame=TRUE)
                  index=cbind(ID=rep(cds[i,2],nrow(index)),index)
                  return(index)
                })
  res=data.frame(ID=unlist(res["ID",]),
                 AA=unlist(res["AA",]),
                 codon=unlist(res["codon",]),
                 eff=unlist(res["eff",]),
                 freq=unlist(res["freq",]),
                 RSCU=unlist(res["RSCU",]))
  write.table(res,out.tsv,sep="\t",row.names=FALSE,quote=FALSE)
}

# Codon usage expressivity measures with ribosomal proteins as reference
# Dependencies: coRdon (R)
expressivity.coRdon=function(cds.fna=cds.fna,
                             ribosomal.lst=ribosomal.lst,
                             out.tsv=out.tsv){
  # cds.fna="/bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/Aaca/Aaca_cds_rep.fna"
  # ribosomal.lst="/bucket/BourguignonU/Cong/termite_pca/geneFunction/ribosomalProtein_bbh/Aaca_ribosomalProtein.lst"
  
  library(coRdon)
  cds=readSet(file=cds.fna)
  codon.table=codonTable(cds)
  ID.lst=getID(codon.table)
  
  ribosomal=readLines(ribosomal.lst)
  ref=ID.lst %in% ribosomal
  
  # ribosomal.codon.table=subset(codon.table,ref)
  # nonRibo.codon.table=subset(codon.table,!ref)
  # CAI(cTobject=nonRibo.codon.table,subsets=list(ribosomal.codon.table))

  res=data.frame(ID=ID.lst)
  
  milc=MILC(cTobject=codon.table,self=FALSE,subsets=list(ref))
  res$milc=milc[,1]
  
  melp=MELP(cTobject=codon.table,subsets=list(ref))
  res$melp=melp[,1]
  
  e=E(cTobject=codon.table,subsets=list(ref))
  res$e=e[,1]
  
  cai=CAI(cTobject=codon.table,subsets=list(ref))
  res$cai=cai[,1]
  
  fop=Fop(cTobject=codon.table,subsets=list(ref))
  res$fop=fop[,1]
  
  write.table(res,out.tsv,sep="\t",row.names=FALSE,quote=FALSE)
}


# tAI: tRAN adaptation index 
# Dependencies: tAI (R), ape (R), stringr (R)
tai=function(cds.fna=cds.fna,
             tRNA.gff3=tRNA.gff3, # tRNAScan-SE + convert_tRNAScanSE_to_gff3.pl (biocode)
             tAI.misc.dir="/home/c/c-liu/Softwares/tai/misc/",
             kingdom="eukaryote", # eukaryote (tAI parameters optimized based on yeast Saccharomyces cerevisiae)
                                  # bacteria (tAI parameters optimized based on E. coli)
             out_prefix=out_prefix){
  cmd=paste(paste(tAI.misc.dir,"/codonM",sep=""),
            cds.fna,
            paste(out_prefix,".m",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  library(ape)
  library(stringr)
  library(Biostrings)
  tRNA=read.gff(tRNA.gff3)
  tRNA=tRNA[tRNA$type=="tRNA","attributes"]
  tRNA=str_extract(tRNA,"anticodon=[ATGC]{3};")
  tRNA=sub("^anticodon=","",tRNA)
  tRNA=sub(";$","",tRNA)
  tRNA=sapply(tRNA,
              function(i){
                bases=c("A","T","G","C")
                names(bases)=c("T","A","C","G")
                i=unlist(strsplit(i,""))
                i=i[c(3,2,1)]
                i=unname(bases[i])
                return(paste(i,collapse=""))
              })
  anticodon_complement.lst=c("TTT","TTC","TTA","TTG",
                  "TCT","TCC","TCA","TCG",
                  "TAT","TAC","TAA","TAG",
                  "TGT","TGC","TGA","TGG",
                  "CTT","CTC","CTA","CTG",
                  "CCT","CCC","CCA","CCG",
                  "CAT","CAC","CAA","CAG",
                  "CGT","CGC","CGA","CGG",
                  "ATT","ATC","ATA","ATG",
                  "ACT","ACC","ACA","ACG",
                  "AAT","AAC","AAA","AAG",
                  "AGT","AGC","AGA","AGG",
                  "GTT","GTC","GTA","GTG",
                  "GCT","GCC","GCA","GCG",
                  "GAT","GAC","GAA","GAG",
                  "GGT","GGC","GGA","GGG")
  tRNA_copy=sapply(anticodon_complement.lst,function(code){return(length(tRNA[tRNA==code]))})
  tRNA_copy=unname(tRNA_copy)
  writeLines(as.character(tRNA_copy),paste(out_prefix,".trna",sep=""))
  
  require("tAI")
  trna=scan(paste(out_prefix,".trna",sep=""))
  if (kingdom=="eukaryote"){sking=0}else{sking=1}
  ws=get.ws(tRNA=trna, sking=sking)
  m=matrix(scan(paste(out_prefix,".m",sep="")), ncol=61, byrow=TRUE)
  m=m[,-33]
  tai=get.tai(m,ws)
  ID.lst=system(paste("grep '>' ",cds.fna," | sed 's/>//'"),intern=TRUE)
  write.table(data.frame(ID=ID.lst,tai=tai),
              paste(out_prefix,"_tai.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
}

# GC123
# Dependencies: coRdon
GC123=function(cds.fna=cds.fna,
               out.tsv=out.tsv){
  # cds.fna="/bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/Aaca/Aaca_cds_rep.fna"
  
  library(coRdon)
  cds=readSet(file=cds.fna)
  codon.table=codonTable(cds)
  codon.count=codonCounts(codon.table)
  ID.lst=getID(codon.table)
  rownames(codon.count)=ID.lst
  
  codon2GC12=sapply(colnames(codon.count),
                    function(i){
                      i=unlist(strsplit(i,""))[c(1,2)]
                      i=i %in% c("G","C")
                      i=which(i)
                      return(length(i))
                    })
  codon2GC3=sapply(colnames(codon.count),
                   function(i){
                     i=unlist(strsplit(i,""))[3]
                     if (i %in% c("G","C")){
                       return(1)
                     }else{
                       return(0)
                     }
                   })
  codon2GC=sapply(colnames(codon.count),
                  function(i){
                    i=unlist(strsplit(i,""))
                    i=i %in% c("G","C")
                    i=which(i)
                    return(length(i))
                  })
  
  res=data.frame(ID=ID.lst)
  res$length.bp=apply(codon.count,1,function(i){return(3*sum(i))})
  res$GC12=apply(codon.count,1,
                 function(i){
                   return( sum(i*codon2GC12[colnames(codon.count)]) / (2*sum(i)) )
                 })
  res$GC3=apply(codon.count,1,
                function(i){
                  return( sum(i*codon2GC3[colnames(codon.count)]) / (sum(i)) )
                })
  res$GC=apply(codon.count,1,
               function(i){
                 return( sum(i*codon2GC[colnames(codon.count)]) / (3*sum(i)) )
               })
  write.table(res,out.tsv,sep="\t",row.names=FALSE,quote=FALSE)
}

# generate random pseudo-CDS from intergenic regions
# Dependencies: seqkit, bedtools
pseudoCDS_randomGenomicRegion=function(chrLen.tsv=chrLen.tsv, # <chromName><TAB><chromSize>
                                       genome.fna=genome.fna,
                                       region.length=1050,
                                       sample_size=1000, # How many regions to be generated
                                       out_prefix=out_prefix){
  cmd=paste("bedtools random",
            "-l",as.character(region.length),
            "-n",as.character(sample_size),
            "-g",chrLen.tsv,">",
            paste(out_prefix,".bed",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("bedtools getfasta -nameOnly -s",
            "-fi",genome.fna,
            "-bed",paste(out_prefix,".bed",sep=""),
            "| seqkit seq -u > ",
            paste(out_prefix,".fna",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# # Angle between vector for codon frequencies of CDS and tRNA pool
# angle=function(tai.m=tai.m,tai.trna=tai.trna,tai_tai.tsv=tai_tai.tsv,# files from tai (R)
#                out.tsv=out.tsv){
#   tai.m=read.table(tai.m,sep="\t",header=FALSE,quote="");tai.m=tai.m[,-62]
#   tai.trna=scan(tai.trna);tai.trna=tai.trna[-c(11,12,15)]
#   tai_tai.tsv=read.table(tai_tai.tsv,sep="\t",header=TRUE,quote="")
# }




# Amino acid 2 nitrogen (20 amino acids)
AA2N=data.frame(AA=c("R","H","K","D","E",
                     "S","T","N","Q",
                     "C","G","P",
                     "A","V","I","L","M","F","Y","W"),
                N=c(4,3,2,1,1,
                    1,1,2,2,
                    1,1,1,
                    1,1,1,1,1,1,1,2))

# codon, nitrogen and amino acid (standard genetic code)
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



