arg=commandArgs(trailingOnly = TRUE)
blast=arg[1]
queries2Length=read.table("../inputPEP.length.tsv",sep="\t",header=FALSE,quote="")
blast=read.table(blast,sep="\t",header=FALSE,quote="")
colnames(blast)=c('#query','chr','ident','overlap','del','insert',
                  'query_start','query_end','chr_start','chr_end',
                  'eval','score','strand')
blast[,"strand"]=sapply(blast[,"strand"],
                        function(i){
                          if (i=="M"){return("-")}else{return("+")}
                        })
queries=blast[,"#query"][!duplicated(blast[,"#query"])]
queries2Length=queries2Length[queries2Length[,1] %in% queries,]
colnames(queries2Length)=c("#query","query_len")
blast=merge(blast,queries2Length,by="#query")
blast=blast[,c('#query','chr','ident','overlap','del','insert',
               'strand','query_start','query_end','query_len','chr_start','chr_end',
               'eval','score')]
blast=blast[blast[,"score"]>50,] # remove hits with score <= 50
blast=blast[blast[,"eval"]<1e-10,] # remove hits with eval >= 1e-10
# blast=blast[blast[,"ident"]>50,] # remove hits with identity <= 50%
blast=blast[blast[,"query_end"]-blast[,"query_start"]>0.75*blast[,"query_len"],] # alignment covers at least 75% of query

Chrs=blast[,"chr"][!duplicated(blast[,"chr"])]
sapply(Chrs,
       function(chr){
         d=blast[blast[,"chr"]==chr,]
         P=d[d[,"strand"]=="+",]
         write.table(P,paste(chr,"_P_blastHits",sep=""),
                     sep="\t",row.names=FALSE,quote=FALSE)
         M=d[d[,"strand"]=="-",]
         write.table(M,paste(chr,"_M_blastHits",sep=""),
                     sep="\t",row.names=FALSE,quote=FALSE)
       })
