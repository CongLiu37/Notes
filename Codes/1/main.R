# Differentail expression analysis using RNA-seq and reference genome.
# Dependencies:
#   Softwares: SRA-toolkit, FastQC, Trimmomatic, Gffread, Hisat2, SAMtools, StringTie, Minpath, eggNOG-mapper (online web). 
#   R packages: DESeq2, ggplot2, ComplexHeatmap, circlize, annotate, stringr, dplyr, pathview.

# SRA-toolkit: Download fastq from NCBI SRA.
DownloadSRA = function(AccessionList=AccessionList, # a tabular table without header. 
                                                    # Its first column is SRA IDs.
                       out_dir=out_dir, # directory in which downloaded fq files are saved.
                       threads=threads){
  threads = as.character(threads)
  
  Accessions = read.table(AccessionList,sep="\t",header=FALSE,quote="")
  for (accession in Accessions[,1]){
      cmd = paste("fasterq-dump","--split-files",accession,"-O",out_dir,"-e",threads,sep=" ")
      print(cmd);system(cmd,wait=TRUE)
  }
  
  return(out_dir)
}

# FastQC: Assess quality of sequencing data.
QualityCheck = function(fq1=fq1, # Input fq file.
                        fq2=fq2, # Input fq file. Set "none" if single-end.
                        out_dir=out_dir, # directory in which quality check reports are saved.
                        threads=threads){
  threads = as.character(threads)
  
  if (fq2!="none"){ # pair-end
      cmd = paste("fastqc","-o",out_dir,"-t",threads,fq1,fq2,sep=" ")
  }else{ # single-end
      cmd = paste("fastqc","-o",out_dir,"-t",threads,fq1,sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  
  return(out_dir)
}

# Trimmomatic: Quality control of sequencing data.
QualityFilter = function(fq1=fq1,fq2=fq2, # Input fq files. Set fq2="none" if single-end.
                         clean_fq1=clean_fq1,clean_fq2=clean_fq2,# Names of clean fq files. Set clean_fq2="none" if single-end.
                         unpaired_fq1=unpaired_fq1,unpaired_fq2=unpaired_fq2, # Names of fq files for unpaired clean reads. Set unpaired_fq2="none" if single-end.
                         QualityFilter=QualityFilter, # String of step options of Trimmomatic.
                                                      # General steps included in quality filter:
                                                      # 1. adaptor; 2. low quality tail; 3. read length; 4. average quality score
                         threads=threads){
  threads = as.character(threads)
    
  if (fq2!="none"){ # pair-end
      cmd = paste("trimmomatic","PE","-threads",threads,"-phred33",fq1,fq2,clean_fq1,unpaired_fq1,clean_fq2,unpaired_fq2,QualityFilter,sep=" ")
  }else{ # single-end
      cmd = paste("trimmomatic","SE","-threads",threads,"-phred33",fq1,clean_fq1,QualityFilter,sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  
  return(0)
}

# gffread: Extract cds from gff.
gffread = function(gff=gff,fna=fna,cds.fna=cds.fna){
  cmd = paste("gffread",gff,"-g",fna,"-x",cds.fna,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(cds.fna)
}

# Hisat2; Build hisat2 index of reference genome.
Hisat2Build = function(fna=fna, # FASTA of reference genome
                       index_prefix=index_prefix){
  cmd = paste("hisat2-build",fna,index_prefix,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(index_prefix)
}

# Hisat2: Map reads to reference genome. 
# SAMtools: Compress SAM to BAM and sort BAM.
Hisat = function(fq1=fq1,fq2=fq2, # Input fq files. Make fq2="none" if single-end.
                 index=index, # Basename of Hisat2 index of reference genome.
                 out_prefix=out_prefix, # Prefix of output BAM file.
                 threads=threads){
  threads = as.character(threads)
  
  bam_filename=paste(out_prefix,".bam",sep="")
  if (fq2!="none"){ # pair-end
      cmd = paste("hisat2","--dta","-x",index,"-p",threads,"-1",fq1,"-2",fq2,"|",
                  "samtools","view","-@",threads,"-bS","|",
                  "samtools","sort","-@",threads,"-o",bam_filename,sep=" ")
  }else{ # single pair
      cmd = paste("hisat2","--dta","-x",index,"-p",threads,"-U",fq1,"|",
                  "samtools","view","-@",threads,"-bS","|",
                  "samtools","sort","-@",threads,"-o",bam_filename,sep=" ")
  }
  print(cmd);system(cmd,wait=TRUE)
  return(bam_filename)
}

# StringTie: Assemble transcripts.
StringTie = function(input_bam=input_bam, # Input sorted BAM.
                     gtf=gtf, # GTF of reference genome.
                     output_prefix=output_prefix,
                     threads=threads){
  threads = as.character(threads)
  
  gtf_filename = paste(output_prefix,".gtf",sep="")
  cmd = paste("stringtie","-p", threads,"-G", gtf, "-o", gtf_filename, "-e -l", output_prefix, input_bam,sep=" ")
  print(cmd);system(cmd,wait = T)
  return(gtf_filename)
}

# prepDE.py3 provided by StringTie: Extract gene/transcript read count matrix from output of stringtie.
run.prepDE.py = function(inputfile=inputfile # A space-separated table without header. 
                                              # Its first column is sample name, and second column is path to corresponding GTF from stringtie.
                         ){
  cmd = paste("python3","prepDE.py3","-i",inputfile,sep=" ")
  print(cmd);system(cmd,wait=T)
  return(0)
}

# DESeq2: Differential expression analysis.
DESeq = function(mat=mat, # Expression matrix generated by stringtie-prepDE.py3.
                 metadata=metadata, # A non-header tabular table whose 
                                    # first column is sample name, 
                                    # second column is a list of paths to fq files seperated by comma, 
                                    # third column is sample description (factor).
                                    # ith column of mat and ith row of metadata must represent the same sample.
                                    # Name experiment group prior to control group.
                 output_prefix=output_prefix){
  library(DESeq2)
    
  countData=read.csv(mat,row.names=1)
  rownames(countData)=sapply(rownames(countData),gsub,pattern="rna-",replacement="")
  colData = read.table(metadata,header=FALSE,sep="\t",row.names=1,quote="",col.names=c("Sample","Fq","Treatment"))
    
  # Normalize read count via variance stabilizing transformation (vst).
  Transformed = as.data.frame(vst(as.matrix(countData),blind=FALSE))
  # This is the expression matrix shall be used for tasks other than differential expression analysis, e.g. PCA and clustering.
  write.table(Transformed,paste(output_prefix,"_NormalizedCount.txt",sep=""),sep="\t",row.names=TRUE,quote=FALSE)
    
  # Differential expression.
  dds = DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~Treatment)
  dds = DESeq(dds)
  res = results(dds)
  resOrdered = res[order(res$padj), ]
  resOrdered = as.data.frame(resOrdered)

  # Whether a gene is induced, inhibited or not influenced.
  GetColor = function(i){
    FoldChange = resOrdered[i,"log2FoldChange"]
    padj = resOrdered[i,"padj"]
    if (is.na(padj)){color="none"}else{
    if (padj<0.05&FoldChange>1){color="up"} 
    if (padj<0.05&FoldChange< -1){color="down"}
    if (padj>=0.05){color="none"}
    if (FoldChange>=-1&FoldChange<=1){color="none"}}
    return(color)
  }
  resOrdered[,"color"] = sapply(1:nrow(resOrdered),GetColor) # Show whether the gene is induced, inhibited or not infulenced by colomn "color".
                                                             # "color" can be "up", "down" or "none".    
  write.table(resOrdered,paste(output_prefix,"DEGs.txt",sep=""),sep="\t",row.names=TRUE,quote=FALSE)
  return(0)
}

# Volcano plot: log2FoldChange against -log10(padj)
Volcano = function(deseq2=deseq2, # Tabular table of differential expression analysis. 
                                  # Columns include "log2FoldChange" (numeric), "padj" (numeric) and "color" (factor, up/down/none)
                   output=output){
    library(ggplot2)
    data = read.table(deseq2,sep="\t",row.names=1,header=TRUE,quote="")
    data = subset(data,!is.na(padj))
    pdf(output)
    p=ggplot()+
        geom_point(aes(x=data$log2FoldChange,y=-log10(data$padj),color=data$color))+
        theme_classic()+
        scale_color_manual(values=c("blue","black","red"))+
        xlab("log2(Treated/Control)")+
        ylab("-log10(adjusted P value)")+
        theme(plot.margin=unit(c(1,1.5,1,1),"lines"),
              axis.title.y=element_text(size=8,face="bold"),
              axis.title.x=element_text(size=8,face="bold"),
              axis.text.y=element_text(size=8,face="bold"),
              axis.text.x=element_text(size=8,face="bold"),
              title=element_text(size=8,face="bold"),
              legend.text=element_text(size=8,face="bold"),
              legend.title = element_blank())
    print(p)
    dev.off()
    return(p)
}

#Principle component analysis
PCA = function(mat=mat, # Expression matrix transformed by vst.
               metadata=metadata, # A non-header tabular table whose 
                                  # first column is sample name, 
                                  # second column is a list of paths to fq files seperated by comma, 
                                  # third column is sample description.
               output=output){
    metadata = read.table(metadata,header=FALSE,sep="\t",row.names=1,quote="",
                          col.names=c("Sample","Fq","Treatment"))
    mat = read.table(mat,sep="\t",row.names=1,header=TRUE,quote="")
    mat = as.matrix(mat)
    mat = t(mat)
    
    pca = prcomp(mat) # pca
    
    library(ggplot2)
    data = as.data.frame(pca$x)
    
    data[,"Treatment"] = metadata$Treatment
    pc1 = summary(pca)$importance[2,"PC1"]
    pc2 = summary(pca)$importance[2,"PC2"]
    pdf(output)
    p=ggplot()+
      geom_point(data=data,
                 aes(x=PC1,y=PC2,color=Treatment))+
      stat_ellipse(data=data,
                   aes(x=PC1,y=PC2,color=Treatment),
                   level=0.95,linetype=2)+
      xlab(paste("PC1 (",as.character(pc1*100),"%)",sep=""))+
      ylab(paste("PC2 (",as.character(pc2*100),"%)",sep=""))+
      theme_classic()+
      theme(plot.margin=unit(c(1,1.5,1,1),"lines"),
            axis.title.y=element_text(size=8,face="bold"),
            axis.title.x=element_text(size=8,face="bold"),
            axis.text.y=element_text(size=8,face="bold"),
            axis.text.x=element_text(size=8,face="bold"),
            title=element_text(size=8,face="bold"),
            legend.text=element_text(size=8,face="bold"),
            legend.title = element_blank())
    print(p)
    dev.off()
    return(p)
}

# Gene expression heatmap
Heatmap = function(mat=mat, # Expression matrix transformed by vst.
                   output=output){
    library(ComplexHeatmap)
    library(circlize)
    mat = read.table(mat,sep="\t",row.names=1,header=TRUE,quote="")
    mat = log10(as.matrix(mat))
    
    pdf(output,onefile = FALSE)
    p = Heatmap(mat,
                heatmap_legend_param = list(title = "log10(VST counts)",
                                            title_position="topleft",
                                            title_gp=gpar(fontsize=5, fontface="bold"),
                                            labels_gp = gpar(fontsize = 5)),
                show_row_names = FALSE,
                column_title="Sample",
                column_names_gp=gpar(fontface ="bold",fontsize=10),
                column_title_gp=gpar(fontface="bold",fontsize=10),
                column_names_side = "top",
                column_names_rot = 0,
                cluster_rows=FALSE,
                cluster_columns=TRUE,
                col=circlize::colorRamp2(c(min(mat),max(mat)),c("white","red"))
    )
    p=draw(p,heatmap_legend_side="left")
    print(p)
    dev.off()
    return(p)
}

# eggNOG-mapper: CDS functional annotation.
eggNOG = function(cds.fna=cds.fna,out_prefix=out_prefix,
                  threads=threads,block_size=block_size){
  threads=as.character(threads);block_size=as.character(block_size)
  
  cmd=paste("emapper.py","-m","diamond","--itype","CDS","-i",cds.fna,"-o",out_prefix,
            "--no_file_comments","--block_size",block_size,"--cpu",threads,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  return(0)
}

# GO enrichment
enrichGO = function(Background=Background, # Tabular table with header and row names.
                                           # Row names are gene names.
                                           # Column "color" indicates whether gene expression is "up", "down" or "none" changed.
                    up=up, # logical. If true, enrich up-regulated genes.
                    Annotation=Annotation, # Output of eggNOG-mapper 
                                           # Fields: query,GOs ("-" if unannotated) 
                                           # Annotations of reference genome
                    output_prefix=output_prefix){
  # Read data
  Background = read.table(Background,sep="\t",header=TRUE,row.names=1,quote="")
  Annotation = read.table(Annotation,sep="\t",header=TRUE,quote="")
  Annotation = data.frame(query=Annotation[,"query"],term=Annotation[,"GOs"])
  
  # DEG names
  if (up){
    DEGs = rownames(subset(Background,color=="up"))
  }else{
    DEGs = rownames(subset(Background,color=="down"))
  }
  
  DEGs2Term = subset(Annotation,Annotation[,"query"]%in%DEGs)
  TermOfDEGs = unname(unlist(sapply(DEGs2Term[,"term"],strsplit,split=",")))
  TermOfDEGs = sort(unname(TermOfDEGs[which(TermOfDEGs!="-")])) # "-" represents unannotated genes
  TermOfDEGs_nr = TermOfDEGs[!duplicated(TermOfDEGs)]
  Term2DEGCount = table(TermOfDEGs)
  
  TermOfAnnotation = unname(unlist(sapply(Annotation[,"term"],strsplit,split=",")))
  Term2AnnotationCount = table(TermOfAnnotation)
  
  # Enrichment by Fisher exact test
  fisher = function(term){ # P value of GO term
    DEG_term = Term2DEGCount[term]
    DEG_NOterm = nrow(DEGs2Term)-DEG_term
    NoDEG_term = Term2AnnotationCount[term]-DEG_term
    NoDEG_NOterm = nrow(Background)-DEG_term-DEG_NOterm-NoDEG_term
    p=fisher.test(matrix(c(DEG_term,NoDEG_term,DEG_NOterm,NoDEG_NOterm),nrow=2),alternative="greater")$p.value
    return(p)
  }
  degCount = function(term){ # DEG count of GO term
    DEG_term = Term2DEGCount[term]
    DEG_NOterm = nrow(DEGs2Term)-DEG_term
    NoDEG_term = Term2AnnotationCount[term]-DEG_term
    NoDEG_NOterm = nrow(Background)-DEG_term-DEG_NOterm-NoDEG_term
    return(DEG_term)
  }
  
  pvalues = sapply(TermOfDEGs_nr,fisher)
  padjs = p.adjust(pvalues,"BH") # False discover rate
  DegCount = sapply(TermOfDEGs_nr,degCount)
  o = data.frame(Term=TermOfDEGs_nr,pvalue=pvalues,padj=padjs,DegCount=DegCount)
  o = subset(o,padj<0.05)
  o = o[order(o[,"padj"]),]
  
  # GO trim. Remove enriched parent if its child is enriched
  library(annotate)
  RawTerms = o$Term
  GetParent = function(ID){
    parents=try(getGOParents(ID),silent=TRUE)[[1]][2][[1]]
    return(unname(parents))
  }
  parentsGO = sapply(RawTerms,GetParent)
  ParentsGO = unlist(parentsGO)
  ParentsGO = ParentsGO[!duplicated(ParentsGO)]
  o = subset(o,o$Term %in% setdiff(o$Term,ParentsGO))
  
  # GO description
  Term2Desc=getGOTerm(o$Term)
  BP=Term2Desc$BP;BP_id=names(BP);BP_desc=unname(BP)
  category=rep("Biological Processes",length(BP))
  BP = data.frame(Term=BP_id,Description=BP_desc,Category=category)
  MF=Term2Desc$MF;MF_id=names(MF);MF_desc=unname(MF)
  category=rep("Molecular Functions",length(MF))
  MF = data.frame(Term=MF_id,Description=MF_desc,Category=category)
  CC=Term2Desc$CC;CC_id=names(CC);CC_desc=unname(CC)
  category=rep("Cellular Components",length(CC))
  CC = data.frame(Term=CC_id,Description=CC_desc,Category=category)
  Term2Desc = rbind(BP,MF,CC)
  p=merge(o,Term2Desc,by.x="Term",by.y="Term")
  p=subset(p,!is.na(p$Description))
  
  # DEG to GO
  library(stringr)
  library(dplyr)
  GOs = str_split(DEGs2Term$term,",")
  time = sapply(GOs,length) 
  DEG2GO = data.frame(DEG=rep(DEGs2Term$query,time),GO=unlist(GOs))
  DEG2GO_des = merge(DEG2GO,p,by.x="GO",by.y="Term",all.y=TRUE)
  DEG2GO_des = DEG2GO_des[order(DEG2GO_des$padj),]
  write.table(data.frame(DEG=DEG2GO_des$DEG,GO=DEG2GO_des$GO,
                         Description=DEG2GO_des$Description,Category=DEG2GO_des$Category,
                         AdjustedP=DEG2GO_des$padj),
              paste(output_prefix,"_GOenriched.txt",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
  return(0)
}

# Bar plot of 10 most significantly enriched BP/MF/CC GO terms.
enrichGOplot = function(df=df, # Tabular table with header.
                               # Columns include "Description" (character), "Category" (factor), "AdjustedP" (numeric).
                        output_prefix=output_prefix){
  df=read.table(df,sep="\t",header=TRUE,quote="")
  df=data.frame(Description=df$Description,Category=df$Category,AdjustedP=df$AdjustedP)
  DEGcount=table(df$Description)
  DEGcount=data.frame(Description=names(DEGcount),count=unname(DEGcount))
  df=df[!duplicated(df),]
  df=merge(df,DEGcount,by="Description",all=TRUE)
  
  # 10 most significant in BP/CC/MF
  BP=subset(df,Category=="Biological Processes");BP=head(BP[order(BP$AdjustedP),],10)
  CC=subset(df,Category=="Cellular Components");CC=head(CC[order(CC$AdjustedP),],10)
  MF=subset(df,Category=="Molecular Functions");MF=head(MF[order(MF$AdjustedP),],10)
  df=rbind(MF,CC,BP)
  
  colnames(df)=c("Description","Category","AdjustedP","buzhidao","DEGcount")
  df$Description=factor(df$Description,levels=df$Description)
  library(ggplot2)
  pdf(paste(output_prefix,".pdf",sep=""),height=12,width=10)
  plot=ggplot(df,aes(y=DEGcount,x=Description))+
    geom_bar(stat="identity",aes(fill=Category))+
    ylab("DEG count")+
    xlab("")+
    labs(fill="Category")+
    scale_y_continuous(expand=c(0,0))+
    theme_classic()+
    theme(plot.margin=unit(c(1,1.5,1,1),"lines"),
          axis.title.x=element_text(size=8,face="bold"),
          axis.title.y=element_text(size=8,face="bold"),
          axis.text.x=element_text(size=11,face="bold",angle=70,hjust=1),
          axis.text.y=element_text(size=8,face="bold"),
          title=element_text(size=8,face="bold"),
          legend.text=element_text(size=8,face="bold"),
          legend.title = element_text(size=8,face="bold"),
          legend.position="left")
    print(plot)
    dev.off()
    return(plot)
}

# Minpath: KEGG pathway inference of whole genome.
PathInf = function(Annotation=Annotation, # Output of eggNOG-mapper 
                                          # Fields: query,KEGG_ko ("-" if unannotated) 
                                          # Annotations of reference genome
                   output_prefix=output_prefix){
  Annotation = read.table(Annotation,sep="\t",header=TRUE,quote="")
  Annotation = data.frame(query=Annotation[,"query"],term=Annotation[,"KEGG_ko"])
  
  # Prepare input of minpath
  library(stringr)
  library(dplyr)
  KOs = str_split(Annotation$term,",")
  time = sapply(KOs,length)
  data = data.frame(query=rep(Annotation$query,time),
                    KO=unlist(KOs))
  data = subset(data,data$KO!="-")
  data$KO = sapply(data$KO,gsub,pattern="ko:",replacement="")
  write.table(data,"temp_query2KO.txt",sep="\t",quote=FALSE,row.names=FALSE) # Temporary table.
  
  # Run MinPath
  cmd = paste("python","MinPath.py","-ko","temp_query2KO.txt","-report",
              paste(output_prefix,".txt",sep=""),
              "-details",
              paste(output_prefix,"_Details.txt",sep=""),
              sep=" ")
  system(cmd,wait=TRUE)
  system(paste("rm","temp_query2KO.txt",sep=" "),wait=TRUE)
  
  cmd = paste("sed","\'s/  /\t/\'",
              paste(output_prefix,".txt",sep=""),
              ">",
              paste(output_prefix,"_tab.txt",sep=""))
  system(cmd,wait=TRUE)
  
  Pathways = read.table(paste(output_prefix,"_tab.txt",sep=""),sep="\t",header=FALSE,quote="")
  system(paste("rm",paste(output_prefix,"_tab.txt",sep=""),sep=" "))
  
  # Wrangle output of MinPath 
  Func = function(str,i){return(strsplit(str," ")[[1]][i])}
  PathwayID = sapply(Pathways$V1,Func,i=2);func=function(str){return(paste("map",str,sep=""))}
  PathwayID = sapply(PathwayID,func)
  Naive = sapply(Pathways$V2,Func,i=2)
  MinPath = sapply(Pathways$V2,Func,i=5)
  TotalGene = sapply(Pathways$V2,Func,i=9)
  GeneFound = sapply(Pathways$V2,Func,i=13)
  Coverage = as.numeric(GeneFound)/as.numeric(TotalGene)
  library(stringr)
  ExtractDes = function(str){
    vector=strsplit(str," ")[[1]];vector=vector[17:length(vector)]
    return(str_c(vector,collapse=" "))
  }
  Description = sapply(Pathways$V2,ExtractDes)
  
  o = data.frame(PathwayID=PathwayID,Naive=Naive,MinPath=MinPath,
                 TotalGene=TotalGene,GeneFound=GeneFound,Coverage=Coverage,
                 Description=Description)
  o = o[order(o$Coverage,decreasing=TRUE),]
  
  write.table(o,paste(output_prefix,"_WholeGenome.txt",sep=""),row.names=FALSE,quote=FALSE)
  return(paste(output_prefix,"_WholeGenome.txt",sep=""))
}

# KEGG pathway enrichment
enrichKEGG = function(Background=Background, # Tabular table with header and row names.
                                             # Row names are gene names.
                                             # Column "color" indicates whether gene expression is "up", "down" or "none" changed.
                      up=up, # logical. If true, enrich up-regulated genes.
                      Annotation=Annotation, # Output of eggNOG-mapper 
                                             # Fields: query,KEGG_Pathway ("-" if unannotated),KEGG_ko ("-" if unannotated)
                                             # Annotations of reference genome
                      WholeGenomePath=WholeGenomePath, # Tabular table with column "PathwayID", which represents KEGG pathways present in reference genome.
                                                       # It can be from "PathInf", or from KEGG Taxonomy.
                      Pathway2KO=Pathway2KO, # From MinPath github: https://github.com/mgtools/MinPath/blob/master/data/KEGG-mapping.txt
                      KO2Des=KO2Des, # From MinPath github: https://github.com/mgtools/MinPath/blob/master/data/KEGG-family.txt
                      output_prefix=output_prefix){
  Background = read.table(Background,sep="\t",header=TRUE,row.names=1,quote="")
  Annotation = read.table(Annotation,sep="\t",header=TRUE,quote="")
  Annotation = data.frame(query=Annotation[,"query"],term=Annotation[,"KEGG_Pathway"],KO=Annotation[,"KEGG_ko"])
  Pathways = read.table(WholeGenomePath,sep="\t",header=TRUE,quote="")
  if ("MinPath"%in%colnames(Pathways)){Pathways = subset(Pathways,Pathways$MinPath!=0)}
  
  if (up){
    DEGs = rownames(subset(Background,color=="up"))
  }else{
    DEGs = rownames(subset(Background,color=="down"))
  }
  
  # Annotated DEGs
  DEGs2Term = subset(Annotation,Annotation[,"query"]%in%DEGs)
  TermOfDEGs = unname(unlist(sapply(DEGs2Term[,"term"],strsplit,split=",")))
  TermOfDEGs = sapply(TermOfDEGs,gsub,pattern="ko[0-9]*",replacement="")
  TermOfDEGs = sort(unname(TermOfDEGs[which(TermOfDEGs!="-"&TermOfDEGs!="")]))
  Term2DEGCount = table(TermOfDEGs)
  TermOfDEGs_nr = TermOfDEGs[!duplicated(TermOfDEGs)]
  
  # Annotations of background genes
  TermOfAnnotation = unname(unlist(sapply(Annotation[,"term"],strsplit,split=",")))
  TermOfAnnotation = sapply(TermOfAnnotation,gsub,pattern="ko[0-9]*",replacement="")
  TermOfAnnotation = sort(unname(TermOfAnnotation[which(TermOfAnnotation!="-"&TermOfAnnotation!="")]))
  Term2AnnotationCount = table(TermOfAnnotation)
  
  # Fisher exact test
  fisher = function(term){
    DEG_term = Term2DEGCount[term]
    DEG_NOterm = nrow(DEGs2Term)-DEG_term
    NoDEG_term = Term2AnnotationCount[term]-DEG_term
    NoDEG_NOterm = nrow(Background)-DEG_term-DEG_NOterm-NoDEG_term
    p=fisher.test(matrix(c(DEG_term,NoDEG_term,DEG_NOterm,NoDEG_NOterm), nrow = 2), alternative = "greater")$p.value
    return(p)
  }
  degCount = function(term){
    DEG_term = Term2DEGCount[term]
    DEG_NOterm = nrow(DEGs2Term)-DEG_term
    NoDEG_term = Term2AnnotationCount[term]-DEG_term
    NoDEG_NOterm = nrow(Background)-DEG_term-DEG_NOterm-NoDEG_term
    return(DEG_term)
  }
  
  pvalues = sapply(TermOfDEGs_nr,fisher)
  padjs = p.adjust(pvalues,"BH") # false discovery rate
  DegCount = sapply(TermOfDEGs_nr,degCount)
  o = data.frame(Term=TermOfDEGs_nr,pvalue=pvalues,padj=padjs,DegCount=DegCount)
  o = subset(o,padj<0.05)
  o = merge(o,Pathways,by.x="Term",by.y="PathwayID") # Only retain KEGG pathways that really present in reference genomes.
  o = o[order(o[,"padj"]),]
  
  # DEGs and enriched pathways
  library(stringr)
  library(dplyr)
  Path = str_split(DEGs2Term$term,",")
  time = sapply(Path,length)
  DEG2KEGG_Pathway = data.frame(DEG=rep(DEGs2Term$query,time),
                                KEGG_Pathway=unlist(Path))
  DEG2KEGG_Pathway$KEGG_Pathway = sapply(DEG2KEGG_Pathway$KEGG_Pathway,gsub,
                                         pattern="ko[0-9]*",replacement="-")
  DEG2KEGG_Pathway = subset(DEG2KEGG_Pathway,DEG2KEGG_Pathway$KEGG_Pathway!="-")
  DEG2KEGG_Pathway_description = merge(DEG2KEGG_Pathway,o,by.x="KEGG_Pathway",by.y="Term",all.y=TRUE)
  DEG2KEGG_Pathway_description = DEG2KEGG_Pathway_description[order(DEG2KEGG_Pathway_description$padj),]
  DEG2Pathway=data.frame(DEG=DEG2KEGG_Pathway_description$DEG,
                         KEGG_Pathway=DEG2KEGG_Pathway_description$KEGG_Pathway,
                         Description=DEG2KEGG_Pathway_description$Description,
                         AdjustedP=DEG2KEGG_Pathway_description$padj)
  
  # DEGs and KOs
  Annotation=data.frame(query=Annotation$query,KO=Annotation$KO)
  Pathway2KO=read.table(Pathway2KO,sep="\t",header=FALSE,quote="",colClasses="character",col.names=c("Pathway","KO"))
  t=function(i){return(paste("map",i,sep=""))}
  Pathway2KO$Pathway=sapply(Pathway2KO$Pathway,t)
  KO2Des=read.table(KO2Des,sep="\t",header=FALSE,quote="",col.names=c("KO","KO_Description"))
  DEGs=DEG2Pathway$DEG;DEGs=DEGs[!duplicated(DEGs)]
  DEGs2KO=subset(Annotation,Annotation$query%in%DEGs)
  KO = str_split(DEGs2KO$KO,",")
  time = sapply(KO,length)
  DEGs2KO = data.frame(DEG=rep(DEGs2KO$query,time),KO=unlist(KO))
  DEGs2KO$KO=sapply(DEGs2KO$KO,gsub,pattern="ko:",replacement="")
  DEGs2KO=merge(DEGs2KO,KO2Des,by="KO",all.x=TRUE)
  
  # DEGs, KOs and pathways
  DEGs2KO2Pathway=merge(DEGs2KO,DEG2Pathway,by="DEG",all=TRUE)
  func=function(i){
    d=DEGs2KO2Pathway[i,]
    TruePathway=subset(Pathway2KO,KO==d[,"KO"])[,"Pathway"]
    if (d[,"KEGG_Pathway"]%in%TruePathway){return(TRUE)}else{return(FALSE)}
  }
  DEGs2KO2Pathway[,"panduan"]=sapply(1:nrow(DEGs2KO2Pathway),func)
  DEGs2KO2Pathway=subset(DEGs2KO2Pathway,panduan) # Remove wrong combinations of KOs and pathways
  
  out=data.frame(DEG=DEGs2KO2Pathway$DEG,
                 KO=DEGs2KO2Pathway$KO,KO_description=DEGs2KO2Pathway$KO_Description,
                 Pathway=DEGs2KO2Pathway$KEGG_Pathway,
                 Pathway_description=DEGs2KO2Pathway$Description,
                 Pathway_padj=DEGs2KO2Pathway$AdjustedP)
  out=out[!duplicated(out),]
  out=out[order(out$Pathway_padj),]
  write.table(out,paste(output_prefix,".txt",sep=""),sep="\t",quote=FALSE,row.names = FALSE)
  return(paste(output_prefix,".txt",sep=""))
}

# Bar plot for 30 most significantly enriched KEGG pathways.
enrichKEGGplot = function(df=df, # Tabular table with header.
                                 # Columns include "Pathway_description" (character) and "AdjustedP" (numeric).
                          output_prefix=output_prefix){
  #df="Results/down_KEGGenriched.txt";output_prefix="Results/down_KEGG"
  df=read.table(df,sep="\t",header=TRUE,quote="")
  df=data.frame(Description=df$Pathway_description,AdjustedP=df$Pathway_padj)
  DEGcount=table(df$Description)
  DEGcount=data.frame(Description=names(DEGcount),count=unname(DEGcount))
  df=df[!duplicated(df),]
  df=merge(df,DEGcount,by="Description",all=TRUE)
  colnames(df)=c("Description","AdjustedP","buzhidao","DEGcount")
  df=df[order(df$AdjustedP),]
  
  df=head(df,30) # 30 most significant
  df$Description=factor(df$Description,levels=df$Description)
  library(ggplot2)
  pdf(paste(output_prefix,".pdf",sep=""),height=12,width=10)
  plot=ggplot(df,aes(y=DEGcount,x=Description))+
    geom_bar(stat="identity",aes(fill=-log10(df$AdjustedP)))+
    ylab("DEG count")+
    xlab("")+
    labs(fill="-log10(Adjusted P)")+
    scale_y_continuous(expand=c(0,0))+
    theme_classic()+
    theme(plot.margin=unit(c(1,1.5,1,1),"lines"),
          axis.title.x=element_text(size=8,face="bold"),
          axis.title.y=element_text(size=8,face="bold"),
          axis.text.x=element_text(size=11,face="bold",angle=70,hjust=1),
          axis.text.y=element_text(size=8,face="bold"),
          title=element_text(size=8,face="bold"),
          legend.text=element_text(size=8,face="bold"),
          legend.title = element_text(size=8,face="bold"),
          legend.position="left")
  print(plot)
  dev.off()
  return(plot)
}

# Visualize pathways enriched by DEGs.
PathwayView=function(up_KO,down_KO, # Two tables of same format (tab, header). 
                                    # Column "Pathway", which represents KEGG pathay IDs (e.g. map00330) enriched by up-/down-regulated genes.
                     FC=FC, # Tabular table with header and row names.
                            # Row names are gene names.
                            # Column "color" indicates whether gene expression is "up", "down" or "none" changed.
                     Genome=Genome # output of eggNOG-mapper
                                   # fields: query,KEGG_ko
                                   # annotations of reference genome
  ){
  # Pathways to be plotted.
  up_KO=read.table(up_KO,sep="\t",header=TRUE,quote="")
  down_KO=read.table(down_KO,sep="\t",header=TRUE,quote="")
  KO=rbind(up_KO,down_KO)
  KEGGs=KO[,"Pathway"];KEGGs=KEGGs[!duplicated(KEGGs)]
  KEGGs=gsub("map","",KEGGs)
  
  # DEGs
  FC=read.table(FC,sep="\t",header=TRUE,quote="")
  FC=data.frame(Gene=rownames(FC),log2fc=FC[,"log2FoldChange"],color=FC[,"color"])
  
  Genome=read.table(Genome,sep="\t",header=TRUE,quote="")
  Genome=data.frame(Gene=Genome[,"query"],KO=Genome[,"KEGG_ko"])
  Genome=Genome[Genome[,"KO"]!="-",]
  Genome=merge(Genome,FC,by="Gene",all=TRUE)
  Genome=Genome[!is.na(Genome[,"KO"]),]
  
  library(stringr);library(dplyr)
  KOs = str_split(Genome[,"KO"],",")
  time = sapply(KOs,length) 
  Genome = data.frame(Gene=rep(Genome[,"Gene"],time),
                      log2fc=rep(Genome[,"log2fc"],time),
                      color=rep(Genome[,"color"],time),
                      KO=unlist(KOs))
  Genome=Genome[!is.na(Genome[,"log2fc"]),]
  
  Genome[Genome[,"color"]=="none","log2fc"]=0
  Genome[Genome[,"color"]=="up","log2fc"]=1
  Genome[Genome[,"color"]=="down","log2fc"]=-1
  func=function(vector){
    f=0
    if (-1 %in% vector){f=-1}
    if (1 %in% vector){f=1}
    return(f)
  }
  gene=tapply(Genome$log2fc,Genome$KO,func)
  IDs=gsub("ko:","",names(gene));gene=as.vector(gene)
  names(gene)=IDs
  
  library(pathview)
  pathview(gene.data=gene,
           pathway.id=KEGGs,
           species="nve",
           node.sum="max",
           limit=list(gene=1),
           low=list(gene="blue"),
           mid=list(gene="grey"),
           high=list(gene="red"))
}

# KO enrichment
enrichKO = function(Background=Background, # Tabular table with header and row names.
                                           # Row names are gene names.
                                           # Column "color" indicates whether gene expression is "up", "down" or "none" changed.
                    up=up,# logical. If true, enrich up-regulated genes.
                    Annotation=Annotation, # output of eggNOG-mapper, 
                                           # annotations of reference genome
                    KOmapper=KOmapper, # From MinPath github: https://github.com/mgtools/MinPath/blob/master/data/KEGG-family.txt
                    output_prefix=output_prefix){
  Background = read.csv(Background,sep="\t",header=TRUE)
  Annotation = read.csv(Annotation,sep="\t",header=TRUE)
  Annotation = data.frame(query=Annotation[,"query"],term=Annotation[,"KEGG_ko"])
  if (up){
    DEGs = rownames(subset(Background,color=="up"))
  }else{
    DEGs = rownames(subset(Background,color=="down"))
  }
  
  DEGs2Term = subset(Annotation,Annotation[,"query"]%in%DEGs)
  
  TermOfDEGs = unname(unlist(sapply(DEGs2Term[,"term"],strsplit,split=",")))
  TermOfDEGs = sort(unname(TermOfDEGs[which(TermOfDEGs!="-")]))
  Term2DEGCount = table(TermOfDEGs)
  TermOfDEGs_nr = TermOfDEGs[!duplicated(TermOfDEGs)]
  
  TermOfAnnotation = unname(unlist(sapply(Annotation[,"term"],strsplit,split=",")))
  Term2AnnotationCount = table(TermOfAnnotation)
  
  # Enrichment by Fisher exact test
  fisher = function(term){
    DEG_term = Term2DEGCount[term]
    DEG_NOterm = nrow(DEGs2Term)-DEG_term
    NoDEG_term = Term2AnnotationCount[term]-DEG_term
    NoDEG_NOterm = nrow(Background)-DEG_term-DEG_NOterm-NoDEG_term
    p=fisher.test(matrix(c(DEG_term,NoDEG_term,DEG_NOterm,NoDEG_NOterm), nrow = 2), alternative = "greater")$p.value
    return(p)
  }
  degCount = function(term){
    DEG_term = Term2DEGCount[term]
    DEG_NOterm = nrow(DEGs2Term)-DEG_term
    NoDEG_term = Term2AnnotationCount[term]-DEG_term
    NoDEG_NOterm = nrow(Background)-DEG_term-DEG_NOterm-NoDEG_term
    #expected=(DEG_term+DEG_NOterm)*(DEG_term+NoDEG_term)/(DEG_term+DEG_NOterm+NoDEG_term+NoDEG_NOterm)
    return(DEG_term)
  }
  pvalues = sapply(TermOfDEGs_nr,fisher)
  padjs = p.adjust(pvalues,"BH")
  DegCount = sapply(TermOfDEGs_nr,degCount)
  o = data.frame(Term=TermOfDEGs_nr,pvalue=pvalues,padj=padjs,DegCount=DegCount)
  o = subset(o,padj<0.05)
  o$Term = sapply(o$Term,gsub,pattern="ko:",replacement="")
  
  # Enriched KO description
  KO2Description=read.csv(KOmapper,sep="\t",quote="",header=FALSE)
  o = merge(o,KO2Description,by.x="Term",by.y="V1",all.x=TRUE)
  o = subset(o,!is.na(o$V2))
  o = o[order(o[,"padj"]),]
  
  # DEG to KO
  library(stringr)
  library(dplyr)
  KOs = str_split(DEGs2Term$term,",")
  time = sapply(KOs,length) 
  DEG2KO = data.frame(DEG=rep(DEGs2Term$query,time),KO=unlist(KOs))
  DEG2KO = subset(DEG2KO,DEG2KO$KO!="-")
  DEG2KO$KO = sapply(DEG2KO$KO,gsub,pattern="ko:",replacement="")
  
  DEG2KO_des = merge(o,DEG2KO,by.x="Term",by.y="KO",all.x=TRUE)
  DEG2KO_des = DEG2KO_des[order(DEG2KO_des$padj),]
  write.table(data.frame(DEG=DEG2KO_des$DEG,KO=DEG2KO_des$Term,
                         Description=DEG2KO_des$V2,
                         AdjustedP=DEG2KO_des$padj),
              paste(output_prefix,"enriched.txt",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
  return(0)
}

# Bar plot of 30 most significantly enriched KOs.
enrichKOplot = function(df=df, # Tabular table with header.
                               # Columns include "Description" (character) and "AdjustedP" (numeric).
                        output_prefix=output_prefix){
  df=read.table(df,sep="\t",header=TRUE,quote="")
  df=data.frame(Description=df$Description,AdjustedP=df$AdjustedP)
  DEGcount=table(df$Description)
  DEGcount=data.frame(Description=names(DEGcount),count=unname(DEGcount))
  df=df[!duplicated(df),]
  df=merge(df,DEGcount,by="Description",all=TRUE)
  colnames(df)=c("Description","AdjustedP","buzhidao","DEGcount")
  df=df[order(df$AdjustedP),]
  
  df=head(df,30) # 30 most significant
  df$Description=factor(df$Description,levels=df$Description)
  library(ggplot2)
  pdf(paste(output_prefix,".pdf",sep=""),height=12,width=10)
  plot=ggplot(df,aes(y=DEGcount,x=Description))+
    geom_bar(stat="identity",aes(fill=-log10(df$AdjustedP)))+
    ylab("DEG count")+
    xlab("")+
    labs(fill="-log10(Adjusted P)")+
    scale_y_continuous(expand=c(0,0))+
    theme_classic()+
    theme(plot.margin=unit(c(1,1.5,1,1),"lines"),
          axis.title.x=element_text(size=8,face="bold"),
          axis.title.y=element_text(size=8,face="bold"),
          axis.text.x=element_text(size=10,face="bold",angle=70,hjust=1),
          axis.text.y=element_text(size=8,face="bold"),
          title=element_text(size=8,face="bold"),
          legend.text=element_text(size=8,face="bold"),
          legend.title = element_text(size=8,face="bold"),
          legend.position="left")
  print(plot)
  dev.off()
  return(plot)
}
