arg=commandArgs(trailingOnly = TRUE)

main.R=arg[1]
main.py=arg[2]
conf=arg[3]

source(main.R)
source(conf)

out_dir=sub("/$","",out_dir)
threads=as.character(threads)
if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}

#####
# Map reads to assembly
#####
if (!file.exists(paste(out_dir,"/minimap2",sep=""))){
  system(paste("mkdir"," ",out_dir,"/minimap2",sep=""))
}
if (!file.exists(paste(out_dir,"/minimap2/",label,sep=""))){
  system(paste("mkdir"," ",out_dir,"/minimap2/",label,sep=""))
}
if (file.exists(paste(out_dir,"/minimap2/",label,"/",label,"_minimap2.finished",sep=""))){
  print("Step 1: minimap2 FINISHED")
}else{
  print("Step 1: minimap2 START")
  minimap2(long_reads=paste(fq1,fq2,sep=" "), # space-separated list for PE
           lr_type="sr", # long read type. 
                    # "map-pb" for PacBio
                    # "map-hifi" for HiFi
                    # "map-ont" for ONT reads.
                    # "sr" for NGS
                    # "asm5" for accurate reads diverging <5% to assembly
            assembly=genome,
            out_prefix=paste(out_dir,"/minimap2/",label,"/",label,sep=""),
            threads=threads)
  system(paste("touch"," ",out_dir,"/minimap2/",label,"/",label,"_minimap2.finished",sep=""))
}
#####
# SprayNPray
#####
if (!file.exists(paste(out_dir,"/SprayNPray",sep=""))){
  system(paste("mkdir"," ",out_dir,"/SprayNPray",sep=""))
}
if (file.exists(paste(out_dir,"/SprayNPray/",label,"/",label,"_SprayNPray.finished",sep=""))){
  print("Step 2: SprayNPray FINISHED")
}else{
  print("Step 2: SprayNPray START")
  SprayNPray(fna=genome, # fna. Input DNA sequences.
             bam=paste(out_dir,"/minimap2/",label,"/",label,".bam",sep=""), # BAM. For coverage computation.
             ref=nr.dmnd, # Diamond database.
             blast="none",
             out_basename=label, # cannot be path.
                      # output files are stored a directory named as out_basename
             out_dir=paste(out_dir,"/SprayNPray",sep=""), # the directory where out_basename are moved to.
             threads=threads)
  system(paste("touch"," ",out_dir,"/SprayNPray/",label,"/",label,"_SprayNPray.finished",sep=""))
}
#####
# diamond-megan
#####
if (!file.exists(paste(out_dir,"/diamond_megan",sep=""))){
  system(paste("mkdir"," ",out_dir,"/diamond_megan",sep=""))
}
if (!file.exists(paste(out_dir,"/diamond_megan/",label,sep=""))){
  system(paste("mkdir"," ",out_dir,"/diamond_megan/",label,sep=""))
}
if (file.exists(paste(out_dir,"/diamond_megan/",label,"/",label,"_diamond_megan.finished",sep=""))){
  print("Step 3: diamond_megan FINISHED")
}else{
  print("Step 3: diamond_megan START")
  Diamond_Megan(fna=genome, # fna. Input DNA sequences.
                out_basename=label,
                blast_dir=paste(out_dir,"/diamond_megan/",label,sep=""), # Directory for diamond output.
                rma_dir=paste(out_dir,"/diamond_megan/",label,sep=""), # Directory for rma output of megan.
                assignment_dir=paste(out_dir,"/diamond_megan/",label,sep=""), # Directory for taxonomy table.
                ref_diamond=nr.dmnd, # Diamond database.
                ref_megan=megan.db, # Megan database.
                threads=threads)
  system(paste("touch"," ",out_dir,"/diamond_megan/",label,"/",label,"_diamond_megan.finished",sep=""))
}
#####
# Classification
#####
library(reticulate)
library(stringr)
source_python(main.py)
if (!file.exists(paste(out_dir,"/sklearn",sep=""))){
  system(paste("mkdir"," ",out_dir,"/sklearn",sep=""))
}
if (file.exists(paste(out_dir,"/sklearn/",label,"_classification.finished",sep=""))){
  print("Step 4: Classification FINISHED")
}else{
  print("Step 4: Classification START")
  spr=read.csv(paste(out_dir,"/SprayNPray/",label,"/",label,".csv",sep=""),header=TRUE,quote="")
  spr=spr[,c("contig","contig_length","hits_per_kb","cov","GC.content")]
  colnames(spr)=c("contig","length","cds_density","coverage","GC")

  dm=read.table(paste(out_dir,"/diamond_megan/",label,"/",label,".tsv",sep=""),sep="\t",quote="")
  colnames(dm)=c("contig","rank","path")
  dm=dm[dm[,"path"]!="",]
  dm[,"assign"]=str_extract(dm[,"path"],"\\[P\\] .*;")
  dm[,"assign"]=sub("; .*$","",dm[,"assign"])
  dm[,"assign"]=sub(";","",dm[,"assign"])
  dm[,"assign"]=sub("\\[P\\] ","",dm[,"assign"])
  dm=dm[,c("contig","assign")]
  dm=dm[!is.na(dm[,"assign"]),]

  df=merge(dm,spr,by="contig",all=TRUE)
  df=df[df[,"length"]>400,] # Remove contigs below 400 bp

  training=df[!is.na(df[,"assign"]),]
  training_contig=training[,"contig"]
  training_sample=training[,c("cds_density","coverage","GC")]
  training_results=training[,"assign"]

  testing=df[is.na(df[,"assign"]),]
  testing_contig=testing[,"contig"]
  testing_sample=testing[,c("cds_density","coverage","GC")]
  # decision tree
  testing_results=decision_tree_classifier(training_sample=training_sample,
                                           training_results=training_results,
                                           testing=testing_sample)
  testing[,"assign"]=testing_results

  sk=rbind(training,testing)  
  write.table(sk,paste(out_dir,"/sklearn/",label,".tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)  
  system(paste("touch"," ",out_dir,"/sklearn/",label,"_classification.finished",sep=""))
}
#####
# Split target contigs and evaluate quality
#####
if (!file.exists(paste(out_dir,"/retrievedGenome",sep=""))){
  system(paste("mkdir"," ",out_dir,"/retrievedGenome",sep=""))
}
cmd=paste("awk -F '\t' -v OFS='\t' '{if ($2==\"",target,"\") print $0}' ",
          out_dir,"/sklearn/",label,".tsv"," > ",
          out_dir,"/retrievedGenome/",label,"_",target,".tsv",sep="")
print(cmd);system(cmd,wait=TRUE)
cmd=paste("awk -F '\t' -v OFS='\t' '{if ($2==\"",target,"\") print $1}' ",
          out_dir,"/sklearn/",label,".tsv"," > ",
          out_dir,"/retrievedGenome/",label,".lst",sep="")
print(cmd);system(cmd,wait=TRUE)
cmd=paste("seqkit grep -f ",out_dir,"/retrievedGenome/",label,".lst ",genome,
          " > ",out_dir,"/retrievedGenome/",label,"_",target,".fna",sep="")
print(cmd);system(cmd,wait=TRUE)
system(paste("rm"," ",out_dir,"/retrievedGenome/",label,".lst",sep=""))

Quast(fna=paste(out_dir,"/retrievedGenome/",label,"_",target,".fna",sep=""), # FASTA of genome assembly.
      out_dir=paste(out_dir,"/retrievedGenome/",label,"_QUAST/",sep=""),
      threads=threads)
BUSCO(fna=paste(out_dir,"/retrievedGenome/",label,"_",target,".fna",sep=""), # Fasta file of nucleotide or protein.
               # Be consistent with Mode.
      Mode="genome", # genome/proteins/transcriptome
      Lineage=busco.db, # Lineage dataset, e.g. insecta_odb10
               # Available datasets: https://busco-data.ezlab.org/v5/data/lineages/
               # BUSCO will download lineage dataset automatically.
      Out_prefix=paste(label,"_BUSCO",sep=""), # Give the analysis run a recognisable short name.
               # Output folders and files will be labelled with this name.
               # Cannot be path
      out_dir=paste(out_dir,"/retrievedGenome/",sep=""),
      Threads=threads)
library(ggplot2);library(ggExtra)
df=read.table(paste(out_dir,"/retrievedGenome/",label,"_",target,".tsv",sep=""),sep="\t",header=TRUE,quote="")
p1=ggplot(df,aes(x=coverage,y=GC))+
  geom_point(size=0.01)+
  theme_classic()+
  labs(title=label,x="coverage",y="GC%")+
  scale_x_continuous(limits=c(0,max(df$cov)),expand=c(0,0))+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))
p=ggMarginal(p1,type="histogram",size=10,
             xparams=list(bins=100),
             yparams=list(bins=100))
pdf(paste(out_dir,"/retrievedGenome/",label,"_",target,".pdf",sep=""));print(p);dev.off()
