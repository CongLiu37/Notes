rm(list=ls())
source("main.R")
args = commandArgs(T)

metadata = args[1] # A non-header tabular table whose 
                   # first column is sample name, 
                   # second column is a list of paths to fq files seperated by comma, 
                   # third column is sample description.
out_dir =  args[2] # Directory for temporary files.
index =    args[3] # Hisat2 index of reference genome.
gtf =      args[4] # GTF of reference genome.
threads =  args[5]


metadata = read.table(metadata,sep="\t",header=FALSE)
for (i in 1:nrow(metadata)){
    sample = metadata[i,1]

    raw_fqs = unlist(strsplit(metadata[i,2],","))
    raw_fq1 = raw_fqs[1]
    if (length(raw_fqs)!=1){
        raw_fq2 = raw_fqs[2] # pair-end
    }else{
        raw_fq2 = "None" # single-end
    }

    # QualityCheck(fq1=raw_fq1,
    #              fq2=raw_fq2,
    #              out_dir=out_dir,
    #              threads=threads)

    QualityFilter(fq1=raw_fq1,
                  fq2=raw_fq2,
                  clean_fq1=paste(out_dir,"/",sample,"_clean.1.fq.gz",sep=""),
                  clean_fq2=paste(out_dir,"/",sample,"_clean.2.fq.gz",sep=""),
                  unpaired_fq1=paste(out_dir,"/",sample,"_unpaired.1.fq.gz",sep=""),
                  unpaired_fq2=paste(out_dir,"/",sample,"_unpaired.2.fq.gz",sep=""),
                  threads=threads)

    #system(paste("rm",paste(out_dir,"/",sample,"_unpaired.1.fq.gz",sep=""),sep=" "),wait=TRUE)
    #system(paste("rm",paste(out_dir,"/",sample,"_unpaired.2.fq.gz",sep=""),sep=" "),wait=TRUE)

    QualityCheck(fq1=paste(out_dir,"/",sample,"_clean.1.fq.gz",sep=""),
                 fq2=paste(out_dir,"/",sample,"_clean.2.fq.gz",sep=""),
                 out_dir=out_dir,
                 threads=threads)

    #system(paste("rm",paste(out_dir,"/",sample,"_clean.1_fastqc.zip",sep=""),sep=" "),wait=TRUE)
    #system(paste("rm",paste(out_dir,"/",sample,"_clean.2_fastqc.zip",sep=""),sep=" "),wait=TRUE)

    Hisat(fq1=paste(out_dir,"/",sample,"_clean.1.fq.gz",sep=""),
          fq2=paste(out_dir,"/",sample,"_clean.2.fq.gz",sep=""),
          index=index,
          out_prefix=paste(out_dir,"/",sample,sep=""),
          threads=threads)

    #system(paste("rm",paste(out_dir,"/",sample,"_clean.1.fq.gz",sep=""),sep=" "),wait=TRUE)
    #system(paste("rm",paste(out_dir,"/",sample,"_clean.2.fq.gz",sep=""),sep=" "),wait=TRUE)

    StringTie(input_bam=paste(out_dir,"/",sample,".bam",sep=""),
              gtf=gtf,
              threads=threads)

    #system(paste("rm",paste(out_dir,"/",sample,".bam",sep=""),sep=" "),wait=TRUE)

    cmd = paste("echo ","\'",sample," ",paste(out_dir,"/",sample,".gtf",sep=""),"\' >>",out_dir,"/SampleList.txt",sep="")
    system(cmd)
}

run.prepDE.py(inputfile=paste(out_dir,"/SampleList.txt",sep=""))
