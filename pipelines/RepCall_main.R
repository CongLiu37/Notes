# MITE-Hunter
# Dependencies: MITE-hunter, seqkit
mite_hunter=function(genome=genome,
                     out_dir=out_dir,
                     out_basename=out_basename,
                     threads=threads){
  threads=as.character(threads)
  wd=getwd()
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  setwd(out_dir)
  
  cmd=paste("seqkit","seq","-u",
            "-j",threads,
            genome,">",
            paste(out_dir,"/genome.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("MITE_Hunter_manager.pl",
            "-i",paste(out_dir,"/genome.fa",sep=""),
            "-g",out_basename,
            "-n",threads,
            "-S 12345678",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("cat",
            paste(out_basename,"_Step8*",sep=""),
            ">",
            paste(out_basename,"_MITE.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm"," ",out_basename,".*",sep=""))
  system(paste("rm"," ",out_basename,"_raw*",sep=""))
  system(paste("rm"," ",out_basename,"_[0-9]*",sep=""))
  system(paste("rm"," ",out_basename,"_Step*",sep=""))
  system("rm genome.*")
  system("rm error.log")
  system("rm formatdb.log")
  setwd(wd)
}

# LTR_finder & LTR_Harvest (gt) & LTR_retriever
ltr=function(fna=fna,
             out_dir=out_dir,
             out_prefix=out_prefix,
             threads=threads){
  threads=as.character(threads)
  wd=getwd()
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  setwd(out_dir)
  
  cmd=paste("seqkit","seq","-u",
            "-j",threads,
            genome,">",
            paste(out_dir,"/genome.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  fna=paste(out_dir,"/genome.fa",sep="")
  
  if (!file.exists(paste(fna,".finder.combine.scn",sep=""))){
    cmd=paste("LTR_FINDER_parallel",
              "-seq",fna,
              "-harvest_out",
              "-threads",threads,
              "-size 1000000 -time 300",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  if (!file.exists(paste(out_prefix,".LTRharvest",sep=""))){
    cmd=paste("gt suffixerator",
              "-db",fna,
              "-indexname",paste(out_prefix,".suffixerator",sep=""),
              "-tis -suf -lcp -des -ssp -sds -dna",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
    cmd=paste("gt ltrharvest",
              "-index",paste(out_prefix,".suffixerator",sep=""),
              "-minlenltr 100 -maxlenltr 7000 -mintsd 4",
              "-maxtsd 6 -motif TGCA -motifmis 1 -similar 85",
              "-vic 10 -seed 20 -seqids yes",
              ">",paste(out_prefix,".LTRharvest",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  cmd=paste("cat",
            paste(fna,".finder.combine.scn",sep=""),
            paste(out_prefix,".LTRharvest",sep=""),
            "> rawLTR.scn")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("LTR_retriever",
            "-genome",fna,
            "-inharvest","rawLTR.scn",
            "-threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("mv genome.fa.LTR.gff3 LTR.gff3")
  system("mv genome.fa.LTRlib.fa LTR.fa")
  
  system(paste("rm",fna,sep=" "),wait=TRUE)
  system(paste("rm"," ",label,".*"))
  system("rm genome.fa.*")
  system("rm rawLTR.scn")
  setwd(wd)
}

# RepeatMasker
repeatmasker=function(fna=fna,# Fasta file of genome.
                      out_dir=out_dir,
                      RepeatLib.fa=RepeatLib.fa, # space-list
                      Threads=Threads){
  Threads=as.character(Threads)
  wd=getwd()
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  setwd(out_dir)
  
  cmd=paste("seqkit","seq","-u",
            "-j",threads,
            fna,">",
            paste(out_dir,"/genome.fa",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  fna=paste(out_dir,"/genome.fa",sep="")
  
  cmd=paste("cat",RepeatLib.fa,"> repeat.fa",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  RepeatLib.fa=paste(out_dir,"/repeat.fa",sep="")
  
  #####################################################################
  # RepeatModeler & RepeatMasker installed in Singularity container
  path="singularity run /home/c/c-liu/Softwares/dfam-tetools-latest.sif"
  #####################################################################
  
  # RepeatMasker: Mask repeats found by RepeatModeler
  cmd=paste(path,
            "RepeatMasker",
            "-xsmall", # soft masking
            "-lib",RepeatLib.fa,
            "-pa",Threads,
            "-gff",
            fna,
            "-dir",out_dir)
  print(cmd);system(cmd,wait=TRUE)
  
  # # RepeatMasker: Mask repeats in RepBase
  # cmd=paste(path,
  #           "RepeatMasker",
  #           "-xsmall", # soft masking
  #           "-pa",Threads,
  #           "-gff",
  #           paste(fna,".masked",sep=""),
  #           "-dir",out_dir,
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  
  # Remove temporaries
  system(paste("rm",fna,sep=" "),wait=TRUE)
  system(paste("rm",RepeatLib.fa,sep=" "),wait=TRUE)
  setwd(wd)
  
  # masked genome
  return(paste(out_dir,"/",fna,".masked.masked",sep=""))
}

# RepeatModeler (no LTR)
# Dependencies: RepeatModeler, Singularity
repeatmodeler_noLTR=function(fna=fna,
                             out_dir=out_dir,
                             out_prefix=out_prefix,
                             Threads=Threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  pwd_begin=getwd();setwd(out_dir)
  Threads=as.character(Threads)
  
  system(paste("cp",fna,out_dir,sep=" "),wait=TRUE)
  fna_name=unlist(strsplit(fna,"/"));fna_name=fna_name[length(fna_name)]
  fna=paste(out_dir,"/",fna_name,sep="")
  #####################################################################
  # RepeatModeler & RepeatMasker installed in Singularity container
  path="singularity run /home/c/c-liu/Softwares/dfam-tetools-latest.sif"
  #####################################################################
  
  # RepeatModeler: de novo repeat library
  cmd=paste(path,
            "BuildDatabase",
            "-name",paste(out_prefix,"_RepeatModeler.db",sep=""),
            "-engine","ncbi",
            fna,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste(path,
            "RepeatModeler",
            "-database",paste(out_prefix,"_RepeatModeler.db",sep=""),
            "-engine","ncbi",
            "-pa",Threads,
            # "-LTRStruct",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",fna,sep=" "),wait=TRUE)
  system("rm -r RM_*")
  
  cmd=paste("seqkit grep -vnrp '#LTR' ",label,"_RepeatModeler.db-families.fa",
            " > noLTR_rep.lib.fa")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(pwd_begin)
  return(paste(out_dir,"/noLTR_rep.lib.fa",sep=""))
}

# Final repeat annotation: repeatmasker out to gff3
rep_elements=function(repeatmasker.out=repeatmasker.out,
                      # Clarify source of repeat elements
                      MITE_Hunter.replib=MITE_Hunter.replib, # RepLib from MITE-Hunter
                      LTR_retriever.replib=LTR_retriever.replib, # RepLib from LTR_retriever
                      RepeatModeler.replib=RepeatModeler.replib, # RepLib from RepeatModeler
                      all.replib=all.replib,
                      out.gff3=out.gff3){
  cmd=paste("cat",MITE_Hunter.replib,LTR_retriever.replib,RepeatModeler.replib,">",all.replib,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("awk 'NR > 3 {if ($11 != \"rRNA\") print $5\"\tRepeatMasker\trepeat_element\t\"$6\"\t\"$7\"\t\"$1\"\t\"$9\"\t\\.\tSequence=\"$10\";Repeat_class=\"$11}'",
            repeatmasker.out,">",out.gff3,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("sed -i 's/Repeat_class=Unspecified/Repeat_class=Unknown/'",out.gff3,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # MITE
  cmd=paste("grep '>' ",MITE_Hunter.replib," | sed 's/>//' | sed 's/ Len:.*$//' | sed 's/Unknow.*/Unknown/'",sep="")
  mite=system(cmd,intern=TRUE)
  for (i in mite){
    i=unlist(strsplit(i," "))
    target=paste("Sequence=",i[1],";Repeat_class=.*$",sep="")
    res=paste("Sequence=",i[1],";Source=MITE-Hunter;Repeat_class=",i[2],sep="")
    cmd=paste("sed -i 's/",target,"/",res,"/'"," ",out.gff3,sep="")
    system(cmd,wait=TRUE)
  }
  # LTR
  cmd=paste("grep '>' ",LTR_retriever.replib," | sed 's/>//'",sep="")
  ltr=system(cmd,intern=TRUE)
  for (i in ltr){
    i=unlist(strsplit(i,"#"))
    target=paste("Sequence=",i[1],";",sep="")
    res=paste("Sequence=",i[1],";Source=LTR_retriever;",sep="")
    cmd=paste("sed -i 's/",target,"/",res,"/'"," ",out.gff3,sep="")
    system(cmd,wait=TRUE)
  }
  # RepeatModeler
  cmd=paste("grep '>' ",RepeatModeler.replib," | sed 's/>//' | sed 's/ .*//'",sep="")
  RM=system(cmd,intern=TRUE)
  for (i in RM){
    i=unlist(strsplit(i,"#"))
    if (!grepl("LTR/.*",i[2])){
      target=paste("Sequence=",i[1],";",sep="")
      res=paste("Sequence=",i[1],";Source=RepeatModeler;",sep="")
      cmd=paste("sed -i 's/",target,"/",res,"/'"," ",out.gff3,sep="")
      system(cmd,wait=TRUE)
    }
  }
  
  cmd=paste("awk -F '\t' -v OFS='\t' '{if ($7==\"C\") $7=\"-\"; print $0}' ",out.gff3,
            " > final.gff3",
            sep="")
  print(cmd);system(cmd,wait=TRUE)
  system(paste("rm",out.gff3,sep=" "))
  system(paste("mv","final.gff3",out.gff3,sep=" "))
}