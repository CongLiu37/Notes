#####
# whole genome alignment
#####
cactus=function(seqFile.cactus="/flash/BourguignonU/Cong/termite_pca/cactus/seqFile.cactus",
                out.hal="/flash/BourguignonU/Cong/termite_pca/cactus/termite_genomes.hal",
                outDir="/flash/BourguignonU/Cong/termite_pca/cactus/outDir",
                jobstore="/flash/BourguignonU/Cong/termite_pca/cactus/jobstore",
                work_dir="/flash/BourguignonU/Cong/termite_pca/cactus/",
                stamps_dir="/flash/BourguignonU/Cong/termite_pca/cactus/stamps"){
  wd=getwd()
  setwd(work_dir)
  
  if (!file.exists(stamps_dir)){dir.create(stamps_dir)}
  
  cmd=paste("cactus-prepare",
            seqFile.cactus,
            "--outDir",outDir,
            "--outSeqFile",paste(outDir,"/seqFile_split.cactus",sep=""),
            "--outHal",out.hal,
            "--jobStore",jobstore,
            sep=" ")
  scripts=system(cmd,wait=TRUE,intern=TRUE)
  
  preprocess=scripts[grepl("^cactus-preprocess",scripts)]
  for (i in 1:length(preprocess)){
    stamp=paste(stamps_dir,"/preprocess_",as.character(i),".finished",sep="")
    if (!file.exists(stamp)){
      print(preprocess[i])
      system(preprocess[i],wait=TRUE)
      system(paste("touch",stamp,sep=" "))
    }
  }
  
  n.round=sum(grepl("^### Round",scripts))
  n.round=1:(n.round)
  for (i in n.round){
    s=which(grepl(paste("^### Round ",as.character(i-1),"$",sep=""),scripts))
    e=which(grepl(paste("^### Round ",as.character(i),"$",sep=""),scripts))
    if (length(e)==0){e=which(grepl("^## HAL merging",scripts))}
    cmds=scripts[s:e]
    cmds=cmds[cmds!="" & !grepl("#",cmds)]
    sapply(1:(length(cmds)/3),
           function(j){
             stamp=paste(stamps_dir,"/Round_",as.character(i),"_",
                         as.character(j),".sh.finished",sep="")
             shell=paste(stamps_dir,"/Round_",as.character(i),"_",
                         as.character(j),".sh",sep="")
             write(cmds[( ((j-1)*3)+1 ):(j*3)],file=shell,append=TRUE,sep="\n")
             write(paste("touch",stamp,sep=" "),file=shell,append=TRUE,sep="\n")
           })
  }
  
  FINISHED=system(paste("ls",stamps_dir,sep=" "),intern=TRUE)
  FINISHED=FINISHED[grepl("Round_[0-9]*_[0-9]*.sh$",FINISHED)]
  FINISHED.bool=file.exists(paste(stamps_dir,"/",FINISHED,".finished",sep=""))
  if (!all(FINISHED.bool)){print("Unfinished:");print(FINISHED[which(!FINISHED.bool)])}
  if (all(FINISHED.bool)){
    cmds=scripts[which(scripts=="## HAL merging"):length(scripts)]
    cmds=cmds[cmds!="" & !grepl("^#",cmds)]
    for (i in 1:length(cmds)){
      stamp=paste(stamps_dir,"/HAL.merge_",as.character(i),".finished",sep="")
      if (!file.exists(stamp)){
        print(cmds[i])
        system(cmds[i],wait=TRUE)
        system(paste("touch",stamp,sep=" "),wait=TRUE)
      }
    }
  }
  
  
  # scripts=scripts[scripts!="" & !grepl("^#",scripts)]
  # 
  # for (i in 1:length(scripts)){
  #   stamp=paste(stamps_dir,"/",as.character(i),".finished",sep="")
  #   if (!file.exists(stamp)){
  #     system(scripts[i],wait=TRUE)
  #     system(paste("touch",stamp,sep=" "))
  #   }
  # }
  
  setwd(wd)
}

cactus_hal2maf=function(hal=hal,
                        out_dir=out_dir,
                        maf=maf,
                        refGenome="Bori",
                        refSequences=NA, # space-lst
                        targetGenomes=NA, # comma-lst, e.g. Aaca,Aban
                        threads=threads){
  if (!dir.exists(out_dir)){dir.create(out_dir)}
  if (!dir.exists(paste(out_dir,"/workDir",sep=""))){dir.create(paste(out_dir,"/workDir",sep=""))}
  if (!dir.exists(paste(out_dir,"/coordinationDir",sep=""))){dir.create(paste(out_dir,"/coordinationDir",sep=""))}
  cmd=paste("cactus-hal2maf",
            "--refGenome",refGenome,
            "--noAncestors",
            "--workDir",paste(out_dir,"/workDir",sep=""),
            "--coordinationDir",paste(out_dir,"/coordinationDir",sep=""),
            "--chunkSize 500000",
            "--batchCores",as.character(threads),
            "--coverage",
            "--filterGapCausingDupes",
            paste(out_dir,"/jobStore",sep=""),hal,maf,
            sep=" ")
  if (!is.na(refSequences)){cmd=paste(cmd,"--refSequence",refSequences,sep=" ")}
  if (!is.na(targetGenomes)){cmd=paste(cmd,"--targetGenomes",targetGenomes,sep=" ")}
  print(cmd);system(cmd,wait=TRUE)
}

mafExtractor=function(in.maf=in.maf,
                      seq="Cges.scaffold_1",
                      Start=0, # inclusive, 0 based
                      End=100, # inclusive, 0 based
                      out.maf=out.maf,){
  cmd=paste("mafExtractor",
            "--maf",in.maf,
            "--seq",seq,
            "--start",as.character(Start),
            "--stop",as.character(End),
            "--first --soft >",
            out.maf,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

mafDuplicateFilter=function(in.maf=in.maf,
                            out.maf=out.maf){
  cmd=paste("mafDuplicateFilter -k",
            "--maf",in.maf,
            ">",out.maf,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

maf_for_phyloP=function(maf=maf,
                        ref_contig2length=ref_contig2length, # seqkit fx2tab -n -l <genome.fna> ## identical with maf
                        out_prefix=ouout_prefix,
                        threads=threads){
  if (!dir.exists(out_prefix)){dir.create(out_prefix)}
  
  contig.lst=read.table(ref_contig2length,header=FALSE,quote="",sep="\t")
  rownames(contig.lst)=contig.lst$V1
  contig.lst=contig.lst[order(-contig.lst$V2),]
  
  library(parallel)
  clus=makeCluster(threads)
  clusterExport(clus,"contig.lst",envir=environment())
  clusterExport(clus,"maf",envir=environment())
  clusterExport(clus,"out_prefix",envir=environment())
  parSapply(clus,contig.lst$V1,
            function(seq){
              cmd=paste("mafExtractor",
                        "--maf",maf,
                        "--seq",seq,
                        "--start 0",
                        "--stop",as.character(contig.lst[seq,2]-1),
                        "--first --soft |",
                        "mafDuplicateFilter -k",
                        "--maf -",
                        ">",paste(out_prefix,"/",seq,"_mafDuplicateFilter.maf",sep=""),
                        sep=" ")
              system(cmd,wait=TRUE)
            })
}



phastCons=function(maf_for_phyloP=maf_for_phyloP, # dir for maf files, e.g. Cges.scaffold_501_mafDuplicateFilter.maf
                   ref_contig2length=ref_contig2length, # seqkit fx2tab -n -l <fna>
                   ref.bed=ref.bed, # seqid e.g.: Cges.scaffold_501. phastCons elements overlapped with ref.bed are removed
                   out_prefix=out_prefix,
                   phyloFit.mod=phyloFit.mod,
                   threads=threads){
  if (!dir.exists(out_prefix)){dir.create(out_prefix)}
  
  ref_contig2length=read.table(ref_contig2length,header=FALSE,sep="\t",quote="")
  ref_contig2length=ref_contig2length[order(-ref_contig2length$V2),]
  
  contig.lst=ref_contig2length$V1
  
  f=function(seq){
    stick=paste(out_prefix,"/",seq,"_phastCons.finished",sep="")
    if (!file.exists(stick)){
      cmd=paste("phastCons",
                "--target-coverage 0.3 --expected-length 45 --rho 0.3", # Cons 124 Insects Track Settings https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=dm6&g=cons124way
                "--most-conserved",paste(out_prefix,"/",seq,"_phastCons.bed",sep=""),
                paste(maf_for_phyloP,"/",seq,"_mafDuplicateFilter.maf",sep=""),
                phyloFit.mod,
                ">",
                paste(out_prefix,"/",seq,"_phastCons.wig",sep=""),
                sep=" ")
      system(cmd,wait=TRUE)
      cmd=paste("touch",stick,sep=" ")
      system(cmd,wait=TRUE)
    }
  }
  
  library(parallel)
  clus=makeCluster(threads)
  clusterExport(clus,"contig.lst",envir=environment())
  clusterExport(clus,"f",envir=environment())
  clusterExport(clus,"out_prefix",envir=environment())
  clusterExport(clus,"phyloFit.mod",envir=environment())
  clusterExport(clus,"maf_for_phyloP",envir=environment())
  parSapply(clus,contig.lst[1:50],f)
  parSapply(clus,contig.lst[51:100],f)
  parSapply(clus,contig.lst[101:150],f)
  parSapply(clus,contig.lst[151:200],f)
  parSapply(clus,contig.lst[201:length(contig.lst)],f)
  
  if (file.exists(paste(out_prefix,"_phastCons.bed",sep=""))){
    system(paste("rm",paste(out_prefix,"_phastCons.bed",sep=""),sep=" "))
  }
  lst=paste(out_prefix,"/",contig.lst,"_phastCons.bed",sep="")
  for (i in lst){
    cmd=paste("cat",i," | awk -F '\t' -v OFS='\t' '{if ($3-$2 >=200) print $0}' >>",
              paste(out_prefix,"_phastCons.bed",sep=""),
              sep=" ")
    system(cmd,wait=TRUE)
  }
  
  phastCons.bed=read.table(paste(out_prefix,"_phastCons.bed",sep=""),
                           sep="\t",header=FALSE,quote="")
  phastCons.bed$V1=sub("_mafDuplicateFilter.*$","",phastCons.bed$V4)
  phastCons.bed$V4=sub("mafDuplicateFilter.","",phastCons.bed$V4)
  
  ref.bed=read.table(ref.bed,sep="\t",header=FALSE,quote="")

  clusterExport(clus,"ref.bed",envir=environment())
  clusterExport(clus,"phastCons.bed",envir=environment())
  phastCons.no_overlap=parSapply(clus,1:nrow(phastCons.bed),
                           function(i){
                             d=ref.bed[ref.bed$V1==phastCons.bed[i,1] &
                                       ref.bed$V2<phastCons.bed[i,3] &
                                       ref.bed$V3>phastCons.bed[i,2] &
                                       ref.bed$V6==phastCons.bed[i,6],]
                             return(nrow(d)==0)
                           })
  phastCons.bed=phastCons.bed[which(phastCons.no_overlap),]
  phastCons.bed$V2=format(phastCons.bed$V2,scientific=FALSE,trim=TRUE)
  phastCons.bed$V3=format(phastCons.bed$V3,scientific=FALSE,trim=TRUE)
  
  write.table(phastCons.bed,
              paste(out_prefix,"_phastCons.bed",sep=""),
              sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}

phastCons2CNE=function(hal="/bucket/BourguignonU/Cong/termite_pca/cactus/termite_genomes.hal",
                       ref.spp=c("Mdar","Hsjo","Ncas","Cges","Ofor","Apac","Munk"), # comma lst
                       phastCons.bed.dir="/bucket/BourguignonU/Cong/termite_pca/phastCons/", # ${ref.spp}_cne.bed from phastCons
                       target.sp=target.sp,
                       target.sp_mask.bed=target.sp_mask.bed, # CNE overlap elements in target.sp_mask.bed are removed
                       threads=threads,
                       out_dir=out_dir){
  if (!dir.exists(out_dir)){dir.create(out_dir)}
  ref.spp=unlist(strsplit(ref.spp,","))
  ref.spp=ref.spp[ref.spp!=target.sp]
  
  for (ref in ref.spp){
    cmd=paste("halLiftover --noDupes --outPSLWithName",
              hal,ref,paste(phastCons.bed.dir,"/",ref,"_cne.bed",sep=""),target.sp,
              paste(out_dir,"/",ref,"2",target.sp,".bedWihPsl",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)

    cmd=paste("awk -F '\t' -v OFS='\t' '{print $15,$17,$18,$1,1,$10}'",
              paste(out_dir,"/",ref,"2",target.sp,".bedWihPsl",sep=""),
              ">",
              paste(out_dir,"/",ref,"2",target.sp,".bed",sep=""),
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)

    cmd=paste("rm",paste(out_dir,"/",ref,"2",target.sp,".bedWihPsl",sep=""),sep=" ")
    print(cmd);system(cmd,wait=TRUE)

    df=read.table(paste(out_dir,"/",ref,"2",target.sp,".bed",sep=""),header=FALSE,sep="\t",quote="")
    df$V6=sub("^.","",df$V6)
    df$V2=format(df$V2,scientific=FALSE,trim=TRUE)
    df$V3=format(df$V3,scientific=FALSE,trim=TRUE)
    write.table(df,paste(out_dir,"/",ref,"2",target.sp,".bed",sep=""),
                sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  
  beds=paste(out_dir,"/",ref.spp,"2",target.sp,".bed",sep="")
  beds=paste(beds,collapse = " ")
  if (file.exists(paste(phastCons.bed.dir,"/",target.sp,"_cne.bed",sep=""))){
    beds=paste(beds,paste(phastCons.bed.dir,"/",target.sp,"_cne.bed",sep=""),sep=" ")
  }
  cmd=paste("cat ",beds," | sort -k1,1 -k2,2n > ",out_dir,"/",target.sp,".bed",sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("bedtools merge",
            "-i",paste(out_dir,"/",target.sp,".bed",sep=""),
            ">",
            paste(out_dir,"/",target.sp,"_merged.bed",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  
  merged.bed=paste(out_dir,"/",target.sp,"_merged.bed",sep="")
  merged.bed=read.table(merged.bed,header=FALSE,sep="\t",quote="")
  target.sp_mask.bed=read.table(target.sp_mask.bed,header=FALSE,sep="\t",quote="")
  
  library(parallel)
  clus=makeCluster(threads)
  clusterExport(clus,"merged.bed",envir=environment())
  clusterExport(clus,"target.sp_mask.bed",envir=environment())
  merged.bed.no_overlap=parSapply(clus,1:nrow(merged.bed),
                                 function(i){
                                   d=target.sp_mask.bed[target.sp_mask.bed$V1==merged.bed[i,1] &
                                                          target.sp_mask.bed$V2<merged.bed[i,3] &
                                                          target.sp_mask.bed$V3>merged.bed[i,2],]
                                   return(nrow(d)==0)
                                 })
  merged.bed=merged.bed[which(merged.bed.no_overlap),]
  merged.bed=merged.bed[merged.bed$V3-merged.bed$V2>=200,]
  merged.bed$V2=format(merged.bed$V2,scientific=FALSE,trim=TRUE)
  merged.bed$V3=format(merged.bed$V3,scientific=FALSE,trim=TRUE)
  merged.bed$V4=paste(target.sp,".",merged.bed$V1,"_",merged.bed$V2,"_",merged.bed$V3,"_CNE",sep="")
  
  write.table(merged.bed,paste(out_dir,"/",target.sp,"_merged_clean.bed",sep=""),
              row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
  return(paste(out_dir,"/",target.sp,"_merged_clean.bed",sep=""))
}

splitBed=function(bed=bed,
                  out_dir=out_dir){
  if (!dir.exists(out_dir)){dir.create(out_dir)}
  contig.lst=paste("awk '{print $1}'",bed,"| sort | uniq",sep=" ")
  contig.lst=system(contig.lst,intern=TRUE)
  split.bed=paste(out_dir,"/",contig.lst,".bed",sep="")
  sapply(1:length(contig.lst),
         function(i){
           cmd=paste("awk -F '\t' -v OFS='\t' ",
                     "'{if ($1==\"",contig.lst[i],"\") print $0}' ",
                     bed," > ",split.bed[i],
                     sep="")
           system(cmd)
         })
}

phyloP_features=function(maf_for_phyloP=maf_for_phyloP, # dir for maf files, e.g. Cges.scaffold_501_mafDuplicateFilter.maf
                         ref_contig2length=ref_contig2length, # seqkit fx2tab -n -l <fna>
                         features.dir=features.dir, # dir for bed, e.g. Cges.scaffold_21.bed
                         out_prefix=out_prefix,
                         test_branches=test_branches, # comma lst
                         test_subtree=NA,
                         phyloFit.mod=phyloFit.mod,
                         threads=threads){
  if (!dir.exists(out_prefix)){dir.create(out_prefix)}
  
  ref_contig2length=read.table(ref_contig2length,header=FALSE,sep="\t",quote="")
  ref_contig2length=ref_contig2length[order(-ref_contig2length$V2),]
  
  contig.lst=ref_contig2length$V1
  
  f=function(seq){
    stick=paste(out_prefix,"/",seq,".mafDuplicateFilter_phyloP.feature.finished",sep="")
    feature=paste(features.dir,"/",seq,".bed",sep="")
    if (!file.exists(stick) & file.exists(feature)){
      cmd=paste("phyloP",
                "-i MAF --method LRT --mode CONACC",
                "--features",feature,
                sep=" ")
      if (!is.na(test_branches)){
        cmd=paste(cmd,"--branch",test_branches,sep=" ")
      }
      if (!is.na(test_subtree)){
        cmd=paste(cmd,"--subtree",test_subtree,sep=" ")
      }
      
      cmd=paste(cmd,
                phyloFit.mod,
                paste(maf_for_phyloP,"/",seq,"_mafDuplicateFilter.maf",sep=""),
                ">",
                paste(out_prefix,"/",seq,".mafDuplicateFilter_phyloP.feature.out",sep=""),
                sep=" ")
      system(cmd,wait=TRUE)
      cmd=paste("touch",stick,sep=" ")
      system(cmd,wait=TRUE)
    }
  }
  
  library(parallel)
  clus=makeCluster(threads)
  clusterExport(clus,"contig.lst",envir=environment())
  clusterExport(clus,"features.dir",envir=environment())
  clusterExport(clus,"f",envir=environment())
  clusterExport(clus,"out_prefix",envir=environment())
  clusterExport(clus,"test_branches",envir=environment())
  clusterExport(clus,"phyloFit.mod",envir=environment())
  clusterExport(clus,"maf_for_phyloP",envir=environment())
  clusterExport(clus,"test_subtree",envir=environment())
  parSapply(clus,contig.lst[1:50],f)
  parSapply(clus,contig.lst[51:100],f)
  parSapply(clus,contig.lst[101:150],f)
  parSapply(clus,contig.lst[151:200],f)
  parSapply(clus,contig.lst[201:length(contig.lst)],f)
  
  out.lst=paste(out_prefix,"/",contig.lst,".mafDuplicateFilter_phyloP.feature.out",sep="")
  out.lst=out.lst[file.exists(out.lst)]
  all.out=paste(out_prefix,"_phyloP.feature.out",sep="")
  for (o in out.lst){
    cmd=paste("cat",o," | awk -F '\t' -v OFS='\t' '{if ($9~/^0/) print $0,\"con\"; else if ($9~/^-0/) print $0,\"acc\"}' >>",
              all.out,sep=" ")
    system(cmd,wait=TRUE)
  }
  
  all.out=read.table(all.out,sep="\t",header=FALSE,quote="")
  all.out[all.out$V10=="acc" & all.out$V9==0 ,"V9"]=-1e-10
  all.out[all.out$V10=="con" & all.out$V9==0 ,"V9"]=1e-10
  write.table(all.out,paste(out_prefix,"_phyloP.feature.out",sep=""),
              sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           "'{if ($9~/^0/ && $9<=0.01) print $0}'",
  #           all.out,
  #           ">",
  #           paste(out_prefix,"_conserved_p0.01.phyloP.feature.out",sep=""),
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  # 
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           "'{if ($9~/^-0/ && $9>=-0.01 ) print $0}'",
  #           all.out,
  #           ">",
  #           paste(out_prefix,"_accelerated_p0.01.phyloP.feature.out",sep=""),
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           "'{if ($9~/^0/ && $9<=0.05) print $0}'",
  #           all.out,
  #           ">",
  #           paste(out_prefix,"_conserved_p0.05.phyloP.feature.out",sep=""),
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  # 
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           "'{if ($9~/^-0/ && $9>=-0.05 ) print $0}'",
  #           all.out,
  #           ">",
  #           paste(out_prefix,"_accelerated_p0.05.phyloP.feature.out",sep=""),
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
}

# phyloP_features(maf_for_phyloP="/bucket/BourguignonU/Cong/termite_pca/phyloP/Cges_maf/", # dir for maf files, e.g. Cges.scaffold_501_mafDuplicateFilter.maf
#                 ref_contig2length="/flash/BourguignonU/Cong/termite_pca/phyloP/Cges_contig2length", # seqkit fx2tab -n -l <fna>
#                 features.dir="/flash/BourguignonU/Cong/termite_pca/phyloP/Cges_anno.splitBed", # dir for bed, e.g. Cges.scaffold_21.bed
#                 out_prefix="/flash/BourguignonU/Cong/termite_pca/phyloP/phyloP.features_Cges_secondary_wood_feeders",
#                 test_branches="Nluj,Hunk,N92,Cpar,Abea,Cwal,Munk", # comma lst
#                 test_subtree=NA,
#                 phyloFit.mod="/bucket/BourguignonU/Cong/termite_pca/phyloFit/models/N48_phyloFit.mod",
#                 threads=100)


phyloP_basewise=function(maf_for_phyloP=maf_for_phyloP, # dir for maf files
                            ref_contig2length=ref_contig2length, # seqkit fx2tab -n -l <fna>
                            out_prefix=out_prefix,
                            test_branches=test_branches, # comma lst
                            test_subtree=NA,
                            phyloFit.mod=phyloFit.mod,
                            threads=threads){
  if (!dir.exists(out_prefix)){dir.create(out_prefix)}
  
  contig.lst=read.table(ref_contig2length,header=FALSE,quote="",sep="\t")
  rownames(contig.lst)=contig.lst$V1
  contig.lst=contig.lst[order(-contig.lst$V2),]
  
  f=function(seq){
    stick=paste(out_prefix,"/",seq,".mafDuplicateFilter_phyloP.finished",sep="")
    if (!file.exists(stick)){
      cmd=paste("phyloP",
                "-i MAF --method LRT --mode CONACC --wig-scores",
                "--chrom",seq,
                sep=" ")
      if (!is.na(test_branches)){
        cmd=paste(cmd,"--branch",test_branches,sep=" ")
      }
      if (!is.na(test_subtree)){
        cmd=paste(cmd,"--subtree",test_subtree,sep=" ")
      }
                
      cmd=paste(cmd,
                phyloFit.mod,
                paste(maf_for_phyloP,"/",seq,"_mafDuplicateFilter.maf",sep=""),
                ">",
                paste(out_prefix,"/",seq,".mafDuplicateFilter_phyloP.wig",sep=""),
                sep=" ")
      system(cmd,wait=TRUE)
      cmd=paste("convert2bed -i wig",
                "<",paste(out_prefix,"/",seq,".mafDuplicateFilter_phyloP.wig",sep=""),">",
                paste(out_prefix,"/",seq,".mafDuplicateFilter_phyloP.bed",sep=""),
                sep=" ")
      system(cmd,wait=TRUE)
      cmd=paste("touch",stick,sep=" ")
      system(cmd,wait=TRUE)
    }
  }
  
  library(parallel)
  clus=makeCluster(threads)
  clusterExport(clus,"contig.lst",envir=environment())
  clusterExport(clus,"f",envir=environment())
  clusterExport(clus,"out_prefix",envir=environment())
  clusterExport(clus,"test_branches",envir=environment())
  clusterExport(clus,"phyloFit.mod",envir=environment())
  clusterExport(clus,"maf_for_phyloP",envir=environment())
  clusterExport(clus,"out_prefix",envir=environment())
  clusterExport(clus,"test_subtree",envir=environment())
  parSapply(clus,contig.lst[1:50,"V1"],f)
  parSapply(clus,contig.lst[51:100,"V1"],f)
  parSapply(clus,contig.lst[101:150,"V1"],f)
  parSapply(clus,contig.lst[151:200,"V1"],f)
  parSapply(clus,contig.lst[201:nrow(contig.lst),"V1"],f)
  
  wig.files=paste(out_prefix,"/",contig.lst$V1,".mafDuplicateFilter_phyloP.wig",sep="")
  wig.files=wig.files[file.exists(wig.files)]
  if (length(wig.files)!=0){sapply(wig.files,function(i){file.remove(i)})}
  
  bed=paste(out_prefix,"/",contig.lst$V1,".mafDuplicateFilter_phyloP.bed",sep="")
  all.bed=paste(out_prefix,".bed",sep="") # # seq start end ID phyloP-score
  if (file.exists(all.bed)){file.remove(all.bed)}
  for (b in bed){
    cmd=paste("cat",b,">>",all.bed,sep=" ")
    system(cmd,wait=TRUE)
  }
  
  # cmd=paste("sort -k5,5g",paste(out_prefix,"_p.tsv",sep=""),"-g --parallel",as.character(threads),
  #           "| awk -F '\t' -v OFS='\t' '{print $0, NR}' >",
  #           paste(out_prefix,"_p_sorted.tsv",sep=""),
  #           sep=" ") # seq start end ID phyloP-score text P P-rank
  # print(cmd);system(cmd,wait=TRUE)
  # file.remove(paste(out_prefix,"_p.tsv",sep=""))
  
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           "'{ if ($5 > 0) print $0, \"conserved\", 10^-$5; else if ($5 < 0) print $0, \"accelerated\", 10^$5; else print $0, \"neutral\", 10^$5 }'",
  #           all.bed,
  #           ">",paste(out_prefix,"_p.tsv",sep=""),
  #           sep=" ") # seq start end ID phyloP-score text P
  # print(cmd);system(cmd,wait=TRUE)
  # file.remove(all.bed)
  # 
  # cmd=paste("sort -k7,7g",paste(out_prefix,"_p.tsv",sep=""),"-g --parallel",as.character(threads),
  #           "| awk -F '\t' -v OFS='\t' '{print $0, NR}' >",
  #           paste(out_prefix,"_p_sorted.tsv",sep=""),
  #           sep=" ") # seq start end ID phyloP-score text P P-rank
  # print(cmd);system(cmd,wait=TRUE)
  # file.remove(paste(out_prefix,"_p.tsv",sep=""))
  # 
  # cmd=paste("tail -n1",paste(out_prefix,"_p_sorted.tsv",sep=""),
  #           " | awk -F '\t' -v OFS='\t' '{print $8}'",
  #           sep=" ")
  # total=system(cmd,intern=TRUE)
  # 
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           paste("-v total=",total,sep=""),
  #           "'{ print $0, ($8 / total)*0.05 }'",
  #           paste(out_prefix,"_p_sorted.tsv",sep=""),
  #           ">",paste(out_prefix,"_p_sorted_BH.tsv",sep=""),
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE) # seq start end ID phyloP-score text P P-rank BH
  # file.remove(paste(out_prefix,"_p_sorted.tsv",sep=""))
  # 
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           "'{if ($7 < $9) print $7}'",
  #           paste(out_prefix,"_p_sorted_BH.tsv",sep=""),
  #           "| tail -n 1 >",
  #           paste(out_prefix,"_pthreshold",sep=""),
  #           sep=" ")
  # print(cmd);system(cmd,intern=TRUE)
  # p.threshold=readLines(paste(out_prefix,"_pthreshold",sep=""))
  # 
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           paste("-v pthreshold=",p.threshold,sep=""),
  #           "'{if ($5 < 0 && $7<= pthreshold) print $0}'",
  #           paste(out_prefix,"_p_sorted_BH.tsv",sep=""),
  #           ">",
  #           paste(out_prefix,"_accelerated_BH.tsv",sep=""),
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  # 
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           paste("-v pthreshold=",p.threshold,sep=""),
  #           "'{if ($5 > 0 && $7<= pthreshold) print $0}'",
  #           paste(out_prefix,"_p_sorted_BH.tsv",sep=""),
  #           ">",
  #           paste(out_prefix,"_conserved_BH.tsv",sep=""),
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  # 
  # cmd=paste("awk -F '\t' -v OFS='\t'",
  #           paste("-v pthreshold=",p.threshold,sep=""),
  #           "'{if ($5 < 0 && $7<= pthreshold) print $0}'",
  #           paste(out_prefix,"_p_sorted_BH.tsv",sep=""),
  #           ">",
  #           paste(out_prefix,"_accelerated_BH.tsv",sep=""),
  #           sep=" ")
  # print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("awk -F '\t' -v OFS='\t'",
            "'{if ($5 >= 2) print $0}'",
            all.bed,
            ">",
            paste(out_prefix,"_conserved_p0.01.bed",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("awk -F '\t' -v OFS='\t'",
            "'{if ($5 <= -2 ) print $0}'",
            all.bed,
            ">",
            paste(out_prefix,"_accelerated_p0.01.bed",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# phyloP_wholeGenome(maf_for_phyloP="/bucket/BourguignonU/Cong/termite_pca/phyloP/Cges_maf/", # dir for maf files
#                    ref_contig2length="/flash/BourguignonU/Cong/termite_pca/phyloP/Cges_contig2length", # seqkit fx2tab -n -l <fna>
#                    out_prefix="/flash/BourguignonU/Cong/termite_pca/phyloP/phyloP_Cges_secondary_wood_feeders",
#                    test_branches="Nluj,Hunk,N92,Cpar,Abea,Cwal,Munk", # comma lst
#                    test_subtree=NA,
#                    phyloFit.mod="/bucket/BourguignonU/Cong/termite_pca/phyloFit/models/N48_phyloFit.mod",
#                    threads=100)
  
phyloPsite_2_element=function(integrated_annotation.bed=integrated_annotation.bed, # 
                              phyloP_p_sorted_BH.tsv=phyloP_p_sorted_BH.tsv, # seq start end ID phyloP-score text P P-rank BH
                              ref_contig2length=ref_contig2length, # seqkit fx2tab -n -l <fna>
                              out.tsv=out.tsv, # same with phyloP_accelerated_BH.tsv with the last column for genomic elements 
                              histogram.tsv=histogram.tsv, #
                              element_pct.tsv=element_pct.tsv,
                              contig_distr.tsv=contig_distr.tsv,
                              threads=threads){
  integrated_annotation=read.table(integrated_annotation.bed,sep="\t",header=FALSE,quote="")
  phyloP_p_sorted_BH=read.table(phyloP_p_sorted_BH.tsv,sep="\t",header=FALSE,quote="")
  
  f=function(i){
    chr=phyloP_p_sorted_BH[i,1]
    site=phyloP_p_sorted_BH[i,2]
    elements=integrated_annotation[integrated_annotation$V1==chr & 
                                     integrated_annotation$V2<=site & 
                                     integrated_annotation$V3>site,]
    if (nrow(elements)!=0){
      elements=elements[elements$V5==max(elements$V5),]
      return(paste(elements$V4,collapse=";;"))
    }else{
      return(NA)
    }
  }
  
  library(parallel)
  clus=makeCluster(threads)
  clusterExport(clus,"f",envir=environment())
  clusterExport(clus,"phyloP_p_sorted_BH",envir=environment())
  clusterExport(clus,"integrated_annotation",envir=environment())
  phyloP_p_sorted_BH$element=parSapply(clus,1:nrow(phyloP_p_sorted_BH),f)
  write.table(phyloP_p_sorted_BH,out.tsv,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
  
  bounds=10*c(min(phyloP_p_sorted_BH$V5),max(phyloP_p_sorted_BH$V5))
  histo=data.frame(phyloP_bin.start=seq(floor(bounds[1])/10,ceiling(bounds[2])/10,0.1),
                   phyloP_bin.end=seq(floor(bounds[1])/10,ceiling(bounds[2])/10,0.1)+0.1)
  clusterExport(clus,"phyloP_p_sorted_BH",envir=environment())
  histo[,c("no.anno","cds","utr5p","utr3p","intron","PoReEl_up2k","PoReEl_down2k",
           "tRNA","rRNA","miRNA","RF")]=t(
             parSapply(clus,histo$phyloP_bin.start,
                       function(i){
                         elements=phyloP_p_sorted_BH[phyloP_p_sorted_BH$V5>=i & 
                                                       phyloP_p_sorted_BH$V5<i+0.1,10]
                         no.anno=sum(is.na(elements))
                         elements=elements[!is.na(elements)]
                    
                         cds=sum(grepl("_cds[0-9]*",elements))
                         utr5p=sum(grepl("_utr5p[0-9]*",elements))
                         utr3p=sum(grepl("_utr3p[0-9]*",elements))
                         intron=sum(grepl("_intron[0-9]*",elements))
                         PoReEl_up2k=sum(grepl("_PoReEl_up2k",elements))
                         PoReEl_down2k=sum(grepl("_PoReEl_down2k",elements))
                         tRNA=sum(grepl("_tRNA[0-9]*",elements))
                         rRNA=sum(grepl("_rRNA[0-9]*",elements))
                         miRNA=sum(grepl("_miRNA[0-9]*",elements))
                         ncRNA=sum(grepl("_RF[0-9]*",elements))
                         return(c(no.anno,cds,utr5p,utr3p,intron,PoReEl_up2k,PoReEl_down2k,
                                  tRNA,rRNA,miRNA,ncRNA))
                       })
           )
  write.table(histo,histogram.tsv,sep="\t",row.names=FALSE,quote=FALSE)
  
  ref_contig2length=read.table(ref_contig2length,header=FALSE,quote="",sep="\t")
  colnames(ref_contig2length)=c("seqid","length")
  ref_contig2length[,basename(phyloP_p_sorted_BH.tsv)]=table(phyloP_p_sorted_BH$V1)[ref_contig2length$seqid]
  ref_contig2length[is.na(ref_contig2length[,basename(phyloP_p_sorted_BH.tsv)]),
                    basename(phyloP_p_sorted_BH.tsv)]=0
  write.table(ref_contig2length,contig_distr.tsv,row.names=FALSE,quote=FALSE,sep="\t")
  
  element_pct=data.frame(Element=c("cds","utr5p","utr3p","intron","PoReEl_up2k","PoReEl_down2k",
                                   "tRNA","rRNA","miRNA","RF"))
  element_pct$total.bp=sapply(element_pct$Element,
                              function(i){
                                p=paste("_",i,"[0-9]*",sep="")
                                if (i=="PoReEl_up2k" | i=="PoReEl_down2k"){p=paste(i,"$",sep="")}
                                d=integrated_annotation[grepl(i,integrated_annotation$V4),]
                                return(sum(d$V3-d$V2))
                              })
  element_pct[,basename(phyloP_p_sorted_BH.tsv)]=apply(histo,2,sum)[element_pct$Element]
  write.table(element_pct,element_pct.tsv,sep="\t",row.names=FALSE,quote=FALSE)
}


taffyCov=function(in.maf=in.maf,
                  out.tsv=out.tsv){
  cmd=paste("taffy coverage",
            "-i",in.maf,
            ">",out.tsv,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

halSynteny=function(hal=hal,
                    out.psl=out.psl,
                    minBlockSize=50000,
                    maxAnchorDistance=5000,
                    queryChromosome=NA,
                    queryGenome="Mnat",
                    targetGenome="Ofor",
                    minChr=5e+6){
  cmd=paste("halSynteny",
            "--minBlockSize",as.character(minBlockSize),
            "--maxAnchorDistance",as.character(maxAnchorDistance),
            sep=" ")
  if (!is.na(queryChromosome)){
    cmd=paste(cmd,
              "--queryChromosome",queryChromosome,
              sep=" ")
  }
  cmd=paste(cmd,
            "--queryGenome",queryGenome,
            "--targetGenome",targetGenome,
            hal,out.psl,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("awk -F '\t' -v OFS='\t'",
            "'{if ($11 > ",as.character(minChr)," && $15 > ",as.character(minChr),") print $9,$10,$11,$12,$13,$14,$15,$16,$17}'",
            out.psl,">",paste(out.psl,".tsv",sep=""),
            sep=" ")
  # 9. strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
  # 10. qName - Query sequence name.
  # 11. qSize - Query sequence size.
  # 12. qStart - Alignment start position in query.
  # 13. qEnd - Alignment end position in query.
  # 14. tName - Target sequence name.
  # 15. tSize - Target sequence size.
  # 16. tStart - Alignment start position in target.
  # 17. tEnd - Alignment end position in target.
  print(cmd);system(cmd,wait=TRUE)
  
}
  
hal2fasta=function(hal=hal,
                   out.fna=out.fna,
                   sp=sp){
  cmd=paste("hal2fasta",
            "--outFaPath",out.fna,
            hal,sp,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}


halLiftover=function(hal=hal,
                     ref.sp=ref.sp,
                     ref.bed=ref.bed,
                     query.sp=query.sp,
                     out.bed=out.bed){
  cmd=paste("halLiftover",
            hal,ref.sp,ref.bed,query.sp,out.bed,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Aban00020367-R0 Apac00012132-R0
hal2ortho=function(hal="/bucket/BourguignonU/Cong/termite_pca/cactus/termite_genomes.hal",
                   bed.dir="/flash/BourguignonU/Cong/termite_pca/hal2ortho/bed", # <bed.dir>/${sp}.bed, ${sp} indicated in name column
                   spp="Cmer,Mdar,Znev", # comma-list of ${sp}
                   out_dir="/flash/BourguignonU/Cong/termite_pca/hal2ortho",
                   threads=6){
  if (!file.exists(out_dir)){dir.create(out_dir)}
  spp.lst=unlist(strsplit(spp,","))
  sp2sp=expand.grid(spp.lst,spp.lst)
  sp2sp=sp2sp[sp2sp$Var1!=sp2sp$Var2,]
  sp2sp$pair=sapply(1:nrow(sp2sp),function(i){return(paste(sort(c(sp2sp[i,1],sp2sp[i,2])),collapse="_"))})
  
  if (!file.exists(paste(out_dir,"/halLiftover",sep=""))){dir.create(paste(out_dir,"/halLiftover",sep=""))}
  library(parallel)
  clus=makeCluster(threads)
  clusterExport(clus,"sp2sp",envir=environment())
  clusterExport(clus,"bed.dir",envir=environment())
  clusterExport(clus,"out_dir",envir=environment())
  clusterExport(clus,"hal",envir=environment())
  stamp=paste(out_dir,"/halLiftover.finished",sep="")
  if (!file.exists(stamp)){
  parSapply(clus,1:nrow(sp2sp),
            function(i){
              reference=sp2sp[i,1]
              reference.bed=paste(bed.dir,"/",reference,".bed",sep="")
              target=sp2sp[i,2]
              target.bed=paste(bed.dir,"/",target,".bed",sep="")
              r.bed=read.table(reference.bed,header=FALSE,sep="\t",quote="")
              rownames(r.bed)=r.bed$V4
              t.bed=read.table(target.bed,header=FALSE,sep="\t",quote="")
              rownames(t.bed)=t.bed$V4
              
              if (!file.exists(paste(out_dir,"/halLiftover/",reference,"2",target,".tsv",sep=""))){
                cmd=paste("halLiftover --noDupes --outPSLWithName",
                          hal,reference,reference.bed,target,
                          paste(out_dir,"/halLiftover/",reference,"2",target,".bed",sep=""),
                          sep=" ")
                system(cmd,wait=TRUE)
                
                cmd=paste("awk -F '\t' -v OFS='\t' '{print $1,$13,$14,$10,$15,$17,$18}'",
                          paste(out_dir,"/halLiftover/",reference,"2",target,".bed",sep=""),
                          ">",
                          paste(out_dir,"/halLiftover/",reference,"2",target,".tmp",sep=""),
                          sep=" ")
                system(cmd,wait=TRUE)
                
                df=read.table(paste(out_dir,"/halLiftover/",reference,"2",target,".tmp",sep=""),
                              sep="\t",header=FALSE,quote="")
                colnames(df)=c("r.element","r.start","r.end","t.strand","t.chr","t.start","t.end")
                df$t.strand=sub("^.","",df$t.strand)
                
                cmd=paste("rm",paste(out_dir,"/halLiftover/",reference,"2",target,".bed",sep=""),sep=" ")
                system(cmd,wait=TRUE)
                cmd=paste("rm",paste(out_dir,"/halLiftover/",reference,"2",target,".tmp",sep=""),sep=" ")
                system(cmd,wait=TRUE)
                
                df[,c("t.element","overlap.bp")]=t(sapply(1:nrow(df),
                                                     function(i){
                                                       overlap=t.bed[t.bed$V2<df[i,"t.end"] & 
                                                                       t.bed$V3>df[i,"t.start"] &
                                                                       t.bed$V6==df[i,"t.strand"] &
                                                                       t.bed$V1==df[i,"t.chr"],]
                                                       t.element=NA;overlap.bp=NA
                                                       if (nrow(overlap)>0){
                                                         for (j in 1:nrow(overlap)){
                                                           tl=max(overlap[j,3],df[i,"t.end"])-min(overlap[j,2],df[i,"t.start"])
                                                           overlap[j,"ol"]=overlap[j,3]-overlap[j,2]+df[i,"t.end"]-df[i,"t.start"]-tl
                                                           #ol=overlap[j,3]-overlap[j,2]+df[i,"t.end"]-df[i,"t.start"]-tl
                                                           #overlap[j,"score"]=ol*2/(overlap[j,3]-overlap[j,2]+df[i,"t.end"]-df[i,"t.start"])
                                                         }
                                                         overlap=overlap[overlap$ol==max(overlap$ol),]
                                                         t.element=overlap[1,4];overlap.bp=overlap[1,"ol"]
                                                       }
                                                       return(c(t.element,overlap.bp))
                                                     }))
                df[is.na(df$overlap.bp),"overlap.bp"]=0
                df$overlap.bp=as.numeric(df$overlap.bp)
                write.table(df,paste(out_dir,"/halLiftover/",reference,"2",target,".tsv",sep=""),
                            sep="\t",row.names=FALSE,quote=FALSE)
              }else{
                  df=read.table(paste(out_dir,"/halLiftover/",reference,"2",target,".tsv",sep=""),
                                sep="\t",header=TRUE,quote="")
              }
              
              df=df[!is.na(df$t.element),]
              if (!file.exists(paste(out_dir,"/halLiftover/",reference,"2",target,"_pairs.tsv",sep=""))){
                pairs=df[!is.na(df$t.element),c("r.element","t.element")]
                pairs=pairs[!duplicated(pairs),]
                pairs$r.len=r.bed[pairs$r.element,3]-r.bed[pairs$r.element,2]
                pairs$t.len=t.bed[pairs$t.element,3]-t.bed[pairs$t.element,2]
                pairs$overlap.len=sapply(1:nrow(pairs),
                                         function(i){
                                           dd=df[df$r.element==pairs[i,"r.element"] & df$t.element==pairs[i,"t.element"],
                                                 c("t.start","t.end")]
                                           ranges=reshape2::melt(dd)
                                           ranges=ranges[order(ranges$value),]
                                           indices=intersect(which(ranges$variable=='t.start')-1, which(ranges$variable=='t.end'))
                                           dd=ranges[c(1, sort(c(indices, indices+1)), nrow(ranges)),]
                                           starts=dd[dd$variable=="t.start","value"]
                                           ends=dd[dd$variable=="t.end","value"]
                                           
                                           t.element.start=t.bed[pairs[i,"t.element"],2]
                                           t.element.end=t.bed[pairs[i,"t.element"],3]
                                           
                                           starts[starts<t.element.start]=t.element.start
                                           ends[ends>t.element.end]=t.element.end
                                           
                                           return(sum(ends-starts))
                                         })
                write.table(pairs,
                            paste(out_dir,"/halLiftover/",reference,"2",target,"_pairs.tsv",sep=""),
                            sep="\t",row.names=FALSE,quote=FALSE)
              }else{
                pairs=read.table(paste(out_dir,"/halLiftover/",reference,"2",target,"_pairs.tsv",sep=""),
                                 header=TRUE,sep="\t",quote="")
              }
              write.table(pairs[pairs$overlap.len/pairs$t.len>0.7,c("r.element","t.element")],
                          paste(out_dir,"/halLiftover/",reference,"2",target,"_0.7.pairs",sep=""),
                          sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
            })
  system(paste("touch",stamp,sep=" "))
  }
    
  if (!file.exists(paste(out_dir,"/bidirectional_pairs",sep=""))){dir.create(paste(out_dir,"/bidirectional_pairs",sep=""))}
  pairs=sp2sp$pair;pairs=pairs[!duplicated(pairs)]
  clusterExport(clus,"out_dir",envir=environment())
  stamp=paste(out_dir,"/bidirectional_pairs.finished",sep="")
  if (!file.exists(stamp)){
  parSapply(clus,pairs,
            function(pair){
              pair=unlist(strsplit(pair,"_"))
              if (!file.exists(paste(out_dir,"/bidirectional_pairs/",pair[1],"_",pair[2],".bidirectional.hits",sep=""))){
                p1=paste(out_dir,"/halLiftover/",pair[1],"2",pair[2],"_0.7.pairs",sep="")
                p1=read.table(p1,header=FALSE,sep="\t",quote="")
                p1$pair=paste(p1$V1,p1$V2,sep="_")
                p2=paste(out_dir,"/halLiftover/",pair[2],"2",pair[1],"_0.7.pairs",sep="")
                p2=read.table(p2,header=FALSE,sep="\t",quote="")
                p2$pair=paste(p2$V2,p2$V1,sep="_")
                
                bidirection=intersect(p1$pair,p2$pair)
                valid=p1[p1$pair %in% bidirection,c(1,2)]
                write.table(valid,
                            paste(out_dir,"/bidirectional_pairs/",pair[1],"_",pair[2],".bidirectional.hits",sep=""),
                            sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
              }
              
            })
  system(paste("touch",stamp,sep=" "))
  }
  
  if (!file.exists(paste(out_dir,"/OGs",sep=""))){dir.create(paste(out_dir,"/OGs",sep=""))}
  edges=paste(out_dir,"/OGs/edges.tsv",sep="")
  if (file.exists(edges)){file.remove(edges)}
  for (pair in pairs){
    pair=unlist(strsplit(pair,"_"))
    bidirectional=paste(out_dir,"/bidirectional_pairs/",pair[1],"_",pair[2],".bidirectional.hits",sep="")
    cmd=paste("cat",bidirectional,">>",edges,sep=" ")
    system(cmd,wait=TRUE)
  }
  stamp=paste(out_dir,"/OGs.finished",sep="")
  if (!file.exists(stamp)){
    df=read.table(edges,header=FALSE,sep="\t",quote="")
    graph=igraph::graph_from_edgelist(as.matrix(df),directed = FALSE)
    graph.clusters=igraph::clusters(graph)
    element2OG=data.frame(element=names(graph.clusters$membership),
                          OG=unname(graph.clusters$membership))
    element2OG$OG=paste("OG",sprintf("%08d",element2OG$OG),sep="")
    element2OG$sp=NA
    for (sp in spp.lst){
      element2OG[grepl(sp,element2OG$element),"sp"]=sp
    }
    write.table(element2OG,paste(out_dir,"/OGs/element2OG.tsv",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
    OG.lst=element2OG$OG;OG.lst=OG.lst[!duplicated(OG.lst)]
    
    OG2sp=paste(element2OG$OG,element2OG$sp,sep=",")
    OG2sp=table(OG2sp)
    OGxSp4copy.longlst=data.frame(OG2sp_=names(OG2sp),
                                  value=unname(OG2sp))
    OGxSp4copy.longlst$OG=sub(",.*$","",OGxSp4copy.longlst$OG2sp_)
    OGxSp4copy.longlst$sp=sub("^.*,","",OGxSp4copy.longlst$OG2sp_)
    OGxSp4copy=reshape2::acast(OGxSp4copy.longlst,OG~sp,value.var = "value.Freq")
    OGxSp4copy[is.na(OGxSp4copy)]=0
    OGxSp4copy=cbind(data.frame(OG=rownames(OGxSp4copy)),OGxSp4copy)
    write.table(OGxSp4copy,paste(out_dir,"/OGs/OGxSp4copy.tsv",sep=""),
                sep="\t",row.names=FALSE,quote=FALSE)
    
    clusterExport(clus,"element2OG",envir=environment())
    clusterExport(clus,"out_dir",envir=environment())
    parSapply(clus,OG.lst,
           function(og){
             writeLines(element2OG[element2OG$OG==og,"element"],
                        paste(out_dir,"/OGs/",og,".lst",sep=""))
           })
    system(paste("touch",stamp,sep=" "))
  }
  #igraph::graph_from_edgelist()
  #g <- graph_from_edgelist(matrix(c("foo", "bar", "bar", "foobar"), nc = 2, byrow = TRUE), directed = FALSE)
  #clusters <- clusters(g)
  
  # if (!file.exists(paste(out_dir,"/OGs",sep=""))){dir.create(paste(out_dir,"/OGs",sep=""))}
  # for (sp in spp.lst){
  #   bed=paste(bed.dir,"/",sp,".bed",sep="")
  #   element=paste("awk -F '\t' -v OFS='\t' '{print$4}' ",bed,sep="")
  #   element=system(element,intern=TRUE)
  #   if (file.exists(paste(out_dir,"/OGs/*.lst",sep=""))){
  #     assigned=system(paste("cat ",out_dir,"/OGs/*.lst",sep=""),intern=TRUE)
  #     element=element[!element %in% assigned]
  #   }
  #   clusterExport(clus,"sp",envir=environment())
  #   clusterExport(clus,"out_dir",envir=environment())
  #   clusterExport(clus,"spp.lst",envir=environment())
  #   parSapply(clus,element,
  #             function(seed.element){
  #               seed.sp=sp
  #               hal2ortho.bidirectional_pairs.dir=paste(out_dir,"/bidirectional_pairs",sep="")
  #               out.lst=paste(out_dir,"/OGs/",seed.element,".lst",sep="")
  #               
  #               df=data.frame(sp=c(spp.lst,seed.sp),element=NA)
  #               df=df[!duplicated(df$sp),]
  #               rownames(df)=df$sp
  #               df$element=sapply(df$sp,
  #                                 function(i){
  #                                   res=NA
  #                                   if (i==seed.sp){res=seed.element}else{
  #                                     edges=paste(hal2ortho.bidirectional_pairs.dir,"/",
  #                                                 seed.sp,"_",i,".bidirectional.hits",sep="")
  #                                     if (!file.exists(edges)){
  #                                       edges=paste(hal2ortho.bidirectional_pairs.dir,"/",
  #                                                   i,"_",seed.sp,".bidirectional.hits",sep="")
  #                                     }
  #                                     cmd=paste("awk -F '\t' -v OFS='\t' ",
  #                                               "'{if ($1==\"",seed.element,"\") print $2; 
  #                                 if ($2==\"",seed.element,"\") print $1}' ",
  #                                               edges,sep="")
  #                                     res=system(cmd,wait=TRUE,intern=TRUE)
  #                                     res=paste(res,collapse = ";")
  #                                   }
  #                                   return(res)
  #                                 })
  #               if (nrow(df[df$element!="",])==1){
  #                 print("No orthologue found")
  #               }else{
  #                 bool=!(nrow(df[df$element=="",])==0)
  #                 while ( bool ){
  #                   new.seed.df=df[df$element!="",]
  #                   new.df=df[df$element=="",]
  #                   for (i in 1:nrow(new.seed.df)){
  #                     new.df[,new.seed.df[i,"element"]]=sapply(new.df$sp,
  #                                                              function(sp){
  #                                                                edges=paste(hal2ortho.bidirectional_pairs.dir,"/",
  #                                                                            new.seed.df[i,"sp"],"_",sp,".bidirectional.hits",sep="")
  #                                                                if (!file.exists(edges)){
  #                                                                  edges=paste(hal2ortho.bidirectional_pairs.dir,"/",
  #                                                                              sp,"_",new.seed.df[i,"sp"],".bidirectional.hits",sep="")
  #                                                                }
  #                                                                cmd=paste("awk -F '\t' -v OFS='\t' ",
  #                                                                          "'{if ($1==\"",new.seed.df[i,"element"],"\") print $2; 
  #                                 if ($2==\"",new.seed.df[i,"element"],"\") print $1}' ",
  #                                                                          edges,sep="")
  #                                                                res=system(cmd,wait=TRUE,intern=TRUE)
  #                                                                res=paste(res,collapse = ";")
  #                                                                
  #                                                                return(res)
  #                                                              })
  #                   }
  #                   df[rownames(new.df),"element"]=apply(new.df[,3:ncol(new.df)],1,function(j){return(paste(j[j!="" & !duplicated(j)],collapse = ";"))})
  #                   bool=!( (nrow(df[df$element=="",])==0) | (nrow(df[df$element=="",])==nrow(new.df)) )
  #                 }
  #                 lst=paste(df$element,collapse=";")
  #                 lst=unlist(strsplit(lst,";"))
  #                 lst=lst[lst!=""]
  #                 writeLines(lst,out.lst)
  #               }
  #             })
  # }
  
  
}

seed2ortho=function(seed.element="Mdar00000233-R0_intron1721",
                    seed.sp="Mdar",
                    spp.lst="Aaca,Aban,Abea,Apac,Aunk,Bori,Cbre,Ccav,Cges,Cmer,Cpar,Csp4,Ctes,Cwal,Dlon,Eunk,Fval,Gfus,Gocu,Hsjo,Hten,Hunk,Isch,Iunk,Kfla,Llab,Lunk,Mdar,Mhub,Mnat,Munk,Ncas,Nluj,Ntar,Ofor,Pada,PAsim,Pred,PRsim,Punk,Rebo,Rfla,Shal,Shey,Ssph,Svic,Znev", # comma-list
                    hal2ortho.bidirectional_pairs.dir="/bucket/BourguignonU/Cong/termite_pca/hal2ortho/bidirectional_pairs/",
                    out.lst=out.lst){
  spp.lst=unlist(strsplit(spp.lst,","))
  df=data.frame(sp=c(spp.lst,seed.sp),element=NA)
  df=df[!duplicated(df$sp),]
  rownames(df)=df$sp
  df$element=sapply(df$sp,
                    function(i){
                      res=NA
                      if (i==seed.sp){res=seed.element}else{
                        edges=paste(hal2ortho.bidirectional_pairs.dir,"/",
                                    seed.sp,"_",i,".bidirectional.hits",sep="")
                        if (!file.exists(edges)){
                          edges=paste(hal2ortho.bidirectional_pairs.dir,"/",
                                      i,"_",seed.sp,".bidirectional.hits",sep="")
                        }
                        cmd=paste("awk -F '\t' -v OFS='\t' ",
                                  "'{if ($1==\"",seed.element,"\") print $2; 
                                  if ($2==\"",seed.element,"\") print $1}' ",
                                  edges,sep="")
                        res=system(cmd,wait=TRUE,intern=TRUE)
                        res=paste(res,collapse = ";")
                      }
                      return(res)
                    })
  if (nrow(df[df$element!="",])==1){
    print("No orthologue found")
  }else{
  bool=!(nrow(df[df$element=="",])==0)
  while ( bool ){
    new.seed.df=df[df$element!="",]
    new.df=df[df$element=="",]
    for (i in 1:nrow(new.seed.df)){
      new.df[,new.seed.df[i,"element"]]=sapply(new.df$sp,
                                               function(sp){
                                                   edges=paste(hal2ortho.bidirectional_pairs.dir,"/",
                                                               new.seed.df[i,"sp"],"_",sp,".bidirectional.hits",sep="")
                                                   if (!file.exists(edges)){
                                                     edges=paste(hal2ortho.bidirectional_pairs.dir,"/",
                                                                 sp,"_",new.seed.df[i,"sp"],".bidirectional.hits",sep="")
                                                   }
                                                   cmd=paste("awk -F '\t' -v OFS='\t' ",
                                                             "'{if ($1==\"",new.seed.df[i,"element"],"\") print $2; 
                                  if ($2==\"",new.seed.df[i,"element"],"\") print $1}' ",
                                                             edges,sep="")
                                                   res=system(cmd,wait=TRUE,intern=TRUE)
                                                   res=paste(res,collapse = ";")
                                                 
                                                 return(res)
                                               })
    }
    df[rownames(new.df),"element"]=apply(new.df[,3:ncol(new.df)],1,function(j){return(paste(j[j!="" & !duplicated(j)],collapse = ";"))})
    bool=!( (nrow(df[df$element=="",])==0) | (nrow(df[df$element=="",])==nrow(new.df)) )
  }
  lst=paste(df$element,collapse=";")
  lst=unlist(strsplit(lst,";"))
  lst=lst[lst!=""]
  writeLines(lst,out.lst)
  }
}
# seed="Mdar00000233-R0_intron1721"
# 
# awk '{if ($1=="") print $2; else if ($2=="Mdar00000233-R0_intron1721") print $1}' ../bidirectional_pairs/* | sort | uniq > round.1.lst
# echo 'Mdar00000233-R0_intron1721' >> round.1.lst
# 
# grep -f round.1.lst ../bidirectional_pairs/* | sed 's/^.*://' | sort | uniq > edges.lst

# genome=read.table("~/quality_genomePeptide.tsv",header=TRUE,sep="\t",quote="")
# hal2ortho(hal="/bucket/BourguignonU/Cong/termite_pca/cactus/termite_genomes.hal",
#                    bed.dir="/flash/BourguignonU/Cong/termite_pca/hal2ortho/bed", # <bed.dir>/${sp}.bed
#                    spp=paste(genome$Label,collapse=","), # comma-list of ${sp}
#                    out_dir="/flash/BourguignonU/Cong/termite_pca/hal2ortho",
#                    threads=100)
  
FastGA=function(genome1.fna=genome1.fna,sp1=sp1,
                genome2.fna=genome2.fna,sp2=sp2,
                tmp_dir=tmp_dir,
                out.paf=out.paf,
                threads=threads){
  if (!dir.exists(tmp_dir)){dir.create(tmp_dir)}
  system(paste("cp",genome1.fna,tmp_dir,sep=" "))
  genome1.fna=paste(tmp_dir,"/",basename(genome1.fna),sep="")
  cmd=paste("sed -i 's/>/>",sp1,"_/' ",genome1.fna,sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("cp",genome2.fna,tmp_dir,sep=" "))
  genome2.fna=paste(tmp_dir,"/",basename(genome2.fna),sep="")
  cmd=paste("sed -i 's/>/>",sp2,"_/' ",genome2.fna,sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("FastGA",
            paste("-T",as.character(threads),sep=""),
            paste("-P",tmp_dir,sep=""),
            sep=" ")
  if (genome1.fna==genome2.fna){
    cmd=paste(cmd,genome1.fna,sep=" ")
  }else{
    cmd=paste(cmd,genome1.fna,genome2.fna,sep=" ")
  }
  cmd=paste(cmd,">",out.paf,sep=" ")
  print(cmd)
  system(cmd,wait=TRUE)
}

ALNplot=function(in.paf=in.paf,
                 threads=threads,
                 out.pdf=out.pdf){
  cmd=paste("ALNplot",
            " -T",as.character(threads),
            " -p",out.pdf,
            " ",in.paf,
            sep="")
  print(cmd);system(cmd,wait=TRUE)
}

make_lastz_chains=function(target.genome.fna=target.genome.fna,
                           sp.target=sp.target, # Dmel, reference,
                           query.genome.fna=query.genome.fna,
                           sp.query=sp.query,
                           tmp_dir=tmp_dir,
                           threads=threads){
  if (!dir.exists(tmp_dir)){dir.create(tmp_dir)}
  wd=getwd()
  setwd(tmp_dir)
  
  # system(paste("cp",target.genome.fna,tmp_dir,sep=" "))
  # target.genome.fna=paste(tmp_dir,"/",basename(target.genome.fna),sep="")
  # cmd=paste("sed -i 's/>/>",sp.target,"_/' ",target.genome.fna,sep="")
  # print(cmd);system(cmd,wait=TRUE)
  
  # system(paste("cp",query.genome.fna,tmp_dir,sep=" "))
  # query.genome.fna=paste(tmp_dir,"/",basename(query.genome.fna),sep="")
  # cmd=paste("sed -i 's/>/>",sp.query,"_/' ",query.genome.fna,sep="")
  # print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("make_chains.py",
            sp.target,sp.query,target.genome.fna,query.genome.fna,
            "--project_dir",paste(tmp_dir,"/make_lastz_chains.projectDir",sep=""),
            "--executor local",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
}

toga=function(target.genome.2bit=target.genome.2bit, # Dmel, reference,
              target.genome.anno.bed12=target.genome.anno.bed12,
              target.gene2transcriptIsoforms.tsv=target.gene2transcriptIsoforms.tsv,
              query.genome.2bit=query.genome.2bit,
              chain=chain,
              tmp_dir=tmp_dir,
              TOGA_dir=TOGA_dir,
              threads=threads){
  wd=getwd()
  setwd(TOGA_dir)
  
  cmd=paste("./toga.py",
            chain,
            target.genome.anno.bed12,
            target.genome.2bit,
            query.genome.2bit,
            "--isoforms",target.gene2transcriptIsoforms.tsv,
            "--project_dir",tmp_dir,
            "--cesar_jobs_num",as.character(threads),
            "--cesar_mem_limit 10",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(wd)
}

nf_LO=function(ref.genome.fna=ref.genome.fna,
               ref.genome.gff=ref.genome.gff,
               target.genome.fna=target.genome.fna,
               out_dir=out_dir,
               threads=64,
               memory=256,
               chain_name="Aaca_Dmel",
               nf_lo.dir="/bucket/BourguignonU/Cong/Softwares/nf-LO"){
  if (!dir.exists(out_dir)){dir.create(out_dir)}
  
  cmd=paste("nextflow run",nf_lo.dir,
            "--source",ref.genome.fna,
            "--target",target.genome.fna,
            "--annotation",ref.genome.gff,
            "--annotation_format gff",
            "--distance far",
            "--aligner lastz",
            "--liftover_algorithm liftover",
            "--outdir",out_dir,
            "--max_cpus",as.character(threads),
            "--max_memory",paste(as.character(memory),".GB",sep=""),
            "--chain_name",chain_name,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  # ${params.outdir}/chainnet/liftover.chain: this is the final chain file that can be used as input for liftover, crossmap and other tools to perform the actual liftover.
  # 
  # ${params.outdir}/chainnet/netfile.net: net file associated with the chain file generate to perform the liftover.
  # 
  # ${params.outdir}/stats/mafCoverage.out: number of bases in the two genomes covered by the chain file, calculated using mafTools (when available)
  # 
  # ${params.outdir}/stats/mafIdentity.out: number of identical bases in the two genomes covered by the chain file, calculated using mafTools (when available)
  # 
  # ${params.outdir}/stats/mafStats.out: generic metrics of the chain file calculated using mafTools (when available)
  # 
  # ${params.outdir}/stats/features.txt: number of features lifted (when a feature set is provided)
  # 
  # ${params.outdir}/lifted/${params.chain_name}.*: set of lifted features using either crossmap or liftover
  
}


#integrated_annotation
# genome=read.table("~/quality_genomePeptide.tsv",
#                   header=TRUE,sep="\t",quote="")
# i=as.numeric( commandArgs(trailingOnly=TRUE) )
# for (sp in genome[i,"Label"]){
# dire="/flash/BourguignonU/Cong/termite_genome_annotation/integrated_annotation/"
# genes.gff3=paste("/bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/",
#                  sp,"/",sp,"_genes.gff3",sep="")
# iso.transcript=paste("/bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/",
#                      sp,"/",sp,"_proteins_iso.faa",sep="")
# iso.transcript.lst=paste(dire,"/",sp,"_iso.lst",sep="")
# cmd=paste("grep '>' ",iso.transcript," | sed 's/>//' > ",iso.transcript.lst,sep="")
# system(cmd)
# rep.transcript.gff3=paste(dire,"/",sp,"_rep.transcript.gff3",sep="")
# cmd=paste("agat_sp_filter_feature_from_kill_list.pl","--gff",genes.gff3,"--kill_list",iso.transcript.lst,
#           "--output",rep.transcript.gff3,sep=" ")
# print(cmd);system(cmd,wait=TRUE)
# rep.transcript_intron.gff3=paste(dire,"/",sp,"_rep.transcript_intron.gff3",sep="")
# cmd=paste("agat_sp_add_introns.pl","--gff",rep.transcript.gff3,"--out",rep.transcript_intron.gff3,sep=" ")
# print(cmd);system(cmd,wait=TRUE)
# genes.gff=ape::read.gff(rep.transcript_intron.gff3)
# genes.bed=data.frame(chrom=genes.gff$seqid,
#                      chromStart=format(genes.gff$start-1,scientific=FALSE,trim=TRUE),
#                      chromEnd=format(genes.gff$end,scientific=FALSE,trim=TRUE))
# genes.bed$name=sapply(1:nrow(genes.bed),
#                       function(i){
#                         type=genes.gff[i,"type"]
#                         attributes=genes.gff[i,"attributes"]
#                         attributes=unlist(strsplit(attributes,";"))
#                         if (type=="gene"){
#                           res=sub("ID=","",attributes[1])
#                         }
#                         if (type=="mRNA"){
#                           res=sub("ID=","",attributes[1])
#                         }
#                         if (type=="exon"){
#                           ID=sub("ID=","",attributes[1]);ID=stringr::str_extract(ID,"exon.*$")
#                           parent=sub("Parent=","",attributes[2])
#                           res=paste(parent,"_",ID,sep="")
#                         }
#                         if (type=="CDS"){
#                           ID=sub("ID=","",attributes[1]);ID=stringr::str_extract(ID,"cds.*$")
#                           ID=sub("\\.","",ID)
#                           parent=sub("Parent=","",attributes[2])
#                           res=paste(parent,"_",ID,sep="")
#                         }
#                         if (type=="intron"){
#                           ID=sub("ID=","",attributes[1]);ID=sub("intron_added-","",ID);ID=paste("intron",ID,sep="")
#                           parent=sub("Parent=","",attributes[2])
#                           res=paste(parent,"_",ID,sep="")
#                         }
#                         if (type=="five_prime_UTR"){
#                           ID=sub("ID=","",attributes[1]);ID=stringr::str_extract(ID,"utr5p.*$")
#                           parent=sub("Parent=","",attributes[2])
#                           if (!is.na(ID)){
#                             res=paste(parent,"_",ID,sep="")
#                           }else{
#                             res=NA
#                           }
# 
#                         }
#                         if (type=="three_prime_UTR"){
#                           ID=sub("ID=","",attributes[1]);ID=stringr::str_extract(ID,"utr3p.*$")
#                           parent=sub("Parent=","",attributes[2])
#                           if (!is.na(ID)){
#                             res=paste(parent,"_",ID,sep="")
#                           }else{
#                             res=NA
#                           }
#                         }
#                         return(res)
#                       })
# genes.bed$score=sapply(1:nrow(genes.bed),
#                        function(i){
#                          type=genes.gff[i,"type"]
#                          if (type=="gene"){
#                            res=1
#                          }
#                          if (type=="mRNA"){
#                            res=2
#                          }
#                          if (type=="exon"){
#                            res=3
#                          }
#                          if (type=="CDS"){
#                            res=4
#                          }
#                          if (type=="intron"){
#                            res=3
#                          }
#                          if (type=="five_prime_UTR"){
#                            res=4
#                          }
#                          if (type=="three_prime_UTR"){
#                            res=4
#                          }
#                          return(res)
#                        })
# genes.bed$strand=genes.gff$strand
# genes.bed=genes.bed[!is.na(genes.bed$name),]
# write.table(genes.bed,paste(dire,"/",sp,"_mRNA.bed",sep=""),
#             sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
# 
# chr2length=paste("/bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/",
#                  sp,"/",sp,"_chrLength.tsv",sep="")
# chr2length=read.table(chr2length,header=FALSE,quote="",sep="\t")
# rownames(chr2length)=chr2length$V1
# df=genes.gff[genes.gff$type=="gene",]
# df[,c("PoReEl_up2k.start","PoReEl_up2k.end",
#       "PoReEl_down2k.start","PoReEl_down2k.end")]=t(sapply(1:nrow(df),
#                             function(i){
#                               g.s=df[i,"start"];g.e=df[i,"end"]
#                               if (df[i,"strand"]=="+"){
#                                 PoReEl_up2k.start=g.s-2001
#                                 PoReEl_up2k.end=g.s-1
#                                 PoReEl_down2k.end=g.e+2001
#                                 PoReEl_down2k.start=g.e+1
#                               }
#                               if (df[i,"strand"]=="-"){
#                                 PoReEl_up2k.start=g.e+1
#                                 PoReEl_up2k.end=g.e+2001
#                                 PoReEl_down2k.end=g.s-1
#                                 PoReEl_down2k.start=g.s-2001
#                               }
#                               return(c(PoReEl_up2k.start,PoReEl_up2k.end,
#                                        PoReEl_down2k.start,PoReEl_down2k.end))
#                             }))
# df$chrLen=chr2length[as.character(df$seqid),2]
# df[df$PoReEl_up2k.start<1,"PoReEl_up2k.start"]=1
# df[df$PoReEl_up2k.start>df$chrLen,"PoReEl_up2k.start"]=df[df$PoReEl_up2k.start>df$chrLen,"chrLen"]
# df[df$PoReEl_up2k.end<1,"PoReEl_up2k.end"]=1
# df[df$PoReEl_up2k.end>df$chrLen,"PoReEl_up2k.end"]=df[df$PoReEl_up2k.end>df$chrLen,"chrLen"]
# df[df$PoReEl_down2k.start<1,"PoReEl_down2k.start"]=1
# df[df$PoReEl_down2k.start>df$chrLen,"PoReEl_down2k.start"]=df[df$PoReEl_down2k.start>df$chrLen,"chrLen"]
# df[df$PoReEl_down2k.end<1,"PoReEl_down2k.end"]=1
# df[df$PoReEl_down2k.end>df$chrLen,"PoReEl_down2k.end"]=df[df$PoReEl_down2k.end>df$chrLen,"chrLen"]
# df$geneID=sapply(df$attributes,
#                  function(i){
#                    i=unlist(strsplit(i,";"))[1]
#                    return(sub("ID=","",i))
#                  })
# PoReEl_up.bed=data.frame(chrom=df$seqid,
#                       chromStart=format(df$PoReEl_up2k.start-1,scientific=FALSE,trim=TRUE),
#                       chromEnd=format(df$PoReEl_up2k.end,scientific=FALSE,trim=TRUE),
#                       name=paste(df$geneID,".PoReEl_up2k",sep=""),
#                       score=2,strand=df$strand)
# write.table(PoReEl_up.bed,paste(dire,"/",sp,"_PoReEl_up2k.bed",sep=""),
#             sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
# PoReEl_down.bed=data.frame(chrom=df$seqid,
#                          chromStart=format(df$PoReEl_down2k.start-1,scientific=FALSE,trim=TRUE),
#                          chromEnd=format(df$PoReEl_down2k.end,scientific=FALSE,trim=TRUE),
#                          name=paste(df$geneID,".PoReEl_down2k",sep=""),
#                          score=2,strand=df$strand)
# write.table(PoReEl_down.bed,paste(dire,"/",sp,"_PoReEl_down2k.bed",sep=""),
#             sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
# 
# tRNA.gff3=paste("/bucket/BourguignonU/Cong/termite_genome_annotation/tRNAscan_filter/",
#                 sp,"/tRNA.gff3",sep="")
# tRNA.gff3=ape::read.gff(tRNA.gff3)
# tRNA.gff3=tRNA.gff3[tRNA.gff3$type=="tRNA",]
# tRNA.gff3[,c("AA","anticodon")]=t(sapply(tRNA.gff3$attributes,
#                                          function(i){
#                                            i=unlist(strsplit(i,";"))
#                                            AA=sub("ID=tRNA-","",i[1])
#                                            anticodon=sub("anticodon=","",i[3])
#                                            return(c(AA,anticodon))
#                                          }))
# tRNA.bed=data.frame(chrom=tRNA.gff3$seqid,
#                     chromStart=format(tRNA.gff3$start-1,scientific=FALSE,trim=TRUE),
#                     chromEnd=format(tRNA.gff3$end,scientific=FALSE,trim=TRUE),
#                     name=paste(sp,"_",sub("[0-9]*_tRNA","",tRNA.gff3$AA),"_anticodon_",tRNA.gff3$anticodon,"_tRNA",
#                                sprintf("%08d",1:nrow(tRNA.gff3)),sep=""),
#                     score=1,strand=tRNA.gff3$strand)
# write.table(tRNA.bed,paste(dire,"/",sp,"_tRNA.bed",sep=""),
#             sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
# rRNA=paste("/bucket/BourguignonU/Cong/termite_genome_annotation/barrnap/",
#            sp,"/filtered_",sp,"_barrnap.gff3",sep="")
# rRNA.gff=ape::read.gff(rRNA)
# rRNA.gff$ID=sapply(1:nrow(rRNA.gff),
#                    function(i){
#                      j=unlist(strsplit(rRNA.gff[i,"attributes"],";"))[1]
#                      j=sub("Name=","",j)
#                      return(paste(j,sprintf("%08d",i),sep=""))
#                    })
# rRNA.bed=data.frame(chrom=rRNA.gff$seqid,
#                     chromStart=format(rRNA.gff$start-1,scientific=FALSE,trim=TRUE),
#                     chromEnd=format(rRNA.gff$end,scientific=FALSE,trim=TRUE),
#                     name=paste(sp,"_",rRNA.gff$ID,sep=""),
#                     score=1,strand=rRNA.gff$strand)
# write.table(rRNA.bed,paste(dire,"/",sp,"_rRNA.bed",sep=""),
#             sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
# miRNA=paste("/bucket/BourguignonU/Cong/termite_genome_annotation/MirMachine/",
#             sp,"/results/predictions/filtered_gff/",sp,".PRE.tsv",sep="")
# miRNA=read.table(miRNA,header=TRUE,sep="\t",quote="")
# miRNA.bed=data.frame(chrom=miRNA$seqid,
#                     chromStart=format(miRNA$start-1,scientific=FALSE,trim=TRUE),
#                     chromEnd=format(miRNA$end,scientific=FALSE,trim=TRUE),
#                     name=paste(miRNA$seqID,"_miRNA",sep=""),
#                     score=1,strand=miRNA$strand)
# write.table(miRNA.bed,paste(dire,"/",sp,"_miRNA.bed",sep=""),
#             sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
#   print(sp)
# infernal.gff=paste("/bucket/BourguignonU/Cong/termite_genome_annotation/infernal/",
#                    sp,"/Infernal.gff3",sep="")
# infernal.gff=ape::read.gff(infernal.gff)
# infernal.gff=infernal.gff[!grepl("microRNA",infernal.gff$attributes) &
#                                    !grepl("rRNA",infernal.gff$attributes),]
# infernal.gff$ID=stringr::str_extract(infernal.gff$attributes,"ID=.*;Class")
# infernal.gff$ID=sub(";Class","",infernal.gff$ID)
# infernal.gff$ID=sub("ID=","",infernal.gff$ID)
# infernal.gff$Rfam=stringr::str_extract(infernal.gff$attributes,"Rfam=.*;Description")
# infernal.gff$Rfam=sub(";Description","",infernal.gff$Rfam)
# infernal.gff$Rfam=sub("Rfam=","",infernal.gff$Rfam)
# ncRNA.bed=data.frame(chrom=infernal.gff$seqid,
#                      chromStart=format(infernal.gff$start-1,scientific=FALSE,trim=TRUE),
#                      chromEnd=format(infernal.gff$end,scientific=FALSE,trim=TRUE),
#                      name=paste(sp,"_",infernal.gff$Rfam,"_",infernal.gff$ID,"_",sprintf("%08d",1:nrow(infernal.gff)),sep=""),
#                      score=1,strand=infernal.gff$strand)
# write.table(ncRNA.bed,paste(dire,"/",sp,"_ncRNA.bed",sep=""),
#             sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
# cmd=paste("cat",
#           paste(dire,"/",sp,"_mRNA.bed",sep=""),
#           paste(dire,"/",sp,"_PoReEl_up2k.bed",sep=""),
#           paste(dire,"/",sp,"_PoReEl_down2k.bed",sep=""),
#           paste(dire,"/",sp,"_tRNA.bed",sep=""),
#           paste(dire,"/",sp,"_rRNA.bed",sep=""),
#           paste(dire,"/",sp,"_miRNA.bed",sep=""),
#           paste(dire,"/",sp,"_ncRNA.bed",sep=""),
#           ">",paste(dire,"/",sp,"_integrated_annotation.bed",sep=""))
# print(cmd);system(cmd,wait=TRUE)
# }
#integrated_annotation


# gff=ape::read.gff("/bucket/BourguignonU/Cong/public_db/model_insects/Drosophila_melanogaster_GCF_000001215.4.gff")
# gff=gff[gff[,3]=="mRNA",]
# mRNA=sub(";.*$","",gff$attributes)
# mRNA=sub("ID=","",mRNA)
# 
# geneID=stringr::str_extract(gff$attributes,"GeneID:[0-9]*,")
# geneID=sub("GeneID:","",geneID)
# geneID=sub(",","",geneID)
# iso=data.frame(geneID,mRNA)
# write.table(iso,"/bucket/BourguignonU/Cong/public_db/Drosophila_melanogaster_genes/Dm_alternative_mRNA.tsv",
#             sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
# 
# bed=read.table("/bucket/BourguignonU/Cong/public_db/Drosophila_melanogaster_genes/Drosophila_melanogaster_GCF_000001215.4.bed",
#                sep="\t",header=FALSE,quote="")
# bed=bed[bed$V4 %in% mRNA,]
# write.table(bed,"/bucket/BourguignonU/Cong/public_db/Drosophila_melanogaster_genes/DmMrna.bed12",
#             sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
# 
# 
# DmGenes=read.table("/bucket/BourguignonU/Cong/public_db/Drosophila_melanogaster_genes/DmGenes.tsv",
#                    sep="\t",header=TRUE,quote="")
# bed=read.table("/bucket/BourguignonU/Cong/public_db/model_insects/Drosophila_melanogaster_GCF_000001215.4.bed",
#                sep="\t",header=FALSE,quote="")
# bed=bed[bed$V4 %in% DmGenes$protein,]


