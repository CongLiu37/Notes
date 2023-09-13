# Ortholog groups, gene functions & metabolic network

#####
# Ortholog groups
#####
# OrthoFinder
# OrthoFinder creates a results directory called ‘OrthoFinder’ inside the input proteome directory and puts the results here.
# Orthofinder with "-M msa" infers phylogeny by mafft & fasttree (too slow)
# Dependencies: Orthofinder, mafft, fasttree
Orthofinder=function(in_dir=in_dir, # Input proteome directory. 
                     # One file per species with extension '.faa'
                     threads=threads){
  threads=as.character(threads)
  
  cmd=paste("orthofinder.py",
            "-f",in_dir,
            "-t",threads,
            "-a",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# Re-run orthofinder with designated phylogenetic tree
Re_orthofinder=function(previous_orthofinder_result_dir=previous_orthofinder_result_dir,
                        tree=tree,
                        threads=threads){
  threads=as.character(threads)
  
  cmd=paste("orthofinder.py",
            "-t",threads,
            "-a",threads,
            "-ft",previous_orthofinder_result_dir,
            "-s",tree,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# OMA
# 不好用。。。

# Best protein-protein hit by blast
best_blastp=function(query.faa=query.faa,
                     reference.faa=reference.faa,
                     out_dir=out_dir,
                     out_basename=out_basename,
                     threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",query.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  query.faa=unlist(strsplit(query.faa,"/"));query.faa=query.faa[length(query.faa)]
  
  cmd=paste("cp",reference.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  reference.faa=unlist(strsplit(reference.faa,"/"));reference.faa=reference.faa[length(reference.faa)]
  
  cmd=paste("makeblastdb",
            "-in",reference.faa,
            "-dbtype","prot",
            "-parse_seqids",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("blastp",
            "-num_threads",threads,
            "-db",reference.faa,
            "-query",query.faa,
            "-outfmt 6",
            "-evalue 1e-5",
            "-num_alignments 1",
            "-out",paste(out_dir,"/",out_basename,".blast",sep=""),
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("rm",query.faa,sep=" "),wait=TRUE)
  system(paste("rm",reference.faa,sep=" "),wait=TRUE)
  system(paste("rm"," ",reference.faa,".*",sep=""),wait=TRUE)
  setwd(wd)
}

#####
# Gene functions
#####
# Interproscan
# Dependencies: Interproscan, seqkit, parallel (R)
# AT LEAST 16 THREADS
interpro=function(proteins.faa=proteins.faa,
                  out_dir=out_dir,
                  out_basename=out_basename,
                  threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",proteins.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  proteins.faa=unlist(strsplit(proteins.faa,"/"));proteins.faa=proteins.faa[length(proteins.faa)]
  cmd=paste("sed -i 's/*/X/g'",proteins.faa,sep=" ")
  print(cmd);system(cmd,wait=TRUE)

  cmd=paste("seqkit split2",proteins.faa,"-s 1000",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  prot.lst=system(paste("ls ",proteins.faa,".split/*",sep=""),intern=TRUE)
  chunks=as.numeric(threads)-1
  
  for (prot in prot.lst){
    cmd=paste("interproscan.sh",
              "-dp",
              "-b",basename(prot),
              "-cpu",1,
              "-f TSV",
              "-goterms",
              "-i",prot,
              "--pathways",
              sep=" ")
    system(paste("echo '",cmd,"' >> jobs",sep=""))
  }
  cmd=paste("split",
            "-d",
            "-n",paste("l/",chunks,sep=""),
            "jobs",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads)-1)
  parSapply(clus,
            1:chunks,
            function(i){
              scr=system("ls x*",wait=TRUE,intern=TRUE)[i]
              
              system(paste("sed -i '1i module load Other\\/interproscan\\/5.60-92.0'",scr,sep=" "))
              system(paste("sed -i '1i module load bioinfo-ugrp-modules'",scr,sep=" "))
              
              cmd=paste("bash"," ",scr,sep="")
              print(cmd);system(cmd,wait=TRUE)})
  stopCluster(clus)
  
  
  for (prot in prot.lst){
    system(paste("cat ",basename(prot),".tsv >> ",out_basename,".interpro.tsv",sep=""))
    system(paste("rm ",basename(prot),".tsv",sep=""))
  }
  
  system("rm x*")
  system(paste("rm",proteins.faa,sep=" "),wait=TRUE)
  system(paste("rm -r ",proteins.faa,".split/",sep=""))
  system("rm -r temp")
  setwd(wd)
}

# eggNOG-mapper
# Dependencies: eggNOG-mapper, seqkit, parallel (R)
eggNOGmapper=function(proteins.faa=proteins.faa,
                      out_dir=out_dir,
                      out_basename=out_basename,
                      threads=threads){
  threads=as.character(threads)
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",proteins.faa,out_dir,sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  proteins.faa=basename(proteins.faa)
  
  cmd=paste("seqkit split2",proteins.faa,"-s 1000",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  prot.lst=system(paste("ls ",proteins.faa,".split/*",sep=""),intern=TRUE)
  chunks=floor( as.numeric(threads)-1 )
  
  for (prot in prot.lst){
    cmd=paste("emapper.py",
              "-i",prot,
              "-o",basename(prot),
              "--cpu",1,
              sep=" ")
    system(paste("echo",cmd,">> jobs",sep=" "))
  }
  cmd=paste("split",
            "-d",
            "-n",paste("l/",chunks,sep=""),
            "jobs",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  library(parallel)
  clus=makeCluster(as.numeric(threads)-1)
  parSapply(clus,
            1:chunks,
            function(i){
              scr=system("ls x*",wait=TRUE,intern=TRUE)[i]
              cmd=paste("bash"," ",scr,sep="")
              print(cmd);system(cmd,wait=TRUE)})
  stopCluster(clus)
  
  system("rm x*")
  system(paste("rm",proteins.faa,sep=" "),wait=TRUE)
  system(paste("rm -r ",proteins.faa,".split/",sep=""))
  
  for (prot in prot.lst){
    system(paste("cat ",basename(prot),".emapper.annotations >> ",out_basename,".emapper.annotations.tsv",sep=""))
    system(paste("rm ",basename(prot),".emapper.annotations",sep=""))
    system(paste("cat ",basename(prot),".emapper.hits >> ",out_basename,".emapper.hits.tsv",sep=""))
    system(paste("rm ",basename(prot),".emapper.hits",sep=""))
    system(paste("cat ",basename(prot),".emapper.seed_orthologs >> ",out_basename,".emapper.seed_orthologs.tsv",sep=""))
    system(paste("rm ",basename(prot),".emapper.seed_orthologs",sep=""))
  }
  setwd(wd)
}
  
# Gene function by best blast hits
# Dependencies: DIAMOND+nr, parallel (R)
geneFun_bbh=function(pep.faa=pep.faa,
                     nr.dmdb="/apps/unit/BioinfoUgrp/DB/diamondDB/ncbi/2022-07/nr.dmnd",
                     out_prefix=out_prefix,
                     threads=threads){
  threads=as.character(threads)
  
  # blast search
  blast=paste(out_prefix,".blast",sep="")
  if (!file.exists(blast)){ # the task here is to find one best hits
    cmd=paste("diamond blastp",
              "--threads",threads,
              "--db",nr.dmdb,
              "--query",pep.faa,
              "--out",blast,
              "--min-score 50",
              "--query-cover 75",
              "--outfmt 6 qseqid sseqid length evalue bitscore stitle",
              sep=" ")
    print(cmd);system(cmd,wait=TRUE)
  }
  
  blast=read.table(blast,sep="\t",header=FALSE,quote="",comment.char="")
  colnames(blast)=c("qseqid","sseqid","length","evalue","bitscore","stitle")
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  clusterExport(clus,list("blast"),envir=environment())
  blast[,"retainHit"]=parSapply(clus,
                                1:nrow(blast),
                                function(i){
                                  query=blast[i,"qseqid"]
                                  evalue=blast[i,"evalue"]
                                  
                                  d=blast[blast[,"qseqid"]==query,]
                                  bestEval=min(d[,"evalue"])
                                  if (evalue!=bestEval){return(FALSE)}else{return(TRUE)}
                                })
  blast=blast[blast[,"retainHit"],]
  clusterExport(clus,list("blast"),envir=environment())
  blast[,"retainHit"]=parSapply(clus,
                                1:nrow(blast),
                                function(i){
                                  query=blast[i,"qseqid"]
                                  length=blast[i,"length"]
                                  
                                  d=blast[blast[,"qseqid"]==query,]
                                  longest=max(d[,"length"])
                                  if (length!=longest){return(FALSE)}else{return(TRUE)}
                                })
  stopCluster(clus)
  
  blast=blast[blast[,"retainHit"],]
  write.table(blast,paste(out_prefix,"_bbh.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
}

# KOfamScan
KOfamScan=function(pep.faa=pep.faa,
                   KOfamScan.ko_list="/bucket/BourguignonU/Cong/public_db/kofamscan/ko_list",
                   KOfamScan.profiles="/bucket/BourguignonU/Cong/public_db/kofamscan/profiles",
                   out_prefix=out_prefix,
                   threads=threads){
  threads=as.character(threads)
  
  cmd=paste("exec_annotation",
            "-o",paste(out_prefix,".tsv",sep=""),
            pep.faa,
            "-p",KOfamScan.profiles,
            "-k",KOfamScan.ko_list,
            paste("--cpu=",threads,sep=""),
            "-f mapper")
  print(cmd);system(cmd,wait=TRUE)
}
  
#####
# Metabolic network based on KEGG
#####
# Assemble KOs into draft metabolic network
# Dependencies: KEGGREST (R), stringr(r)
KO2Network=function(gene2ko=gene2ko, # KOfamScan output, mapper format
                    kegg.compound="/bucket/BourguignonU/Cong/public_db/kofamscan/compound", # https://rest.kegg.jp/list/compound
                    kegg.reaction="/bucket/BourguignonU/Cong/public_db/kofamscan/reaction", # https://rest.kegg.jp/list/reaction
                    out_prefix=out_prefix){
  gene2ko=readLines(gene2ko)
  d=data.frame(geneID=rep(NA,length(gene2ko)),KO=rep(NA,length(gene2ko)))
  d[,c(1,2)]=t(sapply(gene2ko,
                      function(i){
                        i=unlist(strsplit(i,"\t"))
                        return(c(i[1],i[2]))
                      }))
  d=d[!is.na(d[,"KO"]),] 
  gene2ko=d # geneID,KO

  library(KEGGREST)
  ko.lst=gene2ko[,"KO"][!duplicated(gene2ko[,"KO"])]
  reactions=c()
  for (i in seq(0,length(ko.lst),100)){
    KOs=ko.lst[(i+1):(i+100)]
    KOs=KOs[!is.na(KOs)]
    if (length(KOs)!=0){
      reactions=c(reactions,keggLink("reaction",KOs))
    }
  }
  ko2reaction=data.frame(KO=names(reactions),reactionID=unname(reactions))
  ko2reaction[,"KO"]=sub("^ko:","",ko2reaction[,"KO"])
  ko2reaction[,"reactionID"]=sub("^rn:","",ko2reaction[,"reactionID"])
  
  gene.lst=gene2ko[,"geneID"][!duplicated(gene2ko[,"geneID"])]
  ko.lst=unname(sapply(gene.lst,
                       function(gene){return(paste(gene2ko[gene2ko[,"geneID"]==gene,"KO"],collapse=","))}))
  reaction.lst=unname(sapply(ko.lst,
                             function(ko){
                              ko=unlist(strsplit(ko,","))
                              return(paste(ko2reaction[ko2reaction[,"KO"] %in% ko,"reactionID"],collapse=","))
                      }))
  gene2reaction=data.frame(geneID=gene.lst,KO=ko.lst,reactionID=reaction.lst)
  gene2reaction=gene2reaction[gene2reaction[,"reactionID"]!="",]
  write.table(gene2reaction,
              paste(out_prefix,"_gene2reaction.tsv",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
  
  reaction.lst=ko2reaction[,"reactionID"][!duplicated(ko2reaction[,"reactionID"])]
  equations=rep(NA,length(reaction.lst))
  for (i in seq(0,length(reaction.lst),10)){
    reac=reaction.lst[(i+1):(i+10)]
    reac=reac[!is.na(reac)]
    if (length(reac)!=0){
      eq=keggGet(reac)
      for (j in 1:length(eq)){
        equations[j+i]=eq[[j]]$EQUATION
      }
    }
  }
  reactants=unname(sapply(equations,function(equ){return( unlist(strsplit(equ,"<=>"))[1] )}))
  products=unname(sapply(equations,function(equ){return( unlist(strsplit(equ,"<=>"))[2] )}))
  network=data.frame(reactantID=c(),reactantName=c(),
                     productID=c(),productName=c(),
                     reactionID=c(),reactionName=c(),
                     equation=c())
  for (i in 1:length(reaction.lst)){
    r=reaction.lst[i] # R13178
    rname=system(paste("awk -F '\t' -v OFS='\t' '{if ($1==\"",r,"\") print $2}' ",kegg.reaction,sep=""),intern=TRUE)
    if (length(rname)!=0){
      rname=unlist(strsplit(rname,";"))[1]
      e=equations[i]
      # print(c(r,rname,e))
      library(stringr)
      from=unlist(strsplit(reactants[i]," "));from=str_extract(from,"[CG][0-9]{5}");from=from[!is.na(from)]
      to=unlist(strsplit(products[i]," "));to=str_extract(to,"[CG][0-9]{5}");to=to[!is.na(to)]
      f=function(ids){
        id_c=ids[grepl("C[0-9]{5}",ids)]
        id_g=ids[!grepl("C[0-9]{5}",ids)]
        id_g2c=c()
        if (length(id_g)!=0){
          id_g2c=keggLink("compound",id_g);id_g2c=unname(sub("^cpd:","",id_g2c))
        }
        return(c(id_c,id_g2c))
      }
      from=f(from);to=f(to)
      d=merge(from,to);colnames(d)=c("reactantID","productID")
      d[,"reactantName"]=sapply(d[,"reactantID"],
                                function(id){
                                  return( unlist(strsplit(system(paste("awk -F '\t' -v OFS='\t' '{if ($1==\"",id,"\") print $2}' ",kegg.compound,sep=""),
                                                 intern=TRUE),";"))[1] )
                                })
      d[,"productName"]=sapply(d[,"productID"],
                               function(id){
                                 return( unlist(strsplit(system(paste("awk -F '\t' -v OFS='\t' '{if ($1==\"",id,"\") print $2}' ",kegg.compound,sep=""),
                                                                intern=TRUE),";"))[1] )
                               })
      d[,"reactionID"]=rep(r,nrow(d))
      d[,"reactionName"]=rep(rname,nrow(d))
      d[,"equation"]=rep(e,nrow(d))
      network=rbind(network,d)
    }
  }
  network[,"KO"]=sapply(network[,"reactionID"],
                        function(r){
                          return(paste(ko2reaction[ko2reaction[,"reactionID"]==r,"KO"],collapse=","))
                        })
  network[,"geneID"]=sapply(network[,"KO"],
                            function(ko){
                              ko=unlist(strsplit(ko,","))
                              return(paste(gene2ko[gene2ko[,"KO"] %in% ko,"geneID"],collapse = ","))
                            })
  write.table(network,paste(out_prefix,"_network.tsv",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
}

# Curate metabolic network
# Dependencies: KEGGREST (R)
curateNetwork=function(network.tsv=network.tsv, # from KO2Network
                       metaboliteConcentrations="/bucket/BourguignonU/Cong/public_db/Park2016/metConcentrationRanges.tsv",
                       # tsv with header
                       # Fields: compound (Compound ID of KEGG)
                       #         low (lower bound of absolute concentration in mol/m^3 or M)
                       #         high (higher bound of absolute concentration in mol/m^3 or M)
                       standardGibbs="/bucket/BourguignonU/Cong/public_db/eQuilibrator/standardGibbs.tsv",
                       # tsv with header
                       # Fields: reaction (Reaction ID of KEGG)
                       #         Gibbs (Standard Gibbs free energy of reaction in J/mol, temperature=298.15K, all reactants and products are of 1 mol/m^3)
                       out_prefix=out_prefix){
  # network.tsv="/Users/congliu/Desktop/PhD/Results/termite_pca/KO2Network/Nluj_network.tsv"
  # metaboliteConcentrations="~/Desktop/PhD/Results/public_db/Park2016/metConcentrationRanges.tsv"
  # standardGibbs="~/Desktop/PhD/Results/public_db/eQuilibrator/standardGibbs.tsv"
  # out_prefix="/Users/congliu/Desktop/PhD/Results/termite_pca/KO2Network/Nluj"
  # 
  network=read.table(network.tsv,sep="\t",header=TRUE,quote="")
  metaboliteConcentrations=read.table(metaboliteConcentrations,header=TRUE,sep="\t",quote ="")
  rownames(metaboliteConcentrations)=metaboliteConcentrations[,"compound"]
  metaboliteConcentrations=metaboliteConcentrations[metaboliteConcentrations$low!=0 & metaboliteConcentrations$high!=0,]
  standardGibbs=read.table(standardGibbs,header=TRUE,sep="\t",quote="")
  rownames(standardGibbs)=standardGibbs[,"reaction"]
  
  library(KEGGREST)
  compounds=c(network$reactantID,network$productID);compounds=compounds[!duplicated(compounds)]
  formula=rep(NA,length(compounds))
  for (i in seq(0,length(compounds),10)){
    compound=compounds[(i+1):(i+10)]
    compound=compound[!is.na(compound)]
    if (length(compound)!=0){
      f=keggGet(compound)
      for (j in 1:length(f)){
        if (!is.null(f[[j]]$FORMULA)){
          formula[j+i]=f[[j]]$FORMULA
        }
      }
    }
  }
  compound2formula=data.frame(compound=compounds,formula=formula)
  rownames(compound2formula)=compound2formula$compound
  network$reactantFormula=sapply(1:nrow(network),function(i){ return(compound2formula[network[i,"reactantID"],"formula"]) })
  network$productFormula=sapply(1:nrow(network),function(i){ return(compound2formula[network[i,"productID"],"formula"]) })
  
  # I. filter DNA/RNA/Protein reactions
  network=network[!(network[,"reactionID"] %in% 
                     network[network[,"reactantName"]==network[,"productName"],"reactionID"]),]
  network=network[!(network[,"reactionID"] %in% 
                      network[grepl("[DR]NA",network[,"reactantName"]) | 
                                grepl("[DR]NA",network[,"productName"]),"reactionID"]),]
  network=network[!(network[,"reactionID"] %in% 
                     network[grepl("[Pp]rotein",network[,"reactantName"]) | 
                               grepl("[Pp]rotein",network[,"productName"]),"reactionID"]),]
  ambigousCompounds=compound2formula[is.na(compound2formula$formula),"compound"]
  network=network[!(network[,"reactionID"] %in% 
                      network[network[,"reactantID"] %in% ambigousCompounds | 
                                network[,"productID"] %in% ambigousCompounds,"reactionID"]),]
  # II. filter compounds:
  # ambigously defined compounds
  RCompounds=compound2formula[grepl("R",compound2formula$formula),"compound"]
  network=network[!(network[,"reactantID"] %in% RCompounds) &
                    !(network[,"productID"] %in% RCompounds),]
  # water, hydrogen peroxide, oxygen, hydrogen
  network=network[network[,"reactantName"]!="H2O" & network[,"productName"]!="H2O",]
  network=network[network[,"reactantName"]!="Hydrogen peroxide" & network[,"productName"]!="Hydrogen peroxide",]
  network=network[network[,"reactantName"]!="Oxygen" & network[,"productName"]!="Oxygen",]
  network=network[network[,"reactantName"]!="H+" & network[,"productName"]!="H+",]
  # Phosphate
  network=network[network[,"reactantName"]!="Diphosphate" & network[,"productName"]!="Diphosphate",]
  network=network[network[,"reactantName"]!="Orthophosphate" & network[,"productName"]!="Orthophosphate",]
  network=network[network[,"reactantName"]!="Polyphosphate" & network[,"productName"]!="Polyphosphate",]
  # NTP, NMP, dNTP, dNMP
  NTP=c("ATP","TTP","GTP","CTP","UTP","Nucleoside triphosphate");dNTP=paste("d",NTP,sep="")
  NDP=c("ADP","TDP","GDP","CDP","UDP","NDP");dNDP=paste("d",NDP,sep="")
  # NMP=c("AMP","TMP","GMP","CMP","UMP");dNMP=paste("d",NMP,sep="")
  network=network[!(network[,"reactantName"] %in% c(NTP,dNTP,NDP,dNDP)) 
                  & !(network[,"productName"] %in% c(NTP,dNTP,NDP,dNDP)),]
  # H+/e- donors
  network=network[network[,"reactantName"]!="NAD+" & network[,"productName"]!="NAD+",]
  network=network[network[,"reactantName"]!="NADP+" & network[,"productName"]!="NADP+",]
  network=network[network[,"reactantName"]!="NADH" & network[,"productName"]!="NADH",]
  network=network[network[,"reactantName"]!="NADPH" & network[,"productName"]!="NADPH",]
  network=network[network[,"reactantName"]!="FAD" & network[,"productName"]!="FAD",]
  network=network[network[,"reactantName"]!="FADH2" & network[,"productName"]!="FADH2",]
  network=network[network[,"reactantName"]!="FMN" & network[,"productName"]!="FMN",]
  network=network[network[,"reactantName"]!="Reduced FMN" & network[,"productName"]!="Reduced FMN",]
  # CoA
  network=network[network[,"reactantName"]!="CoA" & network[,"productName"]!="CoA",]
  # III. direction of reactions
  reaction2equation=network[,c("reactionID","equation")]
  reaction2equation=reaction2equation[!duplicated(reaction2equation[,"reactionID"]),]
  R=8.3144598 # universal gas constatnt
  Temperature=298.15 # temperature
  reaction2equation[,"direction"]=
    sapply(1:nrow(reaction2equation),
           function(i){
             direction=0
             
             reac=reaction2equation[i,"reactionID"]
             stan.Gibbs=standardGibbs[reac,"Gibbs"]
             
             equa=reaction2equation[i,"equation"]
             equa=unlist(strsplit(equa,"<=>"))
             reactant=equa[1];product=equa[2]
             reactant=unlist(strsplit(reactant,"\\+"));product=unlist(strsplit(product,"\\+"))
             reactant=sub(" $","",sub("^ ","",reactant));product=sub(" $","",sub("^ ","",product))
             d.reactant=data.frame(reactant=rep(NA,length(reactant)),sto=rep(NA,length(reactant)),
                                   low=rep(NA,length(reactant)),high=rep(NA,length(reactant)))
             for (i in 1:nrow(d.reactant)){
               a=unlist(strsplit(reactant[i]," "))
               if (a[length(a)]!="C00001"){ # water
                 d.reactant[i,"reactant"]=a[length(a)]
                 if (a[1]!=a[length(a)]){d.reactant[i,"sto"]=as.numeric(a[1])}else{d.reactant[i,"sto"]=1}
                 l=metaboliteConcentrations[a[length(a)],"low"];h=metaboliteConcentrations[a[length(a)],"high"]
                 if (a[length(a)]=="C00080"){l=1e-7;h=1e-7} # H+: 1e-7 mol/L, pH=7
                 if (is.na(l)){l=1e-5;h=1e+5}else{l=l*1e-3/100;h=h*1e-3*100}
                 d.reactant[i,c("low","high")]=c(l,h)
               }
             }
             d.reactant=d.reactant[!is.na(d.reactant[,"reactant"]),]
             d.product=data.frame(product=rep(NA,length(product)),sto=rep(NA,length(product)),
                                  low=rep(NA,length(product)),high=rep(NA,length(product)))
             for (i in 1:nrow(d.product)){
               a=unlist(strsplit(product[i]," "))
               if (a[length(a)]!="C00001"){ # water
                 d.product[i,"product"]=a[length(a)]
                 if (a[1]!=a[length(a)]){d.product[i,"sto"]=as.numeric(a[1])}else{d.product[i,"sto"]=1}
                 l=metaboliteConcentrations[a[length(a)],"low"];h=metaboliteConcentrations[a[length(a)],"high"]
                 if (a[length(a)]=="C00080"){l=1e-7;h=1e-7} # H+
                 if (is.na(l)){l=1e-5;h=1e+5}else{l=l*1e-3/100;h=h*1e-3*100}
                 d.product[i,c("low","high")]=c(l,h)
               }
             }
             d.product=d.product[!is.na(d.product[,"product"]),]
             
             if (!is.na(stan.Gibbs)){
               Q.min=prod(d.product[,"low"]^d.product[,"sto"])/prod(d.reactant[,"high"]^d.reactant[,"sto"])
               Q.max=prod(d.product[,"high"]^d.product[,"sto"])/prod(d.reactant[,"low"]^d.reactant[,"sto"])
               deltaG.min=stan.Gibbs+R*Temperature*log(Q.min)
               deltaG.max=stan.Gibbs+R*Temperature*log(Q.max)
               if (deltaG.min>0){direction=-1}
               if (deltaG.max<0){direction=1}
               if (deltaG.min<0 & deltaG.max>0){direction=0}
             }
             # CO2 can only be synthesized
             if ("C00011" %in% d.product[,"product"]){direction=1}
             if ("C00011" %in% d.reactant[,"reactant"]){direction=-1}
             # O2 can only be consumed
             if ("C00007" %in% d.product[,"product"]){direction=-1}
             if ("C00007" %in% d.reactant[,"reactant"]){direction=1}
             # ATP can only be consumed
             # C00002, C00131, C00459, C00044, C00286, C00063, C00458, C00201
             if ("C00002" %in% d.product[,"product"]){direction=-1}
             if ("C00002" %in% d.reactant[,"reactant"]){direction=1}
             return(direction)
           })
  # direction:
  #   0: reversible
  #   1: left to right
  #   -1: right to left
  
  # Correct directions
  wrongDir=reaction2equation[reaction2equation[,"direction"]==-1,"reactionID"]
  wrongNetwork=network[network$reactionID %in% wrongDir,]
  wrongNetwork=wrongNetwork[,c("productID","reactantID","productName","reactantName",
                               "reactionID","reactionName","equation","KO","geneID",
                               "productFormula","reactantFormula")]
  colnames(wrongNetwork)=colnames(network)
  correctNetwork=network[!(network$reactionID %in% wrongDir),]
  network=rbind(wrongNetwork,correctNetwork)
  network[,"reversible"]=sapply(network[,"reactionID"],
                                function(ID){
                                  direction=reaction2equation[reaction2equation$reactionID==ID,"direction"]
                                  if (direction==0){return(TRUE)}else{return(FALSE)}
                                })
  write.table(network,paste(out_prefix,"_curateNetwork.tsv",sep=""),
              sep="\t",row.names=FALSE,quote=FALSE)
}


#network=network[!grepl("G[0-9]{5}",network[,"equation"]),]
#network=network[!grepl("C[0-9]{5}\\(n\\)",network[,"equation"]),]
#Alkyl=network[grepl("[Aa]lkyl",network[,"reactantName"]) | grepl("[Aa]lkyl",network[,"productName"]),"reactionID"]
# # 金属
# network=network[network[,"reactantName"]!="Fe2+" & network[,"productName"]!="Fe2+",] # Fe2+
# network=network[network[,"reactantName"]!="Fe3+" & network[,"productName"]!="Fe3+",] # Fe3+
# network=network[network[,"reactantName"]!="Zinc cation" & network[,"productName"]!="Zinc cation",] # Zn2+
# network=network[network[,"reactantName"]!="Calcium cation" & network[,"productName"]!="Calcium cation",] # Ca2+
# network=network[network[,"reactantName"]!="Cobalt ion" & network[,"productName"]!="Cobalt ion",] # Co2+
# network=network[network[,"reactantName"]!="Potassium cation" & network[,"productName"]!="Potassium cation",] # K+
# network=network[network[,"reactantName"]!="Magnesium cation" & network[,"productName"]!="Magnesium cation",] # Mg2+
# network=network[network[,"reactantName"]!="Mercury(2+)" & network[,"productName"]!="Mercury(2+)",] # Hg+
# network=network[network[,"reactantName"]!="Sodium cation" & network[,"productName"]!="Sodium cation",] # Na+
# network=network[network[,"reactantName"]!="Barium cation" & network[,"productName"]!="Barium cation",] # Ba2+
# network=network[network[,"reactantName"]!="Cu2+" & network[,"productName"]!="Cu2+",] # Cu2+
# network=network[network[,"reactantName"]!="Cu+" & network[,"productName"]!="Cu+",] # Cu+
# network=network[network[,"reactantName"]!="Manganese(3+)" & network[,"productName"]!="Manganese(3+)",] # Mn3+
# network=network[network[,"reactantName"]!="Manganese(2+)" & network[,"productName"]!="Manganese(2+)",] # Mn2+
# network=network[network[,"reactantName"]!="Nickel(2+)" & network[,"productName"]!="Nickel(2+)",] # Ni2+
# network=network[network[,"reactantName"]!="Rubidium cation" & network[,"productName"]!="Rubidium cation",] # Rb+
# network=network[network[,"reactantName"]!="Strontium cation" & network[,"productName"]!="Strontium cation",] # Sr2+