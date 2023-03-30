# tRNAscan-SE
# Dependencies: tRNAscan-SE, biocode
tRNAscan=function(genome=genome,
                  mode="eukaryotic", # eukaryotic/bacterial/archaeal/mitochondrial_mammal/mitochondrial_vert/other_organellar/general
                  out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  if (mode=="eukaryotic"){mode="-E"}
  if (mode=="bacterial"){mode="-B"}
  if (mode=="archaeal"){mode="-A"}
  if (mode=="mitochondrial_mammal"){mode="-M mammal"}
  if (mode=="mitochondrial_vert"){mode="-M vert"}
  if (mode=="other_organellar"){mode="-O"}
  if (mode=="general"){mode="G"}
  
  cmd=paste("tRNAscan-SE",
            mode,
            "-o ./final_results",
            "-f ./secondary_structures",
            "-m ./summary",
            genome,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="convert_tRNAScanSE_to_gff3.pl -g --input=final_results > tRNA.gff3"
  print(cmd);system(cmd,wait=TRUE)
  
  setwd(out_dir)
}

# miRNAture: Identify microRNA
miRNAture=function(genome=genome,
                   dataF=dataF, # pre-calculated data directory from miRNAture
                   species=species, # Genus_species
                   out_dir=out_dir,
                   threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  system("source ~/miniconda3/etc/profile.d/conda.sh")
  system("conda init bash")
  system("conda activate mirnature")
  
  cmd=paste("cp -r",dataF,"./dataF",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",genome,"./genome.fa",sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("miRNAture",
            "-stage complete",
            "-dataF","dataF",
            "-speG","./genome.fa",
            "-speN",species,
            "-speT mrna",
            "-w",out_dir,
            "-m hmm,rfam,mirbase,infernal,final",
            "-pe 0",
            "-nbitscore_cut 1.0",
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system("mv ./Final_miRNA_evaluation/miRNA_annotation_mrna_accepted_conf.gff3 miRNAture.gff3")
  system("cat ./Final_miRNA_evaluation/Fasta/Drosophila_mirna_miRNAs_high_confidence.fasta Final_miRNA_evaluation/Fasta/Drosophila_mirna_miRNAs_medium_confidence.fasta > miRNA.fna")
  system("rm -r ./dataF")
  system("rm ./genome.fa")
  system("rm -r LOGS")
  system("rm -r miRNA_prediction")
  system("rm miRNAture_configuration_mrna.yaml")
  system("rm -r miRNA_validation")
  system("rm -r TemporalFiles")
  system("rm -r Final_miRNA_evaluation")
  setwd(wd)
}

# miRanda: microRNA target gene identification
# Dependencies: miRanda, bedtools, seqkit, parallel (R)
miranda=function(gff=gff, # three_prime_UTR feature required
                 genome=genome,
                 miRNA.fna=miRNA.fna,
                 threads=threads,
                 out_dir=out_dir){
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  threads=as.character(threads)
  
  cmd=paste("grep 'three_prime_UTR' ",gff,
            " | awk -F '\t' -v OFS='\t' '{print $1,$4-1,$5,$9,200,$7}' | sed 's/ID=.*;Parent=//' | sed 's/;//' > 3_UTR.bed",
            sep="")
  print(cmd);system(cmd,wait=TRUE)
  system(paste("cp",genome,"./genome.fa",sep=" "))
  cmd=paste("bedtools getfasta -s -name",
            "-fi","./genome.fa",
            "-bed 3_UTR.bed",
            " | sed 's/::.*//' > 3_UTR.fa")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd=paste("seqkit split2",
            "3_UTR.fa",
            "--by-part",as.character(as.numeric(threads)-1),
            "--out-dir","tmp",
            "--threads",threads,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
  
  system(paste("cp",miRNA.fna,"./miRNA.fna",sep=" "))
  library(parallel)
  clus=makeCluster(as.numeric(threads))
  parSapply(clus,
            1:(as.numeric(threads)-1),
            function(i){
              utr=system("ls tmp/*.fa",intern=TRUE)[i]
              cmd=paste("miranda","miRNA.fna",utr,
                        "-sc 150","-en -20",
                        "| grep '>>' | sed 's/>>//' | sort -k 5 -n -r",
                        "| awk -F '\t' -v OFS='\t' '{print $1,$2}' | sort | uniq -u >",
                        paste("tmp/",as.character(i),".tsv",sep=""),
                        sep=" ")
              system(cmd,wait=TRUE)
            })
  stopCluster(clus)
  
  o=system("ls tmp/*.tsv",intern=TRUE)
  cmd=paste("cat",paste(o,collapse=" "),"> miRanda.tsv")
  system(cmd)
  
  system("rm ./miRNA.fna")
  system("rm ./3_UTR.fa")
  system("rm ./3_UTR.bed")
  system("rm ./genome.fa")
  system("rm ./genome.fa.fai")
  system("rm -r tmp")
  setwd(wd)
}

# Infernal: non-coding RNA gene identification
# cmscan treats uppercase and lowercase nucleotides identically
# Dependencies: Infernal, Rfam database, stringr (R)
infernal=function(genome=genome,
                  domain=domain, # eukarya/archaea/bacteria Which domain shall be retained?
                  Rfam.clanin=Rfam.clanin,
                  Rfam.cm=Rfam.cm,
                  out_dir=out_dir,
                  threads=threads){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cmscan",
            "--cut_ga","--rfam","--nohmmonly",
            "--tblout Infernal.tblout",
            "--fmt 2",
            "--clanin",Rfam.clanin,
            "--cpu",threads,
            Rfam.cm,genome,
            "> Infernal.cmscan")
  print(cmd);system(cmd,wait=TRUE)
  
  domains=c("eukarya","archaea","bacteria")
  for (i in 1:3){if (grepl(domain,domains[i])){domains=domains[-i]}}
  cmd=paste("grep -v '=' Infernal.tblout | ",
            "grep -v 'tRNA' | ",
            "grep -v '",domains[1],"' | ",
            "grep -v '",domains[2],"' | ",
            "grep -v 'microRNA'",
            " > final.tblout",
            sep="")
  print(cmd);system(cmd,wait=TRUE)
  
  cmd="awk 'NR>2 {if ($1!=\"#\") print $4\"\tInfernal\t\"$2\"\t\"$10\"\t\"$11\"\t\"$17\"\t\"$12\"\t\\.\tID=\"$2\"_\"$1\";Class=\"$2\";Rfam=\"$3\";Description=\"}' final.tblout  | awk -F '\t' -v OFS='\t' '{if ($7==\"-\") print $1,$2,$3,$5,$4,$6,$7,$8,$9; else print $0}' > Infernal.gff3"
  print(cmd);system(cmd,wait=TRUE)
  Des=system("sed 's/.*- //' final.tblout",intern=TRUE)
  Des=sapply(1:length(Des),
             function(i){
               if (grepl("#",Des[i])){return(NA)}else{return(Des[i])}
             })
  Des=Des[!is.na(Des)]
  Des=Des[-1]
  gff=readLines("Infernal.gff3")
  
  out=sapply(1:length(Des),
             function(i){
               return(paste(gff[i],Des[i],sep=""))
             })
  writeLines(out,"Infernal.gff3")
  setwd(wd)
}

# Integrate ncRNA results from tRNAscan-SE, miRNAture, infernal and miRanda
ncRNA=function(tRNA.gff3=tRNA.gff3,
               miRNAture.gff3=miRNAture.gff3,
               miRNA.fa=miRNA.fa,
               infernal.gff3=infernal.gff3,
               miRanda.tsv=miRanda.tsv,
               species=species,
               out_dir=out_dir){
  out_dir=sub("/$","",out_dir)
  if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "),wait=TRUE)}
  wd=getwd();setwd(out_dir)
  
  cmd=paste("cp",tRNA.gff3,
            paste(species,"_tRNA.gff3",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",miRNAture.gff3,
            paste(species,"_miRNA.gff3",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",miRNA.fa,
            paste(species,"_miRNA.fa",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",infernal.gff3,
            paste(species,"_Rfam.gff3",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  cmd=paste("cp",miRanda.tsv,
            paste(species,"_miRNA_target.tsv",sep=""))
  print(cmd);system(cmd,wait=TRUE)
  setwd(wd)
}