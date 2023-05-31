arg=commandArgs(trailingOnly = TRUE)

main.R=arg[1]
conf=arg[2]

source(main.R)
source(conf)

out_dir=sub("/$","",out_dir)
threads=as.character(threads)
if (!file.exists(out_dir)){system(paste("mkdir",out_dir,sep=" "))}
label=unlist(strsplit(label,"/"));label=label[length(label)]

PseudoPipe(genome=masked_genome, # soft masked
           # protein.fa=paste(
           #   paste(save_dir,"/protein_2/",label,"/",label,"_proteins_iso.faa",sep=""),
           #   paste(save_dir,"/protein_2/",label,"/",label,"_proteins_rep.faa",sep=""),
           #   sep=" "), # space list
           protein.fa=protein.fa,
           gff=gff,
           out_dir=out_dir,
           threads=threads)
PseudoPipe2gff3(Pseudo.out.txt=paste(out_dir,"/",label,"_pgenes.txt",sep=""),
                species=label,
                out.gff3=paste(out_dir,"/",label,"_pseudogene.gff3",sep=""))