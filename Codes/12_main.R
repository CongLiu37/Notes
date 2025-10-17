# Statistics

# Coefficient of determination
coefficient_of_determination=
  function(y.observed=y.observed,y.predicted=y.predicted){
    # mean of observed value
    y.observed_bar=mean(y.observed)
  
    #residual sum of squares
    ssr=sum( (y.observed-y.predicted)^2 )
  
    # total sum of squares
    sst=sum( (y.observed-y.observed_bar)^2 )
  
    # coefficient of determination
    cod=1-ssr/sst
  
    return(cod)
}

# Hierarchical cluster
# Dependencies: ape (R)
hierarchical_cluster=function(df=df,
                              measure="euclidean", # euclidean/maximum/manhattan/canberra/binary/minkowski
                              by.rows=TRUE){
  if (!by.rows){df=t(df)}
  distance.measure=dist(x=df,method=measure)
  clustering=hclust(distance.measure)
  tree=ape::as.phylo(clustering,use.labels=TRUE)
  return(tree)
}

# expand.grid remove redundant
expand.grid.nr=function(vector.1=vector.1,vector.2=vector.2){
  df=expand.grid(vector.1,vector.2,stringsAsFactors=F)
  df=unique(t(apply(df, 1, sort)))
  df=as.data.frame(df)
  df=df[df$V1!=df$V2,]
  return(df)
}

# Correlation test: Pearson's r-square, Kendall's tau, Spearman's rho
correlation_text=function(x=x,y=y){
  pearson.correlation=cor.test(x,y,method="pearson",alternative="two.sided")
  pearson.positive=cor.test(x,y,method="pearson",alternative="greater")
  pearson.negative=cor.test(x,y,method="pearson",alternative="less")
  
  kendall.correlation=cor.test(x,y,method="kendall",alternative="two.sided")
  kendall.positive=cor.test(x,y,method="kendall",alternative="greater")
  kendall.negative=cor.test(x,y,method="kendall",alternative="less")
  
  spearman.correlation=cor.test(x,y,method="spearman",alternative="two.sided")
  spearman.positive=cor.test(x,y,method="spearman",alternative="greater")
  spearman.negative=cor.test(x,y,method="spearman",alternative="less")
  
  res=list(x=x,y=y,
           
           pearson.cor=unname(pearson.correlation$estimate),
           pearson.correlation_p=unname(pearson.correlation$p.value),
           pearson.positive_p=unname(pearson.positive$p.value),
           pearson.negative_p=unname(pearson.negative$p.value),
           
           kendall.tau=unname(kendall.correlation$estimate),
           kendall.correlation_p=unname(kendall.correlation$p.value),
           kendall.positive_p=unname(kendall.positive$p.value),
           kendall.negative_p=unname(kendall.negative$p.value),
           
           spearman.rho=unname(spearman.correlation$estimate),
           spearman.correlation_p=unname(spearman.correlation$p.value),
           spearman.positive_p=unname(spearman.positive$p.value),
           spearman.negative_p=unname(spearman.negative$p.value)
           )
  return(res)
}

# Phylogenetic generalized least squares (PGLS): select correlation structure using AIC
# Depedencies: ape (R), nlme (R)
pgls=function(x=x,y=y,tips=tips,
              phylo=phylo){
  data=data.frame(x=x,y=y)
  rownames(data)=tips
  
  # Brownian correlation structure
  brownian=ape::corBrownian(phy=phylo,form=~tips)
  pgls.brownian=nlme::gls(y~x,correlation=brownian,data=data)
  pgls.brownian.sum=summary(pgls.brownian)
  AIC.brownian=pgls.brownian.sum$AIC
  
  # Pagel's lambda correlation structure
  pagel=ape::corPagel(value=1,phy=phylo,form=~tips,fixed=FALSE)
  pgls.pagel=nlme::gls(y~x,correlation=pagel,data=data)
  pgls.pagel.sum=summary(pgls.pagel)
  AIC.pagel=pgls.pagel.sum$AIC
  
  # Grafen's correlation structure
  grafen=ape::corGrafen(value=1,phy=phylo,form=~tips,fixed=FALSE)
  pgls.grafen=nlme::gls(y~x,correlation=grafen,data=data)
  pgls.grafen.sum=summary(pgls.grafen)
  AIC.grafen=pgls.grafen.sum$AIC
  
  # Martin's correlation structure
  martin=ape::corMartins(value=1,phy=phylo,form=~tips,fixed=FALSE)
  pgls.martin=nlme::gls(y~x,correlation=martin,data=data)
  pgls.martin.sum=summary(pgls.martin)
  AIC.martin=pgls.martin.sum$AIC
  
  # Blomberg's correlation structure
  blomberg=ape::corBlomberg(value=1,phy=phylo,form=~tips,fixed=FALSE)
  pgls.blomberg=nlme::gls(y~x,correlation=blomberg,data=data)
  pgls.blomberg.sum=summary(pgls.blomberg)
  AIC.blomberg=pgls.blomberg.sum$AIC
  
  AICs=c(AIC.brownian,AIC.pagel,AIC.grafen,AIC.martin,AIC.blomberg)
  names(AICs)=c("brownian","pagel","grafen","martin","blomberg")
  best=names(which(AICs==max(AICs)))
  return(get(paste("pgls.",best,".sum",sep="")))
}

# angle of two vectors
angle_of_vectors=function(a,b){
  theta=acos( sum(a*b) / ( sqrt(sum(a*a)) * sqrt(sum(b*b)) ) )
  return(theta)
}

# Shannon diversity
shannon=function(x=x){
  x=x[which(x!=0)]
  x=x/sum(x)
  return( -sum(x*log(x)) )
}

# Hill number
hill=function(x=x,
              order=order){
  # order: 0,1,2
  x=x[which(x!=0)]
  x=x/sum(x)
  if (order==0 | order==2){return( sum(x^order)^(1/(1-order)) )}
  if (order==1){return( exp(-sum(x*log(x))) )}
  
}

# Species richness
richness=function(x){
  x=x[which(x!=0)]
  return(length(x))
}

# Simpson diversity
simpson=function(x){
  x=x[which(x!=0)]
  x=x/sum(x)
  return( sum(x^2) )
}

# reverse Simpson diversity, or Hill number of order 2
reverse_simpson=function(x){
  x=x[which(x!=0)]
  x=x/sum(x)
  return( 1/sum(x^2) )
}

# Gini-Simpson diversity
gini_simpson=function(x){
  x=x[which(x!=0)]
  x=x/sum(x)
  return( 1-sum(x^2) )
}

# Pielou evenness
pielou=function(x){
  x=x[which(x!=0)]
  x=x/sum(x)
  H=-sum(x*log(x)) # Shannon diversity
  S=length(x) # species richness
  return( H/log(S) )
}

# Euclidean distance
euclidean=function(x,y){
  return( sqrt( sum((x-y)^2) ) )
}

# Bray-Curtis dissimilarity
bray_curtis=function(x,y){
  return( 1-2*sum( c(x[which(x<y)],y[which(y<=x)]) ) / sum(x+y) )
}

# Kullback-Leibler divergence of y to x
# Jensen-Shannon divergence 

# likelyhood ratio test
LRT=function(loglik.simple.model=loglik.simple.model,
             loglik.complex.model=loglik.complex.model,
             df=df){
  # df: the difference in #parameters in complex and simple model
  # H0: You should use the simple model.
  # H1: You should use the complex model.
  lr=-2*(loglik.simple.model-loglik.complex.model)
  p=pchisq(lr, df = df, lower.tail = FALSE)
  return(p)
}

# gff->bed->sequences
gff2bed2seq=function(gff=gff,
                     fna=fna,
                     out.fna=out.fna){
  gff=ape::read.gff(gff)
  gff$name=sapply(gff[,9],
                  function(i){
                    i=unlist(strsplit(i,";"))
                    i=i[grepl("^Name=",i)]
                    i=sub("^Name=","",i)
                    i=sub(";$","",i)
                    return(i)
                  })
  bed=gff[,c(1,4,5,10)]
  bed[,2]=bed[,2]-1
  bed[,5]="."
  bed[,6]=gff[,7]
  write.table(bed,paste(out.fna,".bed",sep=""),
              sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  
  cmd=paste("bedtools getfasta",
            "-fi",fna,
            "-bed",paste(out.fna,".bed",sep=""),
            "-nameOnly -s | ",
            "sed 's/([+-])//' > ",
            out.fna,
            sep=" ")
  print(cmd);system(cmd,wait=TRUE)
}

# evolutionary time (eTime, substitutions per site) to calender time
# Dependencies: ape (R), ggtree (R), tidytree (R)
eTime2cTime=function(eTime.nwk="~/Desktop/PhD/Results/termite_pca/molEvol_rate/single_copy_4dtv/4dtv.nwk",
                     cTime.nwk="~/Desktop/PhD/Results/termite_pca/timeTree_label.nwk",
                     sp=sp,eTime=eTime){
  # eTime.nwk="~/Desktop/PhD/Results/termite_pca/molEvol_rate/single_copy_4dtv/4dtv.nwk"
  # cTime.nwk="~/Desktop/PhD/Results/termite_pca/timeTree_label.nwk"
  # eTime=0.1589597
  # sp="Znev"
  
  if (file.exists(eTime.nwk)){
    eTime.tree=ape::read.tree(file=eTime.nwk)
  }else{
    eTime.tree=ape::read.tree(text=eTime.nwk)
  }
  eTime.tree.data=as.data.frame(ggtree::ggtree(eTime.tree)$data)
  
  if (file.exists(cTime.nwk)){
    cTime.tree=ape::read.tree(file=cTime.nwk)
  }else{
    cTime.tree=ape::read.tree(text=cTime.nwk)
  }
  
  cTime.tree.data=as.data.frame(ggtree::ggtree(cTime.tree)$data)
  
  if (!all(eTime.tree$tip.lab==cTime.tree$tip.lab)){
    print("Fatal: input tree topologies are not identical")
    print("Output makes no sense")
  }
  
  cTime.tree.data$mu=eTime.tree.data$branch.length / cTime.tree.data$branch.length
  
  start_node=eTime.tree.data[eTime.tree.data$label==sp,"node"]
  start_node.x=eTime.tree.data[eTime.tree.data$label==sp,"x"]
  
  ancester_nodes=tidytree::ancestor(eTime.tree,start_node)
  ancester_nodes.x=eTime.tree.data[ancester_nodes,"x"]
  
  nodes.lst=c(start_node, ancester_nodes)
  x.lst=c(start_node.x,ancester_nodes.x)
  
  if (eTime>=x.lst[1]){
    print("Events before the root of the tree")
    x.cTime=NULL
    y.cTime=NULL
    x.eTime=NULL
    y.eTime=NULL
  }else{
    x.eTime=x.lst[1]-eTime
    i=which(x.eTime<=x.lst)[length(which(x.eTime<=x.lst))]
    most.recent.node2eTime=nodes.lst[i]
    
    cTime=max(cTime.tree.data[,"x"])-cTime.tree.data[most.recent.node2eTime,"x"]+
      +(eTime.tree.data[most.recent.node2eTime,"x"]-x.eTime)/cTime.tree.data[most.recent.node2eTime,"mu"]
    x.cTime=max(cTime.tree.data[,"x"])-cTime
    y.cTime=cTime.tree.data[most.recent.node2eTime,"y"]
    
    y.eTime=eTime.tree.data[most.recent.node2eTime,"y"]
  }
  
  res=list(sp=sp,eTime=eTime,cTime=cTime,
           coordinate.eTime.tree=c(x.eTime,y.eTime),
           coordinate.cTime.tree=c(x.cTime,y.cTime))
  return(res)
}

plotTree=function(nwk=nwk,outgroup=NA){
  if (file.exists(nwk)){
    UCE=ape::read.tree(file=nwk)
  }else{
    UCE=ape::read.tree(text=nwk)
  }
  
  if (!is.na(outgroup)){
    UCE=ape::root(UCE,outgroup = outgroup,resolve.root = T)
  }
  
  UCE.df=as.data.frame( ggtree(UCE,branch.length = "none")$data )
  nodeColor=UCE.df[!UCE.df$isTip & UCE.df$label!="Root" & UCE.df$label!="",]
  nodeColor$color=sapply(nodeColor$label,
                         function(i){
                           i=as.numeric(i)
                           if (i==100){return("blue")}
                           if (i<100 & i>=90){return("yellow")}
                           if (i<90){return("red")}
                         })
  UCE.f=ggtree(UCE,branch.length = "none")+
    geom_richtext(data=UCE.df[UCE.df$isTip,], 
                  aes(label=label,x=x,y=y),hjust=0,
                  fill = NA,label.colour=NA)+
    geom_point(data=nodeColor,
               aes(x=x,y=y,color=color))+
    scale_color_manual(breaks=c("blue","yellow","red"),
                       values=c("blue","yellow","red"))+
    
    scale_x_continuous(limits = c(0,40))+
    scale_y_continuous(limits = c(0,50))+
    labs(x=t)+
    theme_classic()+
    theme(text=element_text(face="bold",size=15),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          legend.position = "none")
  return(UCE.f)
}

compareP=function(left.tr=left.tr, # phylo
                  right.tr=right.tr){
  left.df=as.data.frame(ggtree(left.tr,branch.length = "none")$data)
  left.nodeColor=left.df[!left.df$isTip & left.df$label!="Root" & left.df$label!="",]
  left.nodeColor$color=sapply(left.nodeColor$label,
                              function(i){
                                i=as.numeric(i)
                                if (i==100){return("blue")}
                                if (i<100 & i>=90){return("yellow")}
                                if (i<90){return("red")}
                              })
  
  or=left.df[left.df$isTip,];or=or[order(or$y),"label"]
  right.tr=ape::rotateConstr(right.tr,constraint = or)
  right.df=as.data.frame(ggtree(right.tr,branch.length = "none")$data)
  right.df$x=max(right.df$x)-right.df$x+max(left.df$x)+80
  right.nodeColor=right.df[!right.df$isTip & right.df$label!="Root" & right.df$label!="",]
  right.nodeColor$color=sapply(right.nodeColor$label,
                               function(i){
                                 i=as.numeric(i)
                                 if (i==100){return("blue")}
                                 if (i<100 & i>=90){return("yellow")}
                                 if (i<90){return("red")}
                               })
  
  dd.lines=left.df
  dd.lines$x=dd.lines$x+40
  dd.lines=rbind(dd.lines,right.df[,-10]);dd.lines=dd.lines[dd.lines$isTip,]
  
  comparephylo.1=ggplot()+
    geom_tree(data=left.df,aes(x=x,y=y))+
    geom_tree(data=right.df,aes(x=x,y=y))+
    geom_line(data=dd.lines,color="red",
              aes(x=x,y=y,group=label),alpha=0.5)+
    geom_richtext(data=left.df[left.df$isTip,], 
                  aes(label=label,
                      x=x+40,y=y),hjust=1,
                  fill = NA,label.colour=NA)+
    
    geom_point(data=left.nodeColor,
               aes(x=x,y=y,color=color))+
    geom_point(data=right.nodeColor,
               aes(x=x,y=y,color=color))+
    scale_color_manual(breaks=c("blue","yellow","red"),
                       values=c("blue","yellow","red"))+
    labs(x=NULL,y=NULL,color="Bootstrap")+
    theme_classic()+
    theme(text=element_text(face="bold",size=15),
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          legend.position = "none")
  return(comparephylo.1)  # 10x8
}

#pgls
{
  tree.phylo=tree.phylo
  df=df
  sp.order=sp.order
  mod="y~x"
  
  fit1=nlme::gls(as.formula(mod),
                 correlation=ape::corBrownian(phy=tree.phylo,form=sp.order),
                 data=df)
  vf=diag(vcv(tree.phylo))
  fit2=nlme::gls(as.formula(mod),
                 correlation=ape::corBrownian(phy=tree.phylo,form=sp.order),
                 data=df,
                 weights=varFixed(~vf))
}

pANOVA=function(x,y,sp.name=sp.name,tree){
    res=phytools::phylANOVA(tree=tree,
                        x=setNames(x,sp.name),
                        y=setNames(y,sp.name),
                        nsim=5000)
    return(res)
}

pgls(tree.phylo=ape::drop.tip(sp.tree,"Blatta_orientalis"),
     df=paranomeSize[paranomeSize$sp!="Bori",],
     sp.order="~sp.name",
     mod= "rawParanomeSize~genome_size+ontogeny")


tree.phylo=ape::drop.tip(sp.tree,"Blatta_orientalis")
ultrametic=F
df=paranomeSize[paranomeSize$sp!="Bori",]
sp.order="sp.name"
weight=nlme::varFixed( ~ diag(ape::vcv(ape::drop.tip(sp.tree,"Blatta_orientalis"))))
mod= "rawParanomeSize~genome_size+ontogeny"
rm(tree.phylo,ultrametic,df,sp.order,mod)
# tree=ape::read.tree(text="(Bori:0.2563775834,Cmer:0.100536575,(Mdar:0.09943275989,((Znev:0.08732374973,Hsjo:0.04542242394):0.05942198095,((Kfla:0.04281909971,(PAsim:0.01774884375,(Gfus:0.1556399208,(Ncas:0.1534238245,((Rebo:0.07456697635,Mhub:0.03698058186):0.006444519778,(Cbre:0.05271272486,Isch:0.01782396019):0.006602280434):0.003024505055):0.005341389972):0.01786051148):0.01220000152):0.0883340121,(Shal:0.06029802325,((Gocu:0.09705845429,Dlon:0.06616272603):0.0138274079,(PRsim:0.05798306263,((Rfla:0.03948792356,(Hten:0.01102499093,(Cges:0.01273356582,Ctes:0.005907411764):0.009193531309):0.02717934078):0.008435693919,((Ssph:0.06263414772,(Aaca:0.08092428146,(Ofor:0.03095486448,Mnat:0.02460853358):1e-08):0.01360062693):0.003741808421,(Fval:0.05110853899,((Aunk:0.08783784523,(Eunk:0.03132652945,(Apac:0.01929907745,Aban:0.01394361536):0.07073307116):0.0140086283):0.02403511752,((Munk:0.03232572391,(Shey:0.02232661998,(Llab:0.01860741821,Cwal:0.2258704564):1e-08):0.01109750977):0.007776511373,(((Punk:0.03767270696,Abea:0.04013948611):0.003712408036,(Pred:0.04628058026,Iunk:0.03525013173):1e-08):0.002405730463,(Cpar:0.07418483895,(Ntar:0.05213310384,(Lunk:0.07921570275,((Nluj:0.03217572383,Hunk:0.02729033179):0.01359651125,(Ccav:0.03958015184,Csp4:0.04221816729):0.002472239263):0.01800280623):0.01744080149):1e-08):1e-08):0.003080828561):0.009307714003):0.01010370701):0.001705247777):0.02490181948):0.00316279518):0.002916154509):0.007133215165):0.06490499687):0.03024594662):0.009268027523):0.009382159741);")
# ggtree::ggtree(ape::unroot(phy = tree))
# ape::write.tree(ape::unroot(phy = tree))


# community=igraph::cluster_fast_greedy(tidygraph::as_tbl_graph(validPairs,directed = FALSE),
#                                       weights=1/validPairs$Ks)
# seq2community=igraph::membership(community)
# seq2community=data.frame(seq=names(seq2community),community=as.numeric(unname(seq2community)))
# rownames(seq2community)=seq2community$seq
# validPairs$community1=seq2community[validPairs$seq1,"community"]
# validPairs$community2=seq2community[validPairs$seq2,"community"]  

# karyotype=read.table("~/Desktop/PhD/Results/termite_pca/chromosomes/karyotype.tsv",
#                      sep="\t",header=TRUE,quote="")
# karyotype$N50_scaffoldAbove1M=sapply(karyotype$sp,
#                                      function(sp){
#                                        d=read.table(paste("~/Desktop/PhD/Results/termite_pca/chromosomes/scaffold_length/",
#                                                           sp,"_chrLength.tsv",sep=""),
#                                                     sep="\t",header=FALSE,quote="")
#                                        return(Biostrings::N50(d[d$V2>1e+6,2]))
#                                      })
# write.table(karyotype,
#             "~/Desktop/PhD/Results/termite_pca/chromosomes/karyotype.tsv",
#             sep="\t",row.names=FALSE,quote=FALSE)


