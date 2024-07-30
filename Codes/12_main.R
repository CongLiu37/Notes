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
angle_of_vectors=function(x,y){
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


  
  
  
  