# Rarefaction analysis
# Dependencies:

SamplingEffortSimulation=function(AbundanceTable=AbundanceTable, # A dataframe with column "Name" and "Abundance".
                                  Time=Time, # How many times a reduced sampling effort is simulated.
                                  Effort=Effort, # Reduced sampling effort. It can be either percent of full effort (e.g. 0.50) or an integer smaller than fuu effort (e.g. 50).
                                  Remove=Remove, # An integer vector indicating which rows of AbundanceTable are removed before computing diversity/dissimilarity. 
                                                 # These rows can represent host contamination or unannonated metagenomic reads.
                                                 # Invalid by set to 0.
                                  Threshold=Threshold # Threshold for removing rare categories.
                                                      # percent or integer.
                                  ){
  True=AbundanceTable[,"Abundance"]
  if(Remove!=0){True=AbundanceTable[,"Abundance"][-Remove]} # Remove rows such as host contamination, etc.
  if (Threshold<1){Threshold=Threshold*Effort} # Threshold for absolute abundance.
  True[True<Threshold]=0 # Remove categories that are too rare to be true
  f1=function(vector){return(vector/sum(vector))} # Absolute abundances to relative abundances
  True=f1(True) # Relative abundance distribution from full effort.
  
  prob=AbundanceTable[,"Abundance"]/sum(AbundanceTable[,"Abundance"]) 
  if (Effort<1){Effort=Effort*sum(AbundanceTable[,"Abundance"])}
  mat=rmultinom(Time,Effort,prob)
  rownames(mat)=AbundanceTable[,"Name"] # Absolute abundance distribution from reduced effort.
  
  if(Remove!=0){mat=mat[-Remove,]} # Remove rows such as host contamination, etc.
  
  mat[mat<Threshold]=0 # Remove categories that are too rare to be true
  mat=mat+1e-3 # add a small constant to all the elements of the abundance matrix.
  mat=apply(mat,2,f1) # # Relative abundance distribution from reduced effort.
  
  o=list(mat=mat,ReducedEffort=Effort,True=True,FullEffort=sum(AbundanceTable[,"Abundance"]))
  return(o)
}

ComputeDiversity=function(input=input, # Value returned by SamplingEffortSimulation
                          Index=Index # Richness(R)/ExponentialShannon(ES)/GiniSimpson(GS)
                          ){
  if (Index=="Richness"|Index=="R"){Index="Richness";index=function(vector){return(length(vector[vector!=0]))}}
  if (Index=="ExponentialShannon"|Index=="ES"){Index="ExponentialShannon";index=function(vector){vector=vector[vector!=0];return(-sum(vector*log(vector)))}}
  if (Index=="GiniSimpson"|Index=="GS"){Index="GiniSimpson";index=function(vector){vector=vector[vector!=0];return(1/sum(vector^2))}}
  diversity=apply(input$mat,2,index)
  o=list(Index=Index,index=diversity,ReducedEffort=input$ReducedEffort,True=index(input$True),FullEffort=input$FullEffort)
  return(o)
}

ComputeDistance=function(input, # Value returned by SamplingEffortSimulation
                         Index=Index # Bray-Curtis(BC)/Kullback-Leibler(KL)/Jensen-Shannon(JS)
                         ){
  if (Index=="Bray-Curtis"|Index=="BC"){Index="Bray-Curtis";index=function(vector,True){return(sum(abs(vector-True))/sum(vector+True))}}
  if (Index=="Kullback-Leibler"|Index=="KL"){Index="Kullback-Leibler";index=function(vector,True){return(sum(True*log(True/vector)) )}}
  if (Index=="Jensen-Shannon"|Index=="JS"){Index="Jensen-Shannon";index=function(vector,True){return(0.5*sum(True*log(True/(0.5*(vector+True))))+0.5*sum(vector*log(vector/(0.5*(vector+True)))) )}}
  distance=apply(input$mat,2,index,True=input$True)
  o=list(Index=Index,index=distance,ReducedEffort=input$ReducedEffort,True=index(input$True,input$True),FullEffort=input$FullEffort)
  return(o)
}

CI=function(input=input # Value returned by ComputeDiversity/ComputeDistance
            ){
    ComputeDelta_star=function(vector,i){
      BootSample = sample(vector,length(vector),replace=TRUE)
      Delta_star = mean(BootSample)-mean(vector)
      return(Delta_star)
    }
    vector=input$index
    Delta_stars = rep(NA,1000) # 1000 bootstrap
    Delta_stars = sapply(Delta_stars,ComputeDelta_star,vector=vector)
    CI_upper = mean(vector)-quantile(Delta_stars,0.025) # .95 conficence
    CI_lower = mean(vector)-quantile(Delta_stars,0.975)
    o=list(Index=input$Index,index=input$index,mean_value=mean(input$index),CI=c(CI_lower,CI_upper),ReducedEffort=input$ReducedEffort,True=input$True,FullEffort=input$FullEffort)
    return(o)
}
