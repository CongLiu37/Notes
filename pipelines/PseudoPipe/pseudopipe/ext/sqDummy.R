arg=commandArgs(trailingOnly = TRUE)
threads=as.numeric(arg[1])
cmd="sed -i 's/\"//g' jobs"
print(cmd);system(cmd,wait=TRUE)
cmd=paste("split -d",
          "-n",paste("l/",as.character(threads-1),sep=""),
          "jobs",
          sep=" ")
print(cmd);system(cmd,wait=TRUE)
library(parallel)
clus=makeCluster(threads)
parSapply(clus,
          1:(threads-1),
          function(i){
            scr=system("ls x*",wait=TRUE,intern=TRUE)[i]
            cmd=paste("bash"," ",scr,sep="")
            print(cmd);system(cmd,wait=TRUE)})
stopCluster(clus)
