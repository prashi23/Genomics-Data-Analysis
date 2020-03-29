#library("preprocessCore")
library("gplots")
options(width=10000);

positions<-c(1:4);
tetras<-read.table("list.tetra");
tetras<-tetras[,1];

param.names<-as.matrix(read.table("list.params"));
param.names<-as.character(param.names[,1]);

databins<-as.matrix(read.table("results/breaks5bins.txt",row.names=1));
databins2<-as.matrix(read.table("results/breaks5bins-mgw.txt",row.names=1));
params<-rownames(databins);
databins<-rbind(databins,t(databins2[,1]));
rownames(databins)<-c(params,"mgw");
show(databins);

### Read all the parameters data from step files

data.all<-c();
myrownames<-c();

for (id in seq(1,length(tetras)))
{
tetra<-tetras[id];
show(paste("Reading ",tetra,"data (",id,"of",length(tetras),")"));
file=paste("steps_data/",tetra,"/bp_step_all.par",sep="");
data.tetra<-as.matrix(read.table(file,skip=3,fill=T,header=F,nrow=-1));
skiprows<-unlist(lapply(seq(15,nrow(data.tetra),by=15),function(x){return(seq(max(x-6,1),min(x+4,nrow(data.tetra))))}));
skiprows<-c(seq(1,4),skiprows,seq(nrow(data.tetra)-3,nrow(data.tetra)));
sequence<-as.character(substr(data.tetra[seq(1,12),1],1,1));
data.tetra<-data.tetra[-skiprows,seq(2,length(param.names)+1)];
colnames(data.tetra)<-param.names;
data.all<-rbind(data.all,data.tetra)
nposes=nrow(data.tetra)/4;
myrownames<-c(myrownames,paste(tetra,rep(seq(1,4),nposes),sep="."));
}

data.mgw<-as.matrix(read.table("groove-data/data-all-MGW.txt",nrow=nrow(data.all)));
data.all<-cbind(myrownames,data.all,data.mgw[,2]);


colnames(data.all)<-c("myrownames",params,"mgw")
write.table(data.all,"results/data-all-all.txt",quote=F);
params<-colnames(data.all[,2:14]);

#### Steps data read

        sdmatrix<-c();
for (param in params)
{
show(param);
        param_hist_all<-c();
        for(tetra in tetras)
        {
                for (pos in positions)
                {
                mypat=paste(tetra,pos,sep=".");
                tmpdata<-data.all[data.all[,1]==mypat,param];
                mybreaks<-databins[param,];
                myhist<-hist(as.numeric(tmpdata),breaks=mybreaks,plot=FALSE);
                myhist<-myhist$counts/length(tmpdata);
                param_hist_all<-rbind(param_hist_all,myhist)
                }
        }
show(dim(param_hist_all));
        param_mysdall<- c();
        for(i in c(1:ncol(param_hist_all)))
        {
        mysd<-sd(param_hist_all[,i])
        param_mysdall<- c(param_mysdall,mysd);
        }
        sdmatrix<-rbind(sdmatrix, param_mysdall);
show(sdmatrix);
}

rownames(sdmatrix)<-params;

write.table(sdmatrix,file="results/sdmatrix.txt",quote=FALSE);
