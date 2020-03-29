#library("preprocessCore")
library("gplots")
options(width=10000);

positions<-c(5:8);
tetras<-read.table("../params_stats/list.tetra");
tetras<-tetras[,1];
databins<-as.matrix(read.table("../params_stats/results/breaks5bins.txt",row.names=1));
params<-rownames(databins);
show(databins);
param.names<-as.matrix(read.table("../params_stats/list.params"));
param.names<-as.character(param.names[,1]);

#### Read params data 
data.all<-c();
myrownames<-c();

#id<-5;
#for (id in seq(1,length(tetras)))
for (id in seq(1,5))
{
tetra<-tetras[id];
show(paste("Reading ",tetra,"data (",id,"of",length(tetras),")"));
file=paste("../params_stats/steps_data/",tetra,"/bp_step_all.par",sep="");
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
data.all<-cbind(myrownames,data.all);

show(data.all[1:20,]);

### Steps data read finished

##Number of bases to be considered in the +/- positions each
nbrsize=4; 

patmatrix<-c();
mypatnames<-c()

#tetra<-"AAAA";        
for(tetra in tetras)
{
show(tetra)
                tetrafull<-paste("CGCG",tetra,"CGCG",sep="");

                #pos<-5;
                for (pos in positions)
                {
                myname<-substr(tetrafull,pos-nbrsize,pos+nbrsize);
                myhistall<-c();
                
                        
                        for (param in params)
                        {
                          name2=paste(tetra,pos-4,sep=".");
                          #tmpdata<-data.matrix(data.all[data.all[,1]==name2,]);
                          tmpdata<-as.numeric(data.all[data.all[,1]==name2,param]);
                          mybreaks<-databins[param,];
                          myhist<-hist(tmpdata,breaks=mybreaks,plot=FALSE);
                          myhist<-myhist$counts/length(tmpdata);
                          myhistall<-c(myhistall,myhist);
                          #show(c(name2,myhistall));
                        }
                mypatnames<-c(mypatnames,myname)
                patmatrix<-rbind(patmatrix,myhistall)
                }
}
rownames(patmatrix)<-mypatnames;
write.table(patmatrix,file="results/patmatrix.txt",quote=FALSE,col.names=FALSE);

