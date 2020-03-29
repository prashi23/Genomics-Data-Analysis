#library("preprocessCore")
library("gplots")
options(width=10000);

positions<-c(5:8);
tetras<-read.table("../params_stats/list.tetra");
tetras<-tetras[,1];
databins<-as.matrix(read.table("../params_stats/results/breaks5bins-mgw.txt",row.names=1));
show(databins);

data.all<-c();
data.all<-as.matrix(read.table("../params_stats/results/data-all-all.txt",header=T,nrow=-1));
show(data.all[1:20,]);

##Number of bases to be considered in the +/- positions each
nbrsize=4; 

patmatrix<-c();
mypatnames<-c()
        
for(tetra in tetras)
{
show(tetra)
                tetrafull<-paste("CGCG",tetra,"CGCG",sep="");
param=c("mgw");
                for (pos in positions)
                {
                myname<-substr(tetrafull,pos-nbrsize,pos+nbrsize);
                myhistall<-c();
                name2=paste(tetra,pos-4,sep=".");
                #tmpdata<-data.matrix(data.all[data.all[,1]==name2,]);
                tmpdata<-as.numeric(data.all[data.all[,1]==name2,param]);
                mybreaks<-databins[,1];
                myhistall<-hist(tmpdata,breaks=mybreaks,plot=FALSE);
                myhistall<-myhistall$counts/length(tmpdata);
                mypatnames<-c(mypatnames,myname)
                patmatrix<-rbind(patmatrix,myhistall)
                }
}
rownames(patmatrix)<-mypatnames;
write.table(patmatrix,file="results/patmatrix-mgw.txt",quote=FALSE,col.names=FALSE);

