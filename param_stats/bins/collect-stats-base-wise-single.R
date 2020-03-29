("preprocessCore")
library("gplots")
options(width=10000);

positions<-c(1:4);
tetras<-read.table("list.tetra");
tetras<-tetras[,1];

databins<-as.matrix(read.table("results/breaks5bins.txt",row.names=1));
params<-rownames(databins);
show(databins);

param.names<-as.matrix(read.table("list.params"));
param.names<-as.character(param.names[,1]);

mycolors=rainbow(n=50,start=0,end=4/6)
mysepcol<-c(rep(c("white"),4),c("black"));
mysepcol<- rep(mysepcol,12);

mycolors=mycolors[seq(50,1,-1)];


        mycolnames<-c();
        mybreaks<-databins;
        rownames(mybreaks)<-rownames(databins);
        show(mybreaks);

        for(param in param.names) 
        { 
                for(j in c("P20","P40","P60","P80","P100")) 
                {
                mycolnames<-c(mycolnames,paste(j,".",param,sep=""))
                }
        }

### Read and process step data

data.all<-c();
myrownames<-c();

id=1;
#for (id in seq(1,length(tetras)))
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
data.all<-cbind(myrownames,data.all);

show(data.all[1:20,]);
#### Steps data read

bases<-c("A","C","G","T");

pdf("results/heatmaps-single.pdf");
summary_matrix<-c();
myrownames<-c();

base="A";
#for (base in bases)
{
        base_param_hist<-c()
        
        param="shear";
        #for(param in params)
        {
                mybasedata<-c();
                
                tetra="AAAA";
                #for(tetra in tetras)
                {
                tetrafull<-paste("CGCG",tetra,"CGCG",sep="")
                
                        pos=1;
                        #for(pos in c(1:4))
                        {
                                
                                if(substr(tetrafull,pos+4,pos+4)==base)
                                {
                                  
                                  #send all those rows when this tetra and base is matching eg. AAAA.1
                                  tmpdata<-data.all[data.all[,1]==paste(tetra,pos,sep="."),]
                                  colnames(tmpdata)<-c("base",params)
                                  show(c(base,param))
                                  tmpdata<-as.numeric(data.all[data.all[,1]==paste(tetra,pos,sep="."),param])
                                  mybasedata<-c(mybasedata,tmpdata);
                                }
                        }
                }
              if(length(mybasedata)>0)
              {
              show(length(mybasedata));
              myhist<-hist(mybasedata,breaks=mybreaks[param,],plot=FALSE);
              myhist<-myhist$counts/length(mybasedata);
              base_param_hist<-c(base_param_hist,myhist)
              }
        }
                if(length(base_param_hist)>0)
                {
        myrownames<-c(myrownames,base);
        summary_matrix<-rbind(summary_matrix,base_param_hist);
                }
}

show(dim(summary_matrix));
show(length(myrownames));
rownames(summary_matrix)<-myrownames;
colnames(summary_matrix)<-mycolnames;

par(cex.main=0.8);
heatmap.2(summary_matrix,main=c("Ensemble population specificity"),col=mycolors,Colv=NA,Rowv=NA,scale=c("none"),density.info=c("none"),trace=c("none"),cexRow=0.8,cexCol=0.5,keysize=c(1.0),lhei=c(1,5),mar=c(4,5),colsep=seq(1,60,1),sepcol=mysepcol,sepwidth=c(0.02,0.02),rowsep=seq(1,nrow(summary_matrix),1));
heatmap.2(summary_matrix,main=c("Population specificity clusters"),col=mycolors,Colv=NA,scale=c("none"),density.info=c("none"),trace=c("none"),cexRow=0.8,cexCol=0.5,keysize=c(1.0),lhei=c(1,5),mar=c(4,5),colsep=seq(1,60,1),sepcol=mysepcol,sepwidth=c(0.02,0.02),rowsep=seq(1,nrow(summary_matrix),1));
write.table(summary_matrix,file=c("results/stats-matrix-single.txt"),quote=FALSE);
dev.off()
