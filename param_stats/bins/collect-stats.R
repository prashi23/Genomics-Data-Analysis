library("gplots")
options(width=10000);

### **** This file was edited by mistake after use. May need corrections if reuse is intended.
#### **********************************************


positions<-c(1:4);
tetras<-read.table("list.tetra");
tetras<-tetras[,1];

param.names<-as.matrix(read.table("list.params"));
param.names<-as.character(param.names[,1]);

### Read and process step data

data.all<-c();
myrownames<-c();
#id=1;

#for (id in seq(1,length(tetras)))
for (id in seq(1,2))
{
tetra<-tetras[id];
show(paste("Reading ",tetra,"data (",id,"of",length(tetras),")"));
file=paste("steps_data/",tetra,"/bp_step_all.par",sep="");
data.tetra<-as.matrix(read.table(file,skip=3,fill=T,header=F,nrow=-1));
#sequence<-as.character(substr(data.tetra[seq(1,12),1],1,1));
#dim(data.tetra)

skiprows<-unlist(lapply(seq(15,nrow(data.tetra),by=15),function(x){return(seq(max(x-6,1),min(x+4,nrow(data.tetra))))}));
skiprows<-c(seq(1,4),skiprows,seq(nrow(data.tetra)-3,nrow(data.tetra)));


sequence<-as.character(substr(data.tetra[seq(1,12),1],1,1));

data.tetra<-data.tetra[-skiprows,seq(2,length(param.names)+1)];  #param.names not found
#data.tetra<-data.tetra[-skiprows,seq(2,13)];

colnames(data.tetra)<-param.names;

data.all<-rbind(data.all,data.tetra)
#nposes=nrow(data.tetra)/4; #surender

#Next line of code surender
nposes=nrow(data.tetra)/4;

myrownames<-c(myrownames,paste(tetra,rep(seq(1,4),nposes),sep="."));
}

data.all<-cbind(myrownames,data.all);

show(data.all[1:20,]);
#### Steps data read

mycolors=rainbow(n=50,start=0,end=5/6)
mysepcol<-c(rep(c("white"),4),c("black"));
mysepcol<- rep(mysepcol,4);

        bin.ids=seq(1,5);
        mycolnames<-c();
        
        for(i in positions) 
        { 
                for(j in bin.ids) 
                {
                mycolnames<-c(mycolnames,paste(i,".",j,sep=""))
                }
        }



pdf("results/heatmaps-tetra.pdf");

#for (param in params)
params=param.names;

param="shear";

#for (param in params)                   #"Shear"   "Stretch" "Stagger" "Buckle" etc
{
        all_tetra_hist<-c();
        
        tetra="AAAA";
        #for(tetra in tetras)            # AAAA AAAC AAAG AAAT AACA etc
        {
        all_pos_hist<-c();
  
                
                pos=1;
                #for (pos in positions)  # 1 2 3 4           
                  
                                        #param tetra pos
                {
                selrow<-paste(tetra,pos,sep=".");
                #tmpdata<-data.all[data.all[,1]%in%selrow,];      # #inserted modifications surender
          
                tmpdata<-data.all[data.all[,1]%in%selrow,];
                
                #surender<-seq(1,10);
                #surenders<-surender[surender%in%5];
                
                #tmpdata<-as.numeric(tmpdata[,param]);
                tmpdata<-as.numeric(tmpdata[,2]);
                
                
                show(paste(tetra,pos,param));
                mybreaks<-databins[param,];
                dim(mybreaks)<-NULL;
                myhist<-hist(tmpdata,breaks=mybreaks,plot=FALSE);
                myhist<-myhist$counts/length(tmpdata);
                all_pos_hist<-c(all_pos_hist,myhist)
                
                
                }
        all_tetra_hist<-rbind(all_tetra_hist,all_pos_hist);
        }
rownames(all_tetra_hist)<-tetras;
#show(dim(all_tetra_hist));
#show(all_tetra_hist);

colnames(all_tetra_hist)<-mycolnames;
myfilename<-paste("results/tetra-heatmap-data-",param,".txt",sep="")
write.table(all_tetra_hist,file=myfilename,quote=FALSE);
heatmap.2(all_tetra_hist,main=param,col=mycolors,Colv=NA,scale=c("none"),density.info=c("none"),trace=c("none"),cexRow=0.3,cexCol=0.5,keysize=c(0.85),lhei=c(1,5),xlab=c("(base position).(ensemble bin ID)"),mar=c(4,15),colsep=seq(1,ncol(all_tetra_hist),1),sepcol=mysepcol,sepwidth=c(0.02,0.02),rowsep=c(1:nrow(all_tetra_hist)));
}

warnings()
dev.off()
