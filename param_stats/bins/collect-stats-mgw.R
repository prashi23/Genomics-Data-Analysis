library("gplots")
options(width=10000);

positions<-c(1:4);
tetras<-read.table("list.tetra");
tetras<-tetras[,1];

databins<-as.matrix(read.table("results/breaks5bins-mgw.txt",row.names=1));
data.all<-as.matrix(read.table("groove-data/data-all-MGW.txt",nrow=-1));

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



pdf("results/heatmaps-tetra-mgw.pdf");
        all_tetra_hist<-c();
        for(tetra in tetras)
        {
        all_pos_hist<-c();
                for (pos in positions)
                {
                selrow<-paste(tetra,pos,sep=".");
                tmpdata<-data.all[data.all[,1]%in%selrow,];
                tmpdata<-as.numeric(tmpdata[,2]);
                show(paste(tetra,pos));
                mybreaks<-databins[,1];
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
myfilename<-paste("results/tetra-heatmap-data-mgw.txt",sep="")
write.table(all_tetra_hist,file=myfilename,quote=FALSE);
heatmap.2(all_tetra_hist,main="MGW",col=mycolors,Colv=NA,scale=c("none"),density.info=c("none"),trace=c("none"),cexRow=0.3,cexCol=0.5,keysize=c(0.85),lhei=c(1,5),xlab=c("(base position).(ensemble bin ID)"),mar=c(4,15),colsep=seq(1,ncol(all_tetra_hist),1),sepcol=mysepcol,sepwidth=c(0.02,0.02),rowsep=c(1:nrow(all_tetra_hist)));

warnings()
dev.off()
