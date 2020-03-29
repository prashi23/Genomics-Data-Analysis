
pdf("results/def5bins-byR.pdf");

data<-read.table("results/breaks5bins.txt",row.names=1);
show(data);
colnames(data)<-seq(1,6);
data<-data.matrix(data);
data2<-read.table("results/breaks5bins-mgw.txt",row.names=1);
dtaa<-rbind(data,t(data2));

barsize<-data[,1:5]*0;
midvalue<-data[,1:5]*0;
labelpos<-data[,1:5]*0;

for(i in seq(1,12))
{
        for(j in seq(1,5))
        {
        barsize[i,j]=abs(data[i,j+1]-data[i,j])
        midvalue[i,j]=c((data[i,c(j+1)]+data[i,j])/2)
        }
barsize[i,]<-barsize[i,]/sum(barsize[i,]);
}

        labelpos[,1]=barsize[,1]/2;

        for(i in seq(2,5))
        {
                for(j in seq(1,12))    
                {
                labelpos[j,i]=sum(barsize[j,seq(1,i-1)])+barsize[j,i]/2;
                }
        }
show(labelpos);
show(barsize);

b<-barplot(t(barsize),horiz=T,beside=F,col=c("lightgreen","lightgrey","yellow","cyan","orange"),las=2,xlab=c("Bin-wise coverage"),cex.lab=1.5,cex.names=0.7)
lapply(seq(1,5),function(x){text(labelpos[,x],b,signif(midvalue[,x],3),cex=0.5,srt=90);return(0)})

show(barsize);

dev.off();
