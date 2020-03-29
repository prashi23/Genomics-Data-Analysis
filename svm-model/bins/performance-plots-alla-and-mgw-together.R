#library("preprocessCore")
library("gplots")
library("graphics")
options(width=10000);

pdf("results/performance-plots-all-and-mgw.pdf");

data<-read.table("results/performance-ws.txt",header=TRUE,row.names=1)
fulldata<-read.table("results/cross-validation-full-ws5.txt",header=TRUE)

data.mgw<-read.table("results/performance-ws-mgw.txt",header=TRUE,row.names=1)
fulldata.mgw<-read.table("results/cross-validation-full-ws5-mgw.txt",header=TRUE)
fulldata.mgw[,1]=fulldata.mgw[,1]+60;
fulldata<-rbind(fulldata,fulldata.mgw)
data<-rbind(data,data.mgw);

params<-c("shear", "stretch", "stagger", "buckle", "prop_tw", "opening", "shift", "slide", "rise", "tilt", "roll", "twist","MGW");
pnames<-c();
for (param  in params)
{
pnames<-c(pnames,c(param,rep("",4)))
}

data<-data.matrix(data);

mycolors=rainbow(n=50,start=0,end=4/6)
mycolors=mycolors[length(mycolors):1]
mysepcol<-c("white");
mysepcol<- rep(mysepcol,13);

par(cex.main=0.8);
#heatmap.2(data[,1:5],main=c("Window size and MAE"),col=mycolors,Colv=NA,Rowv=NA,scale=c("row"),density.info=c("none"),trace=c("none"),cexRow=0.6,cexCol=1.0,keysize=c(1.0),lhei=c(1,5),mar=c(5,10),colsep=seq(1,5,1),sepcol=mysepcol,rowsep=c(1:nrow(data)),sepwidth=c(0.02,0.02));

mymeans<-c();
mysds<-c();

for (i in c(1:ncol(data)))
{
mymeans<-c(mymeans,mean(data[,i]))
mysds<-c(mysds,sd(data[,i]))
}
show(mymeans);
show(mysds);

mynames<-c("SD","1","3","5","7");
par(mfrow=c(2,2));

barplot(mymeans,names=mynames,ylim=c(0.0,0.05),ylab="Cross-validation average MAE",cex.axis=1.2,col=c("red",rep("blue",4)),las=3,main=c("(a) Model input window"),xlab=c("Window size(#bases)"));

fullmae<-abs(fulldata[,3]-fulldata[,4]);
myhist<-hist(fullmae,plot=FALSE,breaks=50,xlim=c(0,0.2));

cumdata<-cumsum(myhist$counts/sum(myhist$counts));

plot( myhist$mids, cumdata, xlab=c("absolute error"),ylab=c("Cumulative frequency"),main=c("(b) Abs error in all cross-validation predictions (WS=5)"),type="b",pch=c(19));
text(0.04+myhist$mids[1:5],cumdata[1:5],labels=signif(cumdata[1:5],2),cex=0.7);




mycor=cor(fulldata[,3],fulldata[,4]);
mytext<-paste("R = ",signif(mycor,2),sep="")
#plot(fulldata[,3],fulldata[,4],xlab=c("Observed bin-wise occupancy"),ylab=c("Predicted occupancy"),main=c("(c) Predicted vs obsereved all individual values"));
#legend("topleft",mytext,bty="n",cex=1.0,text.col=c("red"));

smoothScatter(fulldata[,3],fulldata[,4],xlab=c("Observed bin-wise occupancy"),ylab=c("Predicted occupancy"),main=c("(c) Predicted vs obsereved all individual values (WS=5)"));
legend("topleft",mytext,bty="n",cex=1.0,text.col=c("red"));

mycolors<-rep( c(rep("blue",5),rep("red",5)),6)
rownames(data)<-pnames;
barplot(data[,5],main="(d) Parameter-wise mean absolute errors (WS=5)",las=1,cex.names=0.7,col=mycolors,las=3,ylab="Mean absolute error",xlab=c("Parameter bins"));



dev.off()
