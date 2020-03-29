library(e1071);
library(ROCR);

options(width=1000);
binsmatrix<-as.matrix(read.table("../params_stats/results/breaks5bins-mgw.txt",row.names=1));
midsmatrix<-as.matrix(read.table("../params_stats/results/medians5bins-mgw.txt",row.names=1));
binsmatrix<-as.numeric(binsmatrix[,1]);
midsmatrix<-as.numeric(midsmatrix[,1]);

show(binsmatrix);

for(nbrsize in seq(0,3))
{

##Position of the base where window is located
winpos=5
##Number of patterns which come from the same trajectory (size of a single fold in LOO training)
foldsize<-4;

encoding<-c(0,0,0,1)
encoding<-rbind(encoding,c(0,0,1,0))
encoding<-rbind(encoding,c(0,1,0,0))
encoding<-rbind(encoding,c(1,0,0,0))
rownames(encoding)<-c("A","C","G","T")
#show(encoding);



options(digits=6)
fulldata<-read.table("results/patmatrix-mgw.txt");
seqnames<-as.character(fulldata[,1]);
allnames<-paste(fulldata[,1],"-",c(1:nrow(fulldata)),sep="")
fulldata<-data.matrix(fulldata[,2:ncol(fulldata)]);
allnames<-as.matrix(allnames);
rownames(fulldata)<-allnames;
#show(fulldata);

tetras<-as.matrix(read.table("../params_stats/list.tetra"));
tetras<-tetras[,1];
dim(tetras)<-NULL

allinputs<-c();
performance<-c();
        
        for(i in c(1:nrow(fulldata)))
        {
        mypat<-c();
                for(nbr in seq(-1*nbrsize,nbrsize))
                {
                tmppat<-substr(seqnames[i],winpos+nbr,winpos+nbr);
                #show(c(i,allnames[i],seqnames[i],tmppat));
                tmppat<-encoding[tmppat,];
                mypat<-c(mypat,tmppat);
                }  
        allinputs<-rbind(allinputs,mypat);
        }
        rownames(allinputs)<-allnames;      

fullpred<-c();
dynaseq.mgw<-list();

for(binid in c(1:ncol(fulldata)))
{
show(c("starting bin id ",binid));
        ### Self training results
        trdata<-allinputs;
        trtarget<-fulldata[,binid];

        tstdata<-trdata;
        tsttarget<-trtarget;

show(dim(trdata));
show(length(trtarget));
#show(trdata);
        model<-svm(trdata,as.numeric(trtarget),kernel=c("radial"),type=c("eps"));
        summary(model);
 
        dynaseq.mgw[[binid]]<-model;
        pred<-predict(dynaseq.mgw[[binid]],tstdata);
        pred=as.matrix(pred);
        selfmae<-mean(abs(tsttarget-pred[,1]))
        sd1<-sd(tsttarget);
        sd2<-sd(pred[,1]);
show(c("self-training done. MAE=",selfmae));

### Cross validation results
nfold<-nrow(fulldata)/foldsize
myres<-c();

        for(n in c(1:nfold))
        {
        tetra<-tetras[n];
        show(paste("Window size",nbrsize*2+1,"BinID",binid,"Processing LOO by leaving ",tetra,"data out.",sep=" "));

        tstrange<- seq(foldsize*(n-1)+1,foldsize*(n));
        trrange<- setdiff(c(1:nrow(fulldata)),tstrange);

        trdata<-allinputs[trrange,];
        trtarget<-fulldata[trrange,binid];
        trnames<-allnames[trrange];


        tstdata<-allinputs[tstrange,];
        tsttarget<-fulldata[tstrange,binid];
        tstnames<-allnames[tstrange];

        model<-svm(trdata,trtarget,kernel=c("radial"),type=c("eps"));
        pred<-predict(model,tstdata);
        dim(pred)<-NULL
        pred=as.matrix(pred);

        fullpred<-rbind(fullpred,cbind(binid,allnames[tstrange,1],tsttarget,pred));
        foldmae<-mean(abs(tsttarget-pred))
        myres<-c(myres,foldmae);
        }
        mae<-mean(myres);
        performance<-rbind(performance,c(binid,sd1,sd2,selfmae,mae))
}

dynaseq.mgw[[binid+1]]<-binsmatrix;
dynaseq.mgw[[binid+2]]<-midsmatrix;
colnames(performance)<-c("binID","SD(observed)","SD(predicted)","MAE(self)","MAE(LOO)")
colnames(fullpred)<-c("binID","SeqPattern","Observed","Predicted")

write.table(fullpred,file=paste("results/cross-validation-full-ws",nbrsize*2+1,"-mgw.txt",sep=""),quote=FALSE,row.names=FALSE);
write.table(performance,file=paste("results/performance-loo-ws",nbrsize*2+1,"-mgw.txt",sep=""),quote=FALSE,row.names=FALSE);
save(dynaseq.mgw,file=paste("results/dynaseq-ws",nbrsize*2+1,"-mgw.svm",sep=""));
}
