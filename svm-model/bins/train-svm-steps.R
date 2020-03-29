library(e1071);
library(ROCR);

options(width=1000);
binsmatrix<-as.matrix(read.table("../params_stats/results/breaks5bins.txt",row.names=1));
midsmatrix<-as.matrix(read.table("../params_stats/results/medians5bins.txt",row.names=1));
show(binsmatrix);

nbrsize<-1;
#for(nbrsize in seq(0,3))
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
fulldata<-read.table("results/patmatrix.txt");
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
        
                  #running on 544 rows, rows representing positions, columns are 5X12=60 created by populating 
                  #...the bins of ensemble by pooling data from all frame for each position
        #i<-1;
        for(i in c(1:nrow(fulldata)))          
        {
        mypat<-c();
                #nbr<--1;
                for(nbr in seq(-1*nbrsize,nbrsize))   #loop encoding the base sequence of window to sparse sequence
                {
                #show(seqnames);
                tmppat<-substr(seqnames[i],winpos+nbr,winpos+nbr);
                #show(c(i,allnames[i],seqnames[i],tmppat));
                tmppat<-encoding[tmppat,];
                mypat<-c(mypat,tmppat);
                }  
        allinputs<-rbind(allinputs,mypat);          #allinputs having 544 rows and their sparse sequence of selected window size
        }
        dim(allinputs);
        rownames(allinputs)<-allnames;      

fullpred<-c();
dynaseq<-list();

#binid<-1;                           #fulldata representing the patmatrix
for(binid in c(1:ncol(fulldata)))  #loop running on 5X12=60 columns representing bins of ensemble of 12 parameters except mgw 
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
 
        dynaseq[[binid]]<-model;
        pred<-predict(dynaseq[[binid]],tstdata);
        pred=as.matrix(pred);
        selfmae<-mean(abs(tsttarget-pred[,1]))
        sd1<-sd(tsttarget);
        sd2<-sd(pred[,1]);
show(c("self-training done. MAE=",selfmae));

### Cross validation results
nfold<-nrow(fulldata)/foldsize
#show(nfold);
myres<-c();
        
        #n<-1;
        for(n in c(1:nfold))     #lop running on 544/4 = 136
        {
        tetra<-tetras[n];
        show(paste("Window size",nbrsize*2+1,"BinID",binid,"Processing LOO by leaving ",tetra,"data out.",sep=" "));

        tstrange<- seq(foldsize*(n-1)+1,foldsize*(n)); #range of tetra
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

dynaseq[[binid+1]]<-binsmatrix;
dynaseq[[binid+2]]<-midsmatrix;
colnames(performance)<-c("binID","SD(observed)","SD(predicted)","MAE(self)","MAE(LOO)")
colnames(fullpred)<-c("binID","SeqPattern","Observed","Predicted")

write.table(fullpred,file=paste("results/cross-validation-full-ws",nbrsize*2+1,".txt",sep=""),quote=FALSE,row.names=FALSE);
write.table(performance,file=paste("results/performance-loo-ws",nbrsize*2+1,".txt",sep=""),quote=FALSE,row.names=FALSE);
save(dynaseq,file=paste("results/dynaseq-ws",nbrsize*2+1,".svm",sep=""));
}
