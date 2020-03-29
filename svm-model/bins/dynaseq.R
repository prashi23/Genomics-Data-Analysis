


library(e1071);
options(width=1000);

sequences<-c("AAGTCGTAATG","CGTGAGCAATC","CCCGGAGTCGGTC");

nbrsize=2

encoding<-c(0,0,0,1)
encoding<-rbind(encoding,c(0,0,1,0))
encoding<-rbind(encoding,c(0,1,0,0))
encoding<-rbind(encoding,c(1,0,0,0))
rownames(encoding)<-c("A","C","G","T")

load("ws5/dynaseq-ws5.svm");

allinputs<-c();
mynames<-c();


for (seqid in c(1:length(sequences))) 
{
sequence<-sequences[seqid];

        for(winpos in c(4:6))
        {
        mypat<-c();
        mynames<-c(mynames,paste("[sequence.",seqid,".",sequence,".pos.",winpos,"]",sep=""))
                for(nbr in seq(-1*nbrsize,nbrsize))
                {
                tmppat<-substr(sequence,winpos+nbr,winpos+nbr);
                tmppat<-encoding[tmppat,];
                mypat<-c(mypat,tmppat);
                }
        allinputs<-rbind(allinputs,mypat);
        }
}
show(dim(allinputs));

        params<-rownames(dynaseq[[61]]);
        pheader=c();
        for(i in c(1:nrow(dynaseq[[61]])))
        {
                for(j in c(1:(ncol(dynaseq[[61]])-1)))
                {
        pheader<-c(pheader,paste(params[i],"-bin-",j,sep=""));
                }
        }
#        show(pheader);

        res<-c();
        for(binid in c(1:60))
        {
        pred<-predict(dynaseq[[binid]],allinputs);
        res<-rbind(res,pred);
        }
        rownames(res)<-pheader;
        colnames(res)<-mynames;
        show(t(res));
