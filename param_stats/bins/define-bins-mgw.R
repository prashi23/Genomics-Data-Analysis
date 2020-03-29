

tetras<-as.matrix(read.table("list.tetra"));
tetras<-as.character(tetras[,1]);

##The data-all-MGW was created using an awk script (position 5 to 8 from column 1).
#awk '{if($1>=5 && $1<=8)print substr(FILENAME,13,4)"."$1-4,$4}' groove-data/*/groove.dat  > groove-data/data-all-MGW.txt

data.all<-as.matrix(read.table("groove-data/data-all-MGW.txt"));

#dim(data.all)

data.all<-as.numeric(data.all[,2])

##Now binning the data for each column

breaks.5<-quantile(data.all,probs=seq(0,1,by=0.2));
medians.5<-quantile(data.all,probs=seq(0.1,1,by=0.2));

breaks.3<-quantile(data.all,probs=c(0,0.25,0.75,1));
medians.3<-quantile(data.all,probs=c(0.125,0.5,0.875));

breaks.5<-signif(breaks.5,4);
breaks.3<-signif(breaks.3,4);

show(breaks.5);
show(medians.5);
show(breaks.3);
show(medians.3);

write.table(breaks.5,file="results/breaks5bins-mgw.txt",quote=F,col.names=F);
write.table(medians.5,file="results/medians5bins-mgw.txt",quote=F,col.names=F);

write.table(breaks.3,file="results/breaks3bins-mgw.txt",quote=F,col.names=F);
write.table(medians.3,file="results/medians3bins-mgw.txt",quote=F,col.names=F);
