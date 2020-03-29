d10<-as.matrix(read.table("results/breaks5bins-10ns.txt",row.names=1));
d100<-as.matrix(read.table("results/breaks5bins.txt",header=T));

d10<-data.matrix(d10);
d100<-data.matrix(d100);

show(dim(d10));
show(dim(d100));

diff1<-100*abs(d10-d100)/max(abs(d10),abs(d100));
diff2<-100*abs(d10[,c(1,6)]-d100[,c(1,6)])/max(abs(d10[,c(1,6)]),abs(d100[,c(1,6)]));
diff3<-100*abs(d10[,2:5]-d100[,2:5])/max(abs(d10[,2:5]),abs(d100[,2:5]));

show(paste("difference between non-extreme breakpoints",mean(diff1)));
show(paste("difference between extreme breakpoints",mean(diff2)));
show(paste("difference between all breakpoints",mean(diff3)));
