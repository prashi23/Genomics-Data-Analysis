d10<-as.matrix(read.table("patsmatrix-10ns.txt",header=F));
d100<-as.matrix(read.table("results/patmatrix.txt",header=F));
mynames<-d10[,1];
show(mynames);

d10<-data.matrix(d10[,2:ncol(d10)]);
d100<-data.matrix(d100[,2:ncol(d100)]);

show(dim(d10));
show(dim(d100));
d10.all<-as.numeric(d10);
d100.all<-as.numeric(d100);
ext<-c(seq(1,60,by=5),seq(5,60,by=5));

d10.nonext<-as.numeric(d10[,-ext]);
d100.nonext<-as.numeric(d100[,-ext]);

diff1<-100*abs(d10.all-d100.all)/max(abs(d10.all),abs(d100.all));
diff2<-100*abs(d10.nonext-d100.nonext)/max(abs(d10.nonext),abs(d100.nonext));

show(paste("difference between all pattarens",mean(diff1)));
show(paste("difference between patterns with non-extreme breakpoints",mean(diff2)));
