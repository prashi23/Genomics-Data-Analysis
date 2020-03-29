

tetras<-as.matrix(read.table("list.tetra"));
tetras<-as.character(tetras[,1]);

param.names<-as.matrix(read.table("list.params"));
param.names<-as.character(param.names[,1]);

data.all<-c();#creating an empty vector

for (id in seq(1,length(tetras)))
{
  tetra<-tetras[id];
  show(paste("Reading ",tetra,"data (",id,"of",length(tetras),")"));   
  #displaying this sentence on screen "Reading  AAAA data ( 1 of 136 )"
  
  file=paste("steps_data/",tetra,"/bp_step_all.par",sep="");           
  #it will save the appropriate "path" in file such as "steps_data/AAAT/bp_step_all.par"
  
  # it will 
  data.tetra<-as.matrix(read.table(file,skip=3,fill=T,header=F,nrow=-1));
  
  
  
  #show(data.tetra[1:30,1:5]);
  skiprows<-unlist(lapply(seq(15,nrow(data.tetra),by=15),function(x){return(seq(max(x-6,1),min(x+4,nrow(data.tetra))))}));
  #show(skiprows);
  skiprows<-c(seq(1,4),skiprows,seq(nrow(data.tetra)-3,nrow(data.tetra)));
  #show(skiprows);
  #It will create a list all rows except rows of central tetramers
  
  sequence<-as.character(substr(data.tetra[seq(1,12),1],1,1));
  #show(data.tetra[1:30,1:5]);
  #generating the sequence of the 12mer  e.g CGCGAAAACGCG
  
  #data.tetra<-data.tetra[-skiprows,seq(2,length(param.names)+1)];
  data.tetra<-data.tetra[-skiprows,];
  #deletes all the "skiprows" i.e deletes all the rows except central tetramer
  
  
  
  
  data.tetra2<-c();        #creates an empty vector named data.tetra2
  
  for(ss in seq(2,length(param.names)+1))
  {
    data.tetra2<-cbind(data.tetra2,as.numeric(data.tetra[,ss]))
  }
  
  data.tetra<-data.tetra2;
  
  
  
  
  colnames(data.tetra)<-param.names;
  #mytails<-seq(nrow(data.tetra)-12,nrow(data.tetra));
  #show(data.tetra[mytails,]);
  #show(data.tetra[mytails,1:5]);
  
  
  
  data.all<-rbind(data.all,data.tetra)
  
  #appending all the rows of data.tetra for all combinations of tetras in data.all   
  #rowbind
  
  show(sequence);
  show(dim(data.all));
}


##Now binning the data for each column
#creating empty dataframes 
breaks.all.5<-c();
breaks.all.3<-c();
medians.all.5<-c();
medians.all.3<-c();



for(param in param.names)
{
  show(paste("Processing ",param,"bins"));
  
  tmpdata<-as.numeric(data.all[,param]);
  
  
  #show(tmpdata);
  
  breaks.5<-quantile(tmpdata,probs=seq(0,1,by=0.2));
  medians.5<-quantile(tmpdata,probs=seq(0.1,1,by=0.2));
  
  breaks.3<-quantile(tmpdata,probs=c(0,0.25,0.75,1));
  medians.3<-quantile(tmpdata,probs=c(0.125,0.5,0.875));
  
  breaks.5<-signif(breaks.5,4);
  breaks.3<-signif(breaks.3,4);
  
  
  breaks.all.5<-rbind(breaks.all.5,breaks.5);
  breaks.all.3<-rbind(breaks.all.3,breaks.3);
  
  medians.all.5<-rbind(medians.all.5,medians.5);
  medians.all.3<-rbind(medians.all.3,medians.3);
}

rownames(breaks.all.5)<-param.names;
rownames(breaks.all.3)<-param.names;

rownames(medians.all.5)<-param.names;
rownames(medians.all.3)<-param.names;

show(breaks.all.5);
show(medians.all.5);
show(breaks.all.3);
show(medians.all.3);

write.table(breaks.all.5,file="results/breaks5bins.txt",quote=F);
write.table(medians.all.5,file="results/medians5bins.txt",quote=F);

write.table(breaks.all.3,file="results/breaks3bins.txt",quote=F);
write.table(medians.all.3,file="results/medians3bins.txt",quote=F);
