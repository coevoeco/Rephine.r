#!/usr/bin/env Rscript
library("optparse")
 
option_list = list(
	make_option(c("-d", "--panpath"), type="character", default=NULL, 
              help="path to pangenome directory", metavar="character"),
	make_option(c("-r", "--relbit"), type="double", default=0.25, 
              help="max relative blast threshold [default= %default]", metavar="double"),
	make_option(c("-p", "--maxpercoverlap"),type="double", default=0.25,
			  help="max percent overlap [default= %default]", metavar="double")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#fragclass() returns summary data following a fragment analysis that lists the number of putative fragmented gene calls in each genome. This requires a completed fragfuse.r run with output in the directory clusters_fragfused.
#The only needed input is the path to the main pangenome analysis directory. The function will also use the gene ID information to classify the fragment as due to indels (neighboring) and/or reporting (one of the positions is 0), or insertion (not neighboring)

fragclass<-function(panpath,maxpercoverlap,relbit)
{
	setwd(panpath)
	
	rerun=0																#determines if a second round of fusions has been run
	recheck=grep("allgcfrags_table_wmergers_reFuse.txt",list.files())
	if(length(recheck)>0){rerun=1}
	
	fragmat=as.matrix(read.table("allgcfrags_table_wmergers.txt",header=T,sep='\t'))
	if(rerun==1)														#if there was a second round of fusions, combine info from the two tables before proceeding
	{
		fragmat2=as.matrix(read.table("allgcfrags_table_wmergers_reFuse.txt",header=T,sep='\t'))
		vec1=apply(fragmat[,1:3],1,paste,collapse=',')					#create simple vectors that uniquely identify the possible rows in each matrix based on the first 3 columns (which genes involved)
		vec2=apply(fragmat2[,1:3],1,paste,collapse=',')
		newrows=which(!vec2%in%vec1)									#find any new rows that are in vec2 that weren't in vec1
		if(length(newrows)>0)											#if there are any new rows, add those to the original fragmat. Then remaining code works whether or not the changes were made.
		{
			fragmat=rbind(fragmat,fragmat2[newrows,])
		}
		colcheck1=grep(":",fragmat[,'gid1'])							#need to correct gid's for gdiff later on for any cases that involve a second potential fusion for one gene.
		colcheck2=grep(":",fragmat[,'gid2'])
		if(length(colcheck1)>0)
		{
			for(i in 1:length(colcheck1))
			{
				tempsplit=strsplit(fragmat[colcheck1[i],'gid1'],":")
				fragmat[colcheck1[i],'gid1']=tempsplit[[1]][2]			#later, we need gdiff to help identify indels vs introns, so we need to favor the SMALLEST diff between gid1 and gid2. For + direction genes, that means taking the bigger part of gid1, and for -, the smaller. In either case, that's the second position.
			}	
		}
		
		if(length(colcheck2)>0)
		{
			for(i in 1:length(colcheck2))
			{
				tempsplit=strsplit(fragmat[colcheck2[i],'gid2'],":")
				fragmat[colcheck2[i],'gid2']=tempsplit[[1]][1]			#with gid2 events, the smallest difference to gid1 is always based on the first position, whether + or - oriented.
			}	
		}
	}
	
	maxgene=as.numeric(fragmat[,'ngenes'])-1
	fragmat=cbind(fragmat,maxgene)
	keeps=which(as.numeric(fragmat[,'percoverlap'])<maxpercoverlap&as.numeric(fragmat[,'maxrelblast'])<relbit)
	newmat=fragmat[keeps,c("Genome","gid1","gid2","maxgene")]
	newmat[,2]=as.numeric(newmat[,2])
	newmat[,3]=as.numeric(newmat[,3])
	info=readLines('self.txt')
	infosplit=strsplit(info,'\t')
	if(length(infosplit[[3]])==2){genomelist=strsplit(infosplit[[3]][2],',')[[1]]}			#autodetect if internal or external genomes and act accordingly
	if(length(infosplit[[4]])==2){genomelist=strsplit(infosplit[[4]][2],',')[[1]]}
	
	fragtype=rep('undescribed',length(newmat[,1]))
	edgecase=rep(0,length(newmat[,1]))
	
	gdiff=abs(as.numeric(newmat[,'gid1'])-as.numeric(newmat[,'gid2']))						#this throws an error when using reFuse data. NEED TO CORRECT IT TO ACCOUNT FOR THOSE CASES
	
	for(i in 1:length(newmat[,1]))
	{
		if(gdiff[i]==1){fragtype[i]='indel'}													#Identify all indels
		if(newmat[i,'gid1']=="0"&newmat[i,'gid2']==newmat[i,'maxgene']){edgecase[i]=1}			#Identify all neighboring edgecases (not necessarily indels; could just be a sequence reporting error)
		if(newmat[i,'gid2']=="0"&newmat[i,'gid1']==newmat[i,'maxgene']){edgecase[i]=1}			#Same as above but swapping gid1 and gid2
		if(gdiff[i]>1&edgecase[i]==1){gdiff[i]=1}												#Reset the cases above to a gdiff of 1. This criterion is not true in general, but will be true at this stage in the analysis because edgecases where gdiff>1 have not been defined yet.
		if(gdiff[i]>1){fragtype[i]="SGE"}														#Identify any cases with gdiff > 1 as due to selfish genetic elements ("SGE")
		if(fragtype[i]=="SGE"&abs(as.numeric(newmat[i,'gid1'])-0)<5&abs(as.numeric(newmat[i,'gid2'])-as.numeric(newmat[i,'maxgene']))<5){edgecase[i]=1}			#Identify SGE cases that overlap the beginning/end of the reporter genome (not sure if this is useful). Use a simple heuristic that both gids are within 5 genes of the maxgene and 0
		if(fragtype[i]=="SGE"&abs(as.numeric(newmat[i,'gid2'])-0)<5&abs(as.numeric(newmat[i,'gid1'])-as.numeric(newmat[i,'maxgene']))<5){edgecase[i]=1}			#Same as above but swapping gid1 and gid2
	}
	
	newmat=cbind(newmat,gdiff)
	newmat=cbind(newmat,fragtype)
	newmat=cbind(newmat,edgecase)
	
	misslist=setdiff(genomelist,newmat[,'Genome'])
	nmiss=length(misslist)
	
	if(nmiss>0)
	{
		missmat=array(dim=c(nmiss,7))
		colnames(missmat)=colnames(newmat)
		missmat[,'Genome']=misslist
		for(i in 1:nmiss)
		{
			temp=as.matrix(read.table(paste('./genecalltables/',misslist[i],'_genecalls.txt',sep=''),header=T,sep='\t'))
			missmat[i,'maxgene']=length(temp[,1])-1
		}
	}
	
	newmat=rbind(newmat,missmat)
	
	write.table(newmat,'fragtypemat.txt',row.names=F,quote=F,sep='\t')
	
	#Final step produces a summary file suitable for importing into anvi'o with gids removed. This summary file only includes one line per genome with the number of each type, not a line for each individual fragment case
	anvimat=array(0,dim=c(length(genomelist),6))
	colnames(anvimat)=c("contig","ngenes","nfrags","n_indels","n_SGE","has_edgecase")
	anvimat[,'contig']=genomelist
	for(i in 1:length(genomelist))
	{
		temprows=which(newmat[,'Genome']==genomelist[i])
		
		tempmat=newmat[temprows,]
		if(length(temprows)>1)
		{
			anvimat[i,'ngenes']=as.numeric(tempmat[1,'maxgene'])+1
			anvimat[i,'nfrags']=length(temprows)
			anvimat[i,'n_indels']=length(which(tempmat[,'fragtype']=='indel'))
			anvimat[i,'n_SGE']=length(which(tempmat[,'fragtype']=='SGE'))
			anvimat[i,'has_edgecase']=0+sum(as.numeric(tempmat[,'edgecase']))
		}
		if(length(temprows)==1)
		{
			anvimat[i,'ngenes']=as.numeric(tempmat['maxgene'])+1
			if(!is.na(tempmat['fragtype']))
			{
				anvimat[i,'nfrags']=length(temprows)
				anvimat[i,'n_indels']=length(which(tempmat['fragtype']=='indel'))
				anvimat[i,'n_SGE']=length(which(tempmat['fragtype']=='SGE'))
				anvimat[i,'has_edgecase']=tempmat['edgecase']
			}
		}
	}
	anvimat=cbind(anvimat[,1],anvimat)
	colnames(anvimat)[1]='contig'
	colnames(anvimat)[2]='label'
	
	write.table(anvimat,'anvio_frag_view_data.txt',row.names=F,quote=F,sep='\t')
}
	
fragclass(opt$panpath, opt$maxpercoverlap, opt$relbit)