#!/usr/bin/env Rscript

library("optparse")
library("igraph")
 
option_list = list(
	make_option(c("-d", "--panpath"), type="character", default=NULL, 
              help="path to pangenome directory", metavar="character"),
    make_option(c("-a", "--aligned"), type="character", default=TRUE,
			  help="logical specifying if a concatenated alignment is desired [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#getSCG() returns the single-copy core genome for a pangenome following rephine.r. It takes as input the path to the pangenome directory
# and writes out the SCG for each step of whichever analysis was run (based on what directories are present in the folder). It also includes the option (TRUE by default) to output a concatenated alignment FASTA file.

getSCG<-function(panpath,aligned)
{
	setwd(panpath)
	filelist=list.files()
	fragcheck1=length(which(filelist=="clusters_fragfused"))
	fragcheck2=length(which(filelist=="clusters_reFused"))
	mergecheck=length(which(filelist=="clusters_merged"))
	clustpos=which(filelist=="clusters")
	
	#First, check if any cluster analysis exist before proceeding
	if(length(clustpos)==0){print("No analyses available to work with. Please run a rephine.r analysis first.")}
	
	if(length(clustpos)>0)
	{
		#Build the original SCG without any fusion or merger steps
		scg=clust2core("clusters_aligned")$scg
		nscg=length(scg)
		if(nscg==0){print("No single-copy core genes without merging or fusing")}
		if(nscg>0){write.table(scg,'../scg_orig_gclist.txt',row.names=F,col.names=F,quote=F)}
		
		if(aligned==TRUE & nscg==0)
		{
			print("Alignment not possible. No single-copy core genes")
		}
		if(aligned==TRUE & nscg==1)
		{
			allalign=acc2genome(readLines(scg))
			write.table(allalign,'../scg_orig.fa',row.names=F,col.names=F,quote=F)
		}
		if(aligned==TRUE & nscg>1)
		{
			allalign=alignpaste(acc2genome(readLines(scg[1])),acc2genome(readLines(scg[2])))
			if(nscg>2)
			{
				for(i in 3:nscg)
				{
					allalign=alignpaste(allalign,acc2genome(readLines(scg[i])))
				}
			}
			write.table(allalign,'../scg_orig.fa',row.names=F,col.names=F,quote=F)
		}
		setwd(panpath)
		
		#If a fragfused folder exists, return the SCG based just on those files
		if(fragcheck1>0)
		{
			scgfused=clust2core("clusters_fragfused_aligned")$scg
			nscg=length(scgfused)
			if(nscg==0){print("No single-copy core genes after first fusion")}
			if(nscg>0){write.table(scgfused,'../scg_fused_gclist.txt',row.names=F,col.names=F,quote=F)}
		
			if(aligned==TRUE & nscg==0)
			{
				print("Alignment of SCG with first fusion not possible. No single-copy core genes")
			}
			if(aligned==TRUE & nscg==1)
			{
				allalignfused=acc2genome(readLines(scgfused))
				write.table(allalignfused,'../scg_fused.fa',row.names=F,col.names=F,quote=F)
			}
			if(aligned==TRUE & nscg>1)
			{
				allalignfused=alignpaste(acc2genome(readLines(scgfused[1])),acc2genome(readLines(scgfused[2])))
				if(nscg>2)
				{
					for(i in 3:nscg)
					{
						allalignfused=alignpaste(allalignfused,acc2genome(readLines(scgfused[i])))
					}
				}
				write.table(allalignfused,'../scg_fused.fa',row.names=F,col.names=F,quote=F)
			}
			setwd(panpath)
		}
		
		#If a merged folder exists, return the SCG based on those results (will include effects of prior fusion if those were run as part of "flavor 3")
		if(mergecheck>0)
		{
			scgmerged=clust2core("clusters_merged_aligned")$scg
			nscg=length(scgmerged)
			if(nscg==0){print("No single-copy core genes after first fusion and merging")}
			if(nscg>0){write.table(scgmerged,'../scg_merged_gclist.txt',row.names=F,col.names=F,quote=F)}
		
			if(aligned==TRUE & nscg==0)
			{
				print("Alignment of SCG with first fusion and merging not possible. No single-copy core genes")
			}
			if(aligned==TRUE & nscg==1)
			{
				allalignmerged=acc2genome(readLines(scgmerged))
				write.table(allalignmerged,'../scg_merged.fa',row.names=F,col.names=F,quote=F)
			}
			if(aligned==TRUE & nscg>1)
			{
				allalignmerged=alignpaste(acc2genome(readLines(scgmerged[1])),acc2genome(readLines(scgmerged[2])))
				if(nscg>2)
				{
					for(i in 3:nscg)
					{
						allalignmerged=alignpaste(allalignmerged,acc2genome(readLines(scgmerged[i])))
					}
				}
				write.table(allalignmerged,'../scg_merged.fa',row.names=F,col.names=F,quote=F)
			}
			setwd(panpath)
		}

		#If only mergers and initial fusions exist and only 1 round of fusion was run (basically "flavor 0"), return a concatentated alignment of the two SCG sets (if they would add to each other)
		if(fragcheck1>0 & mergecheck>0 & fragcheck2==0 & aligned==TRUE)
		{
			if(length(scgfused)>0 & length(scgmerged)>0)
			{
				allalign=alignpaste(allalignmerged,allalignfused)
				write.table(allalign,'../scg_merged_fused.fa',row.names=F,col.names=F,quote=F)
			}
		}
		
		#If a second round of fusions was run (essentially only possible if "flavor 3"), return the SCG specific to the new fusions
		if(fragcheck2>0)
		{
			nrefuse=length(list.files('clusters_reFused'))
			if(nrefuse==0){print("No Additional Fusions were added after HMM merging.")}
			if(nrefuse>0)
			{
				scgrefused=clust2core("clusters_reFused_aligned")$scg
				nscg=length(scgrefused)
				
				if(nscg>0){write.table(scgrefused,'../scg_reFused_gclist.txt',row.names=F,col.names=F,quote=F)}
				
				if(aligned==TRUE & nscg==0)
				{
					print("Alignment not possible. No single-copy core genes added by second fusion")
				}
				if(aligned==TRUE & nscg==1)
				{
					allalignrefused=acc2genome(readLines(scgrefused))			
					write.table(allalignrefused,'../scg_reFused.fa',row.names=F,col.names=F,quote=F)
				}
				if(aligned==TRUE & nscg>1)
				{
					allalignrefused=alignpaste(acc2genome(readLines(scgrefused[1])),acc2genome(readLines(scgrefused[2])))
					if(nscg>2)
					{
						for(i in 3:nscg)
						{
							allalignfused=alignpaste(allalignrefused,acc2genome(readLines(scgrefused[i])))
						}
					}
					write.table(allalignrefused,'../scg_reFused.fa',row.names=F,col.names=F,quote=F)
				}
				
				#If new fusions added to the SCG and a previous SCG from merging exists, combine for a final concatenated alignment (if alignment requested)
				if(aligned==TRUE & nscg>0 & length(scgmerged)>0)
				{
					allalign=alignpaste(allalignmerged,allalignrefused)
					write.table(allalign,'scg_merged_reFused.fa',row.names=F,col.names=F,quote=F)
				}
			}
		}	
	}
}

#clust2core() computers a pres/abs and scg based on the gene clusters in fasta format in the specified directory.
clust2core<-function(clustpath)
{
	require(Matrix)
	setwd(clustpath)
	filelist=list.files()
	genomelist=character()
	acclist=list()
	for(i in 1:length(filelist))
	{
		acclist[[i]]=acc2genome(readLines(filelist[i]),out="acc")
	}
	genomelist=unique(unlist(acclist))
	
	pres=Matrix(0,nrow=length(genomelist),ncol=length(filelist))
	rownames(pres)=genomelist
	colnames(pres)=filelist
	
	for(i in 1:length(filelist))
	{
		temptab=table(acclist[[i]])
		pres[names(temptab),i]=as.numeric(temptab)
	}
	
	core=names(which(apply(sign(pres),2,sum)==length(genomelist)))
	scg=intersect(core,names(which(apply(pres,2,sum)==length(genomelist))))
	
	results=list()
	results$pres=pres
	results$scg=scg
	
	return(results)
}

#alignpaste() concatenates two alignment files with the same headers.
alignpaste<-function(file1,file2)
{
	#file1=readLines(align1)		#read in the first file in filelist
	#file2=readLines(align2)		#commented out since better is to give these in already read with a wrapper to iterate through more than 2.
	accpos1=grep(">",file1)		#determine the line numbers where the accessions are found
	accpos2=grep(">",file2)
	seqpos1=cbind(accpos1,append(accpos1[-1]-1,length(file1)))	#finds starts and stops for sequences in file1 (keeping accessions)
	seqpos2=cbind(accpos2+1,append(accpos2[-1]-1,length(file2)))	#finds starts and stop for sequences in file2 (but drops accession line)
	nseq=length(accpos1)			#count number of accessions
	accs1=character()			#prep a variable to accept the accessions
	accs2=character()
	allalign=character()
	resort=numeric()
	for(i in 1:nseq)			#iterate through number of accessions
	{
		accs1[i]=strsplit(file1[accpos1],">")[[i]][2]		#store each accession in accs1. This sets the default order for the others.
		accs2[i]=strsplit(file2[accpos2],">")[[i]][2]
	}

	for(i in 1:nseq)
	{
		resort[i]=which(accs2==accs1[i])		#finds out how the accessions are reordered in the current file relative to file1
	}

	for(i in 1:nseq)
	{
		allalign=append(allalign,append(file1[seqpos1[i,1]:seqpos1[i,2]],file2[seqpos2[resort[i],1]:seqpos2[resort[i],2]]))
	}
				
	return(allalign)
}

#acc2genome() returns just the genome information from the fasta headers
acc2genome<-function(fasta,out="fasta")
{
	accpos=grep(">",fasta)
	acc=fasta[accpos]
	acc=unique(unlist(strsplit(acc,">")))[-1]
	accsplit=strsplit(acc,'_')
	genomelist=character()
	for(i in 1:length(accsplit))
	{
		genomelist[i]=paste(accsplit[[i]][1:(length(accsplit[[i]])-1)],collapse="_")
	}
	
	newacc=paste(">",genomelist,sep='')
	fasta[accpos]=newacc
	
	if(out=="acc"){return(newacc)}
	if(out=="fasta"){return(fasta)}
}

getSCG(opt$panpath, opt$aligned)