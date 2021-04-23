#!/usr/bin/env Rscript

library("optparse")
library("igraph")
 
option_list = list(
	make_option(c("-f", "--flavor"), type="double", default=0,
			  help="Sets the mode (or 'flavor') in which Rephine.r is run. Takes a value in {0,1,2} [default= %default]. \
			  0: HMM and fragments in tandem, followed by a second defragmentation; 1: HMM merging only; 2: fragment fusion only", metavar="double"),
	make_option(c("-d", "--panpath"), type="character", default=NULL, 
              help="path to pangenome directory", metavar="character"),
	make_option(c("-c", "--contigspath"), type="character", default=NULL,
			  help="path to contigs database directory", metavar="character"),
	make_option(c("-b","--bit"), type="double", default=0,
			  help="global minimum HMM bit score threshold [default= %default]",metavar="double"),
    make_option(c("-r", "--relbit"), type="double", default=0.25, 
              help="max relative blast bit threshold [default= %default]", metavar="double"),
	make_option(c("-p", "--maxpercoverlap"), type="double", default=0.25,
			  help="max percent overlap [default= %default]", metavar="double"),
	make_option(c("-T", "--threads"), type="double", default=2,
			  help="number of theads to use when possible [default= %default]", metavar="double")			  
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

rephine<-function(flavor,panpath,contigspath,bit,relbit,maxpercoverlap,threads)
{
	if(flavor==0)
	{
		print("Running Rephine.r with initial fragment identification and HMM merging in tandem, followed by a second round of defragmentation")
		
		print("STEP ONE: RELABELLING ANVIO OUTPUT FILES WITH PREFERRED STYLE AND REBUILDING INITIAL ALIGNMENTS")
		setwd(panpath)
		filelist=list.files()
		panpos=grep("PAN.db",filelist)
		
		clustprep(panpath,threads)					#preps initial files needed after checking to confirm the pangenome db exists
		
		if(length(panpos)>0)	
		{
			print("Doing initial round of fragment fusions")
			fragfuse(panpath,contigspath,flavor,relbit,maxpercoverlap,threads)
			
			print("BUILDING HMM PROFILES OF EACH GENE CLUSTER IDENTIFIED BY ANVIO")
			
			hmmmerge(panpath,flavor,bit,threads)		#checks HMM scan results, identifies clusters to merge, and writes out new merged cluster files in clusters_merged, followed by new alignments in clusters_merged_aligned
			
			print("Doing second round of fragment fusions")
			
			fragfuse(panpath,contigspath,flavor,relbit,maxpercoverlap,threads)		#Identify likely fragmented gene calls and report new clusters with spliced gene calls
			
			print("Rephine.r complete! Use getSCG.r to identify single-copy core genes and produce a concatenated alignment, if needed.")
		}
	}
	
	if(flavor==1)
	{
		print("Running just the HMM merging step of the Rephine.r pipeline.")
		hmmmerge(panpath,flavor,bit,threads)
	}
	
	if(flavor==2)
	{
		print("Running just the fragment fusion step of the Rephine.r pipeline.")
		fragfuse(panpath,contigspath,flavor,relbit,maxpercoverlap,threads)
	}
}

#anvirelabel is a parsing function that relabels anvio outputs that have hash labels with a preferred genome_genecall format and creates a new directory with separate fasta files for each gene cluster
anvirelabel<-function(panpath,external=TRUE)
{
	setwd(panpath)
	system('mkdir clusters')				#creates a new directory for the clusters
	#get the original genomes in same order as hashes
	dbinfo=readLines('self.txt')				#assumes the output of anvi-export-table is in the panpath directory and named "self.txt"
	if(external==TRUE){genomelist=strsplit(dbinfo[4],'\t')[[1]][2]}		#pulls out the comma separated list of genomes from the anvio db list.
	if(external==FALSE){genomelist=strsplit(dbinfo[3],'\t')[[1]][2]}	#uses the appropriate row if from internal genome db and not external.
	genomelist=unlist(strsplit(genomelist,','))
	
	#Extract hash names
	allseqs=readLines("combined-aas.fa")	#read in every sequence (with hash labels)
	allaccpos=grep(">",allseqs)
	allaccs=allseqs[allaccpos]
	allaccs=unique(unlist(strsplit(allaccs,'>')))[-1]
	allseqsvec=allseqs[-allaccpos]
	names(allseqsvec)=allaccs				#now we have the sequences named by their accessions, instead of in fasta format. This will create some indexing efficiency later on and avoid a nested loop.
	hashlist=character()
	accsplit=strsplit(allaccs,'_')
	for(i in 1:length(accsplit)){hashlist[i]=paste(accsplit[[i]][1:(length(accsplit[[i]])-1)],collapse="_")}
	unihashlist=unique(hashlist)
		
	#make a final table matching hashes to genome names
	hashmat=cbind(unihashlist,genomelist)
	colnames(hashmat)=c("hash","genome")
	rownames(hashmat)=hashmat[,1]
	names(genomelist)=unihashlist
	names(unihashlist)=genomelist		#now, we can index the hashes by the genomes and vice versa. This is faster to access than using the whole hash matrix
	#Relabel the sequence file
	newaccs=character()
	for(i in 1:length(accsplit)){newaccs[i]=paste(hashmat[accsplit[[i]][1],2],accsplit[[i]][2],sep="_")}
	#allseqs[accpos]=paste(">",newaccs,sep='')
	newseqvec=allseqsvec
	names(newseqvec)=newaccs
	newseqs=allseqs
	newseqs[allaccpos]=paste(">",newaccs,sep='')
	write.table(newseqs,'combined-aas-relabeled.fa',row.names=F,col.names=F,quote=F)
	
	#Write out MCL clusters as fasta files using the genome info, instead of hashes
	clustlist=strsplit(readLines("mcl-clusters.txt"),"\t")
	newmcl=character()
	for(i in 1:length(clustlist))			#this is slow partly because it's a nested for loop, but really because the allseqs object is just a big character vector. Could speed this up (and avoid that second loop) if restructure allseqs so that the deflines are index names for sequences.
	{
		tempn=2*length(clustlist[[i]])		#total number of positions for a new character vector with the cluster
		tempaccpos=match(clustlist[[i]],allaccs)		#gives index of the accessions in the accpos vector for each cluster gene. Takes advantage of original accpos being in hashes
		tempclust=character()
		tempmcl=names(newseqvec[tempaccpos])
		tempclust[seq(1,tempn,2)]=paste(">",tempmcl,sep='')
		tempclust[seq(2,tempn,2)]=newseqvec[tempaccpos]
		newmcl=c(newmcl,paste(tempmcl,collapse='\t'))
		write.table(tempclust,paste("./clusters/gc_",(i-1),'.fa',sep=''),row.names=F,col.names=F,quote=F)		#use 0 as first index for consistency with gene numbering
	}
	write.table(newmcl,'mcl-clusters-relabeled.txt',row.names=F,col.names=F,quote=F)
}

#hmmmerge is the main function for identifying and merging gene clusters with distant homologs that were missed in the original clustering step.
hmmmerge<-function(panpath,flavor,bit,threads)
{
	setwd(panpath)
	
	hmmpos=grep("allhmmscan_table.txt",list.files())
	if(length(hmmpos)==0)			#Check if HMMs were already produced as part of an earlier run or analysis. If no, then they will be prepared.
	{
		print("HMMs not yet complete. New profiles and queries will be performed now.")
		hmmprep(panpath,threads)
	}
	
	system('mkdir hmmmergers')										#new directory that will contain the merged clusters alone
	
	allhmm=as.matrix(read.table('allhmmscan_table.txt'))			#read in the full hmm results table
	setwd("clusters")
	clustlist=list()
	gclist=paste('gc_',c(0:(length(list.files())-1)),sep='')		#build a list of cluster names
	filelist=paste(gclist,'.fa',sep='')
	
	for(i in 1:length(filelist))
	{
		temp=readLines(filelist[i])
		clustlist[[i]]=unique(unlist(strsplit(temp[grep(">",temp)],">")))[-1]	#create a list object that summarize the cluster presence/absence of each gene call
	}
	
	gcvec=character()
	for(i in 1:length(clustlist)){gcvec=c(gcvec,rep(gclist[i],length(clustlist[[i]])))}	#creates a vector version of the list information
	names(gcvec)=unlist(clustlist)				#Assign the gene call IDs to the associated positions for the GC vectorization
	keeps=which(as.numeric(allhmm[,9])>bit)		#keep hits with minimum bit score (at domain level) above the bit parameter (default 0, i.e. unused)
	hmmsubset=allhmm[keeps,c(1,3,6)]			#record just the cluster, gene call, and full query bit score
	hmmedge=hmmsubset
	hmmedge[,2]=gcvec[hmmedge[,2]]				#create a new edgelist for a future network that relabels gene calls in column 2 with their original cluster identity
	selfmat=hmmedge[which(hmmedge[,1]==hmmedge[,2]),]		#create a smaller matrix of just self hits
	nonselfmat=hmmedge[which(hmmedge[,1]!=hmmedge[,2]),]	#create a second matrix of just the nonself hits
	unigc=unique(selfmat[,1])					#list of unique GC IDs involved in potential mergers
	selfmin=numeric()							#vector that will record the minimum bit score of any hit assigned to a given gene cluster
	bitvec=as.numeric(selfmat[,3])				#make a numeric vector version of the bit score column to speed things up later on
	
	for(i in 1:length(unigc))					#iterate through each GC and record the minimum bit score of a gene call originally assigned to that cluster
	{
		temp=which(selfmat[,1]==unigc[i])
		selfmin[i]=min(bitvec[temp])
	}
	
	names(selfmin)=unigc						#Assign the GC IDs as names to the selfmin vector
	newhits=which(as.numeric(nonselfmat[,3])>selfmin[nonselfmat[,1]])	#find the instances where a gene call outside of an original cluster has a higher hit to the hmm than a gene call within the cluster
	
	hmmnet=graph.edgelist(nonselfmat[newhits,1:2],directed=FALSE)		#create a network object from the connections between clusters with putative new hits
	hmmcomp=components(hmmnet)											#find the connected components. This will be the set of clusters to merge
	complist=character()												#create a character vector that will store comma-separated lists of components with a merger. This will just be used for result summaries later.
	
	if(flavor==0)														#for flavor 4, need the list of gcs with a fusion pre-HMM merging
	{
		fuselist=list.files('../clusters_fragfused')
	}
	
	for(i in 1:length(hmmcomp$csize))									#iterate through every hmm connected component
	{
		tempgc=names(which(hmmcomp$membership==i))
		templist=paste(tempgc,'.fa',sep='')	#build a list of cluster files for the clusters to be merged
		if(flavor==0)													#in flavor 4, update the gcs in each merger with paths to the fused versions when needed
		{
			fusegc=which(paste(tempgc,'_fused.fa',sep='')%in%fuselist)
			if(length(fusegc)>0){templist[fusegc]=paste('../clusters_fragfused/',tempgc[fusegc],'_fused.fa',sep='')}
		}
	
		if(hmmcomp$csize[i]>1){complist=c(complist,paste(templist,collapse=','))}
		newfile=character()
		for(j in 1:length(templist))									#iterate through each cluster and write to a new file. Hard to avoid this nested loop, but in general these are small mergers and not slow to run.
		{
			newfile=c(newfile,readLines(templist[j]))					#this will still work in "flavor 4" because the updated filenames are paths and not just names
		}
		write.table(newfile,paste('../hmmmergers/gcmerged_',i,'.fa',sep=''),row.names=F,col.names=F,quote=F)
	}
	
	unchanged=setdiff(gclist,names(hmmcomp$membership))				#this gives a list of all of the gene clusters not involved in a merger. We will use this in later steps when we bring it all together
	if(flavor==0)													#the unchanged list also has to be udpated (including path) for flavor 4 to account for possible fusions.
	{
		unfusedpos=which(paste(unchanged,'_fused.fa',sep='')%in%fuselist)		#fuselist will exist because of the earlier conditionals.
		if(length(unfusedpos)>0)
		{
			unchanged_fused=paste(unchanged[unfusedpos],"_fused",sep='')	#change from positional info to the actual gc_fused IDs
			unchanged=unchanged[-unfusedpos]					#drop the cases that have a fusion
			write.table(unchanged_fused,'../unmerged_clusters_fused_list.txt',row.names=F,col.names=F,quote=F)	#make a separate list
		}
	}
	write.table(unchanged,'../unmerged_clusters_list.txt',row.names=F,col.names=F,quote=F)
	write.table(complist,'../merged_clusters_list.txt',row.names=F,col.names=F,quote=F)
	
	print("Cluster merging completed!")
		
	print("Realigning clusters that were merged")
	setwd(panpath)
	system('mkdir clusters_merged_aligned; mkdir clusters_merged')
	setwd('./hmmmergers')
	system('cp *.fa ../clusters_merged/')
	if(threads==1){system('for file in $(ls); do muscle -in $file -out ../clusters_merged_aligned/$file; done')}
	if(threads>1){system(paste('ls|parallel -j ',threads,' muscle -in {} -out ../clusters_merged_aligned/{}',sep=''))}
	
	setwd('../clusters_aligned')
	system('while IFS= read -r line; do cp $line.fa ../clusters_merged_aligned/; done<"../unmerged_clusters_list.txt"')
	if(flavor==0)
	{
		if(length(unchanged_fused>0))
		{
			system('while IFS= read -r line; do cp ../clusters_fragfused_aligned/$line.fa ../clusters_merged_aligned/; done<"../unmerged_clusters_fused_list.txt"')
		}
	}
	setwd('../clusters')
	system('while IFS= read -r line; do cp $line.fa ../clusters_merged/; done<"../unmerged_clusters_list.txt"')
	if(flavor==0)
	{
		if(length(unchanged_fused)>0)
		{
			system('while IFS= read -r line; do cp ../clusters_fragfused/$line.fa ../clusters_merged/; done<"../unmerged_clusters_fused_list.txt"')
		}
	}
		
	setwd("..")
}

writefragtables<-function(panpath,rerun,clustdir)
{
	setwd(panpath)
	setwd(paste('./',clustdir,'_aligned',sep=''))			#navigate to the directory with the alignments (post cluster merging step)
		
	filelist=list.files()			#create a list of all of the alignment files
	dupcheck=numeric()				#vector that will count any genomes appearing more than once
	for(i in 1:length(filelist))	#iterate through each alignment
	{
		temp=readLines(filelist[i])
		tempgenomes=acc2genome(temp)	#convert files to a list of genomes
		temptab=table(tempgenomes)		#tabulate the result
		dupcheck=which(temptab==2)	#Check for any cases that are duplicates (ignores cases with more than 2 copies)
		if(length(dupcheck)>0)			#specifically only look for cases with exactly 2 fragments, with reasoning being that we're unlikely to piece together calls with more fragments if they would be "fusable"
		{
			fragcheck=aligncheck(temp)		#use aligncheck() to summarize the information on cases with multiple homologs, to be checked later
			if(is.matrix(fragcheck))
			{
				filename=paste(strsplit(filelist[i],'.fa')[[1]][1],'_fragmat.txt',sep='')
				if(rerun==0){write.table(fragcheck,paste('../genefragtables/',filename,sep=''),row.names=F,col.names=F,quote=F,sep='\t')}			#by default write the output to a file in "genefragtables" directory
				if(rerun==1){write.table(fragcheck,paste('../genefragtables_reFuse/',filename,sep=''),row.names=F,col.names=F,quote=F,sep='\t')}	#change to a _reFuse suffix if genefragtables results already exist
			}
		}
	}
}

#fragfuse is the main function for identifying and correcting fragmented gene calls
fragfuse<-function(panpath,contigspath,flavor,relbit,maxpercoverlap,threads)
{
	setwd(panpath)
	
	genefragpos=grep("genefragtables",list.files())
	rerun=0													#rerun variable tracks if previous fragfuse has been run
	if(length(genefragpos)>0)
	{
		if(length(list.files('./genefragtables'))>0){rerun=1}
	}
	
	clustdir="clusters_merged"						#By default, assume HMM merging was run but switch to original clusters if no HMM merge results are present. Note: this setup still works for "flavor 3" where fragfuse is run before and after HMM merging.
	
	if(!"clusters_merged"%in%list.files()){clustdir="clusters"}
	
	panpos=grep("PAN.db",list.files())
	if(!clustdir%in%list.files()&length(panpos)>0)
	{
		print("No cluster files are currently available in this directory. These will be prepared from existing pangenome results.")
		clustprep(panpath,threads)
	}
	
	print(clustdir)
	if(length(genefragpos)==0)			#Run preparation steps if no previous fragment checking has been performed.
	{
		print("Either there were no gene clusters with a potential fragment, or a previous fragment check has not been run.")
		print("Initial Fragment Check will be run now")
		
		fragprep(panpath,contigspath,rerun,clustdir,threads)
		nfragclust=length(list.files())
		print(paste("Preparation complete. There are ",nfragclust," gene clusters with potential fragmented gene calls.",sep=''))
	}
	
	if(rerun==1)
	{
		print("Preparing second round of fragment checks")
		fragprep(panpath,contigspath,rerun,clustdir,threads)
	}
	
	setwd(panpath)
		
	if(rerun==0){setwd("./genefragtables")}		#navigate to the genefragtables result directory. Note: we could alternatively do this all in one step and not have a separate file for each, but I find it helpful for debugging code and writing new functions to have individual files to work with
	if(rerun==1){setwd("./genefragtables_reFuse")}	#use _reFuse directory if previous results exist
	filelist=list.files()
	fragtable=read.table(filelist[1])
	for(i in 2:length(filelist)){fragtable=rbind(fragtable,read.table(filelist[i]))}	#pool each fragtable into one big matrix. (Could have done this in earlier step as well.)
	fragtable=as.matrix(fragtable)
	colnames(fragtable)=c("Genome","Gene1","Gene2","Overlap","percoverlap")
	gid1=gid2=numeric()				#these will have the gene ID information as numbers, useful for later steps to test adjacency
	tempsplit=strsplit(fragtable[,'Gene1'],"_")
	for(i in 1:length(tempsplit)){gid1[i]=tempsplit[[i]][length(tempsplit[[i]])]}
	tempsplit=strsplit(fragtable[,'Gene2'],"_")
	for(i in 1:length(tempsplit)){gid2[i]=tempsplit[[i]][length(tempsplit[[i]])]}
	fragtable=cbind(fragtable,gid1)
	fragtable=cbind(fragtable,gid2)
	
	setwd("../genecalltables")		#navigate to directory that contains the gene call summaries output by anvio
	
	g1start=g1stop=g2start=g2stop=rep(0,length(fragtable[,1]))		#vectors to track the start and stop positions of each gene call
	g1dir=g2dir=rep("a",length(fragtable[,1]))
	ngenes=numeric()												#vector to record the number of genes in each genome
	fraggenomes=unique(fragtable[,1])								#list of unique genomes that appear with at least one potentially fragmented call
	for(i in 1:length(fraggenomes))									#iterate through each genome
	{
		tempfrags=which(fragtable[,1]==fraggenomes[i])				#find rows in the pooled table with genes from the genome of interest
		tempcalls=as.matrix(read.table(paste(fraggenomes[i],'_genecalls.txt',sep=''),header=T))	#read in the gene call table
		ngenes[i]=length(tempcalls[,1])								#count up the number of gene calls in the genome
		for(j in 1:length(tempfrags))								#iterate through each potential fragment and record information about the start, stop, and read direction
		{
			fusecheck1=length(grep(":",fragtable[tempfrags[j],'gid1']))		#checks to see if any of the gids in fragtable reflect a previous fusion step (should be rare but will otherwise make errors if missed)
			fusecheck2=length(grep(":",fragtable[tempfrags[j],'gid2']))
			if(fusecheck1==0)
			{
				pos1=which(as.numeric(tempcalls[,1])==as.numeric(fragtable[tempfrags[j],'gid1']))	#when running a second round of fusiongs, gid's might contain ":" and won't match any original gene call. Need an if statement to correct this error
				g1start[tempfrags[j]]=tempcalls[pos1,'start']
				g1stop[tempfrags[j]]=tempcalls[pos1,'stop']
				g1dir[tempfrags[j]]=tempcalls[pos1,'direction']
			}
			if(fusecheck1==1)												#if a gid is result of a prior fusion, split the ID and look for both former pieces and report the starts and stops also with colon notation
			{
				fusesplit1=strsplit(fragtable[tempfrags[j],'gid1'],':')[[1]]
				pos1a=which(as.numeric(tempcalls[,1])==fusesplit1[1])
				pos1b=which(as.numeric(tempcalls[,1])==fusesplit1[2])
				g1start[tempfrags[j]]=paste(tempcalls[pos1a,'start'],tempcalls[pos1b,'start'],sep=':')
				g1stop[tempfrags[j]]=paste(tempcalls[pos1a,'stop'],tempcalls[pos1b,'stop'],sep=':')
				g1dir[tempfrags[j]]=tempcalls[pos1a,'direction']
			}
			if(fusecheck2==0)
			{
				pos2=which(as.numeric(tempcalls[,1])==as.numeric(fragtable[tempfrags[j],'gid2']))
				g2start[tempfrags[j]]=tempcalls[pos2,'start']
				g2stop[tempfrags[j]]=tempcalls[pos2,'stop']
				g2dir[tempfrags[j]]=tempcalls[pos2,'direction']
			}
			if(fusecheck2==1)
			{
				fusesplit2=strsplit(fragtable[tempfrags[j],'gid2'],':')[[1]]
				pos2a=which(as.numeric(tempcalls[,1])==fusesplit2[1])
				pos2b=which(as.numeric(tempcalls[,1])==fusesplit2[2])
				g2start[tempfrags[j]]=paste(tempcalls[pos2a,'start'],tempcalls[pos2b,'start'],sep=':')
				g2stop[tempfrags[j]]=paste(tempcalls[pos2a,'stop'],tempcalls[pos2b,'stop'],sep=':')
				g2dir[tempfrags[j]]=tempcalls[pos2a,'direction']
			}
		}
	}
	names(ngenes)=fraggenomes										#index each gene count by the genome name

	fragtable=cbind(fragtable,g1start)
	fragtable=cbind(fragtable,g1stop)
	fragtable=cbind(fragtable,g2start)
	fragtable=cbind(fragtable,g2stop)
	fragtable=cbind(fragtable,g1dir)
	fragtable=cbind(fragtable,g2dir)
	
	if(rerun==0){setwd('../genefragtables/')}
	if(rerun==1){setwd('../genefragtables_reFuse/')}
	
	clustvec=character()
	for(i in 1:length(filelist)){clustvec=c(clustvec,rep(filelist[i],length(read.table(filelist[i])[,1])))}	#clustvec will record the gene cluster containing each of the gene calls of interest
	clustvec=unlist(strsplit(clustvec,"_fragmat.txt"))
	fragtable=cbind(fragtable,clustvec)
	colnames(fragtable)[14]='GeneCluster'
	
	#compare potentially fragmented gene call pieces using blast results calculated outside this script
	if(rerun==0){setwd("../fragblast/")}
	if(rerun==1){setwd("../fragblast_reFuse/")}
	
	blastlist=list.files()
	
	fragblastmat=array(0,dim=c(2,4))													#set up a blast results matrix
	colnames(fragblastmat)=c("Query","Subject","Bit","MeanBit")			
	temp=read.csv(paste(fragtable[1,'GeneCluster'],'.blast',sep=''),header=F)			#read in the first blast file
	q1pos=which(temp[,1]==fragtable[1,'Gene1'])											#identify the positions of gene calls of interest within the query column of the file
	q2pos=which(temp[,1]==fragtable[1,'Gene2'])
	fragblastmat[1,'Query']=fragblastmat[2,'Subject']=fragtable[1,'Gene1']				#fill in the corresponding query and subject info in the new matrix
	fragblastmat[1,'Subject']=fragblastmat[2,'Query']=fragtable[1,'Gene2']
	check1=which(temp[q1pos,2]==fragtable[1,'Gene2'])									#check that the blast file has a hit between the gene calls of interest (it often won't when it's a fragmented gene)
	check2=which(temp[q2pos,2]==fragtable[1,'Gene1'])
	fragblastmat[1,'Bit']=temp[q1pos[check1[1]],12]
	fragblastmat[2,'Bit']=temp[q2pos[check2[1]],12]
	if(length(check1)>0){fragblastmat[1,'MeanBit']=mean(as.numeric(temp[q1pos[-check1],12]))}		#record the mean blast bit score of the queries against the other calls in the blast file (excluding the other gene call of interest, if such a hit exists)
	if(length(check1)==0){fragblastmat[1,'MeanBit']=mean(as.numeric(temp[q1pos,12]))}
	if(length(check2)>0){fragblastmat[2,'MeanBit']=mean(as.numeric(temp[q2pos[-check2],12]))}
	if(length(check2)==0){fragblastmat[2,'MeanBit']=mean(as.numeric(temp[q2pos,12]))}
	
	for(i in 2:length(fragtable[,1]))													#iterate through the remaining potential fragmented genes and repeat the above
	{
		tempmat=array(0,dim=c(2,4))
		colnames(tempmat)=c("Query","Subject","Bit","MeanBit")
		temp=read.csv(paste(fragtable[i,'GeneCluster'],'.blast',sep=''),header=F)
		q1pos=which(temp[,1]==fragtable[i,'Gene1'])
		q2pos=which(temp[,1]==fragtable[i,'Gene2'])
		tempmat[1,'Query']=tempmat[2,'Subject']=fragtable[i,'Gene1']
		tempmat[1,'Subject']=tempmat[2,'Query']=fragtable[i,'Gene2']
		check1=which(temp[q1pos,2]==fragtable[i,'Gene2'])
		check2=which(temp[q2pos,2]==fragtable[i,'Gene1'])
		tempmat[1,'Bit']=temp[q1pos[check1[1]],12]
		tempmat[2,'Bit']=temp[q2pos[check2[1]],12]
		if(length(check1)>0){tempmat[1,'MeanBit']=mean(as.numeric(temp[q1pos[-check1],12]))}
		if(length(check1)==0){tempmat[1,'MeanBit']=mean(as.numeric(temp[q1pos,12]))}
		if(length(check2)>0){tempmat[2,'MeanBit']=mean(as.numeric(temp[q2pos[-check2],12]))}
		if(length(check2)==0){tempmat[2,'MeanBit']=mean(as.numeric(temp[q2pos,12]))}
		fragblastmat=rbind(fragblastmat,tempmat)
	}
	
	relblast=as.numeric(fragblastmat[,'Bit'])/as.numeric(fragblastmat[,'MeanBit'])				#compare the bit scores of the comparison of interest to the mean bit score of the other hits. NAs set to 0
	relblast[which(is.na(relblast))]=0
	fragblastmat=cbind(fragblastmat,relblast)
	maxrelblast=numeric()
	ngenecol=numeric()	#finally, compare the relblast scores for the pair of interest with each as query or subject, and record the maximum value
	
	for(i in 1:length(fragtable[,1]))
	{
		pos1=which(fragblastmat[,'Query']==fragtable[i,'Gene1']&fragblastmat[,'Subject']==fragtable[i,'Gene2'])
		pos2=which(fragblastmat[,'Query']==fragtable[i,'Gene2']&fragblastmat[,'Subject']==fragtable[i,'Gene1'])
		maxrelblast[i]=max(relblast[pos1],relblast[pos2])
		ngenecol[i]=ngenes[fragtable[i,'Genome']]
	}
	
	fragtable=cbind(fragtable,maxrelblast)														#the fragtable object now contains all of the metrics we will use: Overlap, pospvalue, and maxrelblast
	fragtable=cbind(fragtable,ngenecol)															
	colnames(fragtable)[16]='ngenes'
	
	if(rerun==0){write.table(fragtable,'../allgcfrags_table_wmergers.txt',row.names=F,quote=F,sep='\t')}		#writes out a final version of the table that includes the pospvalue and maxrelblast scores
	if(rerun==1){write.table(fragtable,'../allgcfrags_table_wmergers_reFuse.txt',row.names=F,quote=F,sep='\t')}
	
	setwd(paste('../',clustdir,'_aligned/',sep=''))
	fragpos=which(as.numeric(fragtable[,'percoverlap']) < maxpercoverlap & as.numeric(fragtable[,'maxrelblast']) < relbit)			#Note: this pvalue might be worth adjusting for certain datasets, but important to recognize it is not dependent on the size of the data set, it is based on the size of the genes. So, it is very much a heuristic. Might deserve tweaking.
	if(length(fragpos)==0)
	{
		if(rerun==0){print("No potential fragments meet criteria for fusion")}
		if(rerun==1){print("No additional fragments meet criteria for fusion")}
	}
	
	if(length(fragpos)>0)
	{
		subtable=fragtable[fragpos,]								#create a new table with just the subset of needed rows to potentially speed up later steps for bigger data sets
		fraggclist=unique(subtable[,'GeneCluster'])					#list of unique gene clusters that have at least one putative fragmented gene call
		fragfilelist=paste(fraggclist,'.fa',sep='')								#list of the corresponding alignment files
	
		#Make the final fasta files with repaired fragmented gene calls
		for(i in 1:length(fragfilelist))										#iterate through each cluster with a putative fragmented call
		{
			fragfile=readLines(fragfilelist[i])
			fraggcpos=which(subtable[,'GeneCluster']==fraggclist[i])	#find the row positions of all fragmented calls in the current file
			accpos=grep(">",fragfile)
			accs=fragfile[accpos]
			newseq=character()													#make a new fasta file that includes all of the unaligned sequences and the mergers where appropriate.
			for(j in 1:length(accpos))											#to simplify the computation (even if it's a little inefficient) we'll copy everything over. Then when we identify fragments, we'll mark them for deletion, add a new line for the merger, and then clean up all the deletions at the very end.
			{
				if(accpos[j]==max(accpos)){stop=length(fragfile)}				#find end position of each sequence, accounting for end of the file
				if(accpos[j]!=max(accpos)){stop=accpos[j+1]-1}
				tempseq=paste(fragfile[(accpos[j]+1):stop],collapse='')			#concatenate the lines for the sequence
				tempsplit=strsplit(tempseq,'-')									#find all of the gaps
				tempseq=paste(tempsplit[[1]],collapse='')						#paste back together without the gaps
				newseq=c(newseq,fragfile[accpos[j]])
				newseq=c(newseq,tempseq)
			}
			for(j in 1:length(fraggcpos))										#iterate through each fragmented gene in the current gene cluster file
			{
				g1pos=which(accs==paste(">",subtable[fraggcpos[j],'Gene1'],sep=''))			#find the positions of each gene. Note: important to use subtable and subset by fraggcpos to only look at the proposed fusions and not all genes in the cluster.
				g2pos=which(accs==paste(">",subtable[fraggcpos[j],'Gene2'],sep=''))
				g1start=accpos[g1pos]+1
				g2start=accpos[g2pos]+1
				if(g1pos==length(accpos)){g1stop=length(fragfile)}
				if(g1pos!=length(accpos)){g1stop=accpos[g1pos+1]-1}
				if(g2pos==length(accpos)){g2stop=length(fragfile)}
				if(g2pos!=length(accpos)){g2stop=accpos[g2pos+1]-1}
				gene1=paste(fragfile[g1start:g1stop],collapse='')
				gene2=paste(fragfile[g2start:g2stop],collapse='')
				g1nogap=paste(strsplit(gene1,'-')[[1]],collapse='')
				g2nogap=paste(strsplit(gene2,'-')[[1]],collapse='')
				###NOTE: the gene direction and start/stop positions and Gene IDs are NOT sufficient to assume which of the gene fragments actually aligns to the upstream or downstream segments!!
				###Need to check the gaps in their respective alignments to pick who actually goes first! (Not sure if this is an issue with how genes are called by Prodigal, or what.
				seqsplit1=strsplit(gene1,'')[[1]]			#split the sequences to separate each character
				seqsplit2=strsplit(gene2,'')[[1]]
				gapvec1=gapvec2=rep(0,length(seqsplit1))	#each vector has to be the same length since they're in the same alignment file
				gapvec1[seqsplit1!='-']=1					#binary vector that reports a 0 for any gaps, - , and a 1 for any amino acids.
				gapvec2[seqsplit2!='-']=1
				gapcheck=gapvec1-gapvec2					#This is negative when gene2 has an aligned position and positive when gene1 has an aligned position. If gene1 is first, it will be positive early and then decline. But if gene2 is first, this will start out negative before increasing.
				gaphalf1=sum(gapcheck[1:floor(length(gapcheck)/2)])					#one way that is easy and consistent is to compare the sum of gapcheck for first vs second half and assign gene1 to be whichever half is larger. (1st > 2nd = gene1 is actually first; 1st < 2nd = gene1 is actually second)
				gaphalf2=sum(rev(gapcheck)[1:floor(length(gapcheck)/2)])
				if(gaphalf1>gaphalf2)						#merge them with gene1 first and gid1:gid2 in label
				{
					newacc=paste(">",subtable[fraggcpos[j],'Genome'],"_",subtable[fraggcpos[j],'gid1'],":",subtable[fraggcpos[j],'gid2'],sep='')
					newseq=c(newseq,newacc)
					if(fragtable[fragpos[fraggcpos[j]],'Overlap']>0)		#trim the last AA if there's an overlap and gene1 is first and add an X
					{
						tempsplit=strsplit(g1nogap,'')[[1]]
						tempsplit[length(tempsplit)]="X"
						g1nogap=paste(tempsplit,collapse='')
					}
					newseq=c(newseq,paste(g1nogap,"X",g2nogap,sep=''))
				}
				if(gaphalf1<gaphalf2)						#merge them with gene2 first and gid2:gid1 in label
				{
					newacc=paste(">",subtable[fraggcpos[j],'Genome'],"_",subtable[fraggcpos[j],'gid2'],":",subtable[fraggcpos[j],'gid1'],sep='')
					newseq=c(newseq,newacc)
					if(subtable[fraggcpos[j],'Overlap']>0)		#trimming the last AA if there's an overlap and gene1 is first
					{
						tempsplit=strsplit(g2nogap,'')[[1]]
						tempsplit[length(tempsplit)]="X"
						g2nogap=paste(tempsplit,collapse='')
					}
					newseq=c(newseq,paste(g2nogap,"X",g1nogap,sep=''))
				}
				newseq[which(newseq==accs[g1pos])]=paste(accs[g1pos],'_drop',sep='')		#update the original unmerged versions of these genes to have a "_drop" suffix so we can clean it up at the end
				newseq[which(newseq==accs[g2pos])]=paste(accs[g2pos],'_drop',sep='')
			}
			drops=grep("drop",newseq)		#which accs were updated to have a _drop because they were merged
			drops=c(drops,drops+1)			#sequences for the dropped accs
			newseq=newseq[-drops]			#This use of 2 separate loops and then cleaning up could probably be more efficient but it works. Earmark for a future update.
			
			if(rerun==0){write.table(newseq,paste('../clusters_fragfused/',fraggclist[i],'_fused.fa',sep=''),row.names=F,col.names=F,quote=F)}
			if(rerun==1){write.table(newseq,paste('../clusters_reFused/',fraggclist[i],'_fused.fa',sep=''),row.names=F,col.names=F,quote=F)}
		}
		
		setwd(panpath)
		if(rerun==0)
		{
			setwd('clusters_fragfused')
			if(threads==1){system('for file in $(ls); do muscle -in $file -out ../clusters_fragfused_aligned/$file; done')}
			if(threads>1){system(paste('ls|parallel -j ',threads,' muscle -in {} -out ../clusters_fragfused_aligned/{}',sep=''))}
			if(flavor==0)													#if doing the first round of fragments in a mode with 2 rounds planned, then this step helps prep the directory for the HMMs that follow
			{
				origgc=unlist(strsplit(list.files('../clusters'),'.fa'))
				fusegc=unlist(strsplit(list.files(),'_fused.fa'))
				cplist=setdiff(origgc,fusegc)
				write.table(cplist,'../unfused_cluster_list1.txt',row.names=F,col.names=F,quote=F)
				setwd("../clusters")
				system('while IFS= read -r line; do cp $line.fa ../clusters_fragfused/; done<"../unfused_cluster_list1.txt"')
				setwd('../clusters_aligned')
				system('while IFS= read -r line; do cp $line.fa ../clusters_fragfused_aligned/; done<"../unfused_cluster_list1.txt"')
				setwd("..")
			}
		}
		
		if(rerun==1)														#skips the additional copy step if running the 2nd round of fragment fusions
		{
			setwd('clusters_reFused')
			if(threads==1){system('for file in $(ls); do muscle -in $file -out ../clusters_reFused_aligned/$file; done')}
			if(threads>1){system(paste('ls|parallel -j ',threads,' muscle -in {} -out ../clusters_reFused_aligned/{}',sep=''))}
		}
	}
	
	print("Fragment fusion completed!")	
}

#acc2genome strips gene information off of accessions to get just the genome name
acc2genome<-function(fastafile)
{
	accs=fastafile[grep(">",fastafile)]
	accs=unlist(strsplit(accs,'>'))
	accs=accs[-which(accs=='')]
	accsplit=strsplit(accs,"_")
	genomes=character()
	for(i in 1:length(accsplit))
	{
		genomes[i]=paste(accsplit[[i]][1:length(accsplit[[i]])-1],collapse='_')
		#if(length(accsplit[[i]])==2){genomes[i]=accsplit[[i]][1]}
		#if(length(accsplit[[i]])>2){genomes[i]=paste(accsplit[[i]][1],accsplit[[i]][2],sep="_")}
	}
	return(genomes)
}

#aligncheck produces the main fragment table for a sequence alignemnt, with columns Genome, Gene1, Gene2, and Overlap, where Overlap is the size of the intersection in the alignment
aligncheck<-function(alignment)
{
	accpos=grep(">",alignment)
	accs=unique(unlist(strsplit(alignment[accpos],'>')))[-1]
	genomes=character()
	accsplit=strsplit(accs,"_")
	for(i in 1:length(accsplit))
	{
		if(length(accsplit[[i]])==2){genomes[i]=accsplit[[i]][1]}
		if(length(accsplit[[i]])>2){genomes[i]=paste(accsplit[[i]][1:(length(accsplit[[i]])-1)],collapse="_")}		#assumes RefSeq format or similar and that last position is gene information
	}
	gentab=table(genomes)
	
	if(max(gentab)==1)
	{
		results="All genes from unique genomes: No Mergers Needed"
	}
	
	if(max(gentab)>1)
	{
		dupgenomes=names(which(gentab==2))		#for now, only going to worry about genes that are split in 2. More complicated fragmentation is beyond the scope of this for now.
		if(length(dupgenomes)>0)				#possible to have bad mergers or real paralogs where the only duplicates are unfragmented copies with more than 2 present. We specifically just want genes split by early stop codons
		{
			results=array(dim=c(length(dupgenomes),5))
			results[,1]=dupgenomes
			colnames(results)=c("Genome","Gene1","Gene2","Overlap","percoverlap")	#return results in 5 column matrix: column 1 = genome with the merger; column 2 = gene ID that comes "first" (in the direction of the alignment) and column 3 = gene that comes second. column 4 is the absolute overlap and column 5 is the percent of shared aligned positions

			for(i in 1:length(dupgenomes))
			{
				genpos=which(genomes==dupgenomes[i])
				start1=accpos[genpos[1]]+1
				start2=accpos[genpos[2]]+1
				if(genpos[1]==length(accpos)){stop1=length(alignment)}
				else{stop1=accpos[genpos[1]+1]-1}
				if(genpos[2]==length(accpos)){stop2=length(alignment)}
				else{stop2=accpos[genpos[2]+1]-1}
				seq1=paste(alignment[start1:stop1],collapse='')			#do not assume one sequence to each line. 
				seq2=paste(alignment[start2:stop2],collapse='')
				seqsplit1=strsplit(seq1,'')[[1]]			#split the sequences to separate each character
				seqsplit2=strsplit(seq2,'')[[1]]
				gapvec1=gapvec2=rep(0,length(seqsplit1))	#each vector has to be same length since they're in the same alignment file
				gapvec1[seqsplit1!='-']=1					#binary vector that reports a 0 for any gaps, - , and a 1 for any amino acids.
				gapvec2[seqsplit2!='-']=1
				overlap=sum(gapvec1*gapvec2)				#counts the number of positions where each alignment has an amino acid by taking the sum product of the gap vectors
				percoverlap=overlap/(sum(sign(gapvec1+gapvec2)))
				
				results[i,4]=overlap
				results[i,5]=percoverlap
				
				if(max(which(gapvec1==1))<max(which(gapvec2==1)))
				{
					results[i,2]=accs[genpos[1]]
					results[i,3]=accs[genpos[2]]
				}
				else
				{
					results[i,2]=accs[genpos[2]]
					results[i,3]=accs[genpos[1]]
				}
					
			}
		}
	}
	
	return(results)
}

clustprep<-function(panpath,threads)
{
	setwd(panpath)
	filelist=list.files()
	panpos=grep("PAN.db",filelist)
	if(length(panpos)==0)
	{
		print("No pangenome output recognized. Please ensure that your pangenome database file ends with PAN.db")
		break
	}
	if(length(panpos)!=0)
	{
		panfile=filelist[panpos[1]]
		print(paste("pangenome database file found. ",panfile," will be used in downstream analysis. If this is not the correct file, please cancel the script and ensure only one PAN.db file is present in the directory.",sep=''))
			
		print("Exporting anvio information table")
		system(paste('anvi-export-table ',panfile,' --table self -o self.txt',sep=''))
		info=readLines('self.txt')
		infosplit=strsplit(info,'\t')
		external=TRUE											#Assumes external genomes by default
		if(length(infosplit[[3]])==2){external=FALSE}			#If internal genomes, detect this and adjust.
		
		print("Relabeling sequence headers and exporting cluster files")
		anvirelabel(panpath,external)
		print("Relabeling and cluster writing completed!")
		print("Aligning sequences with MUSCLE")
		system('mkdir clusters_aligned')
		setwd('clusters')
		if(threads==1){system('for file in $(ls); do muscle -in $file -out ../clusters_aligned/$file; done')}
		if(threads>1){system(paste('ls|parallel -j ',threads,' muscle -in {} -out ../clusters_aligned/{}',sep=''))}	#alternate version to take advantage of multithreading with "parallel"
		setwd('..')
	}
}

hmmprep<-function(panpath,threads)				#Note: this version of hmmprep() is expanded from the version in rephine_mt.r to include a step that checks for existing cluster directories.
{
		print("Building HMM profiles (warning: this step may be slow)")
		setwd(panpath)
		
		clustdir="clusters"							
	
		panpos=grep("PAN.db",list.files())
		if(length(panpos)==0)
		{
			print("No pangenome analysis results are available in the directory provided.")
			break
		}
		
		if(!clustdir%in%list.files()&length(panpos)>0)
		{
			print("No cluster files are currently available in this directory. These will be prepared from existing pangenome results.")
			clustprep(panpath,threads)
		}
		
		system('mkdir HMMs')
		setwd('clusters_aligned')
		if(threads==1){system('for file in $(ls); do hmmbuild ../HMMs/${file%.fa}.hmm $file; done')}		#debug note: noticed an issue that MAFFT ignores singletons, whereas MUSCLE still prints them out. HMMs can use the singletons, so we might want to change to just use muscle throughout. This would also be easier for anyone with successful anvio pangenome db to start with. They could just export clusters all through anvio.
		if(threads>1){system(paste('ls|parallel -j ',threads,' "hmmbuild ../HMMs/{.}.hmm {}"',sep=''))}
		print("Building the HMM database from the HMM profiles")
		setwd('..')
		system('mkdir HMMdb')
		setwd('HMMs')
		system('cat *.hmm > ../HMMdb/allHMMs')
		setwd('../HMMdb')
		system('hmmpress allHMMs')

		print("Scanning each gene call against the HMM database (warning: this step may be slow)")
		if(threads==1){system('hmmscan -o ../allhmmscan.txt --tblout ../allhmmscan_table.txt allHMMs ../combined-aas-relabeled.fa')}
		if(threads>1){system(paste('hmmscan -o ../allhmmscan.txt --cpu ',threads,' --tblout ../allhmmscan_table.txt allHMMs ../combined-aas-relabeled.fa',sep=''))}
		print("Checking for significant hits between clusters and creating merged fasta files")
		setwd('..')
}

fragprep<-function(panpath,contigspath,rerun,clustdir,threads)
{
	setwd(panpath)
	if(rerun==0)
	{
		system('mkdir genefragtables; mkdir clusters_fragfused; mkdir clusters_fragfused_aligned; mkdir genecalltables; mkdir fragblast')
	
		print("Exporting gene call tables for each genome")
		
		external=TRUE
		if(length(grep(".db",contigspath))>0)			#Identifies if using internal contigs and produces a single concatenated gene call file "allgenecalls.txt." Then parses this out into separate genecalls.txt files for every genome.
		{
			external=FALSE
			system(paste('anvi-export-gene-calls -c ',contigspath,' -o ',panpath,'/allgenecalls.txt --gene-caller prodigal --skip-sequence-reporting',sep=''))
			allgenecalls=as.matrix(read.table('allgenecalls.txt',header=T,sep='\t'))
			genomelist=unique(allgenecalls[,'contig'])
			
			for(i in 1:length(genomelist))				#Iterate through each genome, find the rows corresponding to it in the combined table, create a new file with just those rows, resetting the gene caller IDs to start at 0.
			{
				temprows=which(allgenecalls[,'contig']==genomelist[i])
				nrows=length(temprows)
				maxgene=nrows-1
				tempmat=allgenecalls[temprows,]			#subset the matrix to just the rows for the focal genome
				tempmat[,'gene_callers_id']=seq(0,maxgene)
				write.table(tempmat,paste(panpath,'/genecalltables/',genomelist[i],'_genecalls.txt',sep=''),row.names=F,quote=F,sep='\t')
			}
		}
		
		if(external==TRUE)			#How to export gene calls if using external contigs
		{
			setwd(contigspath)
			if(threads==1){system(paste('for file in $(ls); do anvi-export-gene-calls -c $file -o ',panpath,'/genecalltables/${file%.db}_genecalls.txt --gene-caller "prodigal" --skip-sequence-reporting; done',sep=''))}
			if(threads>1){system(paste('ls|parallel -j ',threads,' "anvi-export-gene-calls -c {} -o ',panpath,'/genecalltables/{.}_genecalls.txt --gene-caller prodigal --skip-sequence-reporting" ',sep=''))}
		}
	}
	
	if(rerun==1){system('mkdir genefragtables_reFuse; mkdir clusters_reFused; mkdir clusters_reFused_aligned; mkdir fragblast_reFuse')}
	
	print("Build initial fragment summary tables")
	setwd(panpath)
	writefragtables(panpath,rerun,clustdir)
	
	print("Blasting the potential fragments against the other members of their respective clusters")
	
	if(rerun==0)
	{
		setwd(paste(panpath,'/genefragtables',sep=''))
		if(threads==1){system(paste('for file in $(ls); do blastp -query ../',clustdir,'/${file%_fragmat.txt}.fa -subject ../',clustdir,'/${file%_fragmat.txt}.fa -outfmt 10 -out ../fragblast/${file%_fragmat.txt}.blast; done',sep=''))}
		if(threads>1){system(paste('ls|parallel -j ',threads,' blastp -query ../',clustdir,'/{= s/_fragmat.txt// =}.fa -subject ../',clustdir,'/{= s/_fragmat.txt// =}.fa -outfmt 10 -out ../fragblast/{= s/_fragmat.txt// =}.blast',sep=''))}
	}
	
	if(rerun==1)
	{
		setwd(paste(panpath,'/genefragtables_reFuse',sep=''))
		if(threads==1){system(paste('for file in $(ls); do blastp -query ../',clustdir,'/${file%_fragmat.txt}.fa -subject ../',clustdir,'/${file%_fragmat.txt}.fa -outfmt 10 -out ../fragblast_reFuse/${file%_fragmat.txt}.blast; done',sep=''))}
		if(threads>1){system(paste('ls|parallel -j ',threads,' blastp -query ../',clustdir,'/{= s/_fragmat.txt// =}.fa -subject ../',clustdir,'/{= s/_fragmat.txt// =}.fa -outfmt 10 -out ../fragblast_reFuse/{= s/_fragmat.txt// =}.blast',sep=''))}
	}
}

rephine(opt$flavor, opt$panpath, opt$contigspath, opt$bit, opt$relbit, opt$maxpercoverlap, opt$threads)