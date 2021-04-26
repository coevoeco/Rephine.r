## Overview

The following is a step-by-step walkthrough for implementing the Rephine.r pipeline. 


This repository also includes the scripts "getSCG.r" and "fragclass.r", which can be used following Rephine.r to build an alignment of single-copy core genes and summarize which original gene calls were 
part of a fragmentation event (described in [Follow-up analyses](https://github.com/coevoeco/Rephine.r/blob/main/README.md#follow-up-analyses)).


### Dependencies:


* [Anvi'o](https://merenlab.org/software/anvio/)
* R packages: igraph, optparse


Rephine.r currently relies on [Anvi'o](https://merenlab.org/software/anvio/) to build initial pangenomes. Future updates will enable using Rephine.r with other sources of gene clusters. An additional benefit of using Anvi'o, however, is it also installs most other
software required by Rephine.r, including MUSCLE, HMMER, and R.


#### Installing Dependencies:


1. Install the latest version of Anvi'o by following the developers' instructions [here](https://merenlab.org/2016/06/26/installation-v2/). **Note: Different versions of Anvi'o come with different versions of HMMER, and this can affect the reproducibility of Rephine.r's HMM merging step. Always make note of which version of Anvi'o and HMMER were used when running your analyses.**

2. Install the required R packages by first activating your Anvi'o conda environment. Then run:

```{bash,eval=FALSE}
R
install.packages('igraph')
install.packages('optparse')
```



## **Running Rephine.r**


### Initial Anvi'o Analysis:
Prior to using Rephine.r, you first need to run an initial pangenome analysis. We recommend following the Anvi'o workflow [here](https://merenlab.org/2016/11/08/pangenomics-v2/). 
Note: You are free to run this analysis with either "internal" or "external" genomes in Anvi'o. Rephine.r will automatically determine which format was used and adjust accordingly.

When working with phage genomes, we recommend:

1. Use minbit 0.35 and mcl-inflation 2.0. This low minbit parameter is used due to the large variation in some phage genes.
2. If you use a very large dataset, it is possible no core gene clusters will be identified by Anvi'o itself, and your anvi-pan-genome run may end with an error. This error will limit your ability to use Anvi'o's
interactive visualization and sequence reporting features, but it will not affect Rephine.r's steps, since Rephine.r relies solely on the initial inputs and mcl-cluster.txt output from Anvi'o.
3. If you plan to do functional analysis, we recommend using anvi-run-pfams in addition to anvi-run-ncbi-cogs and anvi-run-hmms to annotate your phage genomes. Pfams include
functions for many phage genes that are not annotated as NCBI COGs. Note: anvi-run-pfams unpacks a large database into a temporary directory each time. If your default temporary
directory is in a partition with limited space, and you have access to a larger drive, you can change the temporary directory path by updating the environmental variable, $TMPDIR, used by Anvi'o.
You can do this with a command in the form:
 
```{bash,eval=FALSE}
export TMPDIR=/path/to/your/preferred/temp/directory/
```


### A Note on Directory Structure:
Rephine.r and associated R scripts make two important assumptions about directory structure and naming conventions:

1. Your Anvi'o analysis will be in a directory that contains both a pangenome output directory and a directory with your contigs. The full path to the pangenome directory is later used as "panpath" in the script. 
The contigs directory is used as "contigspath". 
2. The pangenome database file within your pangenome directory MUST end with 'PAN.db'.


### Executing the Script:
At minimum, Rephine.r requires only the paths to your pangenome output and contigs to run. By default, it will run both HMM merging and fragment identification steps in tandem, combine the results to create a new set of 
gene clusters, and then run a second step of fragment identification. If desired, the user can choose to run only HMM merging or only fragment identification by changing the "flavor" option (-f) to 1 or 2, respectively.


A default run requires only:

```{bash,eval=FALSE}
Rscript --vanilla rephine.r -d /path/to/pangenome/ -c /path/to/ContigsDB/
```


The complete set of options are:

        -f DOUBLE, --flavor=DOUBLE
                Sets the mode (or 'flavor') in which Rephine.r is run. Takes a value in {0,1,2} [default=%default].
                0: HMM and fragments in tandem, followed by a second defragmentation (default) 
                1: HMM merging only 
                2: fragment fusion only
                
        -d CHARACTER, --panpath=CHARACTER
                path to pangenome directory

        -c CHARACTER, --contigspath=CHARACTER
                path to contigs database directory

        -b DOUBLE, --bit=DOUBLE
                global minimum HMM bit score threshold [default= 0]

        -r DOUBLE, --relbit=DOUBLE
                max relative blast bit threshold [default= 0.25]

        -p DOUBLE, --maxpercoverlap=DOUBLE
                max percent overlap [default= 0.25]
                
        -T DOUBLE, --threads=DOUBLE
                number of threads to use when possible [default= 2]

        -h, --help
                Show this help message and exit




## **Follow-up Analyses**


### Identifying the Single-Copy Core Genes with getSCG.r


The original purpose of Rephine.r was to correct gene calls and clusters so we could obtain more accurate single-copy core genomes and build more accurate phylogenies. 
To facilitate this, we put together a script called getSCG.r which identifies the core and single-copy core genes in a set of clusters and returns a concatenated alignment file (if desired).

Given a completed Rephine.r analysis, getSCG.r can be run as:

```{bash,eval=FALSE}
Rscript --vanilla getSCG.r -d /path/to/pangenome/
```


The full options for getSCG.r are:

        -d CHARACTER, --panpath=CHARACTER
                path to pangenome directory

        -a CHARACTER, --aligned=CHARACTER
                logical specifying if a concatenated alignment is desired [default= TRUE]

For any run of getSCG.r, a list of accessions will be provided named with a prefix in the form "scg_i" and suffix "gclist.txt" with *i* replaced by either "orig", "fused", or "merged". If a concatenated alignment is desired, these will be named analogously in the format scg_i.fa.


getSCG.r will return different sets of SCGs based on each stage in the Rephine.r pipeline. Files with the "scg_orig" prefix reflect the uncorrected gene clusters. The "scg_fused" prefix files depend only on fragment identification and ignore HMM merging, 
whereas files with "scg_merged" include both the outcomes of HMM merging and fragment identification. These are the fully corrected core genes. 
These alignment files are then appropriate for tree inference with your preferred program.


### Predicting the Cause of Fragmentation Events with fragclass.r


Fragmentation events in phage genomes are typically caused by one of 3 events: indels, interruption by a selfish element, such as a homing endonuclease ("SGE"), or being split across the ends of the genome in the original sequence file.


Each of these potential causes can be distinguished by where the gene fragments are located relative to each other in the genome. The script fragclass.r makes these predictions and returns two files:

1. fragtypemat.txt: A table summarizing the type of event that likely caused the fragmentation, including "indel", "SGE", or "undescribed" (for cases that do not resemble either option).
2. anvio_frag_view_data.txt: A data table that is suitable for importing fragment information into Anvi'o genome visualization tools.
