# Little Bootstraps 
Little bootstraps produce accurate bootstrap confidence limits (BCLs) on phylogenies inferred from very long sequene alignmnet. 
<br />

## Directory Structure 
"Codes" directory contains ``lb_sampler`` (lb_sampler.R) and ``aggregator`` (agrregator.R) R functions. <br />
<br />
"Example" directory contains example data (example.fasta), and a candidate tree (ex_candiate_tree.nwk) for performing little bootstrap analyses. <br />
<br />

## Introduction
Little bootstraps analyses have three different steps. 
<br />
#### First step: 
<br />
The first step of the little bootstrap analyses is to create little bootstrap replicates. The lb_sampler function in lb_sampler.R is used.  

```
lb_sampler(data_path,  g,  s,  r)


data_path         : input sequence alignment in fasta format that will be used for little bootstrap analyses. 

g                 : a numeric value within the range (0.6<= g <= 0.9) that specify the little sample size. The little sample size is equal to L^g where L is the sequence length determined from the input alignment.  

s                 : a numeric value that specifies the number of little samples. 

r                 : a numeric value that specifies the number of replicates for each little sample.
```
<br />

#### Second step:

<br />
In this step, the maximum likelihood (ML) tree is inferred for each replicate dataset. The choice of ML inference software and computations are flexible. Users can use any ML tree inference software and compute ML tree for little bootstrap replicates sequentially or parallelly based on their computational architecture. In our analyses, we used IQ-TREE for ML tree inference which can be downloaded from http://www.iqtree.org/. Both Linux and Windows versions of IQ-TREE software are available here. Other ML tree inference program like MEGA (https://www.megasoftware.net/), RAxML (https://cme.h-its.org/exelixis/web/software/raxml/), PHYLIP (https://evolution.genetics.washington.edu/phylip.html)

<br />

#### Final step:

<br />

The last step of little bootstraps analyses is to aggregate results from each subsample. In this step, all inferred trees are aggregated to compute little sample-wise BCLs. The estimated little bootstrap supports are computed using median bagging. The aggregator function in aggregator.R is used for this step. Inputs for the aggregator function are:

```
aggregator(path, tree_format, candiate_tree, s = NULL, r = NULL, output_tree = NULL)


path           : a character vector that specifies the locations where all inferred trees are located. For example, inferred trees for little sample #1 should be stored in a directory named Subsample1 in the input directory.

tree_format    :  a character vector that indicates the tree file format in the directory. Tree file format must be ‘.nwk’,  or ‘.treefile’

candidate_tree :  an object of class "phylo" specifying the candidate tree. The BCLs are placed on this candidate tree. 

s              :  a numeric value input that specifies the number of little samples that will be used. If s = NULL, inferred trees from all little samples in the input directory are used for computing BCLs. 

r              : a numeric value specifying the number of replicate trees will be used from each subsample. If r = NULL, inferred trees from all replicates for a little sample are used. 

output_file    : a character vector specifying the output file name. The output is an object of class `phylo, in ‘.nwk’ format that is the input candidate tree with BCLs. If output_file = NULL, the output file name will be 'output_tree_lb.nwk'.
```
<br />

<br />

## Getting Started:

<br />

To perform the little bootstraps analyses in a local computer, please follow these steps:<br /><br />
1.	Download and install R (https://www.r-project.org/) and Rstudio (https://rstudio.com/products/rstudio/download/).<br />
2.	Download ‘Codes’ directory on the local computer. <br />
3.	In the Rstudio session, type ``setwd(“directory path”)`` to change the working directory to the folder that contains ``lb_sampler`` and ``aggregator`` function<br />
4.	Type ``source(lb_sampler)``, and ``source(aggregator)`` to make available these  function in global environment. <br />
5.	Download and install an ML tree inference software (e.g., IQ-TREE). <br />
6.	Install folowwing R packages if those are not installed. 

```R
install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("stringr")
install.packages("ape")
install.packages("phangorn")
```

<br />

## Little Bootstraps Analyses for Example Dataset:

<br />
To perform the little bootstraps analyses on the example dataset, please follow these steps:<br /><br />
1.	Download the ``Example`` directory on the local computer. <br />
2.	Run the function in the R session by typing 

```R
lb_sampler("~/Example/example.fasta", g= 0.9, s = 3, r = 3)
```

This function will create three directories in the working directories:

```
Subsample1
Subsample2
Subsample2
```

Each subsample directory will contain three little bootstraps replicate datasets. For example, the ``Subsample1`` directory will contain 

```
example_sub1rep1.fasta
example_sub1rep2.fasta
example_sub1rep2.fasta
```
<br />
3.	Infer ML phylogenetic tree for each replicate dataset using the preferred software. Users can flexibly perform this step. Users specify the substitution model and other tree inference settings subjectively for the software. As an example of IQ-TREE analysis for the replicate 1 in the Subsample1:<br />

``` 
iqtree -s ~/Example/Subsample1/example_sub1rep1.fasta -m GTR+G5
```

Trees for each replicate dataset will be stored in each Subsample directory. The tree file name for replicate 1 in Subsample1 will be `` example_sub1rep1.fasta.treefile`` if we use the IQ-TREE. In case of using other software, tree-files need to be converted into ``.treefile``, or ``.nwk`` format. <br /><br />
4.	For the final step, type 

```R
aggregator("~/Example",".treefile", "~/Example/ex_candidate.nwk", s = 3, r = 3, output_file = "example_output")
```

The function will output the candidate tree file with the little bootstraps supports, and the name of the output tree file will be `` example_output.nwk``.<br />

<br />

#### Software and Packages's Version:

<br />

All R codes are tested using R version 3.6.3 in R studio (version 1.2.5033).
<br />  
R packages used:
<br />

```
-BiocManager (version 1.30.10)
-Biostrings  (version 2.54.0)
-stringr     (version 1.4.0)
-ape         (version 5.3)
-phangorn    (version 2.5.5)
```

<br />
<br />
