# Little Bootstraps 
Using the little bootstraps approach to generate bootstrap confidence limits (BCLs) on phylogenies inferred from a long sequence alignment (manuscript in review).
<br />

## Directory Structure 
The "Codes" directory contains files for two R functiosn: ``lb_sampler`` (lb_sampler.R), and ``lb_aggregator`` (lb_agrregator.R), and ``lb_precision`` (lb_precision.R). <br />
<br />
The "Example" directory contains an example data file (example.fasta) and a file containing the phylogenetic tree in the newick format (ex_candiate_tree.nwk) for which BCLs are desired. <br />
<br />

## Introduction
The little bootstraps analyses have three different steps. 
<br />
#### First step: 
<br />
The first step is to create little bootstrap replicates using the ``lb_sampler``  function in lb_sampler.R file. <br /><br /> 

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
In the second step, the maximum likelihood (ML) tree is inferred for each replicate dataset. We used the IQ-TREE software for all our ML tree inference, which can be downloaded from http://www.iqtree.org/. Both Linux and Windows versions of IQ-TREE software are available. 

<br />

#### Final step:

<br />

In the third step, aggregate ML phylogenies from all little samples to compute BCLs using the ``lb_aggregator``  function in the lb_aggregator.R file.  Inputs for the lb_aggregator function are:

```
lb_aggregator(path, tree_format, candiate_tree, s = NULL, r = NULL, output_tree = NULL)


path           : a character vector that specifies locations of the inferred ML trees. For example, inferred trees for little sample #1 should be stored in a directory named Subsample1 in the input directory.

tree_format    : a character vector that indicates the tree file format in the directory. Tree file format must be ‘.nwk’,  or ‘.treefile’

candidate_tree : an object of class "phylo" specifying the phylogeny for which BCLs are desired. 

s              : a numeric value input that specifies the number of little samples to use. If s = NULL, inferred trees from all little samples in the directory are used for computing BCLs. 

r              : a numeric value specifying the number of bootstrap replicate trees for a little sample to use. If r = NULL, inferred trees from all replicates for a little sample are used. 

output_tree    : a character vector specifying the output file name. The output is an object of class "phylo"  in ‘.nwk’ format that contains BCLs. If output_tree = NULL, the output file name will be 'output_tree_lb.nwk'.


```
<br />
If the user wants to output the candidate tree with little bootstrap estimated BCLs and precision (SE), the ``lb_precision.R`` should be used. The lb_precision function aggregates ML phylogenies from all subsamples and outputs two different tree files. One tree file contains little bootstraps BCLs, and another tree file includes the precision of BCLs for each species group. The inputs for lb_precision function are:

```
lb_precision(path, tree_format, candidate_tree, s = NULL, r = NULL, rep = 100, output_tree = NULL)


path           : a character vector that specifies locations of the inferred ML trees. For example, inferred trees for little sample #1 should be stored in a directory named Subsample1 in the input directory.

tree_format    : a character vector that indicates the tree file format in the directory. Tree file format must be ‘.nwk’,  or ‘.treefile’

candidate_tree : an object of class "phylo" specifying the phylogeny for which BCLs are desired. 

s              : a numeric value input that specifies the number of little samples to use. If s = NULL, inferred trees from all little samples in the directory are used for computing BCLs. 

r              : a numeric value specifying the number of bootstrap replicate trees for a little sample to use. If r = NULL, inferred trees from all replicates for a little sample are used. 

rep            : A numeric value indicates the number of bootstrap replicates for calculating the precision (SE) of little bootstrap BLCs. The default number of replicates is equal to 100. A user can change the value. 

output_tree    : a character vector specifying the output file name. The output is an object of class "phylo"  in ‘.nwk’ format that contains BCLs. If output_tree = NULL, the output file name will be 'output_tree_lb.nwk'.


```
<br />

<br />

## Getting Started:

<br />

To perform the little bootstraps analyses on your local computer, please follow these steps:<br /><br />
1.	Download and install R (https://www.r-project.org/) and Rstudio (https://rstudio.com/products/rstudio/download/).<br />
2.	Download ‘Codes’ directory on the local computer. <br />
3.	In the Rstudio session, type ``setwd(“directory path”)`` to change the working directory to the folder that contains ``lb_sampler`` and ``lb_aggregator`` function<br />
4.	Type ``source(lb_sampler)``, and ``source(lb_aggregator)`` to make available these  function in global environment. <br />
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

## Little Bootstraps Analyses for an Example Dataset:

<br />
To perform the little bootstraps analyses on the example dataset, please follow these steps:<br /><br />
1.	Download the ``Example`` directory on the local computer. <br />
2.	Run the function in the R session by typing 

```R
lb_sampler("~/Example/example.fasta", g = 0.9, s = 3, r = 3)
```

This function will create three directories in the working directory:

```
Subsample1
Subsample2
Subsample2
```

Each subsample directory will contain three little bootstrap replicate datasets. For example, the ``Subsample1`` directory will contain 

```
example_sub1rep1.fasta
example_sub1rep2.fasta
example_sub1rep2.fasta
```
<br />
3.	Infer ML phylogenetic tree for each replicate dataset using the IQTREE software. Users specify the substitution model and other tree inference settings subjectively for the software. For example, we used the IQTREE analysis for the replicate 1 dataset in Subsample1:<br />

``` 
iqtree -s ~/Example/Subsample1/example_sub1rep1.fasta -m GTR+G5
```

Trees for replicate datasets will be stored in each Subsample directory. The tree file name for replicate 1 in Subsample1 will be `` example_sub1rep1.fasta.treefile``. <br /><br />
4.	For the final step, type 

```R
lb_aggregator("~/Example",".treefile", "~/Example/ex_candidate.nwk", s = 3, r = 3, precision =TRUE , output_tree = "example_output")
```
The function will output the candidate tree file with BCLs, and the name of the output tree file will be `` example_output.nwk``.<br />
OR
```R
lb_precision("~/Example",".treefile", "~/Example/ex_candidate.nwk", s = 3, r = 3, rep = 100 , output_tree = "example_output")
```

The function will output the candidate tree file with BCLs and another output tree with precision(SE), and the name of output tree files will be `` example_output.nwk`` and `` example_output_precision.nwk``.<br />

<br />

#### Software and Packages' Version:

<br />

All R codes were tested using R version 3.6.3 in R studio (version 1.2.5033).
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
