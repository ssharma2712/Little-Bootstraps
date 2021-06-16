# Little Bootstraps 
Using the little bootstraps approach to generate bootstrap confidence limits (BCLs) on phylogenies inferred from a long sequence alignment (manuscript in review). We provide two R functions for automatic little bootstrap analysis. 
<br />

## Directory Structure 
The "Automatic" directory contains files for two R functions: ``lb_automatic`` (lb_automatic.R), and ``little_bootstraps`` (little_bootstraps.R). This directory also contains ``iqtree.exe`` for inferring ML trees using little bootstrap replicates. The fasta file ``mtCDNA.fas`` and the phylogenetic tree file ``mtCDNA.nwk`` are for the test run of these two R functions. The ``lb_automatic`` (lb_automatic.R) requires only the fasta file and the candidate tree file for inferring the little bootstrap BCLs with their precision (Standard Error, SE). This function automatically determines little bootstraps parameters:

```
g      : Subsample size
s      : Number of subsamples
r      : Number of replicates

```
The ``little_bootstraps`` (little_bootstraps.R) uses the user defined little bootstrap parameters ``(g,s,r)``. 
<br />


## Description of R functions

The R function ``lb_automatic`` can automatically select the little bootstrap parameters and estimate the BCLs, but the ``little_bootstraps`` require user-specific little bootstrap parameters for the analysis.

```
lb_automatic(data_path, candidate_tree, evo_model = NULL, output_tree = NULL, del = 0.001, precision = FALSE)


data_path           : a character vector that specifies locations of the inferred ML trees. For example, fasta file mtCDNA for little bootstrap analysis in the Automatic folder. Therefore, the data_path will be "~/Automatic/mtCDNA.fas"

candidate_tree : an object of class "phylo" specifying the phylogeny for which BCLs are desired. 

evo_model      : a string vector that specifies the substitution model for inferring ML trees from little bootstrap replicates. If NULL, the model will be determined by IQTREE.

output_tree    : a character vector specifying the output file name. The output is an object of class "phylo"  in ‘.nwk’ format that contains BCLs. If output_tree = NULL, the output file name will be 'output_tree_lb.nwk'.

del            : a numeric value to specify the threshold of the change in the average BCLs for selecting the number of subsamples and replicates. The value should be less than 1. For example, if a user allows 1% change in average BCLs, the del = 0.01.

precision      : TRUE/FALSE. If TRUE (T), output files are objects of class "phylo"  in ‘.nwk’ format. If output_file = NULL, the output file name will be 'output_tree_lb.nwk', and 'output_tree_lb_precision.nwk'
```
<br />

<br />

```
little_bootstraps(data_path, candidate_tree, evo_model = NULL, lb_parameter = c(0.8,10,10), output_tree = NULL, del = 0.001, precision = FALSE)


data_path           : a character vector that specifies locations of the dataset. For example, there is a fasta file "mtCDNA.fas" in the Automatic folder. Therefore, the data_path will be "~/Automatic/mtCDNA.fas"

candidate_tree : an object of class "phylo" specifying the phylogeny for which BCLs are desired. 

evo_model      : a string vector that specifies the substitution model for inferring ML trees from little bootstrap replicates. If NULL, the model will be determined by IQTREE.

lb_parameter   : a numeric vector specifying little bootstrap parameters (g,s,r)

output_tree    : a character vector specifying the output file name. The output is an object of class "phylo"  in ‘.nwk’ format that contains BCLs. If output_tree = NULL, the output file name will be 'output_tree_lb.nwk'.

del            : a numeric value to specify the threshold of the change in the average BCLs for selecting the number of subsamples and replicates. The value should be less than 1. For example, if a user allows 1% change in average BCLs, the del = 0.01.

precision      : TRUE/FALSE. If TRUE (T), output files are objects of class "phylo"  in ‘.nwk’ format. If output_file = NULL, the output file name will be 'output_tree_lb.nwk', and 'output_tree_lb_precision.nwk'
```
<br />

<br />

## Getting Started:

<br />

To perform the little bootstraps analyses on your local computer, please follow these steps:<br /><br />
1.	Download and install R (https://www.r-project.org/) and Rstudio (https://rstudio.com/products/rstudio/download/).<br />
2.	Download the ‘Automatic’ directory on the local computer. <br />
3.	In the Rstudio session, type ``setwd(“directory path”)`` to change the working directory to the folder that contains dataset and R codes<br />
4.	Type ``source(lb_automatic.R)``, and ``source(little_bootstraps.R)`` to make available these functions in global environment. <br />
5.	Install the following R packages if those are not installed. 

```R
install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("stringr")
install.packages("ape")
install.packages("phangorn")
install.packages("phyclust")
```

<br />

## Little Bootstraps Analyses for an Example Dataset:

<br />
To perform the little bootstraps analyses on the example dataset, please follow these steps:<br /><br />
1.	Download the ``Automatic`` directory on the local computer. <br />
2.	Run the function in the R session by typing 

```R
lb_automatic("~/Automatic/mtCDNA.fas", candidate_tree = "~/Automatic/mtCDNA.nwk", evo_model = "GTR+G4", output_tree = NULL, del = 0.001, precision = TRUE)

OR

little_bootstraps("~/Automatic/mtCDNA.fas", candidate_tree = "~/Automatic/mtCDNA.nwk", evo_model = "GTR+G4", lb_parameter = c(0.8,10,10), output_tree = NULL, del = 0.001, precision = TRUE)
```

Both functions will output two trees, one tree (output_tree_lb.nwk) with little bootstrap BCLs and another tree (output_tree_lb_precision.nwk) with the precision (SE). If ``precision = FALSE``, there will be one output tree (output_tree_lb.nwk) with BCLs.


#### Software and Packages' Version:

<br />

All R codes were tested using R version 3.6.3 in R studio (version 1.2.5033). We used IQ-TREE (multicore version 1.6.12 for Windows 64-bit built Aug 15 2019) for ML tree inferences.
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
