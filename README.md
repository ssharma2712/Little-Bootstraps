# Little Bootstraps
The little bootstraps analyses have three different steps. <br />
<br />
First step: The lb_sampler function in lb_sampler.R is used for creatin little bootstrap replicates The function inputs are sequence alignment in fasta format, g (0.6<= g <= 0.9, and subsample size is L^0.g), the number of little samples (s), and the number of replicates (r). This function will create directory for little samples (Subsample1, Subsample2, etc.) in the working directory and save the little bootstrap replicate datasets in these subsaple directories.<br /> 
```R
lb_sampler("~/Example/example.fasta", g= 0.9, s = 10, r = 10)
```
<br />
Second step: In this step, the maximum likelihood (ML) tree is inferred for each replicate dataset. In our analyses, we used IQ-TREE for ML tree inference which can be downloaded from http://www.iqtree.org/ - automatic!. Both linux and windows versions of IQ-TREE software are available here. Other ML tree inference program like MEGA (https://www.megasoftware.net/), RAxML (https://cme.h-its.org/exelixis/web/software/raxml/), PHYLIP (https://evolution.genetics.washington.edu/phylip.html) etc., also can be used. <br /> <br />

```md
iqtree -s /Example/Subsample1/example_sub1_rep.fasta -m GTR+G5 
```
<br />
Final step: This is an aggregation step. In this steps, all inferred trees are aggregated from each Subsample directories and compute the final bootstrap confidence limit (BCL) using median bagging. The aggregator function in aggregator.R is used for this this step. The inputs for the aggregator function are the path where inferred trees are located in Subsample directories, tree file format (.nwk, or .treefile), the candiate tree on which the BCLs will be placed, the number of subsample (s, if NULL the function will use all Subsample resultes in the directory), the number of replicates (r, if NULL the function will use all tree file in the Subsample directory), and the output file name (if NULL, the output tree file name will be output_tree_lb.nwk).  
<br /> <br />

```R
aggregator("~/Example",".treefile", "~/Example/ex_candidate.nwk", s = 10, r = 10, output_file = "example_output")
```

<br /> <br />
All R codes are tested using R version 3.6.3 in R studio (version 1.2.5033).
<br />  
R packages used:
<br />
-BiocManager (version 1.30.10)
-Biostrings  (version 2.54.0)
-stringr     (version 1.4.0)
-ape         (version 5.3)
-phangorn    (version 2.5.5) 
