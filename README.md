# Little Bootstraps
The little bootstraps analyses have three different steps. <br />
<br />
<br />
First step: The lb_sampler function in lb_sampler.R is used for creatin little bootstrap replicates The function inputs are sequence alignment in fasta format, g (0.6<= g <= 0.9, and subsample size is L^0.g), the number of little samples (s), and the number of replicates (r). This function will create directory for little samples (Subsample1, Subsample2, etc.) in the working directory and save the little bootstrap replicate datasets in these subsaple directories.<br /> <br />
lb_sampler("~/example.fasta", g= 0.9, s = 2, r = 2)
<br />
<br />
Second step: In this step, the maximum likelihood (ML) tree is inferred for each replicate dataset. In our analyses, we used IQ-TREE for ML tree inference which can be downloaded from http://www.iqtree.org/. Both linux and windows versions of IQ-TREE software are available here. Other ML tree inference program like MEGA (https://www.megasoftware.net/), RAxML (https://cme.h-its.org/exelixis/web/software/raxml/), PHYLIP (https://evolution.genetics.washington.edu/phylip.html) etc., also can be used. <br /> <br />
iqtree -s /Subsample1/example_sub1_rep.fasta -m GTR+G5 
<br />
<br />
The function make_subsample_upsample.R creates little bootstraps replicates. The second step is to infer maximum likelihood phylogenetic tree from these replicates. We used IQ-TREE for this step. The last and final step is to extract little sample-wise bootstrap confidence limit (bcl) and estimate the final bcl for a species group by aggregating little sample-wise bcls.   
