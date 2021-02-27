# Little Bootstraps
The little bootstraps analyses have three different steps. <br />
<br />
<br />
First step:, the lb_sampler function in lb_sampler.R is used for creatin little bootstrap replicates The function inputs are sequence alignment in fasta format, g (0.6<= g <= 0.9, and subsample size is L^0.g), the number of little samples (s), and the number of replicates (r). This function will create directory for little samples (Subsample1, Subsample2, etc.) in the working directory and save the little bootstrap replicate datasets in these subsaple directories.<br /> <br />
lb_sampler("/~example.fasta", g= 0.9, s = 2, r = 2)

<br />
<br />
The function make_subsample_upsample.R creates little bootstraps replicates. The second step is to infer maximum likelihood phylogenetic tree from these replicates. We used IQ-TREE for this step. The last and final step is to extract little sample-wise bootstrap confidence limit (bcl) and estimate the final bcl for a species group by aggregating little sample-wise bcls.   
