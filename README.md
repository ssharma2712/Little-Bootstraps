# Little Bootstraps
The little bootstraps analyses have three different steps. The first step is to create little bootstraps replicates.<br />
<br />
<br />
The function make_subsample_upsample.R creates little bootstraps replicates. The second step is to infer maximum likelihood phylogenetic tree from these replicates. We used IQ-TREE for this step. The last and final step is to extract little sample-wise bootstrap confidence limit (bcl) and estimate the final bcl for a species group by aggregating little sample-wise bcls.   
