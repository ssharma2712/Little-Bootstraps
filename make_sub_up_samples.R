if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
library(Biostrings)

rm(list = ls())
mother_dir <- "C:/Users/Administrator/Downloads/2012PNAS_data/test_rr50/gene162" # Mother data directory
setwd(mother_dir)

f_name <- "MAM2_concatenated.fasta"   # Mother data name
sub_name <- "MAM2_concatenated_sub"   # Replicate file generic name
motherfile <- Biostrings::readAAStringSet(f_name, format = "fasta")  # Reading the mother file 
sln <- as.numeric(fasta.seqlengths(f_name)[1])                       # Getting the sequence length
for(i in 1:20){
  setwd(mother_dir)
  sub_dir <- paste("Subsample", i , sep = "")               
  dir.create(sub_dir)                                                # Create Subsample folder where Replictes will be kept
  setwd(sub_dir)
  s <- sample(1:sln, ceiling(sln^0.8), replace = F)                  # Setting subsample length and sites. The power value changes [0.5, 1]
  subsample <- endoapply(motherfile, function(x) x[s])               # Making subsample
  
  for(j in 1:20){
    s <- sample(1:ceiling(sln^0.8), sln, replace = T)                # Setting the replicate length and sites. 
    upsample <- endoapply(subsample, function(x) x[s])               
    file_name <- paste(sub_name, i, "rep", j, ".fasta", sep = "")
    Biostrings::writeXStringSet(upsample, file_name)                 # Saving Replicates
    print(c(i, j))
  }
}



