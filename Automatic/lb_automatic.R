
aggregator <- function(data_path, candiate_tree, output_tree = NULL){
  
  # data_path     : Directory path where the fasta file located
  # tree_format   : Tree file format (.nwk, .treefile)
  # candidate_tree: The candidate tree that will be assessed
  # s             : Number of subsamples
  # r             : Number of replicates 
  # output_tree   : Output tree file name
  
  ########## Package required ##########
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!requireNamespace("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")
  
  if (!requireNamespace("stringr", quietly = TRUE))
    install.packages("stringr")
  
  if (!requireNamespace("ape", quietly = TRUE))
    install.packages("ape")
  
  if (!requireNamespace("phangorn", quietly = TRUE))
    install.packages("phangorn")
  
  
  ################# Library Required ##################
  
  if (!requireNamespace("ape", quietly = TRUE))
    install.packages("ape")
  
  if (!requireNamespace("phangorn", quietly = TRUE))
    install.packages("phangorn")
  
  if (!library('Biostrings',logical.return = TRUE)){
    stop("'Biostrings' package not found, please install it to run lb_sampler")
  }
  
  if (!library('stringr',logical.return = TRUE)){
    stop("'stringr' package not found, please install it to run lb_sampler")
  }
  ###################################
  
  f_name <- data_path                                                  # Mother data name
  sub_name <- str_replace(basename(data_path), ".fasta", "")                   # Replicate file generic name
  motherfile <- Biostrings::readAAStringSet(f_name, format = "fasta")  # Reading the mother file 
  sln <- as.numeric(fasta.seqlengths(f_name)[1])                       # Getting the sequence length
  directory <- getwd()
  for(i in 1:3){
    setwd(directory)
    sub_dir <- paste("Subsample", i , sep = "")               
    dir.create(sub_dir)                                                # Create Subsample folder where Replictes will be kept
    setwd(sub_dir)
    s <- sample(1:sln, ceiling(sln^g), replace = F)                    # Setting subsample length and sites. The power value changes [0.5, 1]
    subsample <- endoapply(motherfile, function(x) x[s])               # Making subsample
    
    for(j in 1:3){
      s <- sample(1:ceiling(sln^g), sln, replace = T)                  # Setting the replicate length and sites. 
      upsample <- endoapply(subsample, function(x) x[s])               
      file_name <- paste(sub_name, '_sub', i, "rep", j, ".fasta", sep = "")
      Biostrings::writeXStringSet(upsample, file_name)                 # Saving Replicates
      print(c('Subsample=',i, 'Replicates=', j))
    }
  }
  setwd(directory)
  
}
  


setwd("C:/Users/Administrator/Downloads/Manuscript_Bootstrap/Revised/Empirical_data/SBS/lb_auto_test")
shell("iqtree -s Yoneswa.fas -m GTR+G4")
