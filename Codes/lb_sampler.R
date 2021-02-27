lb_sampler <- function(data_path, g, s, r){
  
  # data_path: path for the sequence alignment in fasat format
  # g        : subsample size (L^g) 0.6<= g <= 0.9, L = sequence length
  # s        : Number of subsamples
  # r        : Number of upsamples
  
  ########## Package required ##########
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!requireNamespace("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")
  
  if (!requireNamespace("stringr", quietly = TRUE))
    BiocManager::install("stringr")
  
  ########## Library required ##########
  
  if (!library('Biostrings',logical.return = TRUE)){
    stop("'Biostrings' package not found, please install it to run lb_sampler")
  }
  
  if (!library('stringr',logical.return = TRUE)){
    stop("'stringr' package not found, please install it to run lb_sampler")
  }
  ######################################
  
  f_name <- data_path                                                  # Mother data name
  sub_name <- str_replace(basename(a), ".fasta", "")                   # Replicate file generic name
  motherfile <- Biostrings::readAAStringSet(f_name, format = "fasta")  # Reading the mother file 
  sln <- as.numeric(fasta.seqlengths(f_name)[1])                       # Getting the sequence length
  directory <- getwd()
  for(i in 1:s){
    setwd(directory)
    sub_dir <- paste("Subsample", i , sep = "")               
    dir.create(sub_dir)                                                # Create Subsample folder where Replictes will be kept
    setwd(sub_dir)
    s <- sample(1:sln, ceiling(sln^g), replace = F)                    # Setting subsample length and sites. The power value changes [0.5, 1]
    subsample <- endoapply(motherfile, function(x) x[s])               # Making subsample
    
    for(j in 1:r){
      s <- sample(1:ceiling(sln^g), sln, replace = T)                  # Setting the replicate length and sites. 
      upsample <- endoapply(subsample, function(x) x[s])               
      file_name <- paste(sub_name, '_sub', i, "rep", j, ".fasta", sep = "")
      Biostrings::writeXStringSet(upsample, file_name)                 # Saving Replicates
      print(c('Subsample=',i, 'Replicates=', j))
    }
  }
  setwd(directory)
}
