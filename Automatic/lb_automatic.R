lb_automatic <- function(data_path, candidate_tree, evo_model = NULL, output_tree = NULL, del = 0.001, precision = FALSE){
  
  # data_path     : Directory path where the fasta file located
  # candidate_tree: The candidate tree that will be assessed
  # evo_model     : Substitution model
  # output_tree   : Output tree file name
  # del           : Treshold value for choosing s and r
  # precision     : Precision (SE) for little bootstrap BCL's
  #               : If true there will be an output tree with precision value 
  
  ########## Package required #####################
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  if (!requireNamespace("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")
  
  if (!requireNamespace("stringr", quietly = TRUE))
    install.packages("stringr")
  
  if (!requireNamespace("ape", quietly = TRUE))
    install.packages("ape")
  
  if (!requireNamespace("phyclust", quietly = TRUE))
    install.packages("phyclust")
  
  if (!requireNamespace("phangorn", quietly = TRUE))
    install.packages("phangorn")
  
  
  ################# Library Required ##################
  
  
  if (!library('Biostrings',logical.return = TRUE)){
    stop("'Biostrings' package not found, please install it to run lb_automatic")
  }
  
  if (!library('stringr',logical.return = TRUE)){
    stop("'stringr' package not found, please install it to run lb_automatic")
  }
  
  if (!library('phyclust',logical.return = TRUE)){
    stop("'phyclust' package not found, please install it to run lb_automatic")
  }
  ##############Helper Function ########
  aggregator <- function(path, tree_format, candidate_tree, s = NULL, r = NULL, output_file = NULL){
    
    # path          : Directory path where tree file for each replicate datasets located
    # tree_format   : Tree file format (.nwk, .treefile)
    # candidate_tree: The candidate tree that will be assessed
    # s             : Number of subsamples
    # r             : Number of replicates 
    # output_file   : Output tree file name
    
    ################# Required Package ##################
    
    if (!requireNamespace("ape", quietly = TRUE))
      install.packages("ape")
    
    if (!requireNamespace("phangorn", quietly = TRUE))
      install.packages("phangorn")
    
    ########## Library required ##########
    
    if (!library('ape',logical.return = TRUE)){
      stop("'ape' package not found, please install it to run aggregator")
    }
    
    if (!library('phangorn',logical.return = TRUE)){
      stop("'phangorn' package not found, please install it to run aggregator")
    }
    
    #######################################
    
    sub           <- NULL
    sub_dir       <- dir(path, pattern = "Subsample", full.names = T)
    candidate_tree <- ape::read.tree(candidate_tree)
    
    if(is.null(s) == T) {
      s <- length(sub_dir)
    }
    
    for(k in 1:s){                                                                            # Number of little samples
      di <- sub_dir[k]
      lf <- list.files(di, pattern = tree_format, full.names = TRUE)
      
      if(is.null(r) == T) {
        r <- length(lf)
      }
      
      x<- ape::rmtree(r, candidate_tree$Nnode)
      
      for(l in 1:r){   # number of replicates
        x[l] <- list(ape::read.tree(lf[l]))
        print(c("Subsampe=", k, "Replicate=", l))
      }
      b <- phangorn::plotBS(candidate_tree, x, p =10,  type = "unrooted")
      sub <- cbind(sub, as.numeric(b$node.label)/100)
      
    }
    med_sup <- apply(sub, 1, function(x){median(x, na.rm = T)})
    b$node.label <- med_sup
    if(is.null(output_file) == T) {
      output_tree <- 'output_tree_lb.nwk'
    }else{
      output_tree <- paste(output_file, '.nwk', sep = "")
    }
    
    ape::write.tree(b, file = output_tree)
  }
  lb_avg_bcl <- function(path, tree_format, candidate_tree, s = NULL, r = NULL){
    
    # path          : Directory path where tree file for each replicate datasets located
    # tree_format   : Tree file format (.nwk, .treefile)
    # candidate_tree: The candidate tree that will be assessed
    # s             : Number of subsamples
    # r             : Number of replicates 
    # output_tree   : Output tree file name
    
    ################# Required Package ##################
    
    if (!requireNamespace("ape", quietly = TRUE))
      install.packages("ape")
    
    if (!requireNamespace("phangorn", quietly = TRUE))
      install.packages("phangorn")
    
    ########## Library required ##########
    
    if (!library('ape',logical.return = TRUE)){
      stop("'ape' package not found, please install it to lb_avg_bcl")
    }
    
    if (!library('phangorn',logical.return = TRUE)){
      stop("'phangorn' package not found, please install it to lb_avg_bcl")
    }
    
    #######################################
    
    sub           <- NULL
    sub_dir       <- dir(path, pattern = "Subsample", full.names = T)
    candidate_tree <- ape::read.tree(candidate_tree)
    
    if(is.null(s) == T) {
      s <- length(sub_dir)
    }
    
    for(k in 1:s){                                                                            # Number of little samples
      di <- sub_dir[k]
      lf <- list.files(di, pattern = tree_format, full.names = TRUE)
      
      if(is.null(r) == T) {
        r <- length(lf)
      }
      
      x<- ape::rmtree(r, candidate_tree$Nnode)
      
      for(l in 1:r){   # number of replicates
        x[l] <- list(ape::read.tree(lf[l]))
        print(c("Subsampe=", k, "Replicate=", l))
      }
      b <- phangorn::plotBS(candidate_tree, x, p =10,  type = "unrooted")
      sub <- cbind(sub, as.numeric(b$node.label)/100)
      
    }
    med_sup <- apply(sub, 1, function(x){median(x, na.rm = T)})
    return(mean(med_sup, na.rm = T))
  }
  lb_precision <- function(path, tree_format, candidate_tree, s = NULL, r = NULL, rep = 100, output_tree = NULL){
    
    # path          : Directory path where tree file for each replicate datasets located
    # tree_format   : Tree file format (.nwk, .treefile)
    # candidate_tree: The candidate tree that will be assessed
    # s             : Number of subsamples
    # r             : Number of replicates 
    # rep           : Number of replicate for computing precision of lb BCL's 
    # output_tree   : Output tree file name
    
    ################# Required Package ##################
    
    if (!requireNamespace("ape", quietly = TRUE))
      install.packages("ape")
    
    if (!requireNamespace("phangorn", quietly = TRUE))
      install.packages("phangorn")
    
    ########## Library required ##########
    
    if (!library('ape',logical.return = TRUE)){
      stop("'ape' package not found, please install it to run aggregator")
    }
    
    if (!library('phangorn',logical.return = TRUE)){
      stop("'phangorn' package not found, please install it to run aggregator")
    }
    
    #######################################
    
    sub           <- NULL
    sub_dir       <- dir(path, pattern = "Subsample", full.names = T)
    candidate_tree <- ape::read.tree(candidate_tree)
    
    if(is.null(s) == T) {
      s <- length(sub_dir)
    }
    
    for(k in 1:s){                                                                            # Number of little samples
      di <- sub_dir[k]
      lf <- list.files(di, pattern = tree_format, full.names = TRUE)
      
      if(is.null(r) == T) {
        r <- length(lf)
      }
      
      x<- ape::rmtree(r, candidate_tree$Nnode)
      for(l in 1:r){   # number of replicates
        x[l] <- list(ape::read.tree(lf[l]))
        print(c("Subsampe=", k, "Replicate=", l))
      }
      b <- phangorn::plotBS(candidate_tree, x, p =10,  type = "unrooted")
      sub <- cbind(sub, as.numeric(b$node.label)/100)
      
    }
    med_sup <- apply(sub, 1, function(x){median(x, na.rm = T)})
    b$node.label <- med_sup
    b_tree <- b
    rm(b)
    
    #############Precision calculation##############
    
    sub_all <- NULL
    for(m in 1:rep){
      sub           <- NULL
      sub_dir       <- dir(path, pattern = "Subsample", full.names = T)
      
      if(is.null(s) == T) {
        s <- length(sub_dir)
      }
      
      ran_sub <- sample(1:s, s, replace = T)
      for(k in 1:s){                                                                            # Number of little samples
        di <- sub_dir[ran_sub[k]]
        lf <- list.files(di, pattern = tree_format, full.names = TRUE)
        
        if(is.null(r) == T) {
          r <- length(lf)
        }
        
        x<- ape::rmtree(r, candidate_tree$Nnode)
        ran_rep <- sample(1:r, r, replace = T)
        for(l in 1:r){   # number of replicates
          x[l] <- list(ape::read.tree(lf[ran_rep[l]]))
          
        }
        b <- phangorn::plotBS(candidate_tree, x, p =10,  type = "unrooted")
        sub <- cbind(sub, as.numeric(b$node.label)/100)
        
      }
      med_sup <- apply(sub, 1, function(x){median(x, na.rm = T)})
      sub_all <- cbind(sub_all, med_sup)
      print(c("Replicate = ", m))
    }
    pre <- apply(sub_all, 1, function(x){sd(x, na.rm = T)})
    b$node.label <-  pre
    
    
    if(is.null(output_tree) == T) {
      output_tree <- 'output_tree_lb.nwk'
      output_tree_p <- 'output_tree_lb_precision.nwk'
    }else{
      output_tree <- paste(output_tree, '.nwk', sep = "")
      output_tree_p <- paste(output_tree, '_precision.nwk', sep = "")
    }
    
    ape::write.tree(b_tree, file = output_tree)
    ape::write.tree(b, file = output_tree_p)
    
  }
  
  #########################################
  
  f_name <- data_path                                                  # Mother data name
  sub_name <- str_replace(basename(data_path), ".fasta", "")           # Replicate file generic name
  motherfile <- Biostrings::readAAStringSet(f_name, format = "fasta")  # Reading the mother file 
  sln <- as.numeric(fasta.seqlengths(f_name)[1])                       # Getting the sequence length
  directory <- getwd()
  a <- phyclust::read.fasta(data_path)
  a <- a$org.code
  a <- unique(a, MARGIN = 2)
  if(dim(a)[2] > 100000){
    g <- 0.7
  }else{
    g <- 0.8
  }
  rm(a)
  s1 <- list()
  for(i in 1:4){
    setwd(directory)
    sub_dir <- paste("Subsample", i , sep = "")               
    dir.create(sub_dir)                                                # Create Subsample folder where Replictes will be kept
    setwd(sub_dir)
    s1 <- append(s1, list(sample(1:sln, ceiling(sln^g), replace = F)))                    # Setting subsample length and sites. The power value changes [0.5, 1]
    #subsample[[]] <- endoapply(motherfile, function(x) x[s])               # Making subsample
    
    for(j in 1:4){
      s2 <- sample(s1[[i]], sln, replace = T)                  # Setting the replicate length and sites. 
      upsample <- endoapply(motherfile, function(x) x[s2])               
      file_name <- paste(sub_name, '_sub', i, "rep", j, ".fasta", sep = "")
      Biostrings::writeXStringSet(upsample, file_name)                 # Saving Replicates
      #print(c('Subsample=',i, 'Replicates=', j))
    }
  }
  setwd(directory)
  for(i in 1:4){
    for(j in 1:4){
      dname <- paste(sub_name, '_sub', i, "rep", j, ".fasta", sep = "")
      if(is.null(evo_model) == F){
        shell(paste("iqtree -s", paste("Subsample", i, "/", dname, sep = ""), "-m", evo_model, sep = " "))
      }else{
        shell(paste("iqtree -s", paste("Subsample", i, "/", dname, sep = ""), sep = " "))
      }
    }
  }
  
  ###########
  sub_avg <- c(0,0)
  rep_avg <- c(0,0)
  
  for(j in 3:4){
    a <- lb_avg_bcl(dirname(di), ".treefile", candidate_tree, s = 3, r = j)
    rep_avg <- c(rep_avg, a)
  }
  
  while(abs(rep_avg[j]-rep_avg[j-1])>del){
    j = j +1
    for(m in 1:i){
      for(n in j:j){
        s2 <- sample(s1[[m]], sln, replace = T)                  # Setting the replicate length and sites. 
        upsample <- endoapply(motherfile, function(x) x[s2])               
        file_name <- paste(sub_name, '_sub', m, "rep", n, ".fasta", sep = "")
        setwd(paste("Subsample", m, sep = ""))
        Biostrings::writeXStringSet(upsample, file_name)                 # Saving Replicates
        setwd(directory)
        if(is.null(evo_model) == F){
          shell(paste("iqtree -s", paste("Subsample", m, "/", file_name, sep = ""), "-m", evo_model, sep = " "))
        }else{
          shell(paste("iqtree -s", paste("Subsample", m, "/", file_name, sep = ""), sep = " "))
        }
        #print(c('Subsample=',m, 'Replicates=', n))
      }
    }
    a <- lb_avg_bcl(dirname(di), ".treefile", candidate_tree, s = i, r = j)
    rep_avg <- c(rep_avg, a)
  }
  
  for(i in 3:4){
    a <- lb_avg_bcl(dirname(di), ".treefile", candidate_tree, s = i, r = j)
    sub_avg <- c(sub_avg, a)
  }
  while(abs(sub_avg[i]-sub_avg[i-1])>del){
    i = i +1
    print(i)
    for(m in i:i){
      setwd(directory)
      sub_dir <- paste("Subsample", i , sep = "")               
      dir.create(sub_dir)                                                # Create Subsample folder where Replictes will be kept
      setwd(sub_dir)
      s1 <- append(s1, list(sample(1:sln, ceiling(sln^g), replace = F)))                   # Setting subsample length and sites. The power value changes [0.5, 1]
      #subsample[[]] <- endoapply(motherfile, function(x) x[s])               # Making subsample
      for(n in 1:j){
        s2 <- sample(s1[[m]], sln, replace = T)                  # Setting the replicate length and sites. 
        upsample <- endoapply(motherfile, function(x) x[s2])               
        file_name <- paste(sub_name, '_sub', i, "rep", j, ".fasta", sep = "")
        Biostrings::writeXStringSet(upsample, file_name)                 # Saving Replicates
        setwd(directory)
        if(is.null(evo_model) == F){
          shell(paste("iqtree -s", paste("Subsample", m, "/", file_name, sep = ""), "-m", evo_model, sep = " "))
        }else{
          shell(paste("iqtree -s", paste("Subsample", m, "/", file_name, sep = ""), sep = " "))
        }
        setwd(sub_dir)
        #print(c('Subsample=',i, 'Replicates=', j))
      }
    }
    a <- lb_avg_bcl(dirname(di), ".treefile", candidate_tree, s = i, r = j)
    sub_avg <- c(sub_avg, a)
  }
  setwd(directory)
  
  j = j + 1
  for(m in 1:i){
    for(n in j:j){
      s2 <- sample(s1[[m]], sln, replace = T)                  # Setting the replicate length and sites. 
      upsample <- endoapply(motherfile, function(x) x[s2])               
      file_name <- paste(sub_name, '_sub', m, "rep", n, ".fasta", sep = "")
      setwd(paste("Subsample", m, sep = ""))
      Biostrings::writeXStringSet(upsample, file_name)                 # Saving Replicates
      setwd(directory)
      if(is.null(evo_model) == F){
        shell(paste("iqtree -s", paste("Subsample", m, "/", file_name, sep = ""), "-m", evo_model, sep = " "))
      }else{
        shell(paste("iqtree -s", paste("Subsample", m, "/", file_name, sep = ""), sep = " "))
      }
      #print(c('Subsample=',m, 'Replicates=', n))
    }
  }
  a <- lb_avg_bcl(dirname(di), ".treefile", candidate_tree, s = i, r = j)
  rep_avg <- c(rep_avg, a)
  
  while(abs(rep_avg[j]-rep_avg[j-1])>del){
    j = j +1
    for(m in 1:i){
      for(n in j:j){
        s2 <- sample(s1[[m]], sln, replace = T)                  # Setting the replicate length and sites. 
        upsample <- endoapply(motherfile, function(x) x[s2])               
        file_name <- paste(sub_name, '_sub', m, "rep", n, ".fasta", sep = "")
        setwd(paste("Subsample", m, sep = ""))
        Biostrings::writeXStringSet(upsample, file_name)                 # Saving Replicates
        setwd(directory)
        if(is.null(evo_model) == F){
          shell(paste("iqtree -s", paste("Subsample", m, "/", file_name, sep = ""), "-m", evo_model, sep = " "))
        }else{
          shell(paste("iqtree -s", paste("Subsample", m, "/", file_name, sep = ""), sep = " "))
        }
        #print(c('Subsample=',m, 'Replicates=', n))
      }
    }
    a <- lb_avg_bcl(dirname(di), ".treefile", candidate_tree, s = 3, r = i)
    rep_avg <- c(rep_avg, a)
  }
  print(c("Subsample = ", i, "Replicate = ", j))
  setwd(directory)
  
  if(precision == FALSE){
    aggregator(dirname(di), ".treefile", candidate_tree, s = i, r = j)
  }else{
    lb_precision(dirname(di), ".treefile", candidate_tree, s = i, r = j)
  }
  
}

