#setwd("/working directory")  # Working Directory

lb_aggregator <- function(path, tree_format, candidate_tree, s = NULL, r = NULL, output_tree = NULL){
  
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
  
  if(is.null(output_tree) == T) {
    output_tree <- 'output_tree_lb.nwk'
  }else{
    output_tree <- paste(output_tree, '.nwk', sep = "")
  }
  
  ape::write.tree(b, file = output_tree)
  
}
