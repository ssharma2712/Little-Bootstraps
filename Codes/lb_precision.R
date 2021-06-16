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
    stop("'ape' package not found, please install it to run lb_precision")
  }
  
  if (!library('phangorn',logical.return = TRUE)){
    stop("'phangorn' package not found, please install it to run lb_precision")
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
    ff <- tempfile()
    png(filename = ff)
    b <- phangorn::plotBS(candidate_tree, x, p =10,  type = "unrooted")
    dev.off()
    unlink(ff)
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
      ff <- tempfile()
      png(filename = ff)
      b <- phangorn::plotBS(candidate_tree, x, p =10,  type = "unrooted")
      dev.off()
      unlink(ff)
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

