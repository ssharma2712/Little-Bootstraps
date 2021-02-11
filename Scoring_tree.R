setwd("/working directory")  # Working Directory
true_tree <- ape::read.tree("~/true_tree.nwk")   # True tree

sub <- NULL
for(k in 1:10){     # 10 = number of little sample                                                                       # Number of little samples
  di <- paste("~/Subsample", k, sep = "" )
  lf <- list.files(di, pattern = ".treefile", full.names = TRUE)
  x<- ape::rmtree(10, true_tree$Nnode)
  s <- sample(1:10, 10, replace = F)
  for(l in 1:10){   # 10 = number of replicates
    x[l] <- list(ape::read.tree(lf[s[l]]))
  }
  b <- phangorn::plotBS(true_tree, x, p =10,  "phylogram")
  sub <- cbind(sub, as.numeric(b$node.label)/100)
  print(k)
}
med_sup <- apply(sub, 1, function(x){median(x, na.rm = T)})
b$node.label <- med_sup

ape::write.tree(b, file = "file_name.nwk")