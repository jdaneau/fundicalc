source("convert_tree.R")

p_debug <- function(...) {
  for(p in list(...)) {
    cat(toString(p))
  }
  cat("\n")
}

#####
# Functions from gfmix

LocateRoot <- function(utreec,taxa){
  ##IN: taxa - labels 0,1,.... corresponding to one side of the split
  if(length(taxa)==1) return(taxa)
  ntaxa <- dim(utreec)[1]+1
  
  if(0 %in% taxa){
    spl.of.int <- rep(0,ntaxa)
    spl.of.int[(taxa+1)] <- 1
  }else{
    spl.of.int <- rep(1,ntaxa)
    spl.of.int[(taxa+1)] <- 0
  }
  
  splits <- utreec2splits(utreec)
  spl.label <- -1
  for(i in 1:dim(splits)[1]){
    if(sum(abs(splits[i,]-spl.of.int))==0){
      spl.label <- ntaxa+i-1
      break
    }
  }
  if(spl.label<0) stop("split not found")
  return(spl.label)
}
split.leading1 <- function(spls){
  for(spl in 1:dim(spls)[1])
    if (spls[spl,1] == 0){
      idx <- spls[spl,] == 0
      spls[spl,idx] <- 1
      spls[spl,!idx] <- 0
    }
  return(spls)
}
utreec2splits <- function(utreec){
  # OUT:
  # splits - (ntaxa-3)x ntaxa matrix of 1 and 0; 1 for taxa 1 
  ntaxa <- dim(utreec)[1]+1
  splits <- matrix(0,ncol=ntaxa,nrow=ntaxa-3)
  
  for(split in 1:(ntaxa-3))
    for(k in 1:2){
      b <- utreec[split,k]
      if (b < ntaxa)
        splits[split,(b+1)] <- 1
      else
        splits[split,] <- splits[split,] + splits[(b-ntaxa+1),]
    }
  
  splits <- split.leading1(splits)
  
  return(splits)
}

##-------------------------------------

# root_tree(utreefile,taxa)
# reroots a tree given a list of taxa on one side of the root
root_tree <- function(utreefile,taxa) {
  utc <- utree(utreefile)
  cmdline <- paste("./rert",LocateRoot(utc,taxa),dim(utc)[1]+1,"<",utreefile,"> tmp.root.utreec")
  system(cmdline)
  utc <- matrix(scan("tmp.root.utreec",quiet=TRUE),ncol=4,byrow=TRUE)
  utc <- force_root(utc,file=FALSE) # if rert didn't root the tree, force the root by splitting the last line in half
  did.remove <- file.remove("tmp.root.utreec")
  return(utc)
}

# force_root(utreec) 
# if a utreec matrix ends with an unrooted line (e.g. "x y z 0" or "x y 0 z"), replace "z" and "0" with z/2 to root it.
force_root <- function(utreec,file=TRUE) {
  if(file) {
    f <- file(utreec); l <- readLines(f)
    last_line <- unlist(strsplit(l[length(l)], " "))
  } else last_line <- utreec[nrow(utreec),]
  rewrite <- FALSE
  if(as.numeric(last_line[3]) == 0) {
    n <- as.numeric(last_line[4]); rewrite <- TRUE
  } else if(as.numeric(last_line[4]) == 0) {
    n <- as.numeric(last_line[3]); rewrite <- TRUE
  }
  if(rewrite){
    last_line[3] <- last_line[4] <- n/2
    if(file){
      last_line <- as.character(last_line)
      last_line <- paste(last_line,collapse=" ")
      l[length(l)] <- last_line
      writeLines(l,f); close(f)
    } else {
      utreec[nrow(utreec),] <- last_line
      return(utreec)
    }
  }
  if(file) return() else return(utreec)
}

# removes names for intermediate nodes (numbers, typically) from utreec name files
clean_namesfile <- function(namesfile,out) {
  f <- file(namesfile,"r")
  l <- vector()
  for(line in readLines(f)){
    if(substr(line,2,2) %in% c("0","1","2","3","4","5","6","7","8","9")){
      break
    } else {
      l <- append(l, line)
    }
  }
  close(f)
  f <- file(out,"w+")
  writeLines(l,f); close(f)
}

# make_namesfile(namesfile,remtaxa,outfile)
# generates a new names file for a subtree given the names file of a larger tree and the taxa not present in the subtree
make_namesfile <- function(namesfile,remtaxa,outfile) {
  l <- vector()
  f <- file(namesfile,"r")
  i <- 0
  for(line in readLines(f)){
    if(!(i %in% remtaxa)) {
      l <- append(l, line)
    }
    i <- i+1
  }
  close(f)
  file.create(paste0(outfile,".names"))
  f <- file(paste0(outfile,".names"))
  writeLines(l,f)
  close(f)
}

# split_tree(utreefile,namesfile,outfile,remtaxa)
# routine that splits and roots a subtree from a bigger tree (main function of rtsplit)
split_tree <- function(utreefile, namesfile, outfile, remtaxa){
  ut <- utree(utreefile)
  new.names <- paste0("new_",namesfile); clean_namesfile(namesfile,new.names)
  make_namesfile(new.names, remtaxa, outfile) # make a new names file for the subtree
  
  rem <- del.taxa.utreecc(ut,remtaxa) # remove the taxa
  
  file.create(outfile)
  f <- file(outfile)
  write.table(rem, file=outfile, row.names=FALSE, col.names=FALSE)
  close(f)
  
  # reroot subtree at the split point
  rem <- root_subtree(utreefile, outfile, new.names, paste0(outfile,".names"), remtaxa) 
  
  file.create(outfile)
  f <- file(outfile)
  write.table(rem, file=outfile, row.names=FALSE, col.names=FALSE)
  close(f)
}

# root_subtree(u_ab_file,u_sub_file,ab_names_file,sub_names_file,remtaxa) => u_sub
# subroutine called by split_tree that properly roots a subtree at the split point
root_subtree <- function(u_ab_file, u_sub_file, ab_names_file, sub_names_file, remtaxa) {
  u_ab <- utree(u_ab_file)
  ab_names <- unames(ab_names_file); sub_names <- unames(sub_names_file)
  ntaxa <- length(ab_names)
  keeptaxa <- setdiff(seq(0,ntaxa-1), remtaxa)
  dlist <- dlistf(u_ab)
  split.point <- which(lapply(dlist,function(x){return(setequal(x,keeptaxa))}) == TRUE)
  split.root <- dlist[[u_ab[split.point - ntaxa , 1] + 1]] # add 1 to index since R starts counting at 1 instead of 0
  split.names <- unlist(lapply(split.root, function(x){return(ab_names[x+1])}))
  split.root <- unlist(lapply(which(sub_names %in% split.names), function(x){return(x-1)}))
  
  u_sub <- root_tree(u_sub_file, split.root)
  return(u_sub)
}

# assumes the input tree is rooted, and outputs a file containing taxa that represent the split of that tree's root
make_root_split <- function(tree,outfile,newick=FALSE,seqfile=NULL) {
  if(newick){
    if(is.null(seqfile)) stop("make_root_split: seqfile not provided for newick")
    tmp.tree.name <- paste0("tmp.tmp.",outfile)
    write(tree,tmp.tree.name)
    make.utree(tmp.tree.name,seqfile,"tmp.tmp.tree","tmp.tmp.names")
    tree <- force_root(utree("tmp.tmp.tree"),file=FALSE)
    did.remove <- file.remove(c(tmp.tree.name,"tmp.tmp.tree","tmp.tmp.names"))
  }
  dlist <- dlistf(tree)
  roothalf <- tree[nrow(tree),1] + 1
  split <- dlist[[roothalf]]
  write(split,file=outfile,ncolumns=length(split))
}

# function used by rtsplit to determine whether the root of tree AB is located in subtree A or B 
locate_split_root <- function(utreec,root,split) {
  ntaxa <- nrow(utreec)+1
  all <- seq(0,ntaxa-1)
  root_c <- setdiff(all,root); split_c <- setdiff(all,split)
  split <- sort(split); split_c <- sort(split_c) # ensure the splits are from least -> greatest
  A.contained <- (length(setdiff(split,root)) == 0 | length(setdiff(split,root_c)) == 0)
  B.contained <- (length(setdiff(split_c,root)) == 0 | length(setdiff(split_c,root_c)) == 0)
  
  # determine where the tree's root is located : either inside tree A, tree B, or at the same location as the split
  if(A.contained & !B.contained)
    rootloc <- "B"
  else if(!A.contained & B.contained)
    rootloc <- "A"
  else
    return(list("0",NULL)) # root same as split, no changes needed
  
  # remove all target taxa except for one "ghost" taxon which has no edge length but maintains the structure
  if(rootloc == "A") {
    ghost <- split_c[1]
    subtree <- del.taxa.utreecc(utreec,split_c[-1],root=TRUE)
  }
  else {
    ghost <- split[1]
    subtree <- del.taxa.utreecc(utreec,split[-1],root=TRUE)
  }
  
  # alter edge length for ghost taxon
  subtree[which(subtree == ghost,arr.ind=TRUE)[1],] <- c()
  
  return(list(rootloc,subtree))
}

# dlistf(ut) => dl
# gives a list containing the descendant taxa of each node in a tree
dlistf <- function(ut){
  #OUT: dl - list dl[[b]]: descendents of branch b
  nt <- dim(ut)[1]+1
  dl <- vector("list",(2*nt-1))
  for(b in 1:nt) dl[[b]] <- b
  for(j in 1:(nt-1))
    dl[[(j+nt)]] <- c(dl[[(ut[j,1]+1)]], dl[[(ut[j,2]+1)]])
  for(j in 1:(2*nt-1)) dl[[j]] <- dl[[j]]-1
  return(dl)
}
