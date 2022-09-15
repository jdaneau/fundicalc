######
# convert_tree is a collection of functions for interfacing between newick and utreec formats

# plot_utree_file(utreefile,namesfile)
# creates a visual plot of a utreec from a file
plot_utree_file <- function(utreefile,namesfile){
  plot(read.tree(text=ut2n(utree(utreefile), unames(namesfile))))
}

# plot_utree(utreec,names)
# creates a visual plot of a utreec matrix
plot_utree <- function(utreec,names){
  plot(read.tree(text=ut2n(utreec,names)))
}

# ut2n(utreec,names,blabel) => tstring
# converts a utreec matrix to a newick tree string 
ut2n <- function(utreec, names, blabel){
  # blabel - (optional) labels for 2*ntaxa-3 branches. Can be "" for terminal
  #          and BP
  #          
  ntaxa <- dim(utreec)[1]+1
  tstring <- names
  ## Strange. Below works although the indices of tstring exceed its length
  for(join in 1:(ntaxa-3)){
    r <- utreec[join,1]+1; l <- utreec[join,2]+1
    lr <- utreec[join,3]; ll <- utreec[join,4]
    tstring[(join+ntaxa)] <- paste("(",tstring[r],":",lr,",",tstring[l],":",
                                   ll,")", sep="")
    if (!missing(blabel))
      tstring[(join+ntaxa)] <- paste(tstring[(join+ntaxa)],
                                     blabel[(join+ntaxa)],sep="")
  }
  join <- ntaxa-2
  r <- utreec[join,1]+1; l <- utreec[join,2]+1
  lr <- utreec[join,3]; ll <- utreec[join,4]
  join <- ntaxa-1
  rr <- utreec[join,1]+1; lrr <- utreec[join,3] + utreec[join,4]
  tstring <- paste("(",tstring[rr],":",lrr,",",tstring[r],":",lr,",",
                   tstring[l],":",ll,");",sep="")
  return(tstring)
}

# ut2nrt(utreec,names,blabel) => tstring
# converts a utreec matrix to a *rooted* newick tree string 
ut2nrt <- function(utreec,names,blabel){
  ntaxa <- dim(utreec)[1]+1
  tstring <- names
  for(join in 1:(ntaxa-1)){
    r <- utreec[join,1]+1; l <- utreec[join,2]+1
    lr <- utreec[join,3]; ll <- utreec[join,4]
    tstring <- c(tstring,
                 paste0("(",tstring[r],":",lr,",",tstring[l],":",ll,")"))
    if (!missing(blabel))
      tstring[(join+ntaxa)] <- paste(tstring[(join+ntaxa)],
                                     blabel[(join+ntaxa)],sep="")
    ## print(tstring)
  }
  return(paste0(tstring[(ntaxa-1+ntaxa)],";"))
}

# del.taxa.utreecc(utreec,taxa) => utreec
# removes a set of taxa from a utreec and returns the result
del.taxa.utreecc <- function(utreec, taxa){
  write(t(utreec),ncol=4,file="tmp.del.utreec")
  taxa <- -sort(-taxa)
  for(tc in taxa){
    cmdline <- paste("./del-taxa-utreec -d",tc,"-u tmp.del.utreec > tmp.utreecn")
    system(cmdline)
    system("mv tmp.utreecn tmp.del.utreec")
    utreec <- matrix(scan("tmp.del.utreec",quiet=TRUE),ncol=4,byrow=TRUE)
  }
  utreec <- matrix(scan("tmp.del.utreec",quiet=TRUE),ncol=4,byrow=TRUE)
  system("rm tmp.del.utreec")
  return(utreec)
}

# make.utree(treefile,seqfile,out.utreefile,out.namesfile)
# creates a utreec matrix from a newick tree file (must be a file, since it relies on an external binary call)
make.utree <- function(treefile,seqfile,out.utreefile,out.namesfile) {
  cmd <- paste("./treecns",treefile,out.utreefile,"-1 <",seqfile,">",out.namesfile)
  system(cmd)
}

# utree(utreecfile)
# creates a utreec matrix from a utreec file
utree <- function(utreecfile){
  if("matrix" %in% class(utreecfile)) return(utreecfile) # if provided utree matrix, just return it back
  return(matrix(scan(utreecfile,quiet=TRUE),ncol=4,byrow=TRUE))
}

# unames(namefile)
# creates a utreec names matrix from a names file
unames <- function(namefile){
  if("matrix" %in% class(namefile)) return(namefile) # if provided names matrix, just return it back
  return(scan(namefile,what=character(),quiet=TRUE))
}
