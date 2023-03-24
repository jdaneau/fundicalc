nsite <- 1000
fp <- file("sitellfile","rb")
ll <- readBin(fp,"double",nsite)
close(fp)

file.create("Bsitellfile.txt")
f <- file("Bsitellfile.txt")
writeLines(as.character(ll),f)
close(f)