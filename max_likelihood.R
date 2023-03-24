# Read in files
f <- file("ABsitellfile.txt")
lAB <- scan(f, what=double()); close(f)
f <- file("Asitellfile.txt")
lA <- scan(f, what=double()); close(f)
f <- file("Bsitellfile.txt")
lB <- scan(f, what=double()); close(f)

n <- length(lAB)

# Initialize values
p <- 0.5
pn <- 0
thresh <- 1.0e-6
prev.sum <- -Inf
sums <- vector()
lmax <- vector(length=n)
lf <- vector(length=n)

update_lf <- function(){
  for(i in 1:n) {
    lmax[i] <<- max(lA[i]+lB[i], lAB[i])
    lf[i] <<- lmax[i] + log(p*exp(lA[i]+lB[i]-lmax[i]) + (1-p)*exp(lAB[i]-lmax[i]) )
  }
  s <- sum(lf)
  if(s < prev.sum) stop(paste("lf sum did not increase. something bad happened:",prev.sum,">",s))
  prev.sum <<- s
  sums <<- append(sums, s)
}
update_lf()

mean_i <- function() {
  sum <- 0
  for(i in 1:n) {
    sum <- sum + 
      ( p*exp(lA[i]+lB[i] - lmax[i]) / exp(lf[i] - lmax[i]) )
  }
  return(sum/n)
}

# Main loop
while(abs(pn-p) > thresh) {
  pn <- mean_i()
  p <- pn
  update_lf()
  p <- mean_i()
}

plot(1:length(sums),sums)
print(paste0("p = ",p,"; lf = ",prev.sum))