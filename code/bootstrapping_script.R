

#Setting this up to run remotely
#Need the filtered max SNPs, just the ratios.max object from the following Rdata file
library(boot)
#for parallelizing
library(parallel)
library(MASS)

load('../data/all_v8_GTEx_gene_filtered.skew.est.max.genes.Rdata')


rm(ratios.max, Ns)

p = 50:100/100


# Functions to fold and unfold ratios 
folded <- function(x) { apply( cbind(x,1-x), 1, max)} 
unfold <- function(x) { c(x,1-x) } 

mle_folded <- function(x){ 
  mus = seq(0.5,1, by = 0.001)
  sigmas = seq(0,0.5, by = 0.01)
  
  mles = sapply(mus, function(mu)
    sapply(sigmas, function(sigma)
      -sum( log(dnorm( x,   mean = mu, sd = sigma )
                + dnorm( x,   mean = 1-mu, sd = sigma ) )  )))
  mles[!is.finite(mles)] = NA
  coefs = which(mles==min(mles, na.rm=T), arr.ind=T)
  return (list( mus[coefs[2]] , sigmas[coefs[1]]))
} 


boot_skew = function(data, indices){
  
  x = data[indices]           #For the bootstrap sampling
  x = folded(x)               #Fold for the MLE
  params = mle_folded(x)       #MLE the folded norm
  return(c(params[[1]], params[[2]]))
}


#function for the lapply()
parallel_boot = function(max_ratio){

  sample = max_ratio$C.1 
  
  if(length(sample) == 0 | length(sample) == 1){return(NA) }   #NULL samples, and samples with 1 or 0 SNPs
  #else if(length(sample)==1){return(NA)}  #Only 1 SNP in the sample
  
  else{
    sample = as.numeric(sample)
    sample_boot = boot(data=sample, statistic=boot_skew, R=200)
    return(sample_boot)
  }
}


numCores = 25
print(sprintf('Starting bootstap method, time is: %s', Sys.time()))

boots = mclapply(list.skew.max, parallel_boot, mc.cores = numCores)
print(sprintf('Ending bootstap method, time is: %s', Sys.time()))


save(boots, file = '../data/bootstrap_XCI_skews_all_v8_GTEx.Rdata')


