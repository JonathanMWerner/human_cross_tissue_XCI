

#Setting this up to run remotely
#Need the filtered max SNPs, just the ratios.max object from the following Rdata file
library(boot)
#for parallelizing
library(parallel)
library(MASS)

load('/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/GATK_mod/all_v8_GATK_1_20_2021_update_filt_escape.skew.est.max.genes.Rdata')


rm(ratios.max, Ns)

p = 50:100/100


# Functions to fold and unfold ratios 
folded <- function(x) { apply( cbind(x,1-x), 1, max)} 
unfold <- function(x) { c(x,1-x) } 

#Sara's MLE function
# Maximum likelihood estimate function. Note, x is a global variable! Should fix but this lets me use it in the bootstrap later.  
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



#data is the unfolded reference SNP ratios for a sample
#p is the skews to test
#indicies is needed by the boot() to generate bootstrap samples
boot_skew = function(data, indices){
  
  x = data[indices]           #For the bootstrap sampling
  x = folded(x)               #Fold for the MLE
  params = mle_folded(x)       #MLE the folded norm
  return(c(params[[1]], params[[2]]))
}


#function for the lapply()
parallel_boot = function(max_ratio){

  #high.f = which(max_ratio$A.1 < 3000)      #Filtering out highly expressed SNPs
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


save(boots, file = '/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/GATK_mod/bootstrap_results_all_v8_GATK_1_20_2021_update_filt_escape.Rdata')


