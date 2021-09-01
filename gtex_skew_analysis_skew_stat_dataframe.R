

#Compiling the XCI skew related stats into a single dataframe

#Follows the bootstrapping script
#/home/werner/xchrom_snp_skew/code/bootstrapping_script.R



#For getting the original skew and sigma estimates

get_boot_mu = function(boot){
  if(length(boot)==1){   #If the bootstrap had an NA, length of the object is 1
    return(NA)
  }else{return(boot$t0[1])}
}
get_boot_sig = function(boot){
  if(length(boot)==1){   #If the bootstrap had an NA, length of the object is 1
    return(NA)
  }else{return(boot$t0[2])}
}


mus = sapply(boots, get_boot_mu)
sigmas = sapply(boots, get_boot_sig)



get_CI = function(boot){
  if(length(boot)==1){
    return(c(NA,NA))
  }else{
    mu_distribution = boot$t[ ,1]
    mu_distribution = sort(mu_distribution)
    
    #For 200 bootstraps, 2.5% of the distribution is below the 6th index, 2.5% is above the 195th index
    lower_CI = mu_distribution[6]
    upper_CI = mu_distribution[195]
    return(c(lower_CI, upper_CI))
  }
}


mu_CI = lapply(boots, get_CI)


#Get the lower and upper bounds for the CI
mu_lower_CI = sapply(mu_CI, '[[', 1)
mu_upper_CI = sapply(mu_CI, '[[', 2)
mu_CI_width = mu_upper_CI - mu_lower_CI
#For plotting
mu_lower_CI = mus - mu_lower_CI
mu_upper_CI = mu_upper_CI - mus



index_55 = which(mus < .55)
index_6 = which(mus >= .55 & mus < .6)
index_65 = which(mus >= .6 & mus < .65)
index_7 = which(mus >= .65 & mus < .7)
index_75 = which(mus >= .7 & mus < .75)
index_8 = which(mus >= .75 & mus < .8)
index_85 = which(mus >= .8 & mus < .85)
index_9 = which(mus >= .85 & mus < .9)
index_95 = which(mus >= .9 & mus < .95)

temp_list_skew_max = list.skew.max[index_55]
gene_skews_55 = unlist(mclapply(1:length(temp_list_skew_max), function(i) get_gene_skews(temp_list_skew_max[[i]]), mc.cores = 5))
labels_55 = rep('Tissue skew <.55', length(gene_skews_55))

temp_list_skew_max = list.skew.max[index_6]
gene_skews_6 = unlist(mclapply(1:length(temp_list_skew_max), function(i) get_gene_skews(temp_list_skew_max[[i]]), mc.cores = 5))
labels_6 = rep('.55 >= Tissue skew <.6', length(gene_skews_6))

temp_list_skew_max = list.skew.max[index_65]
gene_skews_65 = unlist(mclapply(1:length(temp_list_skew_max), function(i) get_gene_skews(temp_list_skew_max[[i]]), mc.cores = 5))
labels_65 = rep('.6 >= Tissue skew <.65', length(gene_skews_65))

temp_list_skew_max = list.skew.max[index_7]
gene_skews_7 = unlist(mclapply(1:length(temp_list_skew_max), function(i) get_gene_skews(temp_list_skew_max[[i]]), mc.cores = 5))
labels_7 = rep('.65 >= Tissue skew <.7', length(gene_skews_7))

temp_list_skew_max = list.skew.max[index_75]
gene_skews_75 = unlist(mclapply(1:length(temp_list_skew_max), function(i) get_gene_skews(temp_list_skew_max[[i]]), mc.cores = 5))
labels_75 = rep('.7 >= Tissue skew <.75', length(gene_skews_75))

temp_list_skew_max = list.skew.max[index_8]
gene_skews_8 = unlist(mclapply(1:length(temp_list_skew_max), function(i) get_gene_skews(temp_list_skew_max[[i]]), mc.cores = 5))
labels_8 = rep('.75 >= Tissue skew <.8', length(gene_skews_8))

temp_list_skew_max = list.skew.max[index_85]
gene_skews_85 = unlist(mclapply(1:length(temp_list_skew_max), function(i) get_gene_skews(temp_list_skew_max[[i]]), mc.cores = 5))
labels_85 = rep('.8 >= Tissue skew <.85', length(gene_skews_85))

temp_list_skew_max = list.skew.max[index_9]
gene_skews_9 = unlist(mclapply(1:length(temp_list_skew_max), function(i) get_gene_skews(temp_list_skew_max[[i]]), mc.cores = 5))
labels_9 = rep('.85 >= Tissue skew <.9', length(gene_skews_9))

temp_list_skew_max = list.skew.max[index_95]
gene_skews_95 = unlist(mclapply(1:length(temp_list_skew_max), function(i) get_gene_skews(temp_list_skew_max[[i]]), mc.cores = 5))
labels_95 = rep('.9 >= Tissue skew <.95', length(gene_skews_95))

rm(temp_list_skew_max)


unphased_gene_skew_df = data.frame(gene_skews = c(gene_skews_55, gene_skews_6, gene_skews_65, gene_skews_7, gene_skews_75, gene_skews_8, gene_skews_85, gene_skews_9), 
                                   labels = c(labels_55, labels_6, labels_65, labels_7, labels_75, labels_8, labels_85, labels_9))
unphased_gene_skew_df$labels = factor(unphased_gene_skew_df$labels, levels = rev(c('Tissue skew <.55', '.55 >= Tissue skew <.6', '.6 >= Tissue skew <.65','.65 >= Tissue skew <.7',
                                                                                   '.7 >= Tissue skew <.75','.75 >= Tissue skew <.8', '.8 >= Tissue skew <.85', '.85 >= Tissue skew <.9')))
##############
#Figure 2 panel f
################

ggplot(unphased_gene_skew_df, aes(x = gene_skews, y = labels)) + 
  #stat_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha = 0.7, rel_min_height = 0.005, scale=2 ) +
  geom_density_ridges(scale=2, aes(y = labels), rel_min_height = 0.005) +
  xlim(.5,1) + 
  xlab('Gene skews')+ ylab('Tissue skews') + 
  theme(plot.title=element_text(hjust=.5, size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12),
        legend.position = 'none',
        panel.background = element_rect(fill = 'white'))





#KS-test for the pfold norm fit

library(VGAM)

#Takes in the list of max SNPs, and vectors of the calculated mu and sigma MLE estimates
run_ks_test = function(max_snps, mus, sigmas){
  ks_tests = list()
  
  for(i in 1:length(max_snps)){
    
    if(Ns[i] < 2){ ks_tests[[i]] = NA; next}
    
    
    #high.f = which(max_snps[[i]]$A.1 < 3000)  #High expressed SNPs filter
    skews = folded(max_snps[[i]]$C.1)
    mu = mus[i]
    sig = sigmas[i]
    suppressWarnings( {test = ks.test(skews, y='pfoldnorm', mean=mu, sd=sig)} )
    ks_tests[[i]] = c(test$statistic, test$p.value)
    
  }  
  
  return(ks_tests)
}

ks_tests = run_ks_test(list.skew.max, mus, sigmas)

get_p_vals = function(test){
  if(length(test)==0){return(NA)}
  return(test[2])
}



#pulling together a dataframe with all the combined stats

f = which(!is.na(mus)) #All of the samples that have a mu estimate, they have at least 2 SNPs
skew_and_stats_df = data.frame(sample_index = f, skew = mus[f], skew_sigma = sigmas[f], CI_width = mu_CI_width[f], num_snps=Ns[f], tissue=sample_meta$body_site_s[f], donor=sample_meta$Donor_id[f], germ_layer=sample_meta$germ.layer[f])

skew_and_stats_df = skew_and_stats_df[order(skew_and_stats_df$CI_width), ]
skew_and_stats_df$CI_rank = 1:dim(skew_and_stats_df)[1]
#Reorder on sample index
skew_and_stats_df = skew_and_stats_df[order(skew_and_stats_df$sample_index), ]


#Add ks-test pvals and do BH corrections
p_vals = sapply(ks_tests, get_p_vals)

#Add the p_vals
skew_and_stats_df$ks_pvals = p_vals[f]

#Adjust the pvalues
skew_and_stats_df$ks_pvals = p.adjust(skew_and_stats_df$ks_pvals, method = 'BH')


#Order on the p_vals to assign ranks
skew_and_stats_df = skew_and_stats_df[order(skew_and_stats_df$ks_pvals), ]      #Order on the p values
skew_and_stats_df$p_val_rank = 1:dim(skew_and_stats_df)[1]                      #Assign ranks
#Reorder
skew_and_stats_df = skew_and_stats_df[order(skew_and_stats_df$sample_index), ]

skew_and_stats_df$p_significance = skew_and_stats_df$ks_pvals <= .05

save(skew_and_stats_df, file = "/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/GATK_mod/v8_final_filtering_skew_and_stats_1_20_21_escape_filt_df.Rdata")

