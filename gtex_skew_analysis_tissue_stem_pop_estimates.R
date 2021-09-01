
#Estimating the stem cell population sizes for each tissue using the binomial model of XCI

#Run after gtex_skew_analysis_setting_up.R


#Fitting the distributions really centers on the variance of the distributions
#Making the claim that var = p* q / n   under the binomial model



binomial_fitted_est_n = vector(mode='numeric', length = length(filt_tissues))         #Estimated cell number from the binomial model         
binomial_fitted_error = vector(mode='numeric', length = length(filt_tissues))         #Error associated with that model
binomial_variance_explained = vector(mode='numeric', length = length(filt_tissues))   #percentage of variance explained using that model

normal_fitted_est_n = vector(mode='numeric', length = length(filt_tissues))           #Estimated cell number using the normal model, essentially a smoothened binomial
normal_fitted_error = vector(mode='numeric', length = length(filt_tissues))           #Error associated with that model
normal_variance_explained = vector(mode='numeric', length = length(filt_tissues))     #Variance explained
num_samples_per_tiss= vector(mode='numeric', length = length(filt_tissues))           #Number of samples(donors) per tissue

for(k in 1:length(filt_tissues)){
  
  p = .5    #It's known there is equal chance for either X to inactiavte
  n = 2:50  #Testing these different values of N to find the best fit
  percentiles = seq(.01, .99, 0.01)  #Fitting over the percentiles of the theoretical vs empirical data
  
  bin_error_vector = vector(mode = 'numeric', length = length(n))  #Error for each model based off of N
  norm_error_vector = vector(mode = 'numeric', length = length(n))
  
  #Get the empirical skews for that tissue
  empir_skews = unfold(skew_and_stats_df$skew[skew_and_stats_df$tissue == filt_tissues[k]])
  num_skews = length(empir_skews) 
  num_samples_per_tiss[k] = num_skews / 2 #Number of unfolded skews is twice that of the sample size for that tissue
  
  #Get the quantiles for these skews
  empir_quants = quantile(empir_skews, probs = percentiles, type=1)
  #Only going to be fitting over the confident skews
  quant_filt = empir_quants <= .4 | empir_quants >= .6
  
  for(i in 1:length(n)){
    #Get the normal and binomial percentiles for the n parameter
    norm_sd = sqrt( (p*(1-p)) / n[i] )   
    norm_percentiles = qnorm(percentiles, mean=p, sd = norm_sd)
    bin_percentiles = qbinom(percentiles, size=n[i], prob=p) / n[i]
    #And the error between the empirical and theoretical percentiles
    bin_error_vector[i] = sum((empir_quants[quant_filt] - bin_percentiles[quant_filt])^2)
    norm_error_vector[i] = sum((empir_quants[quant_filt] - norm_percentiles[quant_filt])^2)
  }
  
  #Save the best fitting n
  binomial_fitted_est_n[k] = n[which.min(bin_error_vector)]
  normal_fitted_est_n[k] = n[which.min(norm_error_vector)]
  
  #Save the variance explained by the best model
  norm_sd = sqrt( (p*(1-p)) / normal_fitted_est_n[k] )
  norm_percentiles = qnorm(percentiles, mean=p, sd = norm_sd)
  #Assess variance explained across whole distribution, rather than the confident skews we fitted to
  error = empir_quants - norm_percentiles 
  total_sum_squares = sum( (empir_quants  - mean(empir_quants))^2 )
  sum_squares_error = sum(error^2)
  normal_variance_explained[k] = 1 - (sum_squares_error / total_sum_squares)
  
  #And for the binomial fits
  bin_percentiles = qbinom(percentiles, size=binomial_fitted_est_n[k], prob=p) / binomial_fitted_est_n[k]
  error = empir_quants - bin_percentiles 
  total_sum_squares = sum( (empir_quants  - mean(empir_quants))^2 )
  sum_squares_error = sum(error^2)
  binomial_variance_explained[k] = 1 - (sum_squares_error / total_sum_squares)  
}


#Also run Sara's original estimation method, using the variance formula exactly


cell_pop_est_vector = vector(mode = 'numeric', length = length(filt_tissues))

for( i in 1: length(filt_tissues)){
  
  skews = skew_and_stats_df$skew[skew_and_stats_df$tissue == filt_tissues[i]]
  skew_var = var(unfold(skews))
  estimated_n = .25/skew_var
  cell_pop_est_vector[i] = estimated_n 
}


#Make a dataframe containing the results

#Get the germlayers for the filtered tissues
germ_layers = vector(mode='character', length=length(filt_tissues))

for( i in 1:length(filt_tissues)){
  index = which(skew_and_stats_df$tissue == filt_tissues[i])
  germ_layers[i]= as.character(skew_and_stats_df$germ_layer[index[1]])
  
}

est_stem_cell_pop_df = data.frame(binomial_est_n = binomial_fitted_est_n, normal_est_n = normal_fitted_est_n, var_est_n = cell_pop_est_vector, 
                                  tissues = filt_tissues, num_samples = num_samples_per_tiss, germ_layers = germ_layers)
rownames(est_stem_cell_pop_df) = filt_tissues
#Group by germ layer
est_stem_cell_pop_df = arrange(est_stem_cell_pop_df, germ_layers)
#And reset the factor levels for the tissues
est_stem_cell_pop_df$tissues = factor(est_stem_cell_pop_df$tissues, levels = est_stem_cell_pop_df$tissues)


##########################################
#Plot
##########################################



p_1= ggplot(est_stem_cell_pop_df, aes(x=tissues, y=binomial_est_n, color=germ_layers, size=num_samples)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title=element_text(hjust=.5), plot.margin=margin(.1,.1,.01,.1,'in'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(face="bold", size=12)) + 
  ggtitle('Tissue Stem cell population size estimates - fitting binomial percentiles') + ylab('Estimated Num of cells') +
  xlab('Tissues') + ylim(0,25) + 
  scale_color_manual(breaks = c('Ectoderm','Endoderm','Mesoderm'), values=c(ecto_color, endo_color, meso_color)) + 
  geom_hline(yintercept=16, linetype='dashed', color='blue') +
  geom_hline(yintercept=8, linetype='dashed', color='blue')

p_2 = ggplot(est_stem_cell_pop_df, aes(x=tissues, y=normal_est_n, color=germ_layers, size=num_samples)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title=element_text(hjust=.5), plot.margin=margin(.1,.1,.01,.1,'in'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(face="bold", size=12)) + 
  ggtitle('Tissue Stem cell population size estimates - fitting normal percentiles') + ylab('Estimated Num of cells') +
  xlab('Tissues')+ ylim(0,25) +
  scale_color_manual(breaks = c('Ectoderm','Endoderm','Mesoderm'), values=c(ecto_color, endo_color, meso_color)) + 
  geom_hline(yintercept=16, linetype='dashed', color='blue') +
  geom_hline(yintercept=8, linetype='dashed', color='blue')


p_3 = ggplot(est_stem_cell_pop_df, aes(x=tissues, y=var_est_n, color=germ_layers, size=num_samples)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title=element_text(hjust=.5), plot.margin=margin(.1,.1,.01,.1,'in'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(face="bold", size=12)) + 
  ggtitle('Tissue Stem cell population size estimates - using variance estimate') + ylab('Estimated Num of cells') +
  xlab('Tissues')+ ylim(0,25) + 
  scale_color_manual(breaks = c('Ectoderm','Endoderm','Mesoderm'), values=c(ecto_color, endo_color, meso_color)) + 
  geom_hline(yintercept=16, linetype='dashed', color='blue') +
  geom_hline(yintercept=8, linetype='dashed', color='blue')


p_1
p_2
p_3


#Plotting the best fitted distribution over the empirical data


for(j in 1:length(filt_tissues)){
  
  current_tissue = filt_tissues[j]
  fitted_normal_n = est_stem_cell_pop_df$normal_est_n[est_stem_cell_pop_df$tissues == current_tissue]
  p = .5
  sigma = sqrt(p*(1-p)/fitted_normal_n)
  dnorm_x = 0:1000/1000
  
  fitted_normal = dnorm(dnorm_x, mean=p, sd = sigma)
  
  theoretical_norms = list()
  n = c(4,8,16,32)
  for(i in 1:length(n)){ 
    sigma = sqrt(p*(1-p)/n[i])
    theoretical_norms[[i]] = dnorm(dnorm_x, mean=p, sd = sigma)
  }
  
  
  skews = unfold(skew_and_stats_df$skew[skew_and_stats_df$tissue == current_tissue])
  #Get the densities for each theoretical distribution
  theo_dens = lapply(theoretical_norms, density)
  theo_dens = sapply(theo_dens, '[[', 2)
  theo_maxes = apply(theo_dens,2, max)
  y_max = max( c(theo_maxes, max(density(skews)$y))) #Find the maximum density out of all the distributions
  
  #Fixing tissue name for file name
  tiss_name = current_tissue
  tiss_name = gsub(' ', '_', tiss_name, fixed = TRUE)
  tiss_name = gsub('-', '_', tiss_name, fixed = TRUE)
  
  file_name = sprintf('/home/werner/xchrom_snp_skew/code/graphs/stem_cell_pop_estimates/tissue_specific_fits/%s_normal_fits.pdf', tiss_name)
  #pdf(file = file_name, height = 6, width = 6, useDingbats = FALSE)
  
  hist(skews, main = sprintf('XCI skews: %s', current_tissue), xlim = c(0,1), breaks = 64, freq = FALSE)
  lines(dnorm_x,fitted_normal, lwd = 5, col = 'black')
  
  line_colors = c('tan1','steelblue4','turquoise3','slateblue3')
  for(i in 1:length(n)){ 
    lines(dnorm_x,theoretical_norms[[i]], col = line_colors[i], lwd=2)
  }
  
  legend(x = 'topright', legend=c('4 cells','8 cells', '16 cells', '32 cells', sprintf('Estimated %i cells', fitted_normal_n)), 
         col = c(line_colors, 'black'), lwd = c(2,2,2,2,5))
  
  #dev.off()
}




##########################################
#End plot
##########################################





























