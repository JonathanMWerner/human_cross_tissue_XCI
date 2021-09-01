#Investigating reference bias and gene filtering


load("/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/GATK_mod/all_v8_GATK_mod_stat_filt_1_20_21.skew.est.max.genes.Rdata")



final_filt.list.skew.max = list.skew.max
rm(list.skew, ratios.max.genes)


#Samples that are or are not null

skip = c()


for(i in 1:length(final_filt.list.skew.max)){
  if(is.null(final_filt.list.skew.max[[i]])){skip = c(skip, i)}
}


n = 1:length(final_filt.list.skew.max)
not_skip = n[!n %in% skip]



#look at the entire called SNP distribution before additional filtering

ref_ratios = vector(mode='numeric')
counts = vector(mode='numeric')

for(i in 1:length(final_filt.list.skew.max)){
  
  if(i %in% skip){next}
  
  ref_ratios_1 = final_filt.list.skew.max[[i]]$C.1
  ref_ratios = c(ref_ratios, ref_ratios_1)
  counts_1 = final_filt.list.skew.max[[i]]$A.1
  counts = c(counts, counts_1)
  
}



###################
#Supplemental Figure 2 panel c

mod_df = data.frame(ref_ratios = ref_ratios, tot_counts = counts)
mod_ratios_p = ggplot(mod_df, aes(ref_ratios, log10(tot_counts))) + geom_point(alpha=1/8,size=1) + ylab('log SNP reads') + xlab('SNP ref/total ratio') + ggtitle('SNP ratios vs Reads: after GATK stat filtering') + 
  theme(plot.title=element_text(hjust=.5, size=14), axis.title.x = element_text(size=14), axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) + 
  geom_vline(xintercept=.5, color='red')
ggMarginal(mod_ratios_p, type='density',fill = '#752A8C')


#Filter out the highly expressed genes

for(i in 1:length(final_filt.list.skew.max)){
  
  if(is.null(final_filt.list.skew.max[[i]])){next}
  
  high_filt = final_filt.list.skew.max[[i]]$A.1 < 3000
  
  final_filt.list.skew.max[[i]] = final_filt.list.skew.max[[i]][high_filt, ]
  
}


#Pool together the gene specific ref SNP ration distributions, calculate Pearson's coefficient of skews, rank them and see what I get
#First get a list of all the genes in the dataset


genes = c()

for(i in 1: length(final_filt.list.skew.max)){
  if(is.null(final_filt.list.skew.max[[i]])){next}
  genes = c(genes, as.character(final_filt.list.skew.max[[i]]$name))
  
}

genes = unique(genes)
length(genes)

#Go through the dataset for each gene, building a data frame with all the SNP info for that gene
#Function to get the gene specific data from a sample
gene_filt = function(gene_name, data){
  filt = data$name == gene_name
  return(data[filt, ])
}

#To hold all the gene specific dataframes
compiled_gene_data = list()

for(i in 1:length(genes)){
  print(sprintf('Starting to process: %s %s', genes[i], Sys.time() ))
  compiled_gene_data[[i]] = data.frame(rbindlist(lapply(not_skip, function(x) gene_filt(genes[i], final_filt.list.skew.max[[x]]))))
}




#Get the skew coefficient
ref_skew = sapply(1:length(compiled_gene_data), function(i) skewness(compiled_gene_data[[i]]$C.1, type=3))

names(ref_skew) = genes
ordered_skew_gene_names = names(ref_skew[order(ref_skew)])


#Order genes by their skew coefficients and then
#only look at genes that are present in at least 30 samples


#Get the number of samples a gene is present in
num_samps_per_gene = sapply(1:length(compiled_gene_data), function(i) dim(compiled_gene_data[[i]])[1])


#Filter out the genes with greater than 30 samples
filt = num_samps_per_gene > 30
filt_gene_names = genes[filt]
low_samp_gene_names = genes[!filt]
temp_order_genes = ordered_skew_gene_names[ ordered_skew_gene_names %in% filt_gene_names]


#Order by skew coefficient
index = match(temp_order_genes, genes)
test = compiled_gene_data
test = test[index]


index = match(ordered_skew_gene_names, genes)
skew_ordered_gene_data = compiled_gene_data[index]



#############################
#Perform the diptest, stat test for multi modality. 
#Easiest thing to do is just exclude genes that disprove the null (unimodality)
#############################

library(diptest)


dip_test = mclapply(1:length(compiled_gene_data), function(i) dip.test(compiled_gene_data[[i]]$C.1, simulate.p.value=TRUE), mc.cores = 15)

dip_p_vals = sapply(1:length(dip_test), function(i) dip_test[[i]]$p.value)
names(dip_p_vals) = genes
bimodal_genes =  names(dip_p_vals[dip_p_vals < .05])

########################################
#Or just rank by the gene's skew distribution mean deviance from .5.


gene_mean_skews = vector(mode = 'numeric', length = length(compiled_gene_data))
for(i in 1:length(compiled_gene_data)){
  data = compiled_gene_data[[i]]$C.1
  gene_mean_skews[i] = mean(data)
}

gene_mean_skews_dev = gene_mean_skews - .5
names(gene_mean_skews_dev) = genes

#Only look at the genes that have at least 30 samples
test_mean_skew_devs = gene_mean_skews_dev[names(gene_mean_skews_dev) %in% filt_gene_names]
length(test_mean_skew_devs)
hist(test_mean_skew_devs, main = 'Deviation from .5 for gene specific mean skews')

#Order by the deviances and chop off the top and bottom 5% of genes
test_mean_skew_devs = test_mean_skew_devs[order(test_mean_skew_devs)]
lower_bound = ceiling(length(test_mean_skew_devs)*.05)
upper_bound = ceiling(length(test_mean_skew_devs)*.95)

######################
#Supplemental figure 2 panel a
hist(test_mean_skew_devs, main = ' ', breaks = 50, xlab = 'Deviation from .5 for mean gene-specific reference expression ratios', ylab = 'Number of genes')
abline(v = test_mean_skew_devs[lower_bound], col = 'red', lwd = 2)
abline(v = test_mean_skew_devs[upper_bound], col = 'red', lwd = 2)


filt_mean_skew_devs = test_mean_skew_devs[ceiling(length(test_mean_skew_devs)*.05):ceiling(length(test_mean_skew_devs)*.95)]
keep_genes = c(low_samp_gene_names, names(filt_mean_skew_devs ))


#Get rid of the bimodal genes from the keep_gene list
keep_genes = keep_genes[! keep_genes %in% bimodal_genes]
bad_genes = genes[!genes %in% keep_genes]


#Also exclude known escape genes, going off of the Tukainen classifications, which were determined from the GTEx data

tukiainen_escape_meta = read.table(file = '/home/werner/collabs/entex_Xchrom/code/tukiainen_study.txt', header = TRUE, sep = '\t')
tukiainen_escape_meta = tukiainen_escape_meta[1:683, ] #Random row at the end
escape_genes = as.character(tukiainen_escape_meta$Gene.name[as.character(tukiainen_escape_meta$Reported.XCI.status) == 'Escape'])

#And get rid of those escape genes

keep_genes = keep_genes[! keep_genes %in% escape_genes ]
#But keep XIST
keep_genes = c(keep_genes, 'XIST')
bad_genes = genes[!genes %in% keep_genes]

#Filtering out these SNPs from the dataset

filt_bad_genes = function(dataframe){
  
  gene_names = dataframe$name
  gene_filt = !gene_names %in% bad_genes
  
  dataframe = dataframe[gene_filt, ]
  return(dataframe)
}

final_filt.list.skew.max = lapply(final_filt.list.skew.max, filt_bad_genes)


#Samples that are or are not null



skip = c()


for(i in 1:length(final_filt.list.skew.max)){
  if(is.null(final_filt.list.skew.max[[i]])){skip = c(skip, i)}
}

length(skip)



#Plot the total SNP ref ratio distributions with the ref skewed genes filtered out


ref_ratios = vector(mode='numeric')
counts = vector(mode='numeric')

for(i in 1:length(final_filt.list.skew.max)){
  
  if(i %in% skip){next}
  
  ref_ratios_1 = final_filt.list.skew.max[[i]]$C.1
  ref_ratios = c(ref_ratios, ref_ratios_1)
  counts_1 = final_filt.list.skew.max[[i]]$A.1
  counts = c(counts, counts_1)
  
  
}

################
#Supplemental figure 2 panel d

mod_df = data.frame(ref_ratios = ref_ratios, tot_counts = counts)
mod_ratios_p = ggplot(mod_df, aes(ref_ratios, log10(tot_counts))) + geom_point(alpha=1/8,size=1) + ylab('log SNP reads') + xlab('SNP ref/total ratio') + ggtitle('SNP ratios vs Reads: after filtering ref biased genes') + 
  theme(plot.title=element_text(hjust=.5, size=14), axis.title.x = element_text(size=14), axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) + 
  geom_vline(xintercept=.5, color='red')
ggMarginal(mod_ratios_p, type='density',fill = '#752A8C')







############################################
#Skew modeling
###########################################




# Functions to fold and unfold ratios 
folded <- function(x) { apply( cbind(x,1-x), 1, max)} 
unfold <- function(x) { c(x,1-x) } 

# Maximum likelihood estimate function 
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



Ns = vector(mode='numeric', length = length(final_filt.list.skew.max))
for(i in 1:length(final_filt.list.skew.max)){
  
  if(is.null(final_filt.list.skew.max[[i]])){Ns[i] = 0}
  else{Ns[i] = dim(final_filt.list.skew.max[[i]])[1]}
}

ratios.max = list()
for(i in 1: length(final_filt.list.skew.max)){
  
  if(Ns[i] == 0){next}
  ratios.max[[i]] = final_filt.list.skew.max[[i]][ ,c('C.1', 'G')]  
  
}


list.skew.max = final_filt.list.skew.max


est_skew_func = function(list.skew.max, Ns, ratios.max){
  
  #If there's only 1 or 0 SNPs, don't estimate a skew
  if(Ns <= 0){return(NA)}
  #Fold the reference skews
  folded_ref_skews = folded(ratios.max$C.1)
  #return the mean and sigma estimates
  return(mle_folded(folded_ref_skews))
}

folded_norm_fits = mclapply(1:length(list.skew.max), function(i) est_skew_func(list.skew.max[[i]], Ns[i], ratios.max[[i]]), mc.cores = 20)

#Grab the skew and variance estimates
est_skew = vector(mode = 'numeric', length = length(folded_norm_fits))
est_var = vector(mode = 'numeric', length = length(folded_norm_fits))
for(i in 1:length(folded_norm_fits)){
  
  if(is.na(folded_norm_fits[i])){
    est_skew[i] = NA
    est_var[i] = NA
    next
  }
  
  est_skew[i] = folded_norm_fits[[i]][[1]]
  est_var[i] = folded_norm_fits[[i]][[2]]
}


######
#Save the filtered dataset
save(list.skew.max, Ns, ratios.max, file = "/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/GATK_mod/all_v8_GATK_1_20_2021_update_filt_escape.skew.est.max.genes.Rdata")




#########################################################
#Testing gene skew correlations to tissue estimated skew
#Go Through the list of genes in the data set, 
#pull out the gene skew distribution, 
#generate skew estimates without that gene,
#grab the correlation between the two in folded space

#Seeing which genes follow the skew of the tissue
#########################################################


#permutation test to get a pvalue for the correlation

cor_permute = function(x,y){
  permute_y = sample(y,size=length(y), replace=FALSE)
  corr = suppressWarnings(cor(x,permute_y, method='pearson'))
  return(corr)
}


num_permutes = 10000  #For correlation permutation test

gene_skew_correlations = vector(mode = 'numeric', length = length(keep_genes))
names(gene_skew_correlations) = keep_genes
gene_skew_cor_pvalues = vector(mode = 'numeric', length = length(keep_genes))
names(gene_skew_cor_pvalues) = keep_genes
gene_skew_num_samps = vector(mode = 'numeric', length = length(keep_genes))
names(gene_skew_num_samps) = keep_genes


for(j in 1:length(keep_genes)){
  
  
  current_gene = keep_genes[j]
  
  #Go through the data, and get the index of samples that detect that gene and the reference skew for that gene
  re_gene_ref_skew = c()
  sample_index_for_gene= c()
  for(i in 1:length(list.skew.max)){
    
    data = list.skew.max[[i]]
    present = data$name == current_gene
    if(sum(present) == 0) #If detected in that sample, grab the sample index and the gene skew
    {next} else{
      sample_index_for_gene = c(sample_index_for_gene, i)
      re_gene_ref_skew = c(re_gene_ref_skew ,data[present, 'C.1'])
    }
  }
  
  
  #Just looking at the samples that detect that gene, filter it out and redo the skew estimates
  temp.list.skew.max = list.skew.max[sample_index_for_gene]
  #Filter out the gene
  bad_genes = current_gene  #For the filtering out function, i used a global variable, not smart
  temp.list.skew.max = lapply(temp.list.skew.max, filt_bad_genes)
  
  #Redo the ratios.max and the Ns objects
  temp.Ns = vector(mode='numeric', length = length(temp.list.skew.max))
  for(i in 1:length(temp.list.skew.max)){
    if(is.null(temp.list.skew.max[[i]])){temp.Ns[i] = 0}
    else{temp.Ns[i] = dim(temp.list.skew.max[[i]])[1]}
  }
  
  temp.ratios.max = list()
  for(i in 1: length(temp.list.skew.max)){
    if(temp.Ns[i] == 0){next}
    temp.ratios.max[[i]] = temp.list.skew.max[[i]][ ,c('C.1', 'G')]  
  }
  
  
  #Get the skew estimates
  temp_folded_norm_fits = mclapply(1:length(temp.list.skew.max), function(i) est_skew_func(temp.list.skew.max[[i]], temp.Ns[i], temp.ratios.max[[i]]), mc.cores = 20)
  
  #Grab the skew and variance estimates
  temp_est_skew = vector(mode = 'numeric', length = length(temp_folded_norm_fits))
  for(i in 1:length(temp_folded_norm_fits)){
    if(is.na(temp_folded_norm_fits[i])){
      temp_est_skew[i] = NA
      next
    }
    temp_est_skew[i] = temp_folded_norm_fits[[i]][[1]]
  }
  
  
  #Exclude NA cases, where the sample did not get a skew estimate
  na_index = !is.na(temp_est_skew)
  temp_est_skew = temp_est_skew[na_index]
  re_gene_ref_skew = re_gene_ref_skew[na_index]
  
  
  original_cor = cor(temp_est_skew, folded(re_gene_ref_skew), method = 'pearson')
  permute_cors = unlist(mclapply(1:num_permutes, function(i) cor_permute(temp_est_skew, folded(re_gene_ref_skew)), mc.cores=5))
  corr_pvalue = sum(abs(permute_cors) > abs(original_cor) ) / num_permutes
  
  gene_skew_correlations[j] = original_cor
  gene_skew_cor_pvalues[j] = corr_pvalue
  gene_skew_num_samps[j] = length(temp_est_skew)
  print(sprintf('Finished %s, correlation: %0.3f  pvalue; %0.5f num samps: %i finished %i genes at %s', 
                current_gene, original_cor, corr_pvalue,length(temp_est_skew), j, Sys.time()))
  
}



corrected_correlation_pvalues = p.adjust(gene_skew_cor_pvalues, method = 'BH')
filt_correlations = gene_skew_correlations[corrected_correlation_pvalues <= .05 & gene_skew_num_samps >= 30]
filt_correlations = filt_correlations[order(filt_correlations, decreasing = FALSE)]



#Now look at the expression of those genes, what's most likely the explaining variable in correlation is the expression, higher expressed genes are better correlated with tissue skew
keep_genes_expression_list = list()
for(j in 1:length(keep_genes)){

  
  current_gene = keep_genes[j]
  gene_expression = c()
  #Go through the dataset and get the total expression of each gene in the samples it's detected in
for(i in 1:length(list.skew.max)){
  
  data = list.skew.max[[i]]
  present = data$name == current_gene
  if(sum(present) == 0) #If detected in that sample, grab the sample index and the gene skew
  {next} else{
    gene_expression = c(gene_expression ,data[present, 'A.1'])
  }
}

keep_genes_expression_list[[j]] = gene_expression
}

names(keep_genes_expression_list) = keep_genes





#Get the average expression, can also plot the whole distributions if I want, starting simple
mean_expression_keep_genes = unlist(lapply(keep_genes_expression_list, mean))

gene_skew_corr_df = data.frame(gene_name = keep_genes, correlation_to_tissue_skew = gene_skew_correlations, pvalue =gene_skew_cor_pvalues, 
                               corrected_pvalue = corrected_correlation_pvalues, avg_expression =mean_expression_keep_genes,  num_samples = gene_skew_num_samps)




#save only significant and well powered correlations/genes


index = corrected_correlation_pvalues <= .05 & gene_skew_num_samps >= 30
gene_skew_corr_df = gene_skew_corr_df[index, ]
gene_skew_corr_df = gene_skew_corr_df[order(gene_skew_corr_df$correlation_to_tissue_skew), ]

gene_skew_corr_df$skew_rank = 1:(dim(gene_skew_corr_df)[1])



#Replot above but with better aesthetics
annotation_data = gene_skew_corr_df[gene_skew_corr_df$gene_name %in% c('XIST', 'TSIX', 'AR'), ]

ggplot(gene_skew_corr_df, aes(x = skew_rank, y = correlation_to_tissue_skew)) + geom_jitter(size = 1,alpha = .5, width = 0.01, height = 0.01) + 
  ylim(0,1) + ylab('Similarity between gene and tissue skew (corr)') + xlab('Genes detected in at least 30 samples') +
  ggtitle('Correlation of gene and tissue skews') +
  geom_point(data = annotation_data, aes(x = skew_rank, y = correlation_to_tissue_skew, col = 'red'), size = 3) + 
  geom_label_repel(data = annotation_data, aes(label = gene_name), size = 4, nudge_y = .1) +
  theme(plot.title = element_text(hjust = .5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), 
        legend.position = 'none', panel.background = element_rect(fill = 'white')) 



temp_gene_df = gene_skew_corr_df
temp_gene_df$avg_expression = log10(temp_gene_df$avg_expression)
temp_gene_df = temp_gene_df[order(temp_gene_df$avg_expression), ]
num_samps = dim(temp_gene_df)[1]
num_bins = 5
num_samps_in_bin = num_samps / num_bins 

labels = rep('Binned expression 1', num_samps)
for( i in 1:(num_bins-1)){
  labels[((i*num_samps_in_bin) + 1): ((i+1)*num_samps_in_bin)] = sprintf('Binned expression %i', i+1)
}

temp_gene_df$label = labels


ggplot(temp_gene_df, aes(x = label, y = correlation_to_tissue_skew)) + geom_boxplot()+
  ylab('Similarity between gene and tissue skew (corr)') + ylim(0,1) +
  theme(plot.title=element_text(hjust=.5, size=15), 
        axis.text.y=element_text(size=12),axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),axis.title.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.background = element_rect(fill = 'white'))



#Save the correlation results and the keep_genes and genes 

save(genes, keep_genes, gene_skew_correlations, gene_skew_cor_pvalues, gene_skew_num_samps, bad_gene_skew_df, gene_skew_corr_df,
     file = '/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/GATK_mod/gene_skew_correlations/gene_skew_to_tissue_skew_correlations_v8.Rdata'  )



#############################
#Figure 3 panel a
#############################

#Plot genes and sample size with karyoploteR

library(karyoploteR)

#Get the xchrom starting positions of the genes we kept

kept_gene_meta = attr[attr$name %in% keep_genes, ]

#Match the genes starting position with the number of samples that it was used in for skew estimates
index = match(kept_gene_meta$name,names(gene_skew_num_samps) )
gene_skew_num_samps = gene_skew_num_samps[index]

kp = plotKaryotype(genome="hg38", plot.type=2, chromosomes=c("chrX"))
#White text for the black bands
color_vec = rep('black',40)
color_vec[c(8,10,23,25,33,39)] = 'white'
kpAddCytobandLabels(kp,force.all = T, srt=90, cex = .75, col = color_vec)
#Just a background for the number of samples panel
kpDataBackground(kp, data.panel = 2)
#First do the correlations
#Label a few select genes
gene_skew_cor_pvalues_adjusted = p.adjust(gene_skew_cor_pvalues, method = 'BH')
label_genes = names(which(gene_skew_correlations >= .8 & gene_skew_cor_pvalues_adjusted <= .05))
index = kept_gene_meta$name %in% label_genes
#Random Y axis plaacement
#Now add a panel showing the number of samples each gene is detected in, log10 scale
kpAxis(kp,ymax =4000, tick.pos = c(100,1000,4000), labels = c(2,3,3.6), data.panel = 2)
kpPoints(kp, chr="chrX",data.panel=2, x=kept_gene_meta$start, y= log10(gene_skew_num_samps) / log10(3700), cex=.75, pch=19) 


