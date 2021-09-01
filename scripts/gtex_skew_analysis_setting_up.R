
#Setting up some matrices and general objects that will be used for further downstream analysis for the GTEx skews

#Libraries
library(ComplexHeatmap)
library(circlize)


#Gtex metadata
load("/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/sample_meta.Rdata")
#estimated Skews and stats with escape filtered out
load("/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/GATK_mod/final_filtering_skew_and_stats_11_5_20_escape_filt_df.Rdata")


###############
#Functions
###############
# Functions to fold and unfold ratios 
folded <- function(x) { apply( cbind(x,1-x), 1, max)} 
unfold <- function(x) { c(x,1-x) } 




#Dataframe with the tissues and their germ layers
tissues = names(table(skew_and_stats_df$tissue))
germ_layer = c()

for( i in 1:length(tissues)){
  germ_layer = c(germ_layer, as.character(skew_and_stats_df$germ_layer[skew_and_stats_df$tissue == tissues[i]][1]))
}

germ_layer_df = data.frame(tissues = tissues, germ_layers = germ_layer)




#Grouping tissues into larger groups
sample_tissue_meta = as.character(skew_and_stats_df$tissue)
sample_tissue_meta[sample_tissue_meta == 'Adipose - Subcutaneous' | sample_tissue_meta == 'Adipose - Visceral (Omentum)'] = 'Adipose'
sample_tissue_meta[sample_tissue_meta == 'Artery - Aorta' | sample_tissue_meta == 'Artery - Coronary' | sample_tissue_meta == 'Artery - Tibial' ] = 'Artery'
sample_tissue_meta[sample_tissue_meta == 'Cervix - Ectocervix' | sample_tissue_meta == 'Cervix - Endocervix'] = 'Cervix'
sample_tissue_meta[sample_tissue_meta == 'Colon - Sigmoid' | sample_tissue_meta == 'Colon - Transverse'] = 'Colon'
sample_tissue_meta[sample_tissue_meta == 'Esophagus - Gastroesophageal Junction' | sample_tissue_meta == 'Esophagus - Mucosa' | sample_tissue_meta == 'Esophagus - Muscularis' ] = 'Esophagus'
sample_tissue_meta[sample_tissue_meta == 'Heart - Atrial Appendage' | sample_tissue_meta == 'Heart - Left Ventricle'] = 'Heart'
sample_tissue_meta[sample_tissue_meta == 'Skin - Not Sun Exposed (Suprapubic)' | sample_tissue_meta == 'Skin - Sun Exposed (Lower leg)'] = 'Skin'
sample_tissue_meta[sample_tissue_meta == 'Brain - Amygdala' | sample_tissue_meta == 'Brain - Anterior cingulate cortex (BA24)' | sample_tissue_meta == 'Brain - Caudate (basal ganglia)' | sample_tissue_meta == 'Brain - Cerebellar Hemisphere' | sample_tissue_meta == 'Brain - Cerebellum' | sample_tissue_meta == 'Brain - Cortex' | sample_tissue_meta == 'Brain - Frontal Cortex (BA9)' | sample_tissue_meta == 'Brain - Hippocampus' | sample_tissue_meta == 'Brain - Hypothalamus' | sample_tissue_meta == 'Brain - Nucleus accumbens (basal ganglia)' | sample_tissue_meta == 'Brain - Putamen (basal ganglia)' | sample_tissue_meta == 'Brain - Spinal cord (cervical c-1)' | sample_tissue_meta == 'Brain - Substantia nigra'] = 'Brain'

skew_and_stats_df$grouped.tissue = sample_tissue_meta




#Set up some matrices we'll be reusing and the germ layer annotations/colors.
donors = names(table(skew_and_stats_df$donor))
#Start with all the tissues, exclude the male tissues
tissues =  names(table(skew_and_stats_df$tissue))
tissues = tissues[! tissues %in% c('Prostate', 'Testis','body_site_s')]
#Then look at the grouped tissues
grouped_tissues = names(table(skew_and_stats_df$grouped.tissue))
grouped_tissues = grouped_tissues[! grouped_tissues %in% c('Prostate', 'Testis','body_site_s')]


#Tissues with at least X donors
filt_tissues = names(table(skew_and_stats_df$tissue)[table(skew_and_stats_df$tissue) >= 10])


#Save distinct colors for the germ layers
germ_layer_colors = brewer.pal(n = 4, name = "Paired")
ecto_color = germ_layer_colors[1]
endo_color = germ_layer_colors[2]
meso_color = germ_layer_colors[4]



#Start off looking at a table of pairwise donations, number of donors contributing tissues X,Y
#Just getting a sense of how much pairwise donations there are to work with
#And make a matrix with the donors names that do donate both of those tissues


#Contains umber of donors that contribute to the X.Y tissue
tiss_tiss_donors_df = as.data.frame(matrix(nrow = length(tissues), ncol = length(tissues)))
rownames(tiss_tiss_donors_df) = tissues
colnames(tiss_tiss_donors_df) = tissues
#Contains the donor names that contribute to the X,Y tissue
tiss_tiss_donors_names_matrix = matrix(rep(list(), length(tissues)*length(tissues)), nrow=length(tissues), ncol=length(tissues))
rownames(tiss_tiss_donors_names_matrix) = tissues
colnames(tiss_tiss_donors_names_matrix) = tissues

for(i in 1:length(tissues)){
  
  current_tiss = tissues[i]
  
  for(j in 1: length(tissues)){
    
    next_tiss = tissues[j]
    
    donors_for_tiss = as.character(skew_and_stats_df$donor[skew_and_stats_df$tissue == current_tiss])
    #Just going to ignore the donors that have multiple donations for the same tiss
    duplicates = donors_for_tiss[duplicated(donors_for_tiss)] 
    donors_for_tiss = donors_for_tiss[!donors_for_tiss %in% duplicates]  
    
    donors_for_next_tiss = as.character(skew_and_stats_df$donor[skew_and_stats_df$tissue == next_tiss])
    duplicates = donors_for_next_tiss[duplicated(donors_for_next_tiss)] 
    donors_for_next_tiss = donors_for_next_tiss[!donors_for_next_tiss %in% duplicates]
    
    tiss_tiss_donors_df[i,j] = sum(unique(donors_for_tiss) %in% unique(donors_for_next_tiss))
    tiss_tiss_donors_names_matrix[i,j] = list(unique(donors_for_tiss[donors_for_tiss %in% donors_for_next_tiss]))
  }
}


#binary donor X tissue matrix representing which donors have donated to which tissues
#Contains the donor names that contribute to that tissue
diagonal = diag(tiss_tiss_donors_names_matrix)

binary_tissue_donor_matrix = matrix(nrow = length(names(diagonal)), ncol = length(donors))
rownames(binary_tissue_donor_matrix) = names(diagonal)
colnames(binary_tissue_donor_matrix) = donors

for(i in 1:length(names(diagonal))){
  #Donors for that tissue
  donors_for_tiss = diagonal[[i]]
  index = colnames(binary_tissue_donor_matrix) %in% donors_for_tiss
  binary_tissue_donor_matrix[i, index] = 1 #Donated
  binary_tissue_donor_matrix[i, !index] = 0 #Did not donate
}

binary_tissue_donor_matrix = binary_tissue_donor_matrix[!rownames(binary_tissue_donor_matrix) %in% c('Testis','Prostate','body_site_s'), ]



#####################################
#Plot 
#####################################

#Get the number of tissues per donor
num_tiss_per_donor = colSums(binary_tissue_donor_matrix)
max_tiss_per_donor = max(num_tiss_per_donor)
#Get the number of donors per tissue
num_donors_per_tissue = rowSums(binary_tissue_donor_matrix)
max_donor_per_tiss = max(num_donors_per_tissue)

#Germ layer annotation and number of donors per tissue
row_annot = HeatmapAnnotation(germ_layer = germ_layers, donors_per_tiss = anno_barplot(num_donors_per_tissue, gp = gpar(fill = rep(1,length(tissues)))),
                              which = 'row',
                              col = list(germ_layer=c('Ectoderm'=ecto_color, 'Mesoderm'=meso_color, 'Endoderm' = endo_color)), 
                              gp = gpar(col = "black"), show_annotation_name = c(FALSE,TRUE),
                              annotation_legend_param = list(title = 'Germ layer'))
names(row_annot) = c('Germ layers','Donors per tissue')

#number of tissues per donor column point annotation
tissues_per_donor_annot = HeatmapAnnotation(tissues_per_donor = anno_barplot(num_tiss_per_donor),which = 'column')
names(tissues_per_donor_annot) = c('Tissues per donor')

#pdf(file = '/home/werner/xchrom_snp_skew/code/graphs/summary_graphs/binary_donors_per_tissue_heatmap.pdf', height = 8, width =8, useDingbats = FALSE)

#Colors for the heatmap
colors = structure(c('white','black'), names = c("0", "1"))
Heatmap(binary_tissue_donor_matrix, col=colors, column_title = 'Donor tissue contributions', use_raster = TRUE,
        clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2", 
        show_row_dend=TRUE, row_names_gp = gpar(fontsize = 8), show_column_names = FALSE,
        right_annotation = row_annot,
        top_annotation = tissues_per_donor_annot,
        heatmap_legend_param = list( title = "Donated", at = c(0, 1), labels = c("FALSE", 'TRUE')))
#dev.off()

#####################################
#End plot 
#####################################










