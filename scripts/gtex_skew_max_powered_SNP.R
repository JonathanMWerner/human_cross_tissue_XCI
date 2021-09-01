
#SNP filtering for the GTEx samples
#



load("/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/GATK_mod/all_v8_GATK_mod.Rdata")  #VCF files
load("/data/genomes/human/gene_annotations_v25.Rdata")                                   #Genome annotations

#Grab the metadata for the samples
sample_meta = mapping[f.m, ]
sample_tissues = names(table(sample_meta[ ,'body_site_s']))
sample_tissue_table =table(sample_meta[ ,'body_site_s']) 

sample_ids = as.character(sample_meta[ ,'Sample_Name_s'])
split_sample_ids = lapply(sample_ids, strsplit, split='-', fixed=TRUE)
donor_id = lapply(split_sample_ids, '[[', 1)        #Get rid of first level
donor_id =as.character(lapply(donor_id, '[[', 2) )  #Grab the donor ID
sample_meta$Donor_id = donor_id 



#Grouping tissues into larger groups
sample_tissue_meta = as.character(sample_meta$body_site_s)
sample_tissue_meta[sample_tissue_meta == 'Adipose - Subcutaneous' | sample_tissue_meta == 'Adipose - Visceral (Omentum)'] = 'Adipose'
sample_tissue_meta[sample_tissue_meta == 'Artery - Aorta' | sample_tissue_meta == 'Artery - Coronary' | sample_tissue_meta == 'Artery - Tibial' ] = 'Artery'
sample_tissue_meta[sample_tissue_meta == 'Cervix - Ectocervix' | sample_tissue_meta == 'Cervix - Endocervix'] = 'Cervix'
sample_tissue_meta[sample_tissue_meta == 'Colon - Sigmoid' | sample_tissue_meta == 'Colon - Transverse'] = 'Colon'
sample_tissue_meta[sample_tissue_meta == 'Esophagus - Gastroesophageal Junction' | sample_tissue_meta == 'Esophagus - Mucosa' | sample_tissue_meta == 'Esophagus - Muscularis' ] = 'Esophagus'
sample_tissue_meta[sample_tissue_meta == 'Heart - Atrial Appendage' | sample_tissue_meta == 'Heart - Left Ventricle'] = 'Heart'
sample_tissue_meta[sample_tissue_meta == 'Skin - Not Sun Exposed (Suprapubic)' | sample_tissue_meta == 'Skin - Sun Exposed (Lower leg)'] = 'Skin'
sample_tissue_meta[sample_tissue_meta == 'Brain - Amygdala' | sample_tissue_meta == 'Brain - Anterior cingulate cortex (BA24)' | sample_tissue_meta == 'Brain - Caudate (basal ganglia)' | sample_tissue_meta == 'Brain - Cerebellar Hemisphere' | sample_tissue_meta == 'Brain - Cerebellum' | sample_tissue_meta == 'Brain - Cortex' | sample_tissue_meta == 'Brain - Frontal Cortex (BA9)' | sample_tissue_meta == 'Brain - Hippocampus' | sample_tissue_meta == 'Brain - Hypothalamus' | sample_tissue_meta == 'Brain - Nucleus accumbens (basal ganglia)' | sample_tissue_meta == 'Brain - Putamen (basal ganglia)' | sample_tissue_meta == 'Brain - Spinal cord (cervical c-1)' | sample_tissue_meta == 'Brain - Substantia nigra'] = 'Brain'

sample_meta$grouped.tissue = sample_tissue_meta

#Grouping by germ layer
temp = as.character(sample_meta$grouped.tissue)

temp[temp == 'Adipose' | temp == 'Heart' | temp == 'Whole Blood' | temp == 'Cells - EBV-transformed lymphocytes' | temp == 'Cells - Transformed fibroblasts' | temp == 'Muscle - Skeletal' | temp == 'Adrenal Gland' | temp == 'Artery' | temp == 'Fallopian Tube' | temp == 'Kidney - Cortex' | temp == 'Spleen' | temp == 'Uterus' | temp == 'Cervix' ] = 'Mesoderm'


temp[temp == 'Brain' | temp == 'Skin' | temp == 'Breast - Mammary Tissue' | temp == 'Nerve - Tibial' | temp == 'Minor Salivary Gland' | temp == 'Pituitary' ] = 'Ectoderm'


temp[temp == 'Lung' | temp == 'Thyroid' | temp == 'Pancreas' | temp == 'Liver' | temp == 'Colon' | temp == 'Small Intestine - Terminal Ileum' | temp == 'Vagina' | temp == 'Ovary' | temp == 'Stomach' | temp == 'Esophagus' | temp == 'Bladder' ] = 'Endoderm'

germ_layer = names(table(temp))
sample_meta$germ.layer = temp
save(sample_meta, file = '/home/werner/xchrom_snp_skew/data/GTEx/samples/sample_metadata_with_v8.Rdata')



#####################################
#For the GATK-mod files, create a new list with the filtered SNPs per sample, filtering on the additional statistics
#SOR < 3
#FS < 60
#abs(ReadPosRankSum) < 5

#This will also filter out any remaining indels in the vcf files as the additionaly statistics are only kept for SNPs, not indels

######################################

pos = list.vcf[[2]][ ,1]
present.f.stats = list.vcf.stats[[2]]$POS %in% pos

present.f.vcf = pos %in% list.vcf.stats[[2]]$POS
table(present.f.vcf)

t.list.vcf = list.vcf[[2]][present.f.vcf, ]
t.list.vcf.stats = list.vcf.stats[[2]][present.f.stats, ]

sor.f = t.list.vcf.stats$SOR < 3
sor.f[is.na(sor.f)] = TRUE
fs.f  = t.list.vcf.stats$FS < 60
fs.f[is.na(fs.f)] = TRUE
rpos.f = abs(t.list.vcf.stats$ReadPosRankSum) <= 5
rpos.f[is.na(rpos.f)] = TRUE

all.f  = sor.f & fs.f & rpos.f

t.list.vcf[all.f, ]



list.vcf.snp.filtered = list()

for(i in 1:length(list.vcf)){
  
  if(is.null(list.vcf[[i]])){next}
  
  
  pos = list.vcf[[i]][ ,1]
  present.f.stats = list.vcf.stats[[i]]$POS %in% pos
  present.f.vcf = pos %in% list.vcf.stats[[i]]$POS
  
  #If NA, means that stat is not available for that particular SNP, make NAs = true, really only interested in filtering out SNPs that fail the statistics thresholds
  
  t.list.vcf = list.vcf[[i]][present.f.vcf, ]
  t.list.vcf.stats = list.vcf.stats[[i]][present.f.stats, ]
  
  sor.f = t.list.vcf.stats$SOR < 3
  sor.f[is.na(sor.f)] = TRUE
  fs.f  = t.list.vcf.stats$FS < 60
  fs.f[is.na(fs.f)] = TRUE
  rpos.f = abs(t.list.vcf.stats$ReadPosRankSum) <= 5
  rpos.f[is.na(rpos.f)] = TRUE
  
  all.f  = sor.f & fs.f & rpos.f
  
  list.vcf.snp.filtered[[i]] = t.list.vcf[all.f, ]
  
}


########################################
#Identifying SNPs in annotated genes and choosing the max-powered SNPs per gene, the one with the highest total read count
########################################


list.vcf = list.vcf.snp.filtered

# Create filter of X chromosome genes. 
f.x = attr$chr == "chrX"

# Calculate skew ratios (ref and alt)
list.skew = list()
for( j in 1:length(list.vcf)) {
  if( is.null(list.vcf[[j]]) ) { next }
  f1 = rowSums(list.vcf[[j]][,2:3] > 0 ) == 2  #Only looking at SNPs that have reads > 0 for both the ref and alt allele
  
  skews = list.vcf[[j]][f1,1:6]                      #Grab the SNP position and ref alt read counts, only care about the first three fields
  skews[,4] = rowSums(list.vcf[[j]][f1,2:3])         #Fourth field is the total read counts for that SNP
  skews[,5:6] = list.vcf[[j]][f1,2:3] / skews[ ,4]   #Fifth and sixth fields are the ref and alt ratios (divided by total read counts)
  
  list.skew[[j]] = skews
}


#Get the number of SNPs per sample
Ns = sapply(1:length(list.skew), function(i) dim(list.skew[[i]])[1] ) 
A = sapply(1:length(list.skew), function(i) is.null(Ns[[i]] ) )  * 1     
for( i in 1:length(Ns) ) { if( is.null(Ns[[i]]))  { A[i] = 0 } else { A[i] = Ns[[i]] }}
Ns = A


# Find SNP/gene overlaps
for( j in 1:length(list.skew)) {
  if( is.null(list.skew[[j]]) ) { next }
  
  #For each SNP in the sample, for the attr[f.x, ] (All the xchrom genes) find the SNPs gene index
  genes.snps = sapply( 1:dim(list.skew[[j]])[1], function(i) which(list.skew[[j]][i,1] > attr[f.x,2]  &  list.skew[[j]][i,1] < attr[f.x,3] ))
  #If a SNP is not associated with a gene, put NA
  for(i in 1:length(genes.snps)) { if( length(genes.snps[[i]])==0 ) { genes.snps[[i]] = NA } ;  if( length(genes.snps[[i]]) > 1  ) { genes.snps[[i]] = genes.snps[[i]][1]  }   }
  
  
  if( sum(is.na(genes.snps)*1) == length(genes.snps) ){               #If all the SNPs in the sample are not associated with genes...
    na_genes = attr[f.x,][unlist(genes.snps),][1:length(genes.snps), ]  #An NA matrix with the attr column names, row number = number of SNPs in the sample
    list.skew[[j]] = cbind(list.skew[[j]], unlist(genes.snps), na_genes)
    next
  }
  
  #Add the gene informatino for each SNP in the sample
  list.skew[[j]] = cbind(list.skew[[j]], unlist(genes.snps), attr[f.x,][unlist(genes.snps),]  )
  
}

#Picking the most powered SNPs, the SNP with the most reads in a multiple-SNP gene

## Pick max powered SNPs
list.skew.max = list()
for(j in 1:length(list.skew)){ 
  if(Ns[j] > 0 ){  
    a = list.skew[[j]][,4]   #Get the total read counts for each SNP
    b = list.skew[[j]][,13]  #Get the ensemble ID
    
    abi = tapply(a,b, which.max)   #Find which genes have multiple SNPs, get index of max
    abi = abi[!is.na(abi)]         #Pull those out
    abi = cbind(names(abi), abi )  #Put the ensemble ID and number of SNPs together
    
    if(length(abi)==0){next}       #If there are no SNPs in annotated genes, just skip
    #Grab the SNP data for all the SNPs in the gene...     Grab the max
    maxtest  = lapply(1:dim(abi)[1], function(i) list.skew[[j]][ which(list.skew[[j]][,13]==abi[i,1]),][ as.numeric(abi[i,2] ),] )
    list.skew.max[[j]] = do.call(rbind, maxtest)  #Rind all max SNPs and store it
  }
}
rm(maxtest)



## Store as a separate matrices/lists
ratios.max.genes =  list() 
ratios.max = list()
for(j in 1:length(list.skew.max)){ 
  if(Ns[j] > 0 ){
    
    if(is.null(list.skew.max[[j]])){next}  #If there were no SNPs in annotated genes, just skip    
    
    ratios.max.genes[[j]] = list.skew.max[[j]][,c(13:17,8:12)]
    ratios.max[[j]] = list.skew.max[[j]][,5:6]
  }
}


save(list.skew, list.skew.max,ratios.max.genes,ratios.max, Ns, 
     file="/home/werner/xchrom_snp_skew/data/GTEx/snp_skews/GATK_mod/all_v8_GATK_mod_stat_filt_1_20_21.skew.est.max.genes.Rdata")







