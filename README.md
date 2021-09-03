# human_cross_tissue_XCI

This repository contains code for the manuscript "Cross-tissue analysis of allelic X-chromosome inactivation ratios resolves features of human development"
All associated data generated from these notebooks are available on the FTP server: http://labshare.cshl.edu/shares/gillislab/people/werner/xskew2021_preprint/data


Description of notebooks, 1, 2, and 3 should be run in that order, order does not matter for 4,5, and 6.

1.)GTEx_SNP_gene_filtering_XCI_ratio_estimates.Rmd
Filters identified heterozygous SNPs from vcf files and identifies the max-powered SNP per gene. Filtering genes for reference bias. Calculates inital XCI ratio estimates and computes gene-sample XCI ratio correlations. After running the bootstrapping_script to calculate confidence intervals around the XCI ratio estimates, compiles a dataframe for all XCI ratios and associated metadata for the GTEx dataset. 

2.) bootstrapping_script
Calculates XCI ratio estimates and 95% confidence intervals

3.)tissues_stem_cell_pop_estimates.Rmd
Filters XCI ratio estimates. Models cell numbers from XCI ratio variance in tissue and donor XCI ratio distributions. Computes XCI ratio correlations within and across germ layers

4.)snp_level_skew_correlations.Rmd
Computes correlations of allelic expression ratios for shared SNPs across tissues of individual donors, captures the parental direction of XCI across tissues for individuals 

5.)exploring_escape_genes_gtex.Rmd
Calculates statistics related to escape from XCI.

6.)general_XCI_skew_graphs.Rmd
Miscellaneous graphs 






