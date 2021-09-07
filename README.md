# human_cross_tissue_XCI

This repository contains code for the manuscript "Cross-tissue analysis of allelic X-chromosome inactivation ratios resolves features of human development"
All relevant data generated from these notebooks are available on the FTP server: http://labshare.cshl.edu/shares/gillislab/people/werner/xskew2021_preprint/data
To fully replicate all analysis, you first need access to and download the .fastq data from GTEx, which requires dbGAP permissions. Once you have the fastqs, the VCF and WIG files used as input for estimating XCI ratios can be generated following the snakemake pipeline in human_cross_tissue_XCI/code/snakemake_data_processing_pipeline

Notebooks, 1, 2, and 3 should be run in that order, order does not matter for 4,5, 6, and 7.

1.)GTEx_SNP_gene_filtering_XCI_ratio_estimates.Rmd
Filters identified heterozygous SNPs from vcf files and identifies the max-powered SNP per gene. Filtering genes for reference bias. Calculates inital XCI ratio estimates and computes gene-sample XCI ratio correlations. After running the bootstrapping_script to calculate confidence intervals around the XCI ratio estimates, compiles a dataframe for all XCI ratios and associated metadata for the GTEx dataset. 

2.) bootstrapping_script
Calculates XCI ratio estimates and 95% confidence intervals

3.)GTEx_XCI_ratio_data_analysis_cell_counts_correlations.Rmd
Filters XCI ratio estimates. Models cell numbers from XCI ratio variance in tissue and donor XCI ratio distributions. Computes XCI ratio correlations within and across germ layers

4.) entex_phased_control_folded_normal.Rmd
Compares XCI ratio estimates using phased RNA-sequencing data from the EN-TEx consortium to XCI ratio estimates using the folded-normal model

5.)snp_level_skew_correlations.Rmd
Computes correlations of allelic expression ratios for shared SNPs across tissues of individual donors, captures the parental direction of XCI across tissues for individuals 

6.)exploring_escape_genes_gtex.Rmd
Calculates statistics related to escape from XCI.

7.)general_XCI_skew_graphs.Rmd
Miscellaneous graphs 

Description of data on http://labshare.cshl.edu/shares/gillislab/people/werner/xskew2021_preprint/data

1.)gene_annotations_v25.Rdata 
Genome annotations used in this study

2.)gene_escape_table.txt 

Table containing annotations of inactive or escape status from the escape analysis in this paper

3.)tukiainen_study.txt 
Escape annotations from https://www.nature.com/articles/nature24265


4.)v8_GTEx_skew_and_stats_df.Rdata
All estimated XCI ratios for the female GTEx samples we processed with associated statistics and metadata, generated from GTEx_SNP_gene_filtering_XCI_ratio_estimates.Rmd and bootstrapping_script. Can be used to replicate the main analyses in GTEx_XCI_ratio_data_analysis_cell_counts_correlations.Rmd   

