# 
# Author: Gerard Romero Sola
# Inspiration: Miquel Anglada Girotto
# Github inspiration/copy: https://github.com/MiqG/publication_zadra_ttll11
#
# Script purpose
# --------------
# Prepare correlation of expression and aneuloidy

# libraries
require(readr)
require(dplyr)
require(limma)
require(tibble)
require(doParallel)
require(here)

# variables
ROOT = here()
DATA_DIR = here(ROOT,'data')
TCGA_DIR = here(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = here(DATA_DIR,'prep')
RESULTS_DIR = here(ROOT,'output')

COR_METHOD = 'spearman'

N_CORES = 10

# inputs
aneuploidy_scores_file = here(PREP_DIR,'aneuploidy.tsv')
phenotype_file = here(PREP_DIR,'sample_phenotype.tsv')
phenotype_file_p53 = here(PREP_DIR,'sample_phenotype_p53.tsv')
genexpr_file = here(TCGA_DIR,'rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz')

# outputs
correlations_file_RNF144B = here(RESULTS_DIR,'correlation-genexpr_aneuploidy_RNF144B.tsv')
correlations_file_TP53_mut = here(RESULTS_DIR,'correlation-genexpr_aneuploidy_TP53_mut.tsv')
correlations_file_TP53_WT = here(RESULTS_DIR,'correlation-genexpr_aneuploidy_TP53_WT.tsv')

# load data (NOTE: gene expression data is loaded on the fly)
score = read_tsv(aneuploidy_scores_file) 
metadata = read_tsv(phenotype_file)
metadata_p53 = read_tsv(phenotype_file_p53)

# available samples 
## in gene expression matrix
genexpr_samples = strsplit(readLines(file(genexpr_file,'r'),n=1),split='\t')[[1]]
genexpr_samples = setdiff(genexpr_samples,'sample')
genexpr_genes = read.columns(genexpr_file, required.col = 'sample')[[1]]

metadata_scores_p53 = merge(score, metadata_p53, by = "sample")

## aneuploidy score
score_samples_RNF144B = score$sample
score_samples_TP53_mut = metadata_scores_p53$sample[metadata_scores_p53$TP53_mut == 1]
score_samples_TP53_WT = metadata_scores_p53$sample[metadata_scores_p53$TP53_mut == 0]

# measure correlation between gene expression and aneuplo
cancer_types_RNF144B = metadata %>% pull(cancer_type) %>% unique()
cancer_types_p53 = metadata_p53 %>% pull(cancer_type) %>% unique()

cl =  makeCluster(N_CORES)
registerDoParallel(cl)
correlations_RNF144B = foreach(cancer_type_oi=cancer_types_RNF144B, 
																.combine=cbind,
																.packages=c('dplyr','limma')) %dopar% {
																	# get samples for correlation
																	samples_oi = metadata %>% filter(cancer_type==cancer_type_oi) %>% pull(sample)
																	samples_oi = intersect(samples_oi, intersect(genexpr_samples, score_samples_RNF144B))
																	
																	# read gene expression for these samples
																	genexpr = read.columns(genexpr_file, required.col = c('sample',samples_oi))
																	genexpr[,samples_oi] = genexpr[,samples_oi]
																	
																	# compute correlations
																	x = t(genexpr[,samples_oi])
																	y = score$aneuploidy_score[match(samples_oi, score$sample)]
																	corr = cor(x, y, method=COR_METHOD)    
																	colnames(corr) = cancer_type_oi
																	return(corr)
																}
correlations_TP53_mut = foreach(cancer_type_oi=cancer_types_p53, 
											 .combine=cbind,
											 .packages=c('dplyr','limma')) %dopar% {
											 	# get samples for correlation
											 	samples_oi = metadata_p53 %>% filter(cancer_type==cancer_type_oi) %>% pull(sample)
											 	samples_oi = intersect(samples_oi, intersect(genexpr_samples, score_samples_TP53_mut))
											 	
											 	# read gene expression for these samples
											 	genexpr = read.columns(genexpr_file, required.col = c('sample',samples_oi))
											 	genexpr[,samples_oi] = genexpr[,samples_oi]
											 	
											 	# compute correlations
											 	x = t(genexpr[,samples_oi])
											 	y = score$aneuploidy_score[match(samples_oi, score$sample)]
											 	corr = cor(x, y, method=COR_METHOD)    
											 	colnames(corr) = cancer_type_oi
											 	return(corr)
											 }
correlations_TP53_WT = foreach(cancer_type_oi=cancer_types_p53, 
																.combine=cbind,
																.packages=c('dplyr','limma')) %dopar% {
																	# get samples for correlation
																	samples_oi = metadata_p53 %>% filter(cancer_type==cancer_type_oi) %>% pull(sample)
																	samples_oi = intersect(samples_oi, intersect(genexpr_samples, score_samples_TP53_WT))
																	
																	# read gene expression for these samples
																	genexpr = read.columns(genexpr_file, required.col = c('sample',samples_oi))
																	genexpr[,samples_oi] = genexpr[,samples_oi]
																	
																	# compute correlations
																	x = t(genexpr[,samples_oi])
																	y = score$aneuploidy_score[match(samples_oi, score$sample)]
																	corr = cor(x, y, method=COR_METHOD)    
																	colnames(corr) = cancer_type_oi
																	return(corr)
																}
stopCluster(cl)

# add genes
correlations_RNF144B = data.frame(gene=genexpr_genes, correlations_RNF144B)
correlations_TP53_mut = data.frame(gene=genexpr_genes, correlations_TP53_mut)
correlations_TP53_WT = data.frame(gene=genexpr_genes, correlations_TP53_WT)

# save
write_tsv(correlations_RNF144B, correlations_file_RNF144B)
write_tsv(correlations_TP53_mut, correlations_file_TP53_mut)
write_tsv(correlations_TP53_WT, correlations_file_TP53_WT)

print('Done!')