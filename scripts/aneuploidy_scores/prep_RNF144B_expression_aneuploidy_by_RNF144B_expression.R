# 
# Author: Gerard Romero Sola
# Inspiration: Miquel Anglada Girotto
# Github inspiration/copy: https://github.com/MiqG/publication_zadra_ttll11
#
# Script purpose
# --------------
# Prepare RNF144b and p53 Expression

require(tidyverse)
require(limma)
require(magrittr)
require(reshape2)
require(ggpubr)
require(here)

ROOT = here()
DATA_DIR = here(ROOT,'data')
TCGA_DIR = here(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = here(DATA_DIR,'prep')

# variables
SAMPLE_TYPES_OI = c('Primary Tumor')
GENES_OI = c('RNF144B')

# inputs
phenotype_file = here(PREP_DIR,'sample_phenotype_aneuploidy_RNF144B_expression.tsv')
genexpr_file = here(TCGA_DIR,'rnaseq','AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz')

# outputs
output_file = here(PREP_DIR,'genexpr_RNF144B_aneuploidy_by_RNF144B_expression.tsv')

# load data
metadata = read_tsv(phenotype_file)
samples_oi = metadata %>% pull(sample)
mat = read.columns(genexpr_file, required.col = c('sample',samples_oi))

# preprocess
idx = mat$sample %in% GENES_OI
genexpr = mat[idx,] %>%
	data.frame(row.names=c()) %>%
	column_to_rownames('sample') %>%
	t() %>%
	melt(varnames = c('sample', 'gene'), value.name = 'expression') %>%
	mutate(sample = gsub('\\.','-',sample),
				 expression = as.numeric(expression)) %>%
	left_join(metadata, by='sample') %>%
	drop_na()

# save
write_tsv(genexpr, output_file)

print('Done!')