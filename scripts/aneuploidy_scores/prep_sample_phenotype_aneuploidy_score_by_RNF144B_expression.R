# 
# Author: Gerard Romero Sola
# Inspiration: Miquel Anglada Girotto
# Github inspiration/copy: https://github.com/MiqG/publication_zadra_ttll11
#
# Script purpose
# --------------
# Prepare Sample Phenotype data

### Libraries
require(readr)
require(dplyr)
require(here)

### Variables
ROOT = here()
DATA_DIR = here(ROOT, "data")
TCGA_DIR = here(DATA_DIR, "raw", "UCSCXena", "TCGA")
PREP_DIR = here(DATA_DIR, "prep")

SAMPLE_TYPES_OI = c('Primary Tumor')
GENES_OI = c('TP53')
THRESH = 20

# inputs
raw_survival_file = here(TCGA_DIR, 'phenotype','Survival_SupplementalTable_S1_20171025_xena_sp')
raw_phenotype_file = here(TCGA_DIR, 'phenotype','TCGA_phenotype_denseDataOnlyDownload.tsv.gz')
mutations_file = here(TCGA_DIR, 'mutations','mc3.v0.2.8.PUBLIC.nonsilentGene.xena') 

# outputs
prep_phenotype_file = here(PREP_DIR,'sample_phenotype_aneuploidy_RNF144B_expression.tsv')

# load data
df = read_tsv(raw_phenotype_file)
df_surv = read_tsv(raw_survival_file)
df_mut = read_tsv(mutations_file)

# add cancer type abbreviations to phenotype
abbreviations = merge(df, df_surv, by='sample') %>% distinct(`cancer type abbreviation`,`_primary_disease`)
df = merge(df, abbreviations, by='_primary_disease')
df = df %>%
	select(c(sample,`cancer type abbreviation`, sample_type)) %>%
	filter(sample_type %in% SAMPLE_TYPES_OI)
colnames(df) = c('sample','cancer_type','sample_type')

# Stratify by p53 state
## subset mutations
df_mut_p53 = df_mut %>%
	filter(sample %in% GENES_OI) %>%
	t() %>%
	as.data.frame() %>%
	filter(!row_number() %in% 1)
colnames(df_mut_p53) = c("TP53_mut")
df_mut_p53$sample = rownames(df_mut_p53)

df_all = merge(df, df_mut_p53, by = "sample", all.x = TRUE)

# save
write_tsv(df_all, prep_phenotype_file)

print('Done!')