# 
# Author: Gerard Romero Sola
# Contact: gerard [dot] romero02 [at] estudiant [dot] upf [dot] edu
#
# Script purpose
# --------------
# Subset expression of RNF144B across samples and combine
# with sample metadata.
# 

require(data.table)
require(tidyverse)
require(limma)
require(magrittr)
require(reshape2)
require(ggpubr)
require(here)
require(UCSCXenaTools)
require(data.table)
require(R.utils)
require(readr)
require(TCGAutils)

ROOT = here::here()
DATA_DIR = here(ROOT,'data')
TTG_DIR = here(DATA_DIR,'raw','UCSCXena','TTG')
PREP_DIR = here(DATA_DIR,'prep')

# variable
GENES_OI = c("RNF144B")
SAMPLE_TYPES_OI = c("Normal_Tissue","Primary_Tumor_P53WT","Primary_Tumor_P53mut")

# inputs
phenotype_cancer_type_file = here(PREP_DIR,'TTG','phenotype_cancer_type_TTG.tsv')
phenotype_primary_site_file = here(PREP_DIR,'TTG','phenotype_primary_site_TTG.tsv')

genexpr_RNF144B_file = here(TTG_DIR,'TcgaTargetGtex_RSEM_Hugo_norm_count_RNF144B')

# outputs
output_cancer_type_file = here(PREP_DIR,"TTG",'genexpr_RNF144B_cancer_type_TTG.tsv')
output_primary_site_file = here(PREP_DIR,"TTG",'genexpr_RNF144B_primary_site_TTG.tsv')

# load data
## Phenotype data
metadata_cancer_type = read_tsv(phenotype_cancer_type_file)
metadata_primary_site = read_tsv(phenotype_primary_site_file)
samples_oi_cancer_type = metadata_cancer_type %>% pull(unique(sample))
samples_oi_primary_site = metadata_primary_site %>% pull(unique(sample))

## Expression data
### Since I had to modify the dataset to keep only the row with the expression of RNF144B
### I retrieve the sample names from the database independently
samples_TTG_expr_dataset = c("gene",fetch_dataset_samples(host = "https://toil.xenahubs.net",dataset = "TcgaTargetGtex_RSEM_Hugo_norm_count"))
expr_df_RNF144B = fread(genexpr_RNF144B_file)
colnames(expr_df_RNF144B) = samples_TTG_expr_dataset

# preprocess
df_RNF144B_expr_TTG = expr_df_RNF144B %>%
	t() %>%
	as.data.frame() %>%
	filter(!row_number() %in% 1) %>%
	rename(RNF144B = V1) 

## Cancer type
df_RNF144B_expr_cancer_type_TTG = df_RNF144B_expr_TTG %>%
	filter(row.names(.) %in% samples_oi_cancer_type) %>%
	rownames_to_column("sample") %>%
	left_join(metadata_cancer_type, by = "sample", multiple = "all") %>%
	mutate(RNF144B = as.numeric(RNF144B)) %>%
	mutate(mutation_p53_TTG = factor(mutation_p53_TTG, levels = SAMPLE_TYPES_OI))

df_RNF144B_expr_cancer_type_TTG %>%
	summarise(n = mean(RNF144B))
## Primary site
df_RNF144B_expr_primary_site_TTG = df_RNF144B_expr_TTG %>%
	filter(row.names(.) %in% samples_oi_primary_site) %>%
	rownames_to_column("sample") %>%
	left_join(metadata_primary_site, by = "sample", multiple = "all") %>%
	mutate(RNF144B = as.numeric(RNF144B)) %>%
	mutate(mutation_p53_TTG = factor(mutation_p53_TTG, levels = SAMPLE_TYPES_OI))

# save
write_tsv(df_RNF144B_expr_cancer_type_TTG, output_cancer_type_file)
write_tsv(df_RNF144B_expr_primary_site_TTG, output_primary_site_file)

print('Done!')