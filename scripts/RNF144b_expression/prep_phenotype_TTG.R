# 
# Author: Gerard Romero Sola
# Contact: gerard [dot] romero02 [at] estudiant [dot] upf [dot] edu
#
# Script purpose
# --------------
# Download data, standardize column names and clean table. Create phenotype dataset
# Keep only information for those cancer types that have >20 STN and PT samples
# Within the PT samples, the >20 filter is applied to both p53 WT and p53mut 
#

# References
# https://ucsc-xena.gitbook.io/project/how-do-i/tumor-vs-normal
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8841814

# Libraries
require(UCSCXenaTools)
require(data.table)
require(R.utils)
require(readr)
require(tidyverse)
require(here)
require(TCGAutils)

ROOT = here::here()
DATA_DIR = here(ROOT,'data')
TCGA_DIR = here(DATA_DIR,'raw','UCSCXena','TCGA')
TTG_DIR = here(DATA_DIR,'raw','UCSCXena','TTG')
PREP_DIR = here(DATA_DIR,'prep')
RESULTS_DIR = here(ROOT,'output')

# Variables
CATEGORIES_DISEASE_CODES_NOI = c("Controls","Ffpe Pilot Phase Ii","Miscellaneous")
SAMPLE_TYPE_TUMOR_OI = c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood")
SAMPLE_TYPE_NORMAL_OI = c("Solid Tissue Normal","Normal Tissue")
GENES_OI = c("TP53")
EFFECT_NOI = c("Silent","Intron","5'UTR","3'UTR")

THRESH = 20

# Inputs
mutations_TTG_file = here(TTG_DIR,"mc3.v0.2.8.PUBLIC.toil.xena.gz")
mutations_PANCAN_file = here(TCGA_DIR,"mutations","mc3.v0.2.8.PUBLIC.nonsilentGene.xena")
mutations_PANCAN_public_file = here(TCGA_DIR,"mutations","mc3.v0.2.8.PUBLIC.xena.gz")
phenotype_TTG_file = here(TTG_DIR,"TcgaTargetGTEX_phenotype.txt")

# Output
phenotype_cancer_type_TTG_output = here(PREP_DIR,"TTG","phenotype_cancer_type_TTG.tsv")
sample_size_cancer_type_output =  here(PREP_DIR,"TTG","sample_size_cancer_type_TTG.tsv")
phenotype_primary_site_TTG_output = here(PREP_DIR,"TTG","phenotype_primary_site_TTG.tsv")
sample_size_primary_site_output =  here(PREP_DIR,"TTG","sample_size_primary_site_TTG.tsv")

# Load data
## Mutations file from TCGA TARGET GTEx
df_mutations_TTG = read_tsv(mutations_TTG_file)
## Mutations file from TCGA PANCAN
df_mutations_PANCAN = read_tsv(mutations_PANCAN_file)
df_mutations_PANCAN_public = read_tsv(mutations_PANCAN_public_file)
## Disease Codes Dataset containing the abbreviations for the different Study Cancer Types
data("diseaseCodes")
diseaseCodes = diseaseCodes %>%
	select(Study.Abbreviation, Study.Name) %>%
	rename(Abbreviation = Study.Abbreviation,detailed_category=Study.Name)
## Phenotype file from TCGA TARGET GTEx
df_phenotype_TTG = read_tsv(phenotype_TTG_file)

# Clean and subset datasets
## Disease codes: There are bad annotated cancers, so I will change the annotations of some of them

diseaseCodes = diseaseCodes %>%
	mutate(detailed_category = str_to_title(detailed_category)) %>%
	filter(!detailed_category %in% CATEGORIES_DISEASE_CODES_NOI) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category == "Testicular Germ Cell Tumors", "Testicular Germ Cell Tumor")) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category == "Pheochromocytoma And Paraganglioma", "Pheochromocytoma & Paraganglioma")) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category == "Kidney Renal Clear Cell Carcinoma", "Kidney Clear Cell Carcinoma")) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category == "Kidney Renal Papillary Cell Carcinoma", "Kidney Papillary Cell Carcinoma")) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category == "Head And Neck Squamous Cell Carcinoma", "Head & Neck Squamous Cell Carcinoma")) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category == "Lymphoid Neoplasm Diffuse Large B-Cell Lymphoma", "Diffuse Large B-Cell Lymphoma")) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category == "Cervical Squamous Cell Carcinoma And Endocervical Adenocarcinoma", "Cervical & Endocervical Cancer")) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category == "Adrenocortical Carcinoma", "Adrenocortical Cancer")) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category == "Uterine Corpus Endometrial Carcinoma", "Uterine Corpus Endometrioid Carcinoma"))

## Phenotype_TTG: subset only sample types of interest
df_phenotype_TTG_Tumor = df_phenotype_TTG %>%
	filter(`_sample_type` %in% SAMPLE_TYPE_TUMOR_OI) %>%
	mutate(`_primary_site` = replace(`_primary_site`, `_primary_site` == "Adrenal gland", "Adrenal Gland")) %>%
	mutate(`_primary_site` = replace(`_primary_site`, `_primary_site` == "Lymphatic tissue", "Blood")) %>%
	mutate(`_primary_site` = replace(`_primary_site`, `_primary_site` == "White blood cell", "Blood")) %>%
	filter(!`_primary_site` %in% c("Eye","Lining of body cavities")) %>%
	mutate(detailed_category = replace(detailed_category, detailed_category =="Acute Myeloid Leukemia, Induction Failure Subproject","Acute Myeloid Leukemia"))
df_phenotype_TTG_Tumor = merge(df_phenotype_TTG_Tumor,diseaseCodes,all.x = TRUE, by = "detailed_category")
df_phenotype_TTG_Tumor = df_phenotype_TTG_Tumor %>%
	drop_na(Abbreviation)

## Manual association of _primary_site with type of cancer
cancer_type_site = df_phenotype_TTG_Tumor %>%
	select(Abbreviation,`_primary_site`) %>%
	filter(!duplicated(Abbreviation))

df_phenotype_TTG_Normal = df_phenotype_TTG %>%
	filter(`_sample_type` %in% SAMPLE_TYPE_NORMAL_OI) %>%
	mutate(`_primary_site` = replace(`_primary_site`, `_primary_site` == "Cervix Uteri", "Cervix")) %>%
	mutate(`_primary_site` = replace(`_primary_site`, `_primary_site` == "Fallopian Tube", "Ovary")) %>%
	mutate(`_primary_site` = replace(`_primary_site`, `_primary_site` == "Small Intestine", "Colon")) %>%
	mutate(`_primary_site` = replace(`_primary_site`, `_primary_site` == "Thyroid", "Thyroid Gland")) %>%
	mutate(`_primary_site` = replace(`_primary_site`, `_primary_site` == "Muscle", "Soft tissue,Bone")) %>%
	drop_na(`_primary_site`)
df_phenotype_TTG_Normal = merge(df_phenotype_TTG_Normal,cancer_type_site,all.x = TRUE, by="_primary_site")
df_phenotype_TTG_Normal = df_phenotype_TTG_Normal %>%
	drop_na(Abbreviation)

## Mutations_TTG: subset only P53
p53_mut_TTG_samples = df_mutations_TTG %>%
	filter(gene %in% GENES_OI) %>%
	filter(!effect %in% EFFECT_NOI) %>%
	pull(unique(sample))
p53_WT_TTG_samples = df_mutations_TTG %>%
	pull(unique(sample)) %>%
	setdiff(.,p53_mut_TTG_samples)

patients_with_p53_mut_dataset = df_mutations_TTG %>%
	filter(gene %in% GENES_OI) %>%
	filter(!effect %in% EFFECT_NOI) %>%
	dplyr::select(-c(SIFT,PolyPhen)) %>%
	mutate(Status = "MUT")
patients_with_p53_WT_dataset = df_mutations_TTG %>%
	filter(gene %in% GENES_OI) %>%
	filter(!sample %in% p53_mut_TTG_samples) %>%
	dplyr::select(-c(SIFT,PolyPhen)) %>%
	mutate(Status = "WT")
patients_with_p53_status_dataset = rbind(patients_with_p53_mut_dataset,patients_with_p53_WT_dataset)
write_xlsx(
	x = patients_with_p53_status_dataset,	path = here(PREP_DIR,"TTG","Dataset_of_TP53_MUT_samples.xlsx")
)

## Mutations_PANCAN: 
df_mutations_PANCAN = df_mutations_PANCAN %>%
	filter(sample %in% GENES_OI) %>%
	t() %>%
	as.data.frame() %>%
	filter(!row_number() %in% 1) %>%
	rownames_to_column("sample")
colnames(df_mutations_PANCAN) = c("TP53_mut")
p53_mut_PANCAN_samples = df_mutations_PANCAN$sample[df_mutations_PANCAN$TP53_mut == 1]
p53_WT_PANCAN_samples = df_mutations_PANCAN$sample[df_mutations_PANCAN$TP53_mut == 0]

# Combining datasets
df_phenotype_TTG = full_join(df_phenotype_TTG_Normal,df_phenotype_TTG_Tumor) %>%
	filter(!`_study` == "TARGET") %>%
	mutate(mutation_p53_TTG = if_else(sample %in% p53_mut_TTG_samples, "Primary_Tumor_P53mut",
																		if_else(sample %in% p53_WT_TTG_samples, "Primary_Tumor_P53WT",
																						if_else(sample %in% df_phenotype_TTG_Normal$sample, "Normal_Tissue",NA)))) %>%
	drop_na()

# Filtering
## Filtering by cancer type
sample_size_cancer_types = df_phenotype_TTG %>%
	select(c(sample,Abbreviation,mutation_p53_TTG)) %>%
	group_by(Abbreviation, mutation_p53_TTG) %>%
	summarize(n = n()) %>%
	filter( n >= THRESH )

cancer_types_oi_p53_TTG = sample_size_cancer_types %>%
	filter(duplicated(Abbreviation)) %>%
	filter(duplicated(Abbreviation)) %>%
	pull(Abbreviation)

df_phenotype_cancer_type_TTG = df_phenotype_TTG %>%
	select(c(sample,Abbreviation,mutation_p53_TTG,`_primary_site`)) %>%
	filter(Abbreviation %in% cancer_types_oi_p53_TTG)

## Filtering by primary site
sample_size_primary_sites = df_phenotype_TTG %>%
	select(c(sample,`_primary_site`,mutation_p53_TTG)) %>%
	group_by(`_primary_site`, mutation_p53_TTG) %>%
	summarize(n = n()) %>%
	filter( n >= THRESH )
	
primary_sites_oi_p53_TTG = sample_size_primary_sites %>%
	filter(duplicated(`_primary_site`)) %>%
	filter(duplicated(`_primary_site`)) %>%
	pull(`_primary_site`)

df_phenotype_primary_site_TTG = df_phenotype_TTG %>%
	select(c(sample,Abbreviation,mutation_p53_TTG,`_primary_site`)) %>%
	filter(`_primary_site` %in% primary_sites_oi_p53_TTG)

# save
write_tsv(df_phenotype_cancer_type_TTG, phenotype_cancer_type_TTG_output)
write_tsv(sample_size_cancer_types, sample_size_cancer_type_output)
write_tsv(df_phenotype_primary_site_TTG, phenotype_primary_site_TTG_output)
write_tsv(sample_size_primary_sites, sample_size_primary_site_output)