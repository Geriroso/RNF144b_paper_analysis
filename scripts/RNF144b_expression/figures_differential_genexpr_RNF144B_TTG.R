# 
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Compare gene expression of RNF144B in TCGA.
#

require(tidyverse)
require(limma)
require(magrittr)
require(ggpubr)
require(writexl)
require(latex2exp)
require(here)
require(rstatix)

ROOT = here::here()
DATA_DIR = here(ROOT,'data')
TCGA_DIR = here(DATA_DIR,'raw','UCSCXena','TCGA')
PREP_DIR = here(DATA_DIR,'prep')
RESULTS_DIR = here(ROOT,'output')

# variables
CANCERS_OI = c("BLCA","BRCA","COAD","ESCA","GBM","HNSC","KICH","LGG","LIHC","LUAD","LUSC","PAAD","PRAD","SARC","STAD")
PRIMARY_SITE_OI = c("Bladder","Brain","Breast","Colon","Esophagus","Head&Neck",
										"Kidney","Liver","Lung","Pancreas","Prostate","Soft T.&Bone","Stomach")
SAMPLE_TYPES_OI = c("Normal_Tissue","Primary_Tumor_P53WT","Primary_Tumor_P53mut")

FONT_SIZE = 7 # pt
FONT_FAMILY = 'helvetica'

# inputs
genexpr_cancer_type_file = here(PREP_DIR,"TTG","genexpr_RNF144B_cancer_type_TTG.tsv")
genexpr_primary_site_file = here(PREP_DIR,"TTG","genexpr_RNF144B_primary_site_TTG.tsv")

# outputs
diffexpr_RNF144B_cancer_type_Normal_vs_Tumor_file = here(RESULTS_DIR,'figures','RNF144B-differential_gene_expression-barplot_cancer_types_Normal_Tumor_TTG.pdf')
output_figdata_cancer_type_Normal_vs_Tumor = here(RESULTS_DIR,'files','Normal_Tumor_Cancer_Types_sourcedata_TTG.xlsx')
diffexpr_RNF144B_cancer_type_Normal_WT_mut_file = here(RESULTS_DIR,'figures','RNF144B-differential_gene_expression-barplot_cancer_types_Normal_WT_mut_TTG.pdf')
output_figdata_cancer_type_Normal_WT_mut = here(RESULTS_DIR,'files','Normal_WT_mut_Cancer_Types_sourcedata_TTG.xlsx')
diffexpr_RNF144B_cancer_type_WT_vs_mut_file = here(RESULTS_DIR,'figures','RNF144B-differential_gene_expression-barplot_cancer_types_WT_vs_mut_TTG.pdf')
output_figdata_cancer_type_WT_vs_mut = here(RESULTS_DIR,'files','WT_vs_mut_Cancer_Types_sourcedata_TTG.xlsx')
diffexpr_RNF144B_cancer_type_PANCAN_Normal_vs_Tumor_file = here(RESULTS_DIR,'figures','RNF144B-differential_gene_expression-barplot_cancer_types_PANCAN_Normal_Tumor_TTG.pdf')
output_figdata_cancer_type_PANCAN_Normal_vs_Tumor = here(RESULTS_DIR,'files','Normal_Tumor_Cancer_Types_PANCAN_sourcedata_TTG.xlsx')
diffexpr_RNF144B_cancer_type_PANCAN_WT_vs_mut_file = here(RESULTS_DIR,'figures','RNF144B-differential_gene_expression-barplot_cancer_types_PANCAN_WT_vs_mut_TTG.pdf')
output_figdata_cancer_type_PANCAN_WT_vs_mut = here(RESULTS_DIR,'files','WT_vs_mut_Cancer_Types_PANCAN_sourcedata_TTG.xlsx')
diffexpr_RNF144B_cancer_type_PANCAN_Normal_WT_mut_file = here(RESULTS_DIR,'figures','RNF144B-differential_gene_expression-barplot_cancer_types_PANCAN_Normal_WT_mut_TTG.pdf')
output_figdata_cancer_type_PANCAN_Normal_WT_mut = here(RESULTS_DIR,'files','Normal_WT_mut_Cancer_Types_PANCAN_sourcedata_TTG.xlsx')


diffexpr_RNF144B_primary_site_Normal_vs_Tumor_file = here(RESULTS_DIR,'figures','RNF144B-differential_gene_expression-barplot_primary_sites_Normal_vs_Tumor_TTG.pdf')
output_figdata_primary_site_Normal_vs_Tumor = here(RESULTS_DIR,'files','Normal_vs_Tumor_primary_sites_sourcedata_TTG.xlsx')
diffexpr_RNF144B_primary_site_Normal_WT_mut_file = here(RESULTS_DIR,'figures','RNF144B-differential_gene_expression-barplot_primary_sites_Normal_WT_mut_TTG.pdf')
output_figdata_primary_site_Normal_WT_mut = here(RESULTS_DIR,'files','Normal_WT_mut_primary_sites_sourcedata_TTG.xlsx')
diffexpr_RNF144B_primary_site_WT_vs_mut_file = here(RESULTS_DIR,'figures','RNF144B-differential_gene_expression-barplot_primary_sites_WT_vs_mut_TTG.pdf')
output_figdata_primary_site_WT_vs_mut = here(RESULTS_DIR,'files','WT_vs_mut_primary_sites_sourcedata_TTG.xlsx')

# load data
df_RNF144B_expr_cancer_type = read_tsv(genexpr_cancer_type_file)
df_RNF144B_expr_primary_site = read_tsv(genexpr_primary_site_file)

# change sample_type and cancer_type order
df_RNF144B_expr_cancer_type = df_RNF144B_expr_cancer_type %>%
	filter(RNF144B != 0) %>%
	mutate(mutation_p53_TTG = factor(mutation_p53_TTG, levels = SAMPLE_TYPES_OI),
				 Abbreviation = factor(Abbreviation, levels = CANCERS_OI)) %>%
	mutate(gene = "RNF144B") %>%
	mutate(scaled_expression = scale(RNF144B, scale = FALSE)[,1])
	
df_RNF144B_expr_primary_site = df_RNF144B_expr_primary_site %>%
	filter(RNF144B != 0) %>% 
	mutate(`_primary_site` = replace(`_primary_site`, `_primary_site` == "Head and Neck region", "Head&Neck"),
				 `_primary_site` = replace(`_primary_site`, `_primary_site` == "Soft tissue,Bone", "Soft T.&Bone")) %>%
	mutate(mutation_p53_TTG = factor(mutation_p53_TTG, levels = SAMPLE_TYPES_OI),
				 Primary_site = factor(`_primary_site`, levels = PRIMARY_SITE_OI)) %>%
	dplyr::select(-`_primary_site`) %>%
	mutate(gene = "RNF144B") %>%
	mutate(scaled_expression = scale(RNF144B, scale = FALSE)[,1])

# Add normal or tumor column
df_RNF144B_expr_cancer_type = df_RNF144B_expr_cancer_type %>%
	mutate(condition = if_else(mutation_p53_TTG == "Normal_Tissue", "Normal","Tumor"))
df_RNF144B_expr_primary_site = df_RNF144B_expr_primary_site %>%
	mutate(condition = if_else(mutation_p53_TTG == "Normal_Tissue", "Normal","Tumor"))

comparisons_Normal_WT_mut = list(c("Normal_Tissue", "Primary_Tumor_P53WT"),c("Normal_Tissue", "Primary_Tumor_P53mut"), c("Primary_Tumor_P53WT","Primary_Tumor_P53mut"))
comparisons_WT_mut = list(c("Primary_Tumor_P53WT","Primary_Tumor_P53mut"))
comparisons_Normal_Tumor = list(c("Normal","Tumor"))

df_plot_WT_mut_cancer_type = df_RNF144B_expr_cancer_type %>%
	filter(mutation_p53_TTG != "Normal_Tissue") %>%
	mutate(mutation_p53_TTG = factor(mutation_p53_TTG, levels = c("Primary_Tumor_P53WT","Primary_Tumor_P53mut")))
df_plot_WT_mut_primary_site = df_RNF144B_expr_primary_site %>%
	filter(mutation_p53_TTG != "Normal_Tissue") %>%
	mutate(mutation_p53_TTG = factor(mutation_p53_TTG, levels = c("Primary_Tumor_P53WT","Primary_Tumor_P53mut")))

# for (cancer in CANCERS_OI_NORMAL_VS_TUMOR) {
# 	a = df_RNF144B_expr_cancer_type_Normal_vs_Tumor %>%
# 		filter(Abbreviation == cancer)
# 	print(cancer)
# 	print(var.test(scaled_expression ~ condition, data = a))
# }

stat_test_Normal_Tumor_cancer_type = df_RNF144B_expr_cancer_type %>%
	group_by(Abbreviation) %>%
	pairwise_t_test(scaled_expression ~ condition, p.adjust.method = "BH", ref.group = "Normal", alternative = "two.sided") %>%
	add_xy_position(x = "Abbreviation") %>%
	mutate(y.position = 5) %>%
	mutate(xmax = xmax -0.2)
stat_test_Normal_Tumor_cancer_type_PANCAN = df_RNF144B_expr_cancer_type %>%
	mutate(Abbreviation = factor("PANCAN")) %>%
	group_by(Abbreviation) %>%
	pairwise_t_test(scaled_expression ~ condition, p.adjust.method = "BH", ref.group = "Normal", alternative = "two.sided") %>%
	add_xy_position(x = "Abbreviation") %>%
	mutate(y.position = 5) %>%
	mutate(xmax = xmax -0.2)
stat_test_Normal_WT_mut_cancer_type = df_RNF144B_expr_cancer_type %>%
	group_by(Abbreviation) %>%
	pairwise_t_test(scaled_expression ~ mutation_p53_TTG, p.adjust.method = "BH", alternative = "two.sided") %>%
	add_xy_position(x = "Abbreviation") %>%
	mutate(y.position = rep(c(4.5,5,5.5), 15))
stat_test_Normal_WT_mut_cancer_type_PANCAN = df_RNF144B_expr_cancer_type %>%
	mutate(Abbreviation = factor("PANCAN")) %>%
	group_by(Abbreviation) %>%
	pairwise_t_test(scaled_expression ~ mutation_p53_TTG, p.adjust.method = "BH", alternative = "two.sided") %>%
	add_xy_position(x = "Abbreviation") %>%
	mutate(y.position = c(4.5,5,5.5))
stat_test_WT_mut_cancer_type = df_plot_WT_mut_cancer_type %>%
	group_by(Abbreviation) %>%
	pairwise_t_test(scaled_expression ~ mutation_p53_TTG, p.adjust.method = "BH", ref.group = "Primary_Tumor_P53WT",alternative = "two.sided") %>%
	add_xy_position(x = "Abbreviation") %>%
	mutate(y.position = 5) %>%
	mutate(xmax = xmax -0.2)
stat_test_WT_mut_cancer_type_PANCAN = df_plot_WT_mut_cancer_type %>%
	mutate(Abbreviation = factor("PANCAN")) %>%
	group_by(Abbreviation) %>%
	pairwise_t_test(scaled_expression ~ mutation_p53_TTG, p.adjust.method = "BH", ref.group = "Primary_Tumor_P53WT",alternative = "two.sided") %>%
	add_xy_position(x = "Abbreviation") %>%
	mutate(y.position = 5) %>%
	mutate(xmax = xmax -0.2)


stat_test_Normal_Tumor_primary_site = df_RNF144B_expr_primary_site %>%
	group_by(Primary_site) %>%
	pairwise_t_test(scaled_expression ~ condition, p.adjust.method = "BH", ref.group = "Normal", alternative = "two.sided") %>%
	add_xy_position(x = "Primary_site") %>%
	mutate(y.position = 5) %>%
	mutate(xmax = xmax -0.2)
stat_test_Normal_WT_mut_primary_site = df_RNF144B_expr_primary_site %>%
	group_by(Primary_site) %>%
	pairwise_t_test(scaled_expression ~ mutation_p53_TTG, p.adjust.method = "BH", alternative = "two.sided") %>%
	add_xy_position(x = "Primary_site") %>%
	mutate(y.position = rep(c(4.5,5,5.5), 13))
stat_test_WT_mut_primary_site = df_plot_WT_mut_primary_site %>%
	group_by(Primary_site) %>%
	pairwise_t_test(scaled_expression ~ mutation_p53_TTG, p.adjust.method = "BH", ref.group = "Primary_Tumor_P53WT", alternative = "two.sided") %>%
	add_xy_position(x = "Primary_site") %>%
	mutate(y.position = 5) %>%
	mutate(xmax = xmax -0.2)

## Cancer types
plot_RNF144B_cancer_types_Normal_Tumor = ggboxplot(df_RNF144B_expr_cancer_type,x = "Abbreviation", y = "scaled_expression", fill = "condition") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5, size = 12),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12)) + 
	scale_fill_manual(values = c("#9E9EF2","#824B53"), name = "Condition") +
	ylab("Rnf144b expression") +
	scale_y_continuous(limits = c(-5,5)) +
	stat_pvalue_manual(stat_test_Normal_Tumor_cancer_type, label = "p.signif", size = 10, remove.bracket = TRUE)

plot_RNF144B_cancer_types_Normal_Tumor_PANCAN = ggboxplot(df_RNF144B_expr_cancer_type %>% mutate(Abbreviation = factor("PANCAN")),x = "Abbreviation", y = "scaled_expression", fill = "condition") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5, size = 12),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12)) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "Condition") +
	ylab("Rnf144b expression") +
	scale_y_continuous(limits = c(-5,5)) +
	stat_pvalue_manual(stat_test_Normal_Tumor_cancer_type_PANCAN, label = "p.signif", size = 10, remove.bracket = TRUE)

# save figure
ggsave(plot_RNF144B_cancer_types_Normal_Tumor_PANCAN, filename=diffexpr_RNF144B_cancer_type_PANCAN_Normal_vs_Tumor_file, width = 15, height = 10)


sample_summary_cancer_type_Normal_Tumor = df_RNF144B_expr_cancer_type %>%
	group_by(Abbreviation, condition, gene) %>%
	summarize(n = n(),
						median = median(scaled_expression),
						mean = mean(scaled_expression),
						std = sd(scaled_expression),
						mad = mad(scaled_expression))
stat_test_Normal_Tumor_cancer_type_source_data = stat_test_Normal_Tumor_cancer_type %>%
	dplyr::select(Abbreviation,group1,group2,n1,n2,p,p.adj,p.signif)

write_xlsx(
	x = list(sample_summary_cancer_type = sample_summary_cancer_type_Normal_Tumor, test_result = stat_test_Normal_Tumor_cancer_type_source_data),
	path = output_figdata_cancer_type_Normal_vs_Tumor
)

# save figure
ggsave(plot_RNF144B_cancer_types_Normal_Tumor, filename=diffexpr_RNF144B_cancer_type_Normal_vs_Tumor_file, width = 15, height = 10)

plot_RNF144B_cancer_types_Normal_WT_Mut = ggboxplot(df_RNF144B_expr_cancer_type,x = "Abbreviation", y = "scaled_expression", fill = "mutation_p53_TTG") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5, size = 12),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12)) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2","#824B53"), name = "Condition") +
	ylab("Rnf144b expression") +
	scale_y_continuous(limits = c(-5,5.5), breaks = c(-5,-2.5,0,2.5,5)) +
	stat_pvalue_manual(stat_test_Normal_WT_mut_cancer_type, label = "p.signif", size = 5, remove.bracket = FALSE, bracket.size = 0.3, tip.length = 0.005)

sample_summary_cancer_type_Normal_WT_Mut = df_RNF144B_expr_cancer_type %>%
	group_by(Abbreviation, mutation_p53_TTG, gene) %>%
	summarize(n = n(),
						median = median(scaled_expression),
						mean = mean(scaled_expression),
						std = sd(scaled_expression),
						mad = mad(scaled_expression))
stat_test_Normal_WT_mut_cancer_type_source_data = stat_test_Normal_WT_mut_cancer_type %>%
	dplyr::select(Abbreviation,group1,group2,n1,n2,p,p.adj,p.signif)

write_xlsx(
	x = list(sample_summary_cancer_type = sample_summary_cancer_type_Normal_WT_Mut, test_result = stat_test_Normal_WT_mut_cancer_type_source_data),
	path = output_figdata_cancer_type_Normal_WT_mut
)

# save figure
ggsave(plot_RNF144B_cancer_types_Normal_WT_Mut, filename=diffexpr_RNF144B_cancer_type_Normal_WT_mut_file, width = 15, height = 10)

plot_RNF144B_cancer_types_Normal_WT_Mut_PANCAN = ggboxplot(df_RNF144B_expr_cancer_type %>% mutate(Abbreviation = factor("PANCAN")),x = "Abbreviation", y = "scaled_expression", fill = "mutation_p53_TTG") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5, size = 12),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12)) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2","#824B53"), name = "Condition") +
	ylab("Rnf144b expression") +
	scale_y_continuous(limits = c(-5,5.5), breaks = c(-5,-2.5,0,2.5,5)) +
	stat_pvalue_manual(stat_test_Normal_WT_mut_cancer_type_PANCAN, label = "p.signif", size = 5, remove.bracket = FALSE, bracket.size = 0.3, tip.length = 0.005)

# save figure
ggsave(plot_RNF144B_cancer_types_Normal_WT_Mut_PANCAN, filename=diffexpr_RNF144B_cancer_type_PANCAN_Normal_WT_mut_file, width = 15, height = 10)


plot_RNF144B_cancer_types_WT_mut = ggboxplot(df_plot_WT_mut_cancer_type,x = "Abbreviation", y = "scaled_expression", fill = "mutation_p53_TTG") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5, size = 12),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12)) + 
	scale_fill_manual(values = c("#9E9EF2","#824B53"), name = "Condition") +
	ylab("Rnf144b expression") +
	scale_y_continuous(limits = c(-5,5)) +
	stat_pvalue_manual(stat_test_WT_mut_cancer_type, label = "p.signif", size = 10, remove.bracket = TRUE)

sample_summary_cancer_type_WT_Mut = df_plot_WT_mut_cancer_type %>%
	group_by(Abbreviation, mutation_p53_TTG, gene) %>%
	summarize(n = n(),
						median = median(scaled_expression),
						mean = mean(scaled_expression),
						std = sd(scaled_expression),
						mad = mad(scaled_expression))
stat_test_WT_mut_cancer_type_source_data = stat_test_WT_mut_cancer_type %>%
	dplyr::select(Abbreviation,group1,group2,n1,n2,p,p.adj,p.signif)

write_xlsx(
	x = list(sample_summary_cancer_type = sample_summary_cancer_type_WT_Mut, test_result = stat_test_WT_mut_cancer_type_source_data),
	path = output_figdata_cancer_type_WT_vs_mut
)

# save figure
ggsave(plot_RNF144B_cancer_types_WT_mut, filename=diffexpr_RNF144B_cancer_type_WT_vs_mut_file, width = 15, height = 10)

plot_RNF144B_cancer_types_WT_mut_PANCAN = ggboxplot(df_plot_WT_mut_cancer_type %>% mutate(Abbreviation = factor("PANCAN")),x = "Abbreviation", y = "scaled_expression", fill = "mutation_p53_TTG") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5, size = 12),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12)) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "Condition") +
	ylab("Rnf144b expression") +
	scale_y_continuous(limits = c(-5,5)) +
	stat_pvalue_manual(stat_test_WT_mut_cancer_type_PANCAN, label = "p.signif", size = 10, remove.bracket = TRUE)

sample_summary_cancer_type_WT_Mut_PANCAN = df_plot_WT_mut_cancer_type %>% mutate(Abbreviation = factor("PANCAN")) %>%
	group_by(Abbreviation, mutation_p53_TTG, gene) %>%
	summarize(n = n(),
						median = median(scaled_expression),
						mean = mean(scaled_expression),
						std = sd(scaled_expression),
						mad = mad(scaled_expression))
stat_test_WT_mut_cancer_type_PANCAN_source_data = stat_test_WT_mut_cancer_type_PANCAN %>%
	dplyr::select(Abbreviation,group1,group2,n1,n2,p,p.adj,p.signif)

write_xlsx(
	x = list(sample_summary_cancer_type = sample_summary_cancer_type_WT_Mut_PANCAN, test_result = stat_test_WT_mut_cancer_type_PANCAN_source_data),
	path = output_figdata_cancer_type_PANCAN_WT_vs_mut
)

# save figure
ggsave(plot_RNF144B_cancer_types_WT_mut_PANCAN, filename=diffexpr_RNF144B_cancer_type_PANCAN_WT_vs_mut_file, width = 15, height = 10)


## Primary sites
plot_RNF144B_primary_sites_Normal_Tumor = ggboxplot(df_RNF144B_expr_primary_site,x = "Primary_site", y = "scaled_expression", fill = "condition") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5, size = 12),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12)) + 
	scale_fill_manual(values = c("#9E9EF2","#824B53"), name = "Condition") +
	ylab("Rnf144b expression") +
	scale_y_continuous(limits = c(-5,5)) +
	stat_pvalue_manual(stat_test_Normal_Tumor_primary_site, label = "p.signif", size = 10, remove.bracket = TRUE)

sample_summary_primary_site_Normal_Tumor = df_RNF144B_expr_primary_site %>%
	group_by(Primary_site, condition, gene) %>%
	summarize(n = n(),
						median = median(scaled_expression),
						mean = mean(scaled_expression),
						std = sd(scaled_expression),
						mad = mad(scaled_expression))
stat_test_Normal_Tumor_primary_sites_source_data = stat_test_Normal_Tumor_primary_site %>%
	dplyr::select(Primary_site,group1,group2,n1,n2,p,p.adj,p.signif)

write_xlsx(
	x = list(sample_summary_primary_sites = sample_summary_primary_site_Normal_Tumor, test_result = stat_test_Normal_Tumor_primary_sites_source_data),
	path = output_figdata_primary_site_Normal_vs_Tumor
)

# save figure
ggsave(plot_RNF144B_primary_sites_Normal_Tumor, filename=diffexpr_RNF144B_primary_site_Normal_vs_Tumor_file, width = 15, height = 10)


plot_RNF144B_primary_sites_Normal_WT_Mut = ggboxplot(df_RNF144B_expr_primary_site,x = "Primary_site", y = "scaled_expression", fill = "mutation_p53_TTG") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5, size = 12),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12)) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2","#824B53"), name = "Condition") +
	ylab("Rnf144b expression") +
	scale_y_continuous(limits = c(-5,5.5), breaks = c(-5,-2.5,0,2.5,5)) +
	stat_pvalue_manual(stat_test_Normal_WT_mut_primary_site, label = "p.signif", size = 5, remove.bracket = FALSE, bracket.size = 0.3, tip.length = 0.005)

sample_summary_primary_site_Normal_WT_Mut = df_RNF144B_expr_primary_site %>%
	group_by(Primary_site, mutation_p53_TTG, gene) %>%
	summarize(n = n(),
						median = median(scaled_expression),
						mean = mean(scaled_expression),
						std = sd(scaled_expression),
						mad = mad(scaled_expression))
stat_test_Normal_WT_Mut_primary_sites_source_data = stat_test_Normal_WT_mut_primary_site %>%
	dplyr::select(Primary_site,group1,group2,n1,n2,p,p.adj,p.signif)

write_xlsx(
	x = list(sample_summary_primary_sites = sample_summary_primary_site_Normal_WT_Mut, test_result = stat_test_Normal_WT_Mut_primary_sites_source_data),
	path = output_figdata_primary_site_Normal_WT_mut
)

# save figure
ggsave(plot_RNF144B_primary_sites_Normal_WT_Mut, filename=diffexpr_RNF144B_primary_site_Normal_WT_mut_file, width = 15, height = 10)

plot_RNF144B_primary_sites_WT_Mut = ggboxplot(df_plot_WT_mut_primary_site,x = "Primary_site", y = "RNF144B", fill = "mutation_p53_TTG") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5),
				axis.title.x = element_blank(),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 10, face = "bold")) + 
	scale_fill_manual(values = c("#9E9EF2","#824B53")) +
	ylab("log2(Norm.Counts+1)") +
	scale_y_continuous(limits = c(4,18)) +
	stat_pvalue_manual(stat_test_WT_mut_primary_site, label = "p.adj.signif", size = 6, bracket.size = 0)

plot_RNF144B_primary_sites_WT_Mut = ggboxplot(df_plot_WT_mut_primary_site,x = "Primary_site", y = "scaled_expression", fill = "mutation_p53_TTG") +
	theme_pubr() +
	theme(axis.text.x = element_text(angle = -40, vjust = 0.5, hjust = 0.5, size = 12),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12)) + 
	scale_fill_manual(values = c("#9E9EF2","#824B53"), name = "Condition") +
	ylab("Rnf144b expression") +
	scale_y_continuous(limits = c(-5,5)) +
	stat_pvalue_manual(stat_test_WT_mut_primary_site, label = "p.signif", size = 10, remove.bracket = TRUE)

sample_summary_primary_site_WT_Mut = df_plot_WT_mut_primary_site %>%
	group_by(Primary_site, mutation_p53_TTG, gene) %>%
	summarize(n = n(),
						median = median(scaled_expression),
						mean = mean(scaled_expression),
						std = sd(scaled_expression),
						mad = mad(scaled_expression))
stat_test_WT_Mut_primary_sites_source_data = stat_test_WT_mut_primary_site %>%
	dplyr::select(Primary_site,group1,group2,n1,n2,p,p.adj,p.signif)

write_xlsx(
	x = list(sample_summary_primary_sites = sample_summary_primary_site_WT_Mut, test_result = stat_test_WT_Mut_primary_sites_source_data),
	path = output_figdata_primary_site_WT_vs_mut
)

# save figure
ggsave(plot_RNF144B_primary_sites_WT_Mut, filename=diffexpr_RNF144B_primary_site_WT_vs_mut_file, width = 15, height = 10)

print('Done!')## figure data