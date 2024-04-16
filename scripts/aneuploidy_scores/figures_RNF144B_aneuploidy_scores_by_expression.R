# 
# Author: Gerard Romero Sola
# Inspiration: Miquel Anglada Girotto
# Github inspiration/copy: https://github.com/MiqG/publication_zadra_ttll11
#
# Script purpose
# --------------
# Prepare RNF144b and p53 Expression

require(tidyverse)
library(hrbrthemes)
library(viridis)
require(rstatix)
require(ggpubr)
require(here)

ROOT = here()
DATA_DIR = here(ROOT,'data')
PREP_DIR = here(DATA_DIR,'prep')

# variables
CANCERS_OI = c('LUSC','UCEC','BLCA','BRCA','STAD','LUAD','KIRP','THCA','KICH','COAD','LIHC','HNSC','PRAD','KIRC')

# inputs
RNF144B_expr_file = here(PREP_DIR,'genexpr_RNF144B_aneuploidy_by_RNF144B_expression.tsv')
aneuploidy_scores_file = here(PREP_DIR,"aneuploidy.tsv")

# outputs
figure = here(ROOT, "output", "Violin_plot_RNF144B_expression.pdf")

# load data
RNF144B_expr = read_tsv(RNF144B_expr_file)
aneuploidy_scores = read_tsv(aneuploidy_scores_file)

# merge files
df = merge(RNF144B_expr,aneuploidy_scores, by = "sample")

# Statify by RNF144B expression
df = df %>%
	mutate(RNF144B_expression = if_else(expression > median(df$expression), "High", "Low"),
				 RNF144B_expression = factor(RNF144B_expression, levels = c("High","Low")))

# Select_groups_for_all_tumors
CANCERS_OI = df %>%
	dplyr::select(c(cancer_type, aneuploidy_score, RNF144B_expression)) %>%
	group_by(cancer_type, RNF144B_expression) %>%
	summarise(n = n()) %>%
	filter(n > 20) %>%
	filter(duplicated(cancer_type)) %>%
	pull(cancer_type)

CANCERS_OI_P53_WT = df %>%
	dplyr::select(c(cancer_type, aneuploidy_score, RNF144B_expression, TP53_mut)) %>%
	filter(TP53_mut == 0) %>%
	group_by(cancer_type, RNF144B_expression) %>%
	summarise(n = n()) %>%
	filter(n > 20) %>%
	filter(duplicated(cancer_type)) %>%
	pull(cancer_type)

CANCERS_OI_P53_mut = df %>%
	dplyr::select(c(cancer_type, aneuploidy_score, RNF144B_expression, TP53_mut)) %>%
	filter(TP53_mut == 1) %>%
	group_by(cancer_type, RNF144B_expression) %>%
	summarise(n = n()) %>%
	filter(n > 20) %>%
	filter(duplicated(cancer_type)) %>%
	pull(cancer_type)

stat_test_cancer_type = df %>%
	filter(cancer_type %in% CANCERS_OI) %>%
	dplyr::group_by(cancer_type) %>%
	pairwise_t_test(aneuploidy_score ~ RNF144B_expression, p.adjust.method = "BH", ref.group = "Low", alternative = "two.sided") %>%
	add_xy_position(x = "Abbreviation") %>%
	mutate(y.position = 40) %>%
	mutate(xmax = xmax -0.2)

df %>%
	filter(cancer_type %in% CANCERS_OI) %>%
	ggplot(aes(x = RNF144B_expression , y = aneuploidy_score , fill = RNF144B_expression)) + 
	geom_violin(width=1) +
	geom_boxplot(width=0.1, color="grey", alpha=0.2) +
	facet_wrap(~cancer_type , scales = "free" )

# Plot all tumors
### PANCAN
plot_aneuploidy_by_RNF144b_expression_PANCAN_violin =
	df %>%
	filter(cancer_type %in% CANCERS_OI) %>%
	ggviolin(.,x = "RNF144B_expression", y = "aneuploidy_score", fill = "RNF144B_expression") +
	stat_compare_means(method = "t.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")), comparisons = list(c("High", "Low")), bracket.size = NA, vjust = 0.5, size = 10) +
	stat_summary(fun= "median", mapping = aes(group = RNF144B_expression),geom = "crossbar", width = c(0.663,0.47), color = "black") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(4, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,45))

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_PANCAN_violin, filename=here(ROOT, "output", "Violin_PANCAN_RNF144B_expression.pdf"), device = "pdf", width = 15, height = 10)

stat_test_RNF144b_expression_PANCAN_violin = df %>%
	filter(cancer_type %in% CANCERS_OI) %>%
	pairwise_t_test(aneuploidy_score ~ RNF144B_expression, p.adjust.method = "BH", ref.group = "Low", alternative = "two.sided")

### PANCAN_p53WT
plot_aneuploidy_by_RNF144b_expression_PANCAN_P53_WT_violin = 
	df %>%
	filter(cancer_type %in% CANCERS_OI_P53_WT) %>%
	filter(TP53_mut == 0) %>%
	ggviolin(.,x = "RNF144B_expression", y = "aneuploidy_score", fill = "RNF144B_expression") +
	stat_compare_means(method = "t.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")), comparisons = list(c("High", "Low")), bracket.size = NA, vjust = 0.5, size = 10) +
	stat_summary(fun= "median", mapping = aes(group = RNF144B_expression),geom = "crossbar", width = c(0.657,0.493), color = "black") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(4, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,45))

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_PANCAN_P53_WT_violin, filename=here(ROOT, "output", "Violin_PANCAN_P53WT_RNF144B_expression.pdf"), device = "pdf", width = 15, height = 10)

### PANCAN_p53mut
plot_aneuploidy_by_RNF144b_expression_PANCAN_P53mut_violin =
	df %>%
	filter(cancer_type %in% CANCERS_OI_P53_mut) %>%
	filter(TP53_mut == 1) %>%
	ggviolin(.,x = "RNF144B_expression", y = "aneuploidy_score", fill = "RNF144B_expression") +
	stat_compare_means(method = "t.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")), comparisons = list(c("High", "Low")), bracket.size = NA, vjust = 0.5, size = 10) +
	stat_summary(fun= "median", mapping = aes(group = RNF144B_expression),geom = "crossbar", width = c(0.78,0.97), color = "black") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(4, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,45))

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_PANCAN_P53mut_violin, filename=here(ROOT, "output", "Violin_PANCAN_P53mut_RNF144B_expression.pdf"), device = "pdf", width = 15, height = 10)

### By facets
### Independent of P53 status
plot_aneuploidy_by_RNF144b_expression_all_tumors_facets_violin = 
df %>%
	filter(cancer_type %in% CANCERS_OI) %>%
	ggviolin(.,x = "RNF144B_expression", y = "aneuploidy_score", fill = "RNF144B_expression") +
	stat_compare_means(method = "t.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")), comparisons = list(c("High", "Low")), bracket.size = NA, vjust = 0.5, size = 5) +
	stat_summary(fun= "median", mapping = aes(group = RNF144B_expression),geom = "crossbar", width = c(0.95,0.8,0.65,0.7,0.77,0.9,0.9,0.9,0.9,0.93,0.66,0.35,0.96,0.8,0.97,0.65,0.55,0.93,0.8,0.9,0.92,0.81,0.71,0.35,0.64,0.73,0.67,0.89,0.84,0.54,1,0.15,0.78,0.77,0.93,0.92), color = "black") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(4, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,45)) +
	facet_wrap(~cancer_type, nrow = 3, ncol = 6, scales = "free_x")

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_all_tumors_facets_violin, filename=here(ROOT, "output", "Violin_facets_RNF144B_expression.pdf"), device = "pdf", width = 15, height = 10)

### By facets
### P53WT
plot_aneuploidy_by_RNF144b_expression_all_tumors_P53WT_facets_violin =
	df %>%
	filter(cancer_type %in% CANCERS_OI_P53_WT) %>%
	filter(TP53_mut == 0) %>%
	ggviolin(.,x = "RNF144B_expression", y = "aneuploidy_score", fill = "RNF144B_expression") +
	stat_compare_means(method = "t.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")), comparisons = list(c("High", "Low")), bracket.size = NA, vjust = 0.5, size = 5) +
	stat_summary(fun= "median", mapping = aes(group = RNF144B_expression),geom = "crossbar", width = c(0.7,0.69,0.78,0.5,0.9,0.9,0.87,0.95,0.67,0.4,0.95,0.79,0.96,0.67,0.83,0.8,0.98,0.78,0.57,0.23,0.57,0.83,0.6,0.93,0.91,0.87,1,0.15,0.86,0.78), color = "black") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(4, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,45)) +
	facet_wrap(~cancer_type, nrow = 3, ncol = 5, scales = "free_x")

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_all_tumors_P53WT_facets_violin, filename=here(ROOT, "output", "Violin_facets_RNF144B_expression_P53WT.pdf"), device = "pdf", width = 15, height = 10)

### P53mut
plot_aneuploidy_by_RNF144b_expression_all_tumors_P53mut_facets_violin =
df %>%
	filter(cancer_type %in% CANCERS_OI_P53_mut) %>%
	filter(TP53_mut == 1) %>%
	ggviolin(.,x = "RNF144B_expression", y = "aneuploidy_score", fill = "RNF144B_expression") +
	stat_compare_means(method = "t.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")), comparisons = list(c("High", "Low")), bracket.size = NA, vjust = 0.5, size = 5) +
	stat_summary(fun= "median", mapping = aes(group = RNF144B_expression),geom = "crossbar", width = c(0.965,0.955,0.795,0.93,0.75,0.82,0.8,0.96,0.93,0.62,0.74,0.98,0.77,0.91,0.8,0.915), color = "black") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(4, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,45)) +
	facet_wrap(~cancer_type, nrow = 2, ncol = 4, scales = "free_x")

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_all_tumors_P53mut_facets_violin, filename=here(ROOT, "output", "Violin_facets_RNF144B_expression_P53mut.pdf"), device = "pdf", width = 15, height = 10)

## Boxplot
### All in the same axis
stat_test_cancer_type = df %>%
	filter(cancer_type %in% CANCERS_OI) %>%
	dplyr::group_by(cancer_type) %>%
	pairwise_t_test(aneuploidy_score ~ RNF144B_expression, p.adjust.method = "BH", ref.group = "Low", alternative = "two.sided") %>%
	add_xy_position(x = "cancer_type") %>%
	mutate(y.position = 41) %>%
	mutate(xmax = xmax + 0.2)

plot_aneuploidy_by_RNF144b_expression_all_tumors_boxplot = df %>%
	filter(cancer_type %in% CANCERS_OI) %>%
	ggboxplot(.,x = "cancer_type", y = "aneuploidy_score", fill = "RNF144B_expression") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(5, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,42)) +
	stat_pvalue_manual(stat_test_cancer_type, label = "p.signif", size = 6, remove.bracket = TRUE)

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_all_tumors_boxplot, filename=here(ROOT, "output", "Boxplot_RNF144B_expression.pdf"), device = "pdf", width = 15, height = 10)

## Violin Plot
### All in the same axis
stat_test_cancer_type = df %>%
	filter(cancer_type %in% CANCERS_OI) %>%
	dplyr::group_by(cancer_type) %>%
	pairwise_t_test(aneuploidy_score ~ RNF144B_expression, p.adjust.method = "BH", ref.group = "Low", alternative = "two.sided") %>%
	add_xy_position(x = "cancer_type") %>%
	mutate(y.position = 41) %>%
	mutate(xmax = xmax + 0.2)

plot_aneuploidy_by_RNF144b_expression_all_tumors_violin = df %>%
	filter(cancer_type %in% CANCERS_OI) %>%
	ggviolin(.,x = "cancer_type", y = "aneuploidy_score", fill = "RNF144B_expression", width = 0.5) +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(5, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,42)) +
	stat_pvalue_manual(stat_test_cancer_type, label = "p.signif", size = 6, remove.bracket = TRUE)

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_all_tumors_violin, filename=figure, device = "pdf", width = 20, height = 10)

# Plot LUAD
### PANCAN
plot_aneuploidy_by_RNF144b_expression_LUAD_violin =
	df %>%
	filter(cancer_type == "LUAD") %>%
	ggviolin(.,x = "RNF144B_expression", y = "aneuploidy_score", fill = "RNF144B_expression") +
	stat_compare_means(method = "t.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")), comparisons = list(c("High", "Low")), bracket.size = NA, vjust = 0.5, size = 10) +
	stat_summary(fun= "median", mapping = aes(group = RNF144B_expression),geom = "crossbar", width = c(0.55,0.93), color = "black") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(4, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,45))

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_LUAD_violin, filename=here(ROOT, "output", "Violin_LUAD_RNF144B_expression.pdf"), device = "pdf", width = 15, height = 10)

### PANCAN_p53WT
plot_aneuploidy_by_RNF144b_expression_LUAD_P53_WT_violin =
	df %>%
	filter(cancer_type == "LUAD") %>%
	filter(TP53_mut == 0) %>%
	ggviolin(.,x = "RNF144B_expression", y = "aneuploidy_score", fill = "RNF144B_expression") +
	stat_compare_means(method = "t.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")), comparisons = list(c("High", "Low")), bracket.size = NA, vjust = 0.5, size = 10) +
	stat_summary(fun= "median", mapping = aes(group = RNF144B_expression),geom = "crossbar", width = c(0.8135,0.7935), color = "black") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(4, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,45))

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_LUAD_P53_WT_violin, filename=here(ROOT, "output", "Violin_LUAD_P53WT_RNF144B_expression.pdf"), device = "pdf", width = 15, height = 10)

stat_test_RNF144b_expression_LUAD_violin_WT = df %>%
	filter(cancer_type == "LUAD") %>%
	filter(TP53_mut == 0) %>%
	pairwise_t_test(aneuploidy_score ~ RNF144B_expression, p.adjust.method = "BH", ref.group = "Low", alternative = "two.sided")


### PANCAN_p53mut
plot_aneuploidy_by_RNF144b_expression_LUAD_P53mut_violin =
	df %>%
	filter(cancer_type == "LUAD") %>%
	filter(TP53_mut == 1) %>%
	ggviolin(.,x = "RNF144B_expression", y = "aneuploidy_score", fill = "RNF144B_expression") +
	stat_compare_means(method = "t.test", label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001,0.001, 0.01, 0.05, 1), symbols = c("****","***", "**", "*", "ns")), comparisons = list(c("High", "Low")), bracket.size = NA, vjust = 0.5, size = 10) +
	stat_summary(fun= "median", mapping = aes(group = RNF144B_expression),geom = "crossbar", width = c(0.76,0.99), color = "black") +
	theme_pubr() +
	theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 10),
				axis.title.x = element_blank(),
				axis.title.y = element_text(size = 15),
				axis.text.y = element_text(size = 12),
				legend.position = "right",
				legend.key.width = unit(0.5, "cm"),
				legend.title = element_text(size = 14, face = "bold"),
				legend.text = element_text(size = 12),
				panel.spacing.y = unit(4, "mm")) + 
	scale_fill_manual(values = c("#E2E2FA","#9E9EF2"), name = "RNF144B_expression") +
	ylab("Aneuploidy scores") +
	ylim(c(0,45))

# save figure
ggsave(plot_aneuploidy_by_RNF144b_expression_LUAD_P53mut_violin, filename=here(ROOT, "output", "Violin_LUAD_P53mut_RNF144B_expression.pdf"), device = "pdf", width = 15, height = 10)

stat_test_RNF144b_expression_LUAD_violin_mut = df %>%
	filter(cancer_type == "LUAD") %>%
	filter(TP53_mut == 1) %>%
	pairwise_t_test(aneuploidy_score ~ RNF144B_expression, p.adjust.method = "BH", ref.group = "Low", alternative = "two.sided")

