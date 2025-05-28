
###### load required packages ######
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)
library(tidyr)
library(paletteer)
library(reshape2)



###############################################################

###### differential expression analysis with limma-voom ######
library(limma)
library(edgeR)

# expression matrix and group info

expression_data <- read.table("DLBCL_CH_RNAseq_counts.CH.txt",header=T)
row.names(expression_data) = expression_data$Gene
expression_data <- expression_data[,-1]

group_data = read.table("DLBCL_CH_RNAseq_counts.CH.design.txt",header=T)

group = factor(as.character(group_data[,colnames(expression_data)]))

# limma

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

DGE = DGEList(expression_data)
keep <- filterByExpr(DGE, group=group)

DGE <- DGE[keep, , keep.lib.sizes=FALSE]
DGE = calcNormFactors(DGE,method =c("TMM"))

v = voom(DGE, design, plot = TRUE) 
fit = lmFit(v, design)

# contrast analysis to detect DEGs

contrast.matrix <- makeContrasts(pos_vs_neg = positive - negative, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# DEG table
results <- topTable(fit2, adjust = "fdr", number = Inf)
results$logP <- -log10(results$P.Value)

# normalized expression matrix
normalized_expression = v$E


# volcano plot

results$Significance <- "Not Significant"
results$Significance[results$logFC > log2(1.5) & results$P.Value < 0.05] <- "Upregulated"     # fold change >1.5, P < 0.05
results$Significance[results$logFC < -log2(1.5) & results$P.Value < 0.05] <- "Downregulated"

# select genes to display
target_genes <- c("IL6","CCL8", "C4BPA", "CD163","S100A12","SELE","SELP","C4B_2","C5AR1","ACKR1","SCG2","IL2","LILRA5","LILRB5","NAIP")

ggplot(results, aes(x = logFC, y = logP)) +
  geom_point(aes(color = Significance), alpha = 0.6) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_pubr(base_size=20) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P.value") + xlim(-2,2) +
  theme(axis.text=element_text(colour = "black"), axis.ticks = element_line(colour = "black",linewidth=0.8), legend.position="none") +
  geom_vline(xintercept=c(-log2(1.5),log2(1.5)), linetype="dashed", color = "black") + geom_hline(yintercept=1.3, linetype="dashed", color = "black") +
  geom_text_repel(data = subset(results, rownames(results) %in% target_genes), aes(label = Gene), color = "black", box.padding = 0.65, segment.color = "grey50", max.overlaps = Inf)



###############################################################

###### lasso regression to identify CAPS genes ######


# up-regulated DEGs and their expression data
rownames(results) = results$Gene

DEG_up = rownames(subset(results, logFC>log2(1.5) & P.Value<0.05))

gene_data_up <- normalized_expression[rownames(normalized_expression) %in% DEG_up, ]
gene_data_up = t(gene_data_up)
patient_ids <- rownames(gene_data_up)


# survival data
matched_survival_data <- all_patient[all_patient$RNA_seq_ID %in% patient_ids, c("RNA_seq_ID","OS_time","is_OS")]
matched_survival_data <- matched_survival_data[match(patient_ids, matched_survival_data$RNA_seq_ID), ]

# make sure there is no missing survival data
complete_cases = complete.cases(matched_survival_data$OS_time, matched_survival_data$is_OS)
matched_survival_data <- matched_survival_data[complete_cases,]
gene_data_up <- gene_data_up[complete_cases, ]


## lasso analysis
library(glmnet)

sur <- Surv(matched_survival_data$OS_time, matched_survival_data$is_OS)
exp <- as.matrix(gene_data_up)

cv_lasso_cox <- cv.glmnet(exp, sur, family = "cox", alpha = 1, nfolds = 10)
min_lambda <- cv_lasso_cox$lambda.min
print(min_lambda)

# fitting with lambda
final_model <- glmnet(exp, sur, family = "cox", alpha = 1, lambda = min_lambda)
final_model <- glmnet(exp, sur, family = "cox", alpha = 1, lambda = 0.04)      ## lamda 0.04 for downstream

# genes correlate with OS
coef_matrix <- as.matrix(coef(final_model))
coef_matrix[coef_matrix>0,,drop=F]



###############################################################

###### GSEA with gene sets ######

library(GSVA)
library(GSEABase)
library(AUCell)

# load gene sets
gset_ch <- getGmt("gene_sets.gmt")
gset_ch <- subsetGeneSets(gset_ch, rownames(normalized_expression)) 

# gene set enrichment score with ssGSEA
ssgseaPar <- ssgseaParam(as.matrix(normalized_expression), gset_ch)
gset_ch_ssgsea <- gsva(ssgseaPar, verbose=FALSE)
gset_ch_ssgsea = as.data.frame(gset_ch_ssgsea)

# merge GSEA results with patient information
transposed_data <- as.data.frame(t(gset_ch_ssgsea))
colnames(transposed_data) <- rownames(gset_ch_ssgsea)
rownames(transposed_data) <- colnames(gset_ch_ssgsea)

# add gene expression into the patient info table
transposed_data <- data.frame(Seq_ID = rownames(transposed_data), transposed_data)
merged_data <- merge(all_patient, transposed_data, by.x = "RNA_seq_ID", by.y = "Seq_ID")



###############################################################
library(survival)
library(survminer)

###### KM survival analysis ######

# define groups with cut-off by CAPS median
merged_data$DLBCL_CAPS_group <- ifelse(merged_data$DLBCL_CH_Up_CAPS >= quantile(merged_data$DLBCL_CH_Up_CAPS,0.5), "High", "Low")
merged_data$DLBCL_CAPS_group = factor(merged_data$DLBCL_CAPS_group, levels=c("Low", "High"))

# OS
surv_object <- Surv(time = (merged_data$OS_time)/12, event = merged_data$is_OS)
fit <- survfit(surv_object ~ DLBCL_CAPS_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "OS by CH Signature",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "CH-lasso",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)



# PFS
surv_object <- Surv(time = (merged_data$PFS_time)/12, event = merged_data$is_PFS)
fit <- survfit(surv_object ~ DLBCL_CAPS_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "PFS by CH Signature",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "CH-lasso",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


###### multivariate COX regression analysis ######

# define groups
pt_fisher = merged_data

pt_fisher$Gender_gp <- ifelse(pt_fisher$Gender == "","NA", ifelse(pt_fisher$Gender == "Male", "1", "0"))
pt_fisher$COO_NonvsGCB <- ifelse(pt_fisher$HANS == "","NA", ifelse(pt_fisher$HANS == "nonGCB", "1", "0"))
pt_fisher$IPI_3to5 <- ifelse(pt_fisher$IPI == "","NA", ifelse(pt_fisher$IPI > 2, "1", "0"))

# multivariate cox regression analysis
cox_model <- coxph(Surv(OS_time, is_OS) ~ DLBCL_CAPS_group +Gender_gp +IPI_3to5 +COO_NonvsGCB +Age, data = pt_fisher)
cox_summary = summary(cox_model)
cox_summary

cox_model <- coxph(Surv(PFS_time, is_PFS) ~ DLBCL_CAPS_group +Gender_gp +IPI_3to5 +COO_NonvsGCB +Age, data = pt_fisher)
cox_summary = summary(cox_model)
cox_summary



###############################################################

###### ecotyper results ######

# ecotyper output files, Ecotypes
ecotyper_types_abd = read.table("./ecotyper_output/Lymphoma_Ecotypes/Ecotype_Abundance.txt",header=T,sep="\t")
ecotyper_types = read.table("./ecotyper_output/Lymphoma_Ecotypes/Ecotype_Assignment.txt",header=T,sep="\t")

coldata = merge(merged_data, ecotyper_types_abd, by.x = "RNA_seq_ID", by.y = "ID")
coldata = merge(merged_data, ecotyper_types, by.x = "RNA_seq_ID", by.y = "ID")

chisq.test(table(coldata$DLBCL_CAPS_group, coldata$Lymphoma.Ecotype))


# ecotyper output files, macrophage
ecotyper_mac_abd = read.table("./ecotyper_output/Lymphoma_Cell_States/Monocytes.and.Macrophages/Monocytes.and.Macrophages_Cell_State_Abundance.txt",header=T,sep="\t")
ecotyper_mac = read.table("./ecotyper_output/Lymphoma_Cell_States/Monocytes.and.Macrophages/Monocytes.and.Macrophages_Cell_State_Assignment.txt",header=T,sep="\t")

coldata = merge(merged_data, ecotyper_mac_abd, by.x = "RNA_seq_ID", by.y = "ID")
coldata = merge(merged_data, ecotyper_mac, by.x = "RNA_seq_ID", by.y = "ID")

chisq.test(table(coldata$DLBCL_CAPS_group, coldata$Cell.State))


# ecotyper output files, Fibroblasts
ecotyper_Fib = read.table("./ecotyper_output/Lymphoma_Cell_States/Fibroblasts/Fibroblasts_Cell_State_Assignment.txt",header=T,sep="\t")
ecotyper_Fib_abd = read.table("./ecotyper_output/Lymphoma_Cell_States/Fibroblasts/Fibroblasts_Cell_State_Abundance.txt",header=T,sep="\t")

coldata = merge(merged_data, ecotyper_Fib_abd, by.x = "RNA_seq_ID", by.y = "ID")
coldata = merge(coldata, ecotyper_Fib, by.x = "RNA_seq_ID", by.y = "ID")

chisq.test(table(coldata$DLBCL_CAPS_group, coldata$Cell.State))
