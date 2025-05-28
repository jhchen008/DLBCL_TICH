
###### load required packages ######
library(survival)
library(survminer)
library(paletteer)
library(GSVA)
library(GSEABase)
library(AUCell)

library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)
library(tidyr)

library(viridis)
library(ggpointdensity)
library(RColorBrewer)
library(corrplot)

##############################################################################################################################
# DLBCL cohort in Lenz et al. 2008 NEJM
##############################################################################################################################

###### expression matrix ######

# reading the expression matrix from GEO, with acc# GSE10846
expression_data = read.table("GSE10846_series_matrix.nor.txt",header=T,sep = "\t")
rownames(expression_data) <- expression_data$ID_REF
expression_data <- expression_data[,-1]

# patient information in the series.matrix file
group_data <- read.table("GSE10846_series_matrix.patient.group.txt",header=T,fill = TRUE,sep = "\t")



###### annotation wit gene SYMBOL ######

library("AnnotationDbi")
library("hgu133plus2.db")
library("org.Hs.eg.db")                                          
columns(org.Hs.eg.db)
library(tibble)


id_str <- rownames(expression_data)
annotable = AnnotationDbi::select(hgu133plus2.db,id_str,c("SYMBOL"), keytype="PROBEID")

# add gene symbol to the matrix
expression_data_df <- expression_data %>%  tibble::rownames_to_column(var = "PROBEID")
expression_data_df <- expression_data_df %>% left_join(annotable, by = "PROBEID")

# remove rows without annotation results
df_filtered <- expression_data_df %>% filter(!is.na(SYMBOL))

# collapse probe sets to gene symbol with mean values
df_final <- df_filtered %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), 
   ~ if(n() > 1) mean(.x, na.rm = TRUE) else .x[1])) %>% ungroup()

df_final <- df_final %>% column_to_rownames(var = "SYMBOL")



###### gene set analysis ######

gset_ch <- getGmt("gene_sets.gmt")
gset_ch <- subsetGeneSets(gset_ch, rownames(df_final)) 

ssgseaPar <- ssgseaParam(as.matrix(df_final), gset_ch, )
gset_ch_ssgsea <- gsva(ssgseaPar, verbose=FALSE)
gset_ch_ssgsea = as.data.frame(gset_ch_ssgsea)

# merge with table containing patient information
transposed_data <- as.data.frame(t(gset_ch_ssgsea))
transposed_data <- data.frame(ACC = rownames(transposed_data), transposed_data)
merged_data <- merge(group_data, transposed_data, by.x = "GSE_acc", by.y = "ACC")



###### survival analysis ######

merged_data = subset(merged_data, OSS != "NA")

merged_data$DLBCL_CH_group <- ifelse(merged_data$DLBCL_CH_Up_CAPS >= quantile(merged_data$DLBCL_CH_Up_CAPS,0.5), "High", "Low")
merged_data$DLBCL_CH_group <- factor(merged_data$DLBCL_CH_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data$FU_time, event = merged_data$OSS)
fit <- survfit(surv_object ~ DLBCL_CH_group, data = merged_data)

plot.new()
p = ggsurvplot(fit, data = merged_data, pval = TRUE, risk.table = TRUE,
           title = "OS by CAPS score",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "CAPS score",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)



###### multivariate COX regression ######

merged_data$Gender_gp <- ifelse(merged_data$Gender == "","NA", ifelse(merged_data$Gender == "female","1", "0"))
merged_data$HANS <- ifelse(merged_data$HANS == "","NA", ifelse(merged_data$HANS != "ABC","1", "0"))

cox_model <- coxph(Surv(FU_time, OSS) ~ DLBCL_CH_group +Gender_gp +nonGCB +IPI_gp, data = merged_data)
cox_summary = summary(cox_model)
cox_summary

cox_model <- coxph(Surv(FU_time, OSS) ~ DLBCL_CH_group +Gender_gp +HANS +IPI +Age, data = merged_data)
cox_summary = summary(cox_model)
cox_summary



###### ecotyper results ######

# ecotyper_types
ecotyper_types = read.table("./ecotyper_output/Lymphoma_Ecotypes/Ecotype_Assignment.txt",header=T,sep="\t")
coldata = merge(merged_data, ecotyper_types, by.x = "GSE_acc", by.y = "ID")

chisq.test(table(coldata$DLBCL_CH_group, coldata$Lymphoma.Ecotype))


# monocytes and macrophage
ecotyper_mac = read.table("./ecotyper_output/Lymphoma_Cell_States/Monocytes.and.Macrophages/Monocytes.and.Macrophages_Cell_State_Assignment.txt",header=T,sep="\t")
coldata = merge(merged_data, ecotyper_mac, by.x = "GSE_acc", by.y = "ID")

chisq.test(table(coldata$DLBCL_CH_group, coldata$Cell.State))

#### Fibroblasts
ecotyper_Fib = read.table("./ecotyper_output/Lymphoma_Cell_States/Fibroblasts/Fibroblasts_Cell_State_Assignment.txt",header=T,sep="\t")
coldata = merge(merged_data, ecotyper_Fib, by.x = "sampleID", by.y = "ID")

chisq.test(table(coldata$DLBCL_CH_group, coldata$Cell.State))




##############################################################################################################################
# DLBCL cohort in Schmitz et al. 2018 NEJM
##############################################################################################################################

###### expression matrix ######

# expression matrix and patient data were download from: 
# https://gdc.cancer.gov/about-data/publications/DLBCL-2018

expression_data = read.table("RNAseq_gene_expression_562.txt",header=T,sep = "\t")
rownames(expression_data) <- expression_data$Gene
expression_data <- expression_data[,c(-1,-2,-3)]

group_data <- read.table("Supplementary_Appendix_2.TableS9.txt",header=T,fill = TRUE,sep = "\t")
group_data$sampleID = gsub("-", ".", group_data$dbGaP.submitted.subject.ID)



###### gene set analysis ######

gset_ch <- getGmt("gene_sets.gmt")
gset_ch <- subsetGeneSets(gset_ch, rownames(expression_data)) 

ssgseaPar <- ssgseaParam(as.matrix(expression_data), gset_ch, )
gset_ch_ssgsea <- gsva(ssgseaPar, verbose=FALSE)
gset_ch_ssgsea = as.data.frame(gset_ch_ssgsea)

# merge with table containing patient information
transposed_data <- as.data.frame(t(gset_ch_ssgsea))
transposed_data <- data.frame(ACC = rownames(transposed_data), transposed_data)
merged_data <- merge(group_data, transposed_data, by.x = "sampleID", by.y = "ACC")



###### survival analysis ######

# subset patients with OS/PFS data
merged_data_ch = subset(merged_data_all, Included.in.Survival.Analysis == "Yes")

merged_data_ch$DLBCL_CH_group <- ifelse(merged_data_ch$DLBCL_CH_Up_CAPS >= quantile(merged_data_ch$DLBCL_CH_Up_CAPS,0.5), "High", "Low")
merged_data_ch$DLBCL_CH_group <- factor(merged_data_ch$DLBCL_CH_group, levels = c("Low", "High"))


## PFS
surv_object <- Surv(time = merged_data_ch$Progression_Free.Survival._PFS_.Time._yrs, event = merged_data_ch$Progression_Free.Survival._PFS_.Status_.0.No.Progressoin_.1.Progression)
fit <- survfit(surv_object ~ DLBCL_CH_group, data = merged_data_ch)

plot.new()
p = ggsurvplot(fit, data = merged_data_ch, pval = TRUE, risk.table = TRUE,
           title = "PFS by CAPS score",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "CAPS score",
           legend.labs = c("Low", "High"),
           legend = c(0.8, 0.2), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)

## OS
surv_object <- Surv(time = merged_data_ch$Follow_up.Time._yrs, event = merged_data_ch$Status.at.Follow_up_.0.Alive_.1.Dead)
fit <- survfit(surv_object ~ DLBCL_CH_group, data = merged_data_ch)

plot.new()
p = ggsurvplot(fit, data = merged_data_ch, pval = TRUE, risk.table = TRUE,
           title = "OS by CAPS score",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "CAPS score",
           legend.labs = c("Low", "High"),
           legend = c(0.8, 0.2), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)



###### multivariate COX regression ######

merged_data_ch$nonGCB <- ifelse(merged_data_ch$Gene.Expression.Subgroup == "","NA", ifelse(merged_data_ch$Gene.Expression.Subgroup == "ABC","1", "0"))
merged_data_ch$Gender_gp <- ifelse(merged_data_ch$Gender == "","NA", ifelse(merged_data_ch$Gender == "F","1", "0"))
merged_data_ch$IPI_3to5 <- ifelse(merged_data_ch$IPI_Score == "","NA", ifelse(merged_data_ch$IPI_Score >2,"1", "0"))

## PFS
cox_model <- coxph(Surv(Progression_Free.Survival._PFS_.Time._yrs, Progression_Free.Survival._PFS_.Status_.0.No.Progressoin_.1.Progression) ~ DLBCL_CH_group +Gender_gp +IPI_3to5 +nonGCB +Age, data = pt_fihser)
cox_summary = summary(cox_model)
cox_summary

## OS
cox_model <- coxph(Surv(Follow_up.Time._yrs, Status.at.Follow_up_.0.Alive_.1.Dead) ~ DLBCL_CH_group +Gender_gp +IPI_3to5 +nonGCB +Age, data = merged_data_ch)
cox_summary = summary(cox_model)
cox_summary



###### ecotyper results ######

# ecotyper_types
ecotyper_types = read.table("./ecotyper_output/Lymphoma_Ecotypes/Ecotype_Assignment.txt",header=T,sep="\t")
coldata = merge(merged_data_ch, ecotyper_types, by.x = "sampleID", by.y = "ID")

chisq.test(table(coldata$DLBCL_CH_group, coldata$Lymphoma.Ecotype))

# monocytes and macrophage
ecotyper_mac = read.table("./ecotyper_output/Lymphoma_Cell_States/Monocytes.and.Macrophages/Monocytes.and.Macrophages_Cell_State_Assignment.txt",header=T,sep="\t")
coldata = merge(merged_data_ch, ecotyper_mac, by.x = "sampleID", by.y = "ID")

chisq.test(table(coldata$DLBCL_CH_group, coldata$Cell.State))

#### Fibroblasts
ecotyper_Fib = read.table("./ecotyper_output/Lymphoma_Cell_States/Fibroblasts/Fibroblasts_Cell_State_Assignment.txt",header=T,sep="\t")
coldata = merge(merged_data_ch, ecotyper_Fib, by.x = "sampleID", by.y = "ID")

chisq.test(table(coldata$DLBCL_CH_group, coldata$Cell.State))





##############################################################################################################################
# DLBCL cohort in Lacy et al. 2020 Blood
##############################################################################################################################


###### expression matrix ######

# reading the expression matrix from GEO, with acc# GSE181063
expression_data = read.table("GSE181063_series_matrix.nor.txt",header=T,sep = "\t")
rownames(expression_data) <- expression_data$ID_REF
expression_data <- expression_data[,-1]

# patient information in the series.matrix file
group_data <- read.table("GSE181063_series_matrix.sample.info.txt",header=T,fill = TRUE,sep = "\t")

# PFS data from subset of samples were available in the supplement table S4
group_data_pfs <- read.table("Blood_PMID32187361_TableS4.txt",header=T,fill = TRUE,sep = "\t")


###### annotation wit gene SYMBOL ######

library("AnnotationDbi")
library("illuminaHumanv4.db")
library("org.Hs.eg.db")
columns(illuminaHumanv4.db)
library(tibble)

id_str <- rownames(expression_data)

annotable = AnnotationDbi::select(illuminaHumanv4.db,id_str,c("SYMBOL"), keytype="PROBEID")

# add gene symbol to the matrix
expression_data_df <- expression_data %>%  tibble::rownames_to_column(var = "PROBEID")
expression_data_df <- expression_data_df %>% left_join(annotable, by = "PROBEID")

df_final <- expression_data_df %>% filter(!is.na(SYMBOL))

# collapse probe sets to gene symbol with mean values
df_final <- df_final %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), 
   ~ if(n() > 1) mean(.x, na.rm = TRUE) else .x[1])) %>% ungroup()

df_final <- df_final %>% column_to_rownames(var = "SYMBOL")



###### gene set analysis ######

gset_ch <- getGmt("gene_sets.gmt")
gset_ch <- subsetGeneSets(gset_ch, rownames(expression_data)) 

ssgseaPar <- ssgseaParam(as.matrix(expression_data), gset_ch, )
gset_ch_ssgsea <- gsva(ssgseaPar, verbose=FALSE)
gset_ch_ssgsea = as.data.frame(gset_ch_ssgsea)

# merge with table containing patient information
transposed_data <- as.data.frame(t(gset_ch_ssgsea))
transposed_data <- data.frame(ACC = rownames(transposed_data), transposed_data)

# OS and PFS tables
merged_data_all <- merge(group_data, transposed_data, by.x = "sampleID", by.y = "ACC")
merged_data_pfs_all <- merge(group_data_pfs, transposed_data, by.x = "sampleID", by.y = "ACC")

# only CHOP/R-CHOP based
merged_data_ch = subset(merged_data_all, disease == "DLBCL" & FL %in% c("CHOP-R","CHOP"))
merged_data_ch_pfs = subset(merged_data_pfs_all, disease == "DLBCL" & FL %in% c("CHOP-R","CHOP"))



###### survival analysis ######

# OS
merged_data_ch$DLBCL_CH_group <- ifelse(merged_data_ch$DLBCL_CH_Up_CAPS >= quantile(merged_data_ch$DLBCL_CH_Up_CAPS,0.5), "High", "Low")
merged_data_ch$DLBCL_CH_group <- factor(merged_data_ch$DLBCL_CH_group, levels = c("Low", "High"))

surv_object <- Surv(time = merged_data_ch$OS_time, event = merged_data_ch$OSS)
fit <- survfit(surv_object ~ DLBCL_CH_group, data = merged_data_ch)

plot.new()
p = ggsurvplot(fit, data = merged_data_ch, pval = TRUE, risk.table = TRUE,
           title = "OS by CAPS score",
           xlab = "Time (Years)", ylab = "OS (%)",
           legend.title = "CAPS score",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)


# PFS
merged_data_ch_pfs = subset(merged_data_ch_pfs, !is.na(merged_data_ch_pfs$PFS_status))
merged_data_ch_pfs$DLBCL_CH_group <- ifelse(merged_data_ch_pfs$DLBCL_CH_Up_CAPS  >= quantile(merged_data_ch_pfs$DLBCL_CH_Up_CAPS,0.5), "High", "Low")
merged_data_ch_pfs$DLBCL_CH_group <- factor(merged_data_ch_pfs$DLBCL_CH_group, levels = c("Low", "High"))

surv_object <- Surv(time = (merged_data_ch_pfs$PFS_time)/365, event = merged_data_ch_pfs$PFS_status)
fit <- survfit(surv_object ~ DLBCL_CH_group, data = merged_data_ch_pfs)

plot.new()
p = ggsurvplot(fit, data = merged_data_ch_pfs, pval = TRUE, risk.table = TRUE,
           title = "PFS by CAPS score",
           xlab = "Time (Years)", ylab = "PFS (%)",
           legend.title = "CAPS score",
           legend.labs = c("Low", "High"),
           legend = c(0.75, 0.8), font.legend=12,
           risk.table.height = 0.28,
           palette = c("#2E9FDF", "#FF2400"))
print(p, newpage = FALSE)



###### multivariate COX regression ######

# OS
merged_data_ch$nonGCB <- ifelse(merged_data_ch$Pred_combine == "","NA", ifelse(merged_data_ch$Pred_combine == "ABC","1", "0"))
merged_data_ch$Gender_gp <- ifelse(merged_data_ch$Sex == "","NA", ifelse(merged_data_ch$Sex == "F","1", "0"))
merged_data_ch$IPI_3to5 <- ifelse(merged_data_ch$Ipi_score == "","NA", ifelse(merged_data_ch$Ipi_score >2,"1", "0"))

cox_model <- coxph(Surv(Os_followup_y, Os_status) ~ DLBCL_CH_group +Gender_gp +IPI_3to5 +nonGCB+Age_at_diagnosis, data = merged_data_ch)
cox_summary = summary(cox_model)
cox_summary

# PFS
merged_data_ch_pfs$nonGCB <- ifelse(merged_data_ch_pfs$cell_of_origin == "","NA", ifelse(merged_data_ch_pfs$cell_of_origin == "ABC","1", "0"))
merged_data_ch_pfs$Gender_gp <- ifelse(merged_data_ch_pfs$Sex == "","NA", ifelse(merged_data_ch_pfs$Sex == "F","1", "0"))
merged_data_ch_pfs$IPI_3to5 <- ifelse(merged_data_ch_pfs$Ipi_score == "","NA", ifelse(merged_data_ch_pfs$Ipi_score >2,"1", "0"))

cox_model <- coxph(Surv(PFS_time, PFS_status) ~ DLBCL_CH_group +Gender_gp +IPI_3to5 +nonGCB+Age_at_diagnosis, data = merged_data_ch_pfs)
cox_summary = summary(cox_model)
cox_summary



###### ecotyper results ######

# ecotyper_types
ecotyper_types = read.table("./ecotyper_output/Lymphoma_Ecotypes/Ecotype_Assignment.txt",header=T,sep="\t")
coldata = merge(merged_data_ch, ecotyper_types, by.x = "GSE_acc", by.y = "ID")

chisq.test(table(coldata$DLBCL_CH_group, coldata$Lymphoma.Ecotype))


# monocytes and macrophage
ecotyper_mac = read.table("./ecotyper_output/Lymphoma_Cell_States/Monocytes.and.Macrophages/Monocytes.and.Macrophages_Cell_State_Assignment.txt",header=T,sep="\t")
coldata = merge(merged_data_ch, ecotyper_mac, by.x = "GSE_acc", by.y = "ID")

chisq.test(table(coldata$DLBCL_CH_group, coldata$Cell.State))


# Fibroblasts
ecotyper_Fib = read.table("./ecotyper_output/Lymphoma_Cell_States/Fibroblasts/Fibroblasts_Cell_State_Assignment.txt",header=T,sep="\t")
coldata = merge(merged_data_ch, ecotyper_Fib, by.x = "GSE_acc", by.y = "ID")

chisq.test(table(coldata$DLBCL_CH_group, coldata$Cell.State))