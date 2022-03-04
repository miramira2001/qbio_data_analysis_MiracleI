if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
install.packages(c("survival", "survminer"))
library(TCGAbiolinks)
library(SummarizedExperiment)
library("DESeq2")
library(survival)
library(survminer)
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)
rowData(sum_exp)$external_gene_name
# Wanted to see what genes I could choose from 
"MYC" %in% rowData(sum_exp)$external_gene_name
# After some research online, the gene MYC looked interesting to me and I decdided that this is the 
# gene I want to explore. This line is to ensure that it is present in this dataset. 
colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_index < 50, "Young","Old")
# This is a column in colData(sum_exp) that categorizes patients as young and old 

bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)
num_na

age_cat_no_NAs = colData(sum_exp)$age_category[!bool_age_na]
sum(is.na(age_cat_no_NAs))
length(age_cat_no_NAs)

geneA_id_mask = (rowData(sum_exp)$external_gene_name == "MYC")
sum(geneA_id_mask)
ensembl_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_id_mask]


gene_counts = assays(sum_exp)$"HTSeq - Counts"[geneA_id_mask, !bool_age_na]

barplot(table(age_cat_no_NAs), col = rainbow(2), main = "Number of Old and Young")
boxplot(gene_counts ~ age_cat_no_NAs, xlab = "Age category", ylab = "Counts of MYC")

sum_exp$days_to_death = ifelse(is.na(sum_exp$days_to_death), sum_exp$days_to_last_follow_up, sum_exp$days_to_death)
sum_exp$death_event = as.integer(sum_exp$vital_status == "Dead")

surv_object = Surv(time = sum_exp$days_to_death, event = sum_exp$death_event)

gender_fit = survfit(surv_object ~ sum_exp$gender, data = sum_exp)
survfit(surv_object ~ sum_exp$gender, data = sum_exp)

survplot = ggsurvplot(gender_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")
p = survplot$plot + theme_bw()
p + theme(legend.position = "bottom", legend.direction = "vertical")
