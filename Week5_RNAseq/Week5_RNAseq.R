if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
library(TCGAbiolinks)
library(SummarizedExperiment)
library("DESeq2")
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
GDCdownload(query) # only need to download the data once! Comment this out once you have completed it once
sum_exp <- GDCprepare(query)
assays(sum_exp)$"HTSeq - Counts"[1:5, 1:5]
dim(colData(sum_exp))
dim(rowData(sum_exp))
str(colData(sum_exp))
head(colData(sum_exp))
colnames(colData(sum_exp))
colData(sum_exp)[1:5, 25:29]
"days"["age_at_diagnosis[1:10]"]
colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis 
/365 


colData(sum_exp)$age_category == (colData(sum_exp)$age_category) %%  365
# Tried to use mod to find remainder which would give values in years

colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_diagnosis<50, "Young", "Old")

head(rowData(sum_exp))
dim(rowData(sum_exp))
# Summary of the experiment 

"MSH2" %in%rowData(sum_exp)$external_gene_name 
"MSH6" %in%rowData(sum_exp)$external_gene_name 

assays(sum_exp)$"HTSeq - Counts"[20:25, 30:35]
# Rows are clinical variables 

geneA_id_mask = rowData(sum_exp)$external_gene_name == "MSH2"
sum(geneA_id_mask)
ensembl_geneA = rowData(sum_exp)$external_gene_name["HTSeq - Counts"]

geneB_id_mask = rowData(sum_exp)$external_gene_name == "MSH6"
sum(geneB_id_mask)
ensembl_geneB = rowData(sum_exp)$external_gene_name["HTSeq - Counts"]
# Ensembl gene ID should be a column

min(assays(sum_exp)$"HTSeq - Counts"[ ,ensembl_geneA]  )
# Goes on right side since it is a column 

max(assays(sum_exp)$"HTSeq - Counts"[ ,ensembl_geneA]  )

summary(assays(sum_exp)$"HTSeq - Counts")

plot(geneA_id_mask, 
     geneB_id_mask,
     xlab = "Gene MSH2", 
     ylab = "Gene MSH6"
)
# Genes do not seem correlated 

bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)
num_na
# Age category colum is in colData dataframe

age_cat_no_NAs = !(is.na(colData(sum_exp)$age_category))
length(age_cat_no_NAs)  
num_na + age_cat_no_NAs == dim( colData(sum_exp) )[1]

dim( colData(sum_exp) )[1]
# Number of patients is 521

identical( rownames(colData(sum_exp)), colnames(assays(sum_exp)$"HTSeq - Counts")  )

length(age_cat_no_NAs) == length(assays(sum_exp)$"HTSeq - Counts")
# Returns false value so they are not the same length but if they were, the boxplot would be made as shown below

boxplot(assays(sum_exp)$"HTSeq - Counts" ~ age_cat_no_NAs,
        xlab = "Data without NAs", 
        ylab = "Gene Counts")

# To access HTSeq - Counts - assays(sum_exp)$"HTSeq - Counts"
# Rows are clinical variables adn columns are category names
# Use boolean indexing to access certain information in counts dataframe [rows,] for the rows or [,col] for column

