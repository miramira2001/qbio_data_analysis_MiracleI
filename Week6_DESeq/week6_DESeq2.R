if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
BiocManager::install("DESeq2")
library(TCGAbiolinks)
library(SummarizedExperiment)
library("DESeq2")
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)
head(colData(sum_exp))
is.na(colData(sum_exp))
# colData has the information about patients with NA values for their ages
patient_data = colData(sum_exp)
counts = assays(sum_exp)$"HTSeq - Counts"
# Instructed to skip 1.2 part 3
patient_data_na_mask = is.na(colData(sum_exp)$age_at_index)
counts = counts[ , !(patient_data_na_mask) ]
patient_data = patient_data[ !(patient_data_na_mask) , ]
patient_data$age_category = ifelse(patient_data$age_at_index < 50, "young", "old") 
patient_data$age_category = factor(patient_data$age_category, levels = c("young", "old"))


if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}
counts_row_sums = rowSums(counts)
low_counts_mask = ifelse(counts_row_sums < 10, FALSE, TRUE)
sum(low_counts_mask) 
counts = counts[ low_counts_mask, ]
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = patient_data, 
                             design = ~age_category)

dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

head(results)
str(results)

my_df = data.frame(x = c('b', 'd', 'c', 'e', 'a'),
                   y = c(2,4,3,5,1))

order_indices = order(my_df$y)
# we expect c(5, 1, 3, 2, 4) because:
# 1 is index 5
# 2 is index 1
# 3 is index 3
# 4 is index 2
# 5 is index 4
order_indices  # note the order!

my_df = my_df[order_indices, ]
my_df

row_order = order(results)
head(results, 20)

# Gene chosen is CFTR and this stands for CF Transmembrane Conductance Regulator and this gene 
# encodes a member of the ATP-binding cassette (ABC) transporter superfamily and the encoded
# proteins then act as chloride channels. Gene is more highly expressed among old patients
# (above age of 50) since the log2FoldChange value is negative.

log2FoldChange_threshold =  ifelse(results$log2FoldChange_threshold < 1, "young", "old") 
padj_threshold = ifelse(results$padj_threshold < 0.05, "young", "old") 
results$log2FoldChange_threshold > results$log2FoldChange 
# Tried indexing to see only the values where log2FoldChange_threshold is greater than 
# log2FoldChange
results$padj_threshold > results$padj 
results = results [ results$padj_threshold > results$padj  , results$log2FoldChange_threshold > results$log2FoldChange ]



fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

plot(x = log2FoldChange_threshold,
     y = -log10(padj_threshold),
     xlab = "-log10(p value)", # be sure the specify that it's young over old!
     ylab = "Log2 Fold Change",
     pch = 20) # smaller solid circles

abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")


fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05


library(ggplot2)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange_threshold, y = -log10(padj_threshold))) + 
  geom_point(aes(color = ifelse(log2FoldChange_threshold < -1 & padj_threshold < 0.05, "lower in young",
                                ifelse(log2FoldChange_threshold > 1 & padj_threshold < 0.05, "higher in young", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Young/Old)",
       y = "-log10 Adjusted p-value")


volcano_plot

