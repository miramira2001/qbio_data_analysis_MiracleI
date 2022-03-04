BiocManager::install("maftools")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(maftools)

## Question 1
## options(warn=-1)

setwd("/Users/miracle/qbio_data_analysis_MiracleI")
clinic <- data.table::fread("/Users/miracle/qbio_data_analysis_MiracleI/week4_clinical/Week_4_HW.R", data.table = F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] = "Tumor_Sample_Barcode"
length(colnames(clinic))
## colnames(clinic) has a length of 2. 

length(clinic[ "Tumor_Sample_Barcode" , ])

## It has a length of 524 
sum(colnames(clinic) == "bcr_patient_barcode")
## There are 0 trues 
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)
maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)
getwd()
setwd("/Users/miracle/qbio_data_analysis_MiracleI/week7_MAF/GDCdata/")
list.files()

maf_dataframe = data.table::fread("TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv", data.table=F)
clinic = data.table::fread("/Users/miracle/qbio_data_analysis_MiracleI/week4_clinical/Week_4_HW.R", data.table = F)
colnames(clinic)[colnames(clinic) == "bcr_patient_barcode"] = "Tumor_Sample_Barcode"
maf_object = read.maf(maf = maf_dataframe, clinicalData = clinic, isTCGA = T)

## Question 2
maf_object
str(maf_object)
maf_object@data
str(maf_object@data)
maf_object@clinical.data
str(maf_object@clinical.data)

colnames(maf_object@data)
colnames(maf_object@clinical.data)
## Tumor_Sample_Barcode is shared between the 2 data frames and this makes sense because the mutation information must relate 
## to the patients accordingly

## Question 3
oncoplot(maf = maf_object,
         top = 19) 

ggsave("/Users/miracle/qbio_data_analysis_MiracleI/week7_MAF/oncoplot.png")

clinic = maf_object@clinical.data
young_patient_ids = clinic$Tumor_Sample_Barcode[clinic$age_category == "young"]
young_maf = subsetMaf(maf = maf_object, young_patient_ids)

old_maf = subsetMaf(maf = maf_object, clinic$Tumor_Sample_Barcode[clinic$age_category == "old"])

coOncoplot(m1 = young_maf,
           m2 = old_maf,
           m1Name = "Young Patients         ",
           m2Name = "Old Patients")

ggsave("/Users/miracle/qbio_data_analysis_MiracleI/week7_MAF/coOncoplot.png")

## Here I noticed that all the genes aren't more mutated in one group more than another. Specifically 4 have greater mutation in 
## old and 2 have greater mutation in young patients. This was surprising to me becuase I thought old patients would have more 
## mutations in all cateogries 

lollipopPlot(maf_object, gene = "TP53")
ggsave("/Users/miracle/qbio_data_analysis_MiracleI/week7_MAF/lollipopTP53.png")

lollipopPlot2(m1 = young_maf,
              m2 = old_maf,
              m1_name = "Young Patients         ",
              m2_name = "Old Patients",
              gene = "TP53")
## More mutation for old patients, especially in P53 and most commonly mis-sense mutations.

ggsave("/Users/miracle/qbio_data_analysis_MiracleI/week7_MAF/lollipopPlot2TP53.png")

## 10 don't have mutations for Gene A or B, which amounts to 10% and mutations are not independent of each other
## as a mutation in A shows likliness of mutation in B as well

## b=7; c=2; d=35; e=37; f=42;
## 6   7   13
## 2   35  37
## 8   42  50

geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")

geneB_maf <- subsetMaf(maf = maf_object, 
                       genes = "KRAS")

geneA_maf
geneB_maf
## subsetMaf() makes dataframes of a smaller size and are useful becuase they contain the necessary variables  
maf_object@data
maf_object@clinical.data
## I think while possible, it is unlikely that there is only one mutation per sample in Gene A since there are 
## so many possible mutations 
geneA_maf@clinical.data 
geneB_maf@clinical.data 
geneA_maf@data 
geneB_maf@data 
## The number of samples in data is not the same as in clinical.data since a few patients have a greater amount of mutations
## so they are included in different data

mut_bc_geneA = geneA_maf@clinical.data$Tumor_Sample_Barcode
mut_bc_geneB = geneB_maf@clinical.data$Tumor_Sample_Barcode

num_mut_geneA = length(mut_bc_geneA)
num_mut_geneA
num_mut_geneB = length(mut_bc_geneB)
num_mut_geneB

mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB)
num_mut_geneAB

num_mut_geneA_only = 135
num_mut_geneB_only = 85

num_neither_mutation = 78

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_neither_mutation), 
                       nrow=2)

contig_table

fe_results <- fisher.test(contig_table)
fe_results
# There is a low value for p so the null (mutations are independent) is rejected  
