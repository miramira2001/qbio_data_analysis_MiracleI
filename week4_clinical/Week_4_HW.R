if (!requireNamespace("BiocManager", quietly = TRUE, force = TRUE))
  install.packages("BiocManager")
    BiocManager::install(version = "3.13")
    
if(!require(TCGAbiolinks)) {
  BiocManager::install("TCGAbiolinks", force = TRUE)
}

library(TCGAbiolinks)
    clin_query <- GDCquery(project = "TCGA-COAD", 
                           data.category = "Clinical",
                           file.type = "xml")
    GDCdownload(clin_query) 
    
    clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
    
    names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"
    head(clinic)
str(clinic)    
clinic$race_list = as.character(clinic$race_list)
names(clinic)
plot(clinic$age_at_initial_pathologic_diagnosis,clinic$weight, xlab="Age", ylab="Weight")
boxplot(clinic$age ~ clinic$race_list)
unique(clinic$race_list)
boxplot(clinic$age ~ clinic$race_list,  las = 2, 
        cex.axis = 0.5, pars=list(par(mar=c(10,1,1,1))))
#lengthened margins
sum(is.na(clinic$race_list))
clinic$race_list[clinic$race_list == ""] <- "No data" 

min(clinic$age)
max(clinic$age)
mean(clinic$age)
median(clinic$age)
summary(clinic$age)

sum(clinic$age > 50)
sum(clinic$age < 50)
sum(clinic$age ==50)
summary(clinic)
# 524 rows in clinic; sum of < 50, > 50 and == 50 equals 524 
names(clinic)
x = clinic$age_at_initial_pathologic_diagnosis
young_patient_ids = (clinic$bcr_patient_barcode < x)
old_patient_ids = (clinic$bcr_patient_barcode < x)
# I tried to make a separate column called x to represent the boolean vector and then use this vector to separate the 
# barcodes of the old and young
clinic$age_category <- ifelse(is.na ==young_patient_ids, "young", 
                                       ifelse(is.na == old_patient_ids, "old"))
  
# I tried a couple things, this one still doesn't work but I saw something similar on stack overflow to 
# use set is.na equal to the informnation you want to put in the new column. 

clinic[1,1] 
# Top left in data frame
clinic[1,]
# Leaving the blank shows everything associated with the top left entry and since it is blank it is not specifying 
# any one single column value so it shows all
clinic[2:5,]
# Associated with the second row all the way to the 5th row and it gives all associated columns for those 
clinic[,3]
# Doesnt specify a row but tells what the value at the 3rd column gives. Some are blank as opposed to "Colon" showing 
# there is no data 

young_clinic = c("young_patient_ids")
old_clinic = c("old_patient_ids")
# Age category should be a column with patients being rows

young_clinic_one_line = clinic$age_at_initial_pathologic_diagnosis
identical(dim(young_clinic), dim(young_clinic_one_line))
# yields true result

install.packages("survival")
install.packages("survminer")
library(survival)
library(survminer) 
sum(is.na(clinic$days_to_death))
clinic$days_to_death[clinic$days_to_death == clinic$days_to_last_follow_up]  
# I tried setting them equal to each other to replace the data although I'm not 100% sure this is correct

clinic$vital_status
death_event = if((clinic$vital_status == Alive death_event = 0), (if (clinic$vital_status == Dead death_event = 1))

# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# We then create a fit object
race_fit <- surv_fit( surv_object ~ clinic$race_list, data = clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

#Plots can be improved with axis titles or adding lines for better clarity 

write.csv(clinic, "/Users/miracle/sp22_qbio_public_data_analysis_resources/week4_clinical", row.names = F)

