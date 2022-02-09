# 1. 
# Categorical - variables that categorize observations into groups - ex. favourite colours
# Discrete - data whose values are whole numbers or other discrete values - ex  population
# Continuous - measurable quantities; values that can take on any value within an interval so can be decimal values - ex. height

# 2.
names(clinic)
is.na(clinic$age_at_initial_pathologic_diagnosis)
# Since there are no true values, which represent NA values, this is a good category 

# 3.
# values are obtained by collecting the data. This is an example of continuous data however interestingly enough, in our specific 
# case the data is considered discrete as there are only integers invovled and no fractions or decimals.

# 4.
# 1. "Age and Cancer Risk" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544764/ - While it was once commonly believed that cancer 
# cannot be prevented among older patients, this talks about what can be done to help improve the prevalence of cancer among older
# individuals and the importance of certain factors like a healthy environment and a proper transition into old age can be to preventing 
# cancer.
# 2. The Biology of Aging and Cancer: A Brief Overview of Shared and Divergent Molecular Hallmarks" 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5614326/ - This article talks more about ageing markers and how they affect
# ones susceptibility to cancer. 

# 5.
# Second variable we chose is height (clinic$height). I think this is a very interesting category because it tells us about the size 
# of a patient lengthwise. We were quite curious to see if there was any possible relation to our first variable age, especially when
# it comes to cancer. When we are talking about height, this is specifically the person's length from the top of their head to the 
# bottom of the feet and this is measured with a measuring instrument such as a tape measurer or ruler on a scale. This is also a 
# good variable to look at is taking a look for na values results in 0, showing that each patient has data filled out and we don't have
# to clean to account for the na values.

sum(is.na(clinic$height))
 
# 6.
# It is common for age and height to be related. When it comes to older age, as a person gets older they tend to become shorter in height. 
# Our first variable, age, is negatively correlated with survival rate in cancer patients. The older the patient (higher age), the lower 
# the survival rate. Relating the second variable, height, to survival rate of cancer patients, the shorter the height of the patient,
# the more likely they have a decreased survival rate. 

# 7. 
# The first plot we had showed us that the lower the age, the higher the height of the patient. However, we also saw that the older people
# had a lot of varying heights. The survival plot showed that as a whole, the younger individuals had a higher survival rate. When comparing 
# height to survival, we were surprised to find that the shorter individials lived longer than the taller individuals after the initial 
# onset of colorectal cancer. 

# CODING 

# 1.
# Turned age into a categorical variable to make for the data useful for our plots. Here we have the continuous vaariable being turned into 
# categorical data.

mask=!is.na(clinic$height)
# This boolean mask has NA values in the heigh data as "false" while everything else gives "true" values 

cleaned_clinic = clinic[mask,]
# Must clean the data so that we can see the rows that contain the height values 


boxplot(cleaned_clinic$height~cleaned_clinic$age_at_initial_pathologic_diagnosis,
        xlab = "Age",
        ylab ='Height', main ="Age vs Heigh Grapht", ylim =c(120,200))

# Now we decided to do a boxplot as we felt this was a better representation of the relationship between our two variables chosen 
# Here we have age categorized into young (< 50) and old (> = 50) 

library(survival)
library(survminer)
# Installing packages 

# 2.
surv_object = Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)
# We wanted a specific time frame so we used this to initialize the survival object 

age_fit = surv_fit( surv_object ~ clinic$age_at_initial_pathologic_diagnosis, data = clinic)
# Now initialized object for age variable 


survplot = ggsurvplot(age_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# Made the survival x age plot with some changes to default to make better to read as well as a legend for clarity


age_surv_plot = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=10))
# Making some visual changes in the plot

age_surv_plot
# Now to show the survival x age plot 


# 3.

median(cleaned_clinic$height)
# To find median of our new height data so that we can see what is above average and below average height 

cleaned_clinic$height_category = ifelse(cleaned_clinic$height <= 170, "short", "tall")
# Now we needed height to be categorical for the survival plot 

surv_object <- Surv(time = cleaned_clinic$days_to_death, 
                    event = cleaned_clinic$death_event)
# Initialized a survival object for the time frame we chose

height_fit <- surv_fit( surv_object ~ cleaned_clinic$height_category, data = cleaned_clinic )
# Same as above step for fit object 

survplot = ggsurvplot(height_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
# Made the survival x age plot with some changes to default to make better to read as well as a legend for clarity

height_surv_plot = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
# Making some visual changes in the plot

height_surv_plot
# Now to show the survival x age plot 





