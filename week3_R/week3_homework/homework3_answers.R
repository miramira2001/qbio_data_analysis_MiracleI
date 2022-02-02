attenu
is.na(x)
is.na(attenu)
is.na(data)


#Excerise 1.1
attenu
x <- attenu
is.na(x)
na.omit(x)
attenu.cleaned <- na.omit(x)
head(attenu.cleaned, n=6)

#Exercise 1.2
Theoph
Theoph_2 <- Theoph
str(Theoph_2)
# Median treatment dose is 4.02 
within(Theoph_2,Theoph_2$Dose_Class <- (ifelse(a>4.02, high,ifelse(a<4.02, low)))
       head(Theoph_2, n=6)
       
#Exercise 1.3
print(getwd())
setwd("/Users/miracle")
starbucks <- read.csv("starbucks.csv")
print(starbucks)
is.na(starbucks)
is_row_empty <- c(rowSums=6)
nrow(starbucks)
y <- starbucks
is.na(y)
starbucks_cleaned <- na.omit(y)  
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories,xlab="Carbs (grams)", ylab="Carlories")
# The plot has a positive correlation between calories and amount of carbs in the drinks
max(starbucks_cleaned$Calories)

within(starbucks_cleaned,starbucks_cleaned$is_highest_fat<- (ifelse(a>9.0, TRUE,ifelse(a<11.0, FALSE)))
       
plot(starbucks_cleaned$is_highest_fat, starbucks_cleaned$Calories,xlab="fat", ylab="Carlories")

#Exercise 1.4
print(getwd())
setwd("/Users/miracle")
Batting <- read.csv("Batting.csv")
print(Batting)
Batting$HR
plot(Batting$HR, Batting$yearID, xlab = "Homeruns", ylab = "Year")
data.frame(LAA)

#Exercise 1.5
easy_plot = function(x,y color_data){
 return (x)
}
levels = ifelse(levels < median, "low", "high")
levels = factor(levels)

plot(x, y, color, pch = 20)
cor.test(x,y)

#Exercise 2.1
library(datasets)
data(iris)
summary(iris)
# Iris data set has information on 4 different attributes/variables of flowers from 3 species (setosa, versicolor, virginica)
# Contains 4 features - sepal length, sepal width, petal length and petal width 

#Exercise 2.2
#Sepal length, sepal width, petal length, and petal width are continuous because they have numerical values that can be 
#added, subtracted, etc. and do not have any jumps. Species is categorical as it does not have numeric values.

#Exercise 2.3
hist(iris$Petal.Length, prob = TRUE, lines(density(iris$Petal.Length))
hist(iris$Sepal.Length, prob = TRUE, lines(density(iris$Sepal.Length))
hist(iris$Sepal.Width)
hist(iris$Petal.Width)

#I found it interesting that the rectangles of a histogram touch the adjacent one which shows this is continuous data.

#Exercise 2.4

results_mean = mean(iris$Sepal.Width)
iris
iris_copy <- iris
numeric_values = ifelse(iris$Sepal.Width < results_mean, "small", "large")
iris$numeric_values = ifelse(iris$Sepal.Width < results_mean, "small", "large")

boxplot(iris[,1],xlab="Sepal.Width",ylab="Width(in centemeters)",
        main="Sepal.Width (Iris Data)")


#Exercise 2.5
# I think setosa looks the most unique while versicolor and virginica look more similar to each other.

pairs
pairs (iris[1:4])

#Exercise 3.1
install.packages("TCGAbiolinks")
library(TCGAbiolinks)
