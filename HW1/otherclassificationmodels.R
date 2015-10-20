setwd("/Users/ellen/BMI215/HW1")
getwd()
library(LiblineaR)
library(randomForest)

# Load raw data from .csv file
trainingData <- read.csv("./training.csv", header=TRUE, as.is=TRUE)
testData <- read.csv("./test.csv", header=TRUE, as.is=TRUE)
diseaseName <- "Obesity"

# Convert our features of interest and output of interest to a matrix
X <- do.call(rbind, trainingData)
X <- t(X)
X <- subset(X, select=c("Obesity","obese","obesity","morbid","sleep","apnea","albuterol","obstructive","wheezing","chronic","complications","morbidly","allergy","bowel","study"))
dfX <- data.frame(X)

# Convert test data to a dataframe
newX <- do.call(rbind, testData)
newX <- t(newX)
newX <- subset(newX, select=c("Obesity","obese","obesity","morbid","sleep","apnea","albuterol","obstructive","wheezing","chronic","complications","morbidly","allergy","bowel","study"))
dfNewX <- data.frame(newX)

# Mean absolute error
mae <- function(errors) { mean(abs(errors)) }

# Logistic regression
glmModel <- glm(formula = Obesity ~ obese + obesity + morbid + sleep + apnea + albuterol + obstructive + wheezing + chronic + complications + morbidly + allergy + bowel + study, data = dfX, family="binomial")
summary(glmModel)
# Training error
glmTrainingPredictions <- predict(glmModel, se.fit=T)
# Convert to Y/N
glmTrainingPredictions$fit[glmTrainingPredictions$fit < 0] <- "N"
glmTrainingPredictions$fit[glmTrainingPredictions$fit != "N"] <- "Y"
table(glmTrainingPredictions$fit, trainingData$Obesity)
# Test error
glmTestingPredictions <- predict(glmModel, newdata = dfNewX, se.fit=T)
glmTestingPredictions$fit[glmTestingPredictions$fit < 0] <- "N"
glmTestingPredictions$fit[glmTestingPredictions$fit != "N"] <- "Y"
table(glmTestingPredictions$fit, testData$Obesity)

# Support vector machine
svmModel <- LiblineaR(dfX[,2:length(dfX)], dfX[,1], type=1)
summary(svmModel)
svmTrainingPredictions <- predict(svmModel, dfX[,2:length(dfX)])
svmTrainingPredictions <- table(svmTrainingPredictions$predictions,trainingData$Obesity)
svmTestPredictions <- predict(svmModel, dfNewX[,2:length(dfNewX)])
svmTestPredictions <- table(svmTestPredictions$predictions,testData$Obesity)
svmTrainingPredictions
svmTestPredictions

# Random forest
rfModel <- randomForest(dfX[,2:length(dfX)], dfX[,1], ntree=500, do.trace=T)
rfTrainingPredictions <- predict(rfModel)
length(rfTrainingPredictions)
table(rfTrainingPredictions,trainingData$Obesity)
rfTestingPredictions <- predict(rfModel, dfNewX[,2:length(dfX)])
table(rfTestingPredictions,testData$Obesity)
rfModel$importance


# Create vector y with answers
y <- trainingData[,diseaseName]
y[y=="Y"] <- 1
y[y=="N"] <- 0
class(y)<- "numeric"
y <- data.frame(y)
y$ID <- 1:length(y)