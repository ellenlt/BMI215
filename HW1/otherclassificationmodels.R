setwd("/Users/ellen/BMI215/HW1")
getwd()
library(LiblineaR)
library(randomForest)

# Load raw data from .csv file
trainingData <- read.csv("./training.csv", header=TRUE, as.is=TRUE)
testData <- read.csv("./test.csv", header=TRUE, as.is=TRUE)
diseaseName <- "Obesity"

# Convert our features of interest to a matrix
X <- do.call(rbind, trainingData[1:wordsInVocab])
X <- t(X)
X <- subset(X, select=c("obese","obesity","morbid","sleep","apnea","albuterol","obstructive","wheezing","chronic","complications","morbidly","allergy","bowel","study"))
dfX <- data.frame(X)

# Convert test data to a dataframe
newX <- do.call(rbind, testData[1:wordsInVocab])
newX <- t(newX)
newX <- subset(newX, select=c("obese","obesity","morbid","sleep","apnea","albuterol","obstructive","wheezing","chronic","complications","morbidly","allergy","bowel","study"))
newX <- data.frame(matrix(unlist(testData[1:wordsInVocab]), nrow=length(testData[,1]), byrow=F),stringsAsFactors=F)
newX <- data.frame(X)

# Convert our disease of interest into a binary vector
y <- trainingData[,diseaseName]
y[y=="Y"] <- 1
y[y=="N"] <- 0
class(y)<- "numeric"

# Mean absolute error
mae <- function(errors) { mean(abs(errors)) }

# Logistic regression
glmModel <- glm(formula = y ~ X, family="binomial")
summary(glmModel)
# Training error
glmTrainingPredictions <- predict(glmModel, se.fit=T)
glmTrainingError <- mae(glmTrainingPredictions$se.fit)
glmTrainingError
# Test error
glmTestingPredictions <- predict(glmModel, newX, se.fit=T)
glmTestingError <- mae(glmTestingPredictions$se.fit)
glmTestingError

# Support vector machine
# Create matrix testData and trainData with data
# Create matrix y with answers
svmModel <- LiblineaR(trainData, y, type=1)
svmTrainingPredictions <- predict(svmModel, trainData)
svmTestPredictions <- predict(svmModel, testData)

# Random forest
rfModel <- randomForest(trainData, y, ntree=500, do.trace=T)
rfTrainingPredictions <- predict(rfModel)
rfTestingPredictions <- predict(rfModel, testData)
rfImportantFeatures <- rfModel$importance