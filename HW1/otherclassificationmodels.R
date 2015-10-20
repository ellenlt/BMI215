setwd("/Users/ellen/BMI215/HW1")
getwd()
library(LiblineaR)
library(randomForest)

trainingData <- read.csv("./training.csv", header=TRUE, as.is=TRUE)
testData <- read.csv("./test.csv", header=TRUE, as.is=TRUE)
diseaseName <- "Obesity"
wordsInVocab <<- 10582

# Logistic regression
glmModel <- glm(y ~ X, family="binomial")
glmTrainingPredictions <- predict(glmModel)
glmTestingPredictions <- predict(glmModel, testdata)

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

# filterFeatures
# Input: Data as a list
# Output: Data frame containing columns for all the words as well as the disease of interest.
# Returned data frame doesn't include columns for other diseases
filterFeatures <- function(data, disease) {
  words <- data.frame(matrix(unlist(data[1:wordsInVocab]), nrow=length(data[,1]), byrow=F),stringsAsFactors=F)
  names(words) <- names(data[1:wordsInVocab])
  words$ID <- 1:nrow(words)
  
  disease <- data.frame(unlist(data[,diseaseName]),stringsAsFactors=F)
  names(disease) <- diseaseName
  disease$ID <- 1:nrow(disease)
  
  dataTable <- merge(disease,words,by="ID")
  dataTable$ID <- NULL
  dataTable
}