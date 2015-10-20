# evaluateModel
# Input: Vector of predictions and vector of answers
# Output: Error rate
evaluateModel <- function(predictions, answers) {
  result <- table(predictions, answers)
}

# useModel
# Input: Model of log conditional probabilities (Data frame where 1st and 2nd rows contain log(P(w|D)P(D))
#       for disease D = No and D = Yes, respectively. One column per word w.)
#       Disease string name, and 
#       data as data frame
# Output: A vector of predictions of whether each document is associated with given disease
useModel <- function(disease, data, model) {
  data[,disease] <- NULL
  result <- c(rep(0,length(data[,1])))
  for(i in 1:length(data[,1])) {
    logProbDisease <- 0
    logProbNoDisease <- 0
    for(w in 1:wordsInVocab) {
      if(data[i,w]==1) {
        logProbNoDisease <- logProbNoDisease + model[1,w]
        logProbDisease <- logProbDisease + model[2,w]
      }
    }
    if(logProbDisease > logProbNoDisease) {
      result[i] = "Y"
    } else {
      result[i] = "N"
    }
  }
  result
}

#trainModel
# Input: Data frame containing vectors for each word and disease of interest
#         Disease of interest, as a string
# Output: Data frame where 1st and 2nd rows contain log(P(w|D)P(D))
#         for disease D = No and D = Yes, respectively. One column per word w.
trainModel <- function(disease, data) {
  counts <- generateCounts(disease, data)
  model <- generateCondProbs(counts, data)
  model
}

# csvToDf
# Input: Data as a list
# Output: Data frame containing columns for all the words as well as the disease of interest.
# Returned data frame doesn't include columns for other diseases
csvToDf <- function(data, diseaseName) {
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

# generateCounts
# Input: Disease name (as a string), data (in a dataframe)
# Output: Data frame where 1st and 2nd rows are how many times each word appears
# when disease = N and Y, respectively. 3rd and 4th rows are how many times each word
# does not appear when disease = N and Y, respectively.
generateCounts <- function(disease, data) {
  # Count when word appears (w=1)
  counts1 <- ddply(data,disease,colwise(sum))
}

# generateCondProbs
# Input: Data for a given disease and all words, in correctly processed data frame
# Output: Data frame where 1st and 2nd rows contain log(P(w=1|D)P(D))
# for disease D = No and Yes, respectively, and the 3rd and 4th rows contain
# log(P(w=0|D)P(D)) for D = No and Yes, respectively. One column per word w.
generateCondProbs <- function(counts, data) {
  wordsPerClass <- rowSums(counts[2:length(counts)])
  probDisease <- sum(data[,1]=="Y")/length(data[,1])
  probNoDisease <- sum(data[,1]=="N")/length(data[,1])
  
  condProbs <- counts
  condProbs[,diseaseName] <- NULL
  if(LaPlace==TRUE) {
    condProbs[1,] <- log10((condProbs[1,]+1)/(wordsPerClass[1]+length(condProbs[1,]))) + log10(probNoDisease)
    condProbs[2,] <- log10((condProbs[2,]+1)/(wordsPerClass[2]+length(condProbs[1,]))) + log10(probDisease)    
  } else {
    condProbs[1,] <- log10((condProbs[1,])/(wordsPerClass[1])) + log10(probNoDisease)
    condProbs[2,] <- log10((condProbs[2,])/(wordsPerClass[2])) + log10(probDisease)
  }
  condProbs
}

sample <- trainingDataFrame[1:20,1:20]
sampleCounts <- generateCounts(diseaseName, sample)
#sampleProbs <- generateCondProbs(sampleCounts, sample)

zero <- function(x) sum(x == 0)

setwd("/Users/ellen/BMI215/HW1")
getwd()
library(plyr)

trainingData <- read.csv("./training.csv", header=TRUE, as.is=TRUE)
testData <- read.csv("./test.csv", header=TRUE, as.is=TRUE)
diseaseName <- "Obesity"
wordsInVocab <<- 10582
# If set to TRUE, will use LaPlace smoothing to calculate probs
LaPlace <<- F

trainingDataFrame <- csvToDf(trainingData, diseaseName)
model <- trainModel(diseaseName, trainingDataFrame)
trainingDataPredictions <- useModel(diseaseName, trainingDataFrame[,2:length(trainingDataFrame)], model)
trainingDataAnswers <- trainingDataFrame[,diseaseName]
trainingError <- evaluateModel(trainingDataPredictions, trainingDataAnswers)
trainingError

testDataFrame <- csvToDf(testData, diseaseName)
testDataPredictions <- useModel(diseaseName, testDataFrame, model)
testDataAnswers <- testDataFrame[,diseaseName]
testError <- evaluateModel(testDataPredictions, testDataAnswers)
testError