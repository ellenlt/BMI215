# Function takes in two vectors and returns 2x2 table
# for disease ("Y" or "N") vs word occurrence ("1" vs "0")
getTable <- function(wordOccurrences, diseasePresence) {
  result <- table(factor(wordOccurrences, levels=c("0", "1")), factor(diseasePresence, levels=c("N", "Y")))
}

# Function which takes in a 2x2 table for a given word/disease pair
# and returns the term frequency difference
termFreqDiff <- function(table) {
  probWordGivenDisease <- table["1","Y"]
  probWordGivenNoDisease <- table["1","N"]
  result <- probWordGivenDisease - probWordGivenNoDisease
}



# Takes in a 2x2 table for a given word/disease pair
# and returns the information gain
infoGain <- function(table) {
  result <- 0
  probabilities <- prop.table(table)
  probWord <- prop.table(margin.table(table,1))
  probDisease <- prop.table(margin.table(table,2))
  for (w in 1:2) {
    for (d in 1:2) {
      if (probabilities[w, d] > 0) {
        result <- result + probabilities[w, d]*log2(probabilities[w, d]/(probWord[[w]]*probDisease[[d]]))  
      }
    }
  }
  result
}

# Takes in a 2x2 table for a given word/disease pair
# and returns the chi-squared value
chiSquared <- function(table) {
  result <- 0
  probabilities <- prop.table(table)
  probWord <- prop.table(margin.table(table,1))
  probDisease <- prop.table(margin.table(table,2))
  for (w in 1:2) {
    for (d in 1:2) {
      if(probWord[[w]]*probDisease[[d]] == 0) break;
      num1 <- probabilities[w, d]
      num2 <- probWord[[w]]*probDisease[[d]]
      numdiff <- probabilities[w, d] - probWord[[w]]*probDisease[[d]]
      numerator <- (probabilities[w, d] - probWord[[w]]*probDisease[[d]])^2
      denominator <- probWord[[w]]*probDisease[[d]]
      result <- result + (probabilities[w, d] - probWord[[w]]*probDisease[[d]])^2 / (probWord[[w]]*probDisease[[d]])
    }
  }
  result
}

# For a given disease, calculates the term frequency difference, information gain,
# and chi-squared value for each possible word and returns them in a data frame
calcStats = function(disease, data) {
  numFeatures <- dim(data)[2]
  words <- names(data[1:(numFeatures-16)])
  termFreqDiff <- c(rep(0,(numFeatures-16)))
  infoGain <- c(rep(0,(numFeatures-16)))
  chiSquared <- c(rep(0,(numFeatures-16)))
  result <- data.frame(termFreqDiff, infoGain, chiSquared)
  row.names(result) <- words
  
  for(w in 1:length(words)) {
    wordVec <- data[,w]
    diseaseVec <- data[,disease]
    table <- getTable(wordVec, diseaseVec)
    result[w,"termFreqDiff"] <- termFreqDiff(table)
    result[w,"infoGain"] <- infoGain(table)
    result[w,"chiSquared"] <- chiSquared(table)
  }
  result
}

# For a given disease, returns a data table with the ten most highly-selected words,
# for each of the three feature selection techniques
selectTopTenWords = function(disease, data) {
  statsVals <- calcStats(disease, data)
  
  topTfd <- rownames(statsVals[order(-statsVals$termFreqDiff),][1:10,])
  topIg <- rownames(statsVals[order(-statsVals$infoGain),][1:10,])
  topX2 <- rownames(statsVals[order(-statsVals$chiSquared),][1:10,])
  result <- data.frame(topTfd, topIg, topX2)
}

setwd("/Users/ellen/BMI215/HW1")
getwd()
trainingData <- read.csv("./training.csv", header=TRUE, as.is=TRUE)
testData <- read.csv("./test.csv", header=TRUE, as.is=TRUE)

asthma <- selectTopTenWords("Asthma", trainingData)
depression <- selectTopTenWords("Depression", trainingData)
diabetes <- selectTopTenWords("Diabetes", trainingData)
gallstones <- selectTopTenWords("Gallstones", trainingData)
OA <- selectTopTenWords("OA", trainingData)
obesity <- selectTopTenWords("Obesity", trainingData)
pvd <- selectTopTenWords("PVD", trainingData)