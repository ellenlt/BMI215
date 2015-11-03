setwd("/Users/ellen/BMI215/HW2")
library(plyr)
library(ggplot2)

#REMOVE
drug=c("chol1","chol2","chol3","chol4","chol5","nchol1","nchol2","nchol3","nchol4","nchol5","chol1","nchol1","nchol1","nchol2","nchol3","nchol4","nchol5")
event=c("backache","backache","backache","backache","backache","backache","backache","backache","backache","backache","headache","stomachache","stomachache","stomachache","stomachache","stomachache","stomachache")
freq=c(0.1,0.2,0.3,0.1,0.2,0.3,0.1,0.2,0.3,0.1,0.2,0.3,0.2,0.3,0.1,0.2,0.3)
aeFreqs = data.frame(drug,event,freq)
cholDrugs=c("chol1","chol2","chol3","chol4","chol5")

# Dataframe with 3 columns:
# drug - drug names (string)
# event - name of adverse event (AE) a drug was reported with (string)
# freq - fraction of reports for the drug that list the event (float)
aeFreqs <- read.csv("./single_drug_event_frequencies.csv", header=T, as.is=TRUE)
names(aeFreqs) <- c("drug", "event", "freq")
# List/vector of cholesterol drug names (strings)
cholDrugs <- scan("./cholesterol_drugs.txt", what="ch", sep="\n")

# Dataframe where 1st column contains all adverse events (AEs) reported with a cholesterol drug
# and second column is number of different drugs with which that event co-occurred
cholAEs <- unique(aeFreqs$event[aeFreqs$drug %in% cholDrugs])
#------------------------------------------------------------------------------------------
# 1.1 - FILTER ADVERSE EVENTS
#------------------------------------------------------------------------------------------
# Dataframe where 1st column: all AEs reported with a cholesterol drug
# and 2nd column: # of different cholesterol drugs co-occuring with that AE
cholAECounts <- count(aeFreqs$event[aeFreqs$drug %in% cholDrugs])
names(cholAECounts) <- c("event", "numdrugs")
# String of AEs that co-occur with at least 5 different cholesterol drugs
commonCholAEs <- cholAECounts$event[cholAECounts$numdrugs>=5]
commonCholAEs

# Dataframe where 1st column: all AEs reported with a non-cholesterol drug
# and 2nd column: # of different non-cholesterol drugs co-occuring with that AE
nonCholAECounts <- count(aeFreqs$event[!aeFreqs$drug %in% cholDrugs])
names(nonCholAECounts) <- c("event", "numdrugs")
# String of AEs that co-occur with at least 5 different non-cholesterol drugs
commonNonCholAEs <- nonCholAECounts$event[nonCholAECounts$numdrugs>=5]
commonNonCholAEs

# String of AEs that co-occur with at least 5 different cholesterol
# and 5 different non-cholesterol drugs
filteredAEs <- intersect(commonCholAEs, commonNonCholAEs)
filteredAEs
#------------------------------------------------------------------------------------------
# 1.1 - PERFORM STUDENT'S T-TEST, MANN-WHITNEY TEST, and FISHER'S EXACT TEST
#------------------------------------------------------------------------------------------
# For each adverse event co-occuring with at least 5 cholesterol drugs and 5 other drugs,
# Perform three statistical test to determine which adverse events are most enriched/diminished
# in the response variable (whether or not a cholesterol drug or other drug is known to cause
# the adverse event)

# Inputs: index i of the event for which we would like to print a contingency table, and the
#         vector in which the event name can be found
# Output: 2x2 contingency table for the frequency with which a cholesterol or other drug
#         occurs with the event
contingencyTable <- function(i, AEs) {
  chol <- subset(aeFreqs, event %in% AEs[i] & drug %in% cholDrugs)$freq
  other <- subset(aeFreqs, event %in% AEs[i] & !drug %in% cholDrugs)$freq
  # Convert raw event freqs to bins (<=0.01, >0.01)
  chol[chol<=0.01] <- "freq<=0.01"; chol[chol!="freq<=0.01"] <- "freq>0.01"
  other[other<=0.01] <- "freq<=0.01"; other[other!="freq<=0.01"] <- "freq>0.01"
  # Create 2x2 contingency table
  drugs <- c(rep("cholesterol drugs", length(chol)),rep("other drugs", length(other)))
  table <- table(factor(drugs, levels=c("cholesterol drugs", "other drugs")), factor(c(chol, other), levels=c("freq>0.01", "freq<=0.01"))) 
  table
}

tTest<-data.frame()
mannWhitney<-data.frame()
fisher<-data.frame()

for(i in 1:length(filteredAEs)){
  # Generate two vectors containing frequencies that the current AE was reported with each drug
  cholesterol <- subset(aeFreqs, event %in% filteredAEs[i] & drug %in% cholDrugs)$freq
  notCholesterol <- subset(aeFreqs, event %in% filteredAEs[i] & !drug %in% cholDrugs)$freq
  
  # Perform t-test and store results
  result <- t.test(cholesterol, notCholesterol)
  tTest <- rbind(tTest, c(result$p.value, result$estimate[1], result$estimate[2]))
  
  # Perform Mann-Whitney test with Bonferroni Correction
  result <- wilcox.test(cholesterol, notCholesterol, correct=T)
  mannWhitney <- rbind(mannWhitney, c(result$p.value, median(cholesterol), median(notCholesterol)))
  
  # Perform Fisher's exact test
  result <- fisher.test(contingencyTable(i, filteredAEs))
  fisher <- rbind(fisher, c(result$p.value, result$estimate))
}

names(tTest) = c("pVal","meanOccurrenceInCholDrugs","meanOccurrenceInNonCholDrugs")
names(mannWhitney) = c("pVal","medianOccurrenceInCholDrugs","medianOccurrenceInNonCholDrugs")
names(fisher) = c("pVal","OR")

#------------------------------------------------------------------------------------------
# 1.1 - REPORT EVENTS MOST ENRICHED and DIMINISHED FOR CHOLESTEROL DRUGS
#------------------------------------------------------------------------------------------
# For each of the three tests performed, report the 10 events most enriched for cholesterol drugs
# as well as the 10 events most diminished for cholesterol drugs, based on lowest p-value

# for student's t-test (assuming unequal variances)
tTest <- tTest[order(tTest$pVal),]
enriched_tTest <- tTest[tTest$meanOccurrenceInCholDrugs > tTest$meanOccurrenceInNonCholDrugs,][1:10,]
enriched_tTest$event <- filteredAEs[strtoi(row.names(enriched_tTest))]
diminished_tTest <- tTest[tTest$meanOccurrenceInCholDrugs < tTest$meanOccurrenceInNonCholDrugs,][1:10,]
diminished_tTest$event <- filteredAEs[strtoi(row.names(diminished_tTest))]

# for Mann-Whitney test (with Bonferroni Correction)
mannWhitney <- mannWhitney[order(mannWhitney$pVal),]
enriched_mannWhitney <- mannWhitney[mannWhitney$medianOccurrenceInCholDrugs > mannWhitney$medianOccurrenceInNonCholDrugs,][1:10,]
enriched_mannWhitney$event <- filteredAEs[strtoi(row.names(enriched_mannWhitney))]
diminished_mannWhitney <- mannWhitney[mannWhitney$medianOccurrenceInCholDrugs < mannWhitney$medianOccurrenceInNonCholDrugs,][1:10,]
diminished_mannWhitney$event <- filteredAEs[strtoi(row.names(diminished_mannWhitney))]

# for Fisher's exact test
fisher <- fisher[order(fisher$pVal),]
fisher <- fisher[is.finite(fisher$OR),]
enriched_fisher <- fisher[fisher$OR>1,][1:10,]
enriched_fisher$event <- filteredAEs[strtoi(row.names(enriched_fisher))]
diminished_fisher <- fisher[fisher$OR<1,][1:10,]
diminished_fisher$event <- filteredAEs[strtoi(row.names(diminished_fisher))]

#------------------------------------------------------------------------------------------
# 1.2 - PLOT EVENT FREQUENCY DISTRIBUTIONS
#------------------------------------------------------------------------------------------
# For the top 2 events enriched in cholesterol drugs and the top 2 events diminished in cholesterol drugs
# (based on lowest Mann-Whitney p-value), plot event frequency distributions for cholesterol drugs and other drugs.

# Inputs: the index i of the adverse event in aeFreqs that is to be plotted,
#         and a boolean indicating whether or not the event is enriched
# Output: frequency distribution for an adverse event and its co-occurrence with cholesterol and other drugs
plotMannWhitneyDist <- function(i, enriched) {
  if(enriched) {
    # Event enriched in cholesterol drugs
    chol <- subset(aeFreqs, event %in% enriched_mannWhitney$event[i] & drug %in% cholDrugs)$freq
    other <- subset(aeFreqs, event %in% enriched_mannWhitney$event[i] & !drug %in% cholDrugs)$freq
    title <- paste("Event Frequency Distribution for\n",enriched_mannWhitney$event[i],"\n#",i,"Event Enriched in Cholesterol Drugs")
  } else {
    # Event diminished in cholesterol drugs
    chol <- subset(aeFreqs, event %in% diminished_mannWhitney$event[i] & drug %in% cholDrugs)$freq
    other <- subset(aeFreqs, event %in% diminished_mannWhitney$event[i] & !drug %in% cholDrugs)$freq
    title <- paste("Event Frequency Distribution for\n",diminished_mannWhitney$event[i],"\n#",i,"Event Diminished in Cholesterol Drugs")
  }
  data <- data.frame(drugtype=factor(c(rep("Cholesterol Drugs",length(chol)),rep("Other Drugs",length(other)))),frequency=c(chol,other))
  plot <- ggplot(data, aes(x=frequency, colour=drugtype))+geom_density()
  plot + ggtitle(title)
}

# Plot freq distributions
plotMannWhitneyDist(1, T) # for most enriched
plotMannWhitneyDist(2, T) # for 2nd most enriched
plotMannWhitneyDist(1, F) # for most diminished
plotMannWhitneyDist(2, F) # for second most diminished

#------------------------------------------------------------------------------------------
# 1.3 - 2X2 CONTINGENCY TABLES
#------------------------------------------------------------------------------------------
# Print out contingency tables for top two enriched and top two diminished events
# used in Fisher's test

paste("Contingency Table for",enriched_fisher$event[1])
print(contingencyTable(1, enriched_fisher$event))
paste("Contingency Table for",enriched_fisher$event[2])
print(contingencyTable(2, enriched_fisher$event))
paste("Contingency Table for",diminished_fisher$event[1])
print(contingencyTable(1, diminished_fisher$event))
paste("Contingency Table for",diminished_fisher$event[2])
print(contingencyTable(2, diminished_fisher$event))

#------------------------------------------------------------------------------------------
# 1.4 - PLOT P-VALUE VS. ODDS RATIO
#------------------------------------------------------------------------------------------
# Plot the Fisher's exact test p-values vs odds ratios for each adverse event
ggplot(data=fisher, aes(x=OR, y=pVal, group=1)) + geom_point() + ggtitle("Fisher's Exact Test:\np-value vs odds ratio")

#------------------------------------------------------------------------------------------
# 1.5 - ANALYZING ENRICHED EVENTS FOR PRESCRIPTION BIAS AND DESIRED SIDE EFFECTS
#------------------------------------------------------------------------------------------
# Print out top 30 most enriched events for cholesterol drugs from Mann-Whitney data
top30enriched_mannWhitney <- mannWhitney[mannWhitney$medianOccurrenceInCholDrugs > mannWhitney$medianOccurrenceInNonCholDrugs,][1:30,]
top30enrichedAEs <- filteredAEs[strtoi(row.names(top30enriched_mannWhitney))]
for(i in 1:length(top30enrichedAEs)){
  print(top30enrichedAEs[i])
}

#------------------------------------------------------------------------------------------
# 2.1 - ASSEMBLING SIDE-EFFECT PROFILE FOR CHOLESTEROL DRUGS
#------------------------------------------------------------------------------------------
# We filter features down to the top 5 events most associated with cholesterol drugs 
# (according to Fisher's test), and we fit a logistic regression model to the data,
# using these 5 events as predictors for the outcome variable, whether the drug is
# a cholesterol drug or not.

# Vector of top 5 events most associated with cholesterol drugs (in either direction) based
# on Fisher's test
top5AEsByFisher <- filteredAEs[strtoi(row.names(fisher[1:5,]))]

# Input: 1) Dataframe with at least the following 3 columns:
#             "drug": drug names (string)
#             "event": adverse event names (string)
#             "freq": frequency of co-occurence between a drug and adverse event (float)
#       2) Vector containing all adverse events (strings)
#       3) Dataframe column containing all drugs
#       4) For labeled data, a vector containing all sample names (strings) that should be labeled as positive
#           If this is an empty list, does not label the data and instead retains a column of drug names.
# Output: Dataframe with drugs as rows and events as columns,
#         with a binary output "Y" column
sampleFeatureFormat <- function(data, events, Y, positives) {
  result = data.frame(Y)
  # For each feature, merge a column of frequencies into the data frame
  for(i in 1:length(events)){
    column <- data[data$event %in% events[i],]
    column <- merge(Y, column, all=T)
    result[,i+1] <- c(column$freq)
  }
  names(result)=c("drug", make.names(events))
  result[is.na(result)] <- 0  # Add zeros for missing values
  # If positive label names are provided, label
  # data. Otherwise, retain a column with the drug names
  if(length(positives)>0) {
    result$Y[result$drug %in% positives] <- 1
    result$Y[!result$drug %in% positives] <- 0
    result$drug <- NULL    
  }
  result
}

Y <- data.frame(drug=unique(aeFreqs$drug)) # Dataframe column of all possible drugs
trainingData <- sampleFeatureFormat(aeFreqs, top5AEsByFisher, Y, cholDrugs)
# Build model
glmModel <- glm(formula = Y ~ RHABDOMYOLYSIS + MUSCULAR.WEAKNESS  + BLOOD.CREATINE.PHOSPHOKINASE.INCREASED +  BLOOD.TRIGLYCERIDES.INCREASED + MUSCLE.SPASMS, data = trainingData, family=binomial(link="logit"))
summary(glmModel)

#------------------------------------------------------------------------------------------
# 2.2 - VALIDATING SIDE-EFFECT PROFILE FOR CHOLESTEROL DRUGS AND AUC CURVE
#------------------------------------------------------------------------------------------
# Dataframe with 3 columns:
# drug - drug names (string)
# event - name of adverse event (AE) a drug was reported with (string)
# freq - fraction of reports for the drug that list the event (float)
aeFreqsRecent <- read.csv("./single_drug_event_frequencies_validation.csv", header=T, as.is=TRUE)
names(aeFreqsRecent) <- c("drug", "event", "freq")
Y <- data.frame(drug=unique(aeFreqsRecent$drug)) # Dataframe column of all possible drugs
testData <- sampleFeatureFormat(aeFreqsRecent, top5AEsByFisher, Y, cholDrugs)

# Vector of probability thresholds above which to classify a sample as positive
threshold <- c(seq(0,1,by=0.001))

# Predictions
glmPredictions <- predict(glmModel, newdata = testData, type="response")

# Dataframe where 1st column is TPR and 2nd column is FPR
rocData <- data.frame()

for(i in 1:length(threshold)) {
  # Convert fitted probabilities to binary depending on the current threshold
  predictions <- glmPredictions
  predictions[predictions >= threshold[i]] <- 1
  predictions[predictions < threshold[i]] <- 0
  # Compute true positive rate and false positive rate
  table <- table(factor(predictions, levels=c("0", "1")), factor(testData$Y, levels=c("0", "1")), dnn=c("predicted", "actual")) 
  tpr <- table["1","1"]/(table["1","1"] + table["0","1"])
  fpr <- 1 - table["0","0"]/(table["0","0"] + table["1","0"])
  rocData <- rbind(rocData, c(tpr, fpr))
}
names(rocData) <- c("TPR","FPR")
rocData <- cbind(rocData, threshold)

# Plot receiver operating characteristic curve
ggplot(data=rocData, aes(x=FPR, y=TPR, group=1)) + geom_line() + geom_point() + ggtitle("ROC Curve")

#------------------------------------------------------------------------------------------
# 2.3 - AUC
#------------------------------------------------------------------------------------------

# Inputs: vector of x coordinates, vector of y coordinates
# Output: area under curve
auc <- function(x, y) {
  result <- 0
  for(p in 1:(length(x)-1)) {
    y0 <- y[p]
    x0 <- x[p]
    y1 <- y[p+1]
    x1 <- x[p+1]
    square <- abs(x0-x1)*min(y0,y1)
    triangle <- (abs(x0-x1)*abs(y0-y1))/2
    result <- result + square + triangle
  }  
  result
}

# Area under the ROC curve
aucValue <- auc(rocData$FPR, rocData$TPR)

#------------------------------------------------------------------------------------------
# 3.1 - APPLY MODEL TO NEW DATA: PROCESSING THE NEW DATA
#------------------------------------------------------------------------------------------
# Read in pair drug data and convert to a samples x features matrix
pairDrugAeFreqs <- read.csv("./pair_drug_event_frequencies.csv", header=T, as.is=TRUE)
names(pairDrugAeFreqs) <- c("drug","drug1","drug2","event","freq")
Y <- data.frame(drug=unique(pairDrugAeFreqs$drug)) # Dataframe column of all possible drugs
pairDrugData <- sampleFeatureFormat(pairDrugAeFreqs, top5AEsByFisher, Y, c())

# Number of unique drug pairs in the dataset
nrow(pairDrugData)

# Filter out all non-cholesterol drugs
pairDrugData <- cbind(ID=rownames(pairDrugData),pairDrugData) # Split drug pair into separate columns with name of each drug
pairDrugData <- cbind(pairDrugData, data.frame(do.call('rbind', strsplit(as.character(pairDrugData$drug),',',fixed=TRUE))))
pairDrugData <- pairDrugData[!as.character(pairDrugData$X1) %in% cholDrugs & !as.character(pairDrugData$X2) %in% cholDrugs,]  #Filter

# Number of unique non-cholesterol drug pairs
nrow(pairDrugData)

# Clean up dataframe
pairDrugData$ID <- NULL
pairDrugData$X1 <- NULL
pairDrugData$X2 <- NULL
row.names(pairDrugData) <- NULL

#------------------------------------------------------------------------------------------
# 3.3 - APPLYING MODEL TO NEW DATA:
#       FIND DRUG PAIRS THAT MATCH THE SIDE-EFFECT PROFILE OF A CHOLESTEROL DRUG
#------------------------------------------------------------------------------------------
# At or above this threshold, classify drug pair as matching side effect profile for cholesterol drugs
optimalProbThreshold <- 0.009

# Apply logistic model we built in part 2.1 to the new data
glmPredictions <- predict(glmModel, newdata = pairDrugData, type="response")

# Restrict list to those whose score is greater than the threshold and rank pairs by profile score
glmPredictions <- glmPredictions[glmPredictions >= optimalProbThreshold]
glmPredictions <- sort(glmPredictions[glmPredictions >= optimalProbThreshold], decreasing=T)
length(glmPredictions)  # Sanity check - should be ~1700

# Submit list of drug-pairs as tab-delimited file called "ps2_problem3.tsv"
drugPairMatches <- data.frame(pairDrugData[names(glmPredictions),]$drug, glmPredictions)
rownames(drugPairMatches) <- NULL
names(drugPairMatches) <- c("drug", "score")
write.table(drugPairMatches, file="ps2_problem3.tsv", sep="\t")

#------------------------------------------------------------------------------------------
# 3.4 - COMPARE OUR RESULTS TO VA'S LIST
#------------------------------------------------------------------------------------------
# List of 3086 known drug-drug interactions; obtained from Veteran's Association Hospital in Tucson, AZ
knownInteractions <- read.csv("./va_drug_drug_interactions.csv", header=T, as.is=TRUE)

# Merge two drugs into single comma-delimited string for easy comparison
knownInteractions <- within(knownInteractions, drug <- paste(drug1,drug2,sep=','))

# Find drug pairs we predicted as interacting to cause cholesterol-drug-like side effects,
# which also appear in the VA's list of known drug interactions
results <- drugPairMatches[as.character(drugPairMatches$drug) %in% knownInteractions$drug,]
indices <- c(which(knownInteractions$drug %in% results$drug))
# List of drug pairs and their severities
results <- knownInteractions[indices,c("drug","type")]




