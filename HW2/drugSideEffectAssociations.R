setwd("/Users/ellen/BMI215/HW2")
library(plyr)
library(ggplot2)

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
tTest<-data.frame()
mannWhitney<-data.frame()
fisher<-data.frame()

for(i in 1:length(filteredAEs)){
  # Generate two vectors containing frequencies that the current AE was reported with each drug
  cholesterol <- subset(aeFreqs, event %in% filteredAEs[i] & drug %in% cholDrugs)$freq
  notCholesterol <- subset(aeFreqs, event %in% filteredAEs[i] & !drug %in% cholDrugs)$freq
  
  # Perform t-test and store results
  result <- t.test(cholesterol, notCholesterol, var.equal=T)
  tTest <- rbind(tTest, c(result$p.value, result$estimate[1], result$estimate[2]))
  
  # Perform Mann-Whitney test
  result <- wilcox.test(cholesterol, notCholesterol, correct=F)
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

# for student's t-test
tTest <- tTest[order(tTest$pVal),]
enriched_tTest <- tTest[tTest$meanOccurrenceInCholDrugs > tTest$meanOccurrenceInNonCholDrugs,][1:10,]
enriched_tTest$event <- filteredAEs[strtoi(row.names(enriched_tTest))]
diminished_tTest <- tTest[tTest$meanOccurrenceInCholDrugs < tTest$meanOccurrenceInNonCholDrugs,][1:10,]
diminished_tTest$event <- filteredAEs[strtoi(row.names(diminished_tTest))]

# for Mann-Whitney test
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
top5AEs <- filteredAEs[strtoi(row.names(fisher[1:5,]))]


# Dataframe with 6 rows (top 5 events + response variable)
# and 2574 rows (number of drugs)
Y <- data.frame(unique(aeFreqs$drug)) # Dataframe of drugs
glmData = data.frame(Y)

for(i in 1:length(top5AEs)){
  column <- aeFreqs[aeFreqs$event == top5AEs[i],]
  column <- merge(Y, column, all=T)
  glmData[,i+1] <- c(column$freq)
}
names(glmData)=c("drug", make.names(top5AEs))
glmData[is.na(glmData)] <- 0
glmData$Y[glmData$drug %in% cholDrugs] <- 1
glmData$Y[!glmData$drug %in% cholDrugs] <- 0
glmData$drug <- NULL

glmModel <- glm(formula = Y ~ RHABDOMYOLYSIS + MUSCULAR.WEAKNESS	+ BLOOD.CREATINE.PHOSPHOKINASE.INCREASED +	BLOOD.TRIGLYCERIDES.INCREASED + MUSCLE.SPASMS, data = glmData, family=binomial(link="logit"))
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
