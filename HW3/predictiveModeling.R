setwd("/Users/ellen/BMI215/HW3/")
library(pROC)
library(ResourceSelection)
library(plyr)
library(matrixStats)
library(MASS)

# Read in data
sapsScores <- read.table("dataFiles/saps_scores.txt", header=T)
akiData <- read.table("dataFiles/aki_data.txt", header=T)
akiData <- as.data.frame(akiData)
cvIndices <- read.table("dataFiles/cvIndices.txt", header=T)

#------------------------------------------------------------------------------------------
# 2.1 - EVALUATING PREDICTIVE ACCURACY OF SAPS SCORES
#------------------------------------------------------------------------------------------
# Convert SAPS scores to predicted mortality
predictedMortality <- sapsScores
logit <- -7.7631+0.0737*predictedMortality$saps+0.9971*log(predictedMortality$saps+1)
predictedMortality$saps <- exp(logit)/(1+exp(logit))
names(predictedMortality) = c("subject_id", "predicted_mortality")

# Plot histogram of predicted mortality
hist(predictedMortality$predicted_mortality,main="Histogram of Predicted Mortality",xlab="Predicted Mortality")
# Calculate AUC using predicted mortality
response <- as.numeric(akiData$expire_flg)-1
roc(response, predictedMortality$predicted_mortality)
# Plot ROC
plot.roc(response, predictedMortality$predicted_mortality,main="ROC for Predicted Mortality")
# Calculate Hosmer-Lemeshow p-value using predicted mortality
hoslem.test(response, predictedMortality$predicted_mortality, g = 10)

#------------------------------------------------------------------------------------------
# 2.2 - MISSING VALUE IMPUTATION
#------------------------------------------------------------------------------------------
# Perform missing value imputation
akiDataImputed <- akiData
for(i in seq(4,length(akiData),3)){
  # For variables with missing values for all days (this only occurs with Bilirubin), replace with 0.7
  rowWithMissingVal <- akiDataImputed[is.na(akiDataImputed[,i]) & is.na(akiDataImputed[,i+1]) & is.na(akiDataImputed[,i+2]), i:(i+2)]
  akiDataImputed[row.names(rowWithMissingVal),i:(i+2)] <- 0.7
  
  # For variables missing only the Day 2 value, replace with average of Day 1 and 3 values
  rowWithMissingVal <- akiDataImputed[!is.na(akiDataImputed[,i]) & is.na(akiDataImputed[,i+1]) & !is.na(akiDataImputed[,i+2]), i:(i+2)]
  akiDataImputed[row.names(rowWithMissingVal),(i+1)] <- unname(rowMeans(rowWithMissingVal[,c(1,3)]))
  
  # For variables missing only the Day 1 value, replace with Day 2 value
  rowWithMissingVal <- akiDataImputed[is.na(akiDataImputed[,i]) & !is.na(akiDataImputed[,i+1]) & !is.na(akiDataImputed[,i+2]), i:(i+2)]
  akiDataImputed[row.names(rowWithMissingVal),i] <- unname(rowWithMissingVal[2])
  
  # For variables missing only the Day 1 value, replace with Day 2 value
  rowWithMissingVal <- akiDataImputed[is.na(akiDataImputed[,i]) & !is.na(akiDataImputed[,i+1]) & !is.na(akiDataImputed[,i+2]), i:(i+2)]
  akiDataImputed[row.names(rowWithMissingVal),i] <- unname(rowWithMissingVal[2])
  
  # For variables missing only the Day 3 value, replace with Day 2 value
  rowWithMissingVal <- akiDataImputed[!is.na(akiDataImputed[,i]) & !is.na(akiDataImputed[,i+1]) & is.na(akiDataImputed[,i+2]), i:(i+2)]
  akiDataImputed[row.names(rowWithMissingVal),i+2] <- unname(rowWithMissingVal[2])
  
  # For variables missing values from two days, replace with the value that is present
  rowWithMissingVal <- akiDataImputed[is.na(akiDataImputed[,i]) & is.na(akiDataImputed[,i+1]) & !is.na(akiDataImputed[,i+2]), i:(i+2)]
  akiDataImputed[row.names(rowWithMissingVal),i:(i+1)] <- unname(rowWithMissingVal[3])
  rowWithMissingVal <- akiDataImputed[is.na(akiDataImputed[,i]) & !is.na(akiDataImputed[,i+1]) & is.na(akiDataImputed[,i+2]), i:(i+2)]
  akiDataImputed[row.names(rowWithMissingVal),i:(i+2)] <- unname(rowWithMissingVal[2])
  rowWithMissingVal <- akiDataImputed[!is.na(akiDataImputed[,i]) & is.na(akiDataImputed[,i+1]) & is.na(akiDataImputed[,i+2]), i:(i+2)]
  akiDataImputed[row.names(rowWithMissingVal),i:(i+2)] <- unname(rowWithMissingVal[1])
}

# Report means and standard deviations for each completed feature
mean <- colMeans(akiDataImputed[,3:length(akiDataImputed)])
stdev <- colSds(data.matrix(akiDataImputed[,3:length(akiDataImputed)], rownames.force = NA))
featureStats <- data.frame(mean, stdev)

#------------------------------------------------------------------------------------------
# 2.3 - BUILDING MODELS WITH 3 CROSS-VALIDATION DATASETS
#------------------------------------------------------------------------------------------
# Separate out the 3 cross-validation sets and format for glm call
cvSet1 <- akiDataImputed[akiDataImputed$subject_id==cvIndices$subject_id & cvIndices$cv_index==1,]
cvSet1$subject_id <- NULL
cvSet1$age <- NULL
row.names(cvSet1) <- NULL
cvSet1$expire_flg <- as.numeric(cvSet1$expire_flg)-1
cvSet2 <- akiDataImputed[akiDataImputed$subject_id==cvIndices$subject_id & cvIndices$cv_index==2,]
cvSet2$subject_id <- NULL
cvSet2$age <- NULL
row.names(cvSet2) <- NULL
cvSet2$expire_flg <- as.numeric(cvSet2$expire_flg)-1
cvSet3 <- akiDataImputed[akiDataImputed$subject_id==cvIndices$subject_id & cvIndices$cv_index==3,]
cvSet3$subject_id <- NULL
cvSet3$age <- NULL
row.names(cvSet3) <- NULL
cvSet3$expire_flg <- as.numeric(cvSet3$expire_flg)-1

# Build 3 models, holding out three different sets
fullModel1 <- glm(formula = expire_flg ~ ., data = rbind(cvSet2,cvSet3), family=binomial(link="logit"))
summary(fullModel1)
fullModel2 <- glm(formula = expire_flg ~ ., data = rbind(cvSet1,cvSet3), family=binomial(link="logit"))
summary(fullModel2)
fullModel3 <- glm(formula = expire_flg ~ ., data = rbind(cvSet1,cvSet2), family=binomial(link="logit"))
summary(fullModel3)

# Make predictions
predict1 <- predict(fullModel1, newdata = cvSet1, type="response")
predict2 <- predict(fullModel2, newdata = cvSet2, type="response")
predict3 <- predict(fullModel3, newdata = cvSet3, type="response")

# Calculate average AUC and Hosmer-Lemeshow p-value for each cross validation model
mean(c(roc(cvSet1$expire_flg, predict1)$auc, roc(cvSet2$expire_flg, predict2)$auc, roc(cvSet3$expire_flg, predict3)$auc))
sd(c(roc(cvSet1$expire_flg, predict1)$auc, roc(cvSet2$expire_flg, predict2)$auc, roc(cvSet3$expire_flg, predict3)$auc))
hoslem.test(cvSet1$expire_flg, predict1, g = 10)
hoslem.test(cvSet2$expire_flg, predict2, g = 10)
hoslem.test(cvSet3$expire_flg, predict3, g = 10)

#------------------------------------------------------------------------------------------
# 2.4 - FORWARD AND BACKWARD STEPWISE FEATURE SELECTION
#------------------------------------------------------------------------------------------
# Perform Forward Stepwise Selection on each of the 3 cross-validation training sets
emptyModel1 <- glm(formula = expire_flg ~ 1, data = rbind(cvSet2,cvSet3), family=binomial(link="logit"))
fwd1model <- stepAIC(emptyModel1, scope=list(lower=emptyModel1, upper=fullModel1), direction="forward", trace=0)
  length(fwd1model$coefficients)-1
emptyModel2 <- glm(formula = expire_flg ~ 1, data = rbind(cvSet1,cvSet3), family=binomial(link="logit"))
fwd2model <- stepAIC(emptyModel2, scope=list(lower=emptyModel2, upper=fullModel2), direction="forward", trace=0)
  length(fwd2model$coefficients)-1
emptyModel3 <- glm(formula = expire_flg ~ 1, data = rbind(cvSet1,cvSet2), family=binomial(link="logit"))
fwd3model <- stepAIC(emptyModel3, scope=list(lower=emptyModel3, upper=fullModel3), direction="forward", trace=0)
  length(fwd3model$coefficients)-1
# Perform Backward Stepwise Selection on each of the 3 cross-validation training sets
back1model <- stepAIC(fullModel1, trace=0)
  length(back1model$coefficients)-1
back2model <- stepAIC(fullModel2, trace=0)
  length(back2model$coefficients)-1
back3model <- stepAIC(fullModel3, trace=0)
  length(back3model$coefficients)-1

# Finds the majority vote features in multiple cross-validation models
# I: Features - list of all feature names
#    model1, model2, model3 - three cross-validation models
# O: List of features that appear in at least 2 cross-validation models
findMajorityVoteFeatures <- function(features, model1, model2, model3) {
  featureCounts <- data.frame(feature = features,
                                 model1 = c(rep(0, length(features))),
                                 model2 = c(rep(0, length(features))),
                                 model3 = c(rep(0, length(features))))
  featureCounts$model1[featureCounts$feature %in% names(model1$coefficients)] <- 1    
  featureCounts$model2[featureCounts$feature %in% names(model2$coefficients)] <- 1
  featureCounts$model3[featureCounts$feature %in% names(model3$coefficients)] <- 1
  result <- featureCounts$feature[rowSums(featureCounts[,2:4])>=2]
}

# Find majority vote features by forward and backward selection
fwdSelectionFeatures <-findMajorityVoteFeatures(names(cvSet1), fwd1model,fwd2model,fwd3model)
length(fwdSelectionFeatures)
backSelectionFeatures <-findMajorityVoteFeatures(names(cvSet1), back1model,back2model,back3model)
length(backSelectionFeatures)

# Build models using forward selection majority vote features
fml <- as.formula(paste("expire_flg ~", paste(fwdSelectionFeatures, collapse="+")))
fwd1modelMajority <- glm(formula = fml, data = rbind(cvSet2,cvSet3), family=binomial(link="logit"))
summary(fwd1modelMajority)
fwd2modelMajority <- glm(formula = fml, data = rbind(cvSet1,cvSet3), family=binomial(link="logit"))
summary(fwd2modelMajority)
fwd3modelMajority <- glm(formula = fml, data = rbind(cvSet1,cvSet2), family=binomial(link="logit"))
summary(fwd3modelMajority)

# Build models using backward selection majority vote features
fml <- as.formula(paste("expire_flg ~", paste(backSelectionFeatures, collapse="+")))
back1modelMajority <- glm(formula = fml, data = rbind(cvSet2,cvSet3), family=binomial(link="logit"))
summary(back1modelMajority)
back2modelMajority <- glm(formula = fml, data = rbind(cvSet1,cvSet3), family=binomial(link="logit"))
summary(back2modelMajority)
back3modelMajority <- glm(formula = fml, data = rbind(cvSet1,cvSet2), family=binomial(link="logit"))
summary(back3modelMajority)

# Make predictions for forward selection models
predict1 <- predict(fwd1modelMajority, newdata = cvSet1, type="response")
predict2 <- predict(fwd2modelMajority, newdata = cvSet2, type="response")
predict3 <- predict(fwd3modelMajority, newdata = cvSet3, type="response")
# Calculate average AUC and Hosmer-Lemeshow p-value for each cross validation model
mean(c(roc(cvSet1$expire_flg, predict1)$auc, roc(cvSet2$expire_flg, predict2)$auc, roc(cvSet3$expire_flg, predict3)$auc))
sd(c(roc(cvSet1$expire_flg, predict1)$auc, roc(cvSet2$expire_flg, predict2)$auc, roc(cvSet3$expire_flg, predict3)$auc))
hoslem.test(cvSet1$expire_flg, predict1, g = 10)
hoslem.test(cvSet2$expire_flg, predict2, g = 10)
hoslem.test(cvSet3$expire_flg, predict3, g = 10)

# Make predictions for backward selection models
predict1 <- predict(back1modelMajority, newdata = cvSet1, type="response")
predict2 <- predict(back2modelMajority, newdata = cvSet2, type="response")
predict3 <- predict(back3modelMajority, newdata = cvSet3, type="response")
# Calculate average AUC and Hosmer-Lemeshow p-value for each cross validation model
mean(c(roc(cvSet1$expire_flg, predict1)$auc, roc(cvSet2$expire_flg, predict2)$auc, roc(cvSet3$expire_flg, predict3)$auc))
sd(c(roc(cvSet1$expire_flg, predict1)$auc, roc(cvSet2$expire_flg, predict2)$auc, roc(cvSet3$expire_flg, predict3)$auc))
hoslem.test(cvSet1$expire_flg, predict1, g = 10)
hoslem.test(cvSet2$expire_flg, predict2, g = 10)
hoslem.test(cvSet3$expire_flg, predict3, g = 10)

#------------------------------------------------------------------------------------------
# 2.5 - TRYING SOMETHING NEW: INTERACTION TERMS
#------------------------------------------------------------------------------------------
# Build models with interaction terms, using backward selection majority vote features
fml <- as.formula("expire_flg ~ max_bili_1st_day * max_bili_3rd_day +
                            max_sodium_2nd_day + 
                           min_sysbp_1st_day * min_sysbp_3rd_day +
                          min_wbc_2nd_day * min_wbc_3rd_day +
                          max_bun_1st_day * max_bun_3rd_day + 
                          max_temp_2nd_day * max_temp_3rd_day +
                          min_bicarbonate_3rd_day + 
                          out_2nd_day + out_3rd_day")

interactionModel1 <- glm(formula = fml, data = rbind(cvSet2,cvSet3), family=binomial(link="logit"))
summary(interactionModel1)
interactionModel2 <- glm(formula = fml, data = rbind(cvSet1,cvSet3), family=binomial(link="logit"))
summary(interactionModel2)
interactionModel3 <- glm(formula = fml, data = rbind(cvSet1,cvSet2), family=binomial(link="logit"))
summary(interactionModel3)

# Make predictions for forward selection models
predict1 <- predict(interactionModel1, newdata = cvSet1, type="response")
predict2 <- predict(interactionModel2, newdata = cvSet2, type="response")
predict3 <- predict(interactionModel3, newdata = cvSet3, type="response")

cvSet1$expire_flg <- as.numeric(cvSet1$expire_flg)-1
cvSet2$expire_flg <- as.numeric(cvSet2$expire_flg)-1
cvSet3$expire_flg <- as.numeric(cvSet3$expire_flg)-1

# Calculate average AUC and Hosmer-Lemeshow p-value for each cross validation model
mean(c(roc(cvSet1$expire_flg, predict1)$auc, roc(cvSet2$expire_flg, predict2)$auc, roc(cvSet3$expire_flg, predict3)$auc))
sd(c(roc(cvSet1$expire_flg, predict1)$auc, roc(cvSet2$expire_flg, predict2)$auc, roc(cvSet3$expire_flg, predict3)$auc))
hoslem.test(cvSet1$expire_flg, predict1, g = 10)
hoslem.test(cvSet2$expire_flg, predict2, g = 10)
hoslem.test(cvSet3$expire_flg, predict3, g = 10)
