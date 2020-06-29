# ============================================================================= #
# INSTALL & LOAD REQUIRED LIBRARIES
# 01. grDevices: Support for Colours and Fonts
# 02. foreign: Reading SPSS Source
# 03. caret: nearZeroVar
# 04. DMwR: knnImputation
# 05. mice
# 06. Boruta
# 07. randomForest
# 08. ggplot2
# 09. pacman
# 10. caTools
# 11. party
# 12. earth
# 13. naivebayes
# ============================================================================= #

load_lib <- c('grDevices', 'foreign', 'caret', 'DMwR', 'mice', 'Boruta', 'randomForest', 'ggplot2', 'pacman', 'missMDA', 'caTools', 'party', 'earth', 'ROSE', 'kernlab', 'naivebayes', 'RANN')
install_lib <- load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies=TRUE)
sapply(load_lib, require, character=TRUE)


# ============================================================================= #
# DEFINE REQUIRED VARIABLES
# 01. SRC: Source file location with file name
# 02. TRAIN_TEST_RATIO: Ratio to split the sample into train and test data set
# 03. NULL_THRESHOLD: Remove Features with more than a certain percentage of NA
# ============================================================================= #

SRC <- 'C:/Users/unnat/Desktop/YDS/YS_SEMCON/Source/secom_mod.SAV'
TGT <- "C:/Users/unnat/Desktop/YDS/YS_SEMCON/Target/"
TRAIN_TEST_RATIO <- 0.8
NULL_THRESHOLD <- 0.6


# ============================================================================= #
# FUNCTION DEFINITIONS
# 01. write_csv
# 02. read_sav
# 03. analyse_data
# 04. class_distribution
# 05. missing_value_analysis
# 06. outliers_detection
# 07. outlier_analysis
# 08. split_data
# 09. null_feature_removal
# 10. variance_feature_removal
# 11. impute_outlier_3s
# 12. impute_outlier_NA
# 13. mean_imputation
# 14. KNN_imputation
# 15. mice_imputation
# 16. PCA_Analysis
# 17. selection_boruta
# ============================================================================= #

write_csv <- function(df, filename){
  # Function to write df as filename at the below directory
  file <- paste(dir, filename, sep = "/")
  write.csv(df, file)
}

read_sav <- function(file) {
  # Function to read data from source where file is path and filename
  read.spss(file, to.data.frame = TRUE)
}

analyse_data <- function(df) {
  # structure and dimension of df
  message('Structure of data:')
  str(df)
  message('Dimension of data:')
  dim(df)
}

class_distribution <- function(df, framename) {
  # Function to create barplot and visualize the distribution of class in the df
  message('Percentage of distribution of class variable:')
  tN <- table(df$class) * 100 / nrow(df)
  print(tN)
  t = paste('Class Distribution in', framename)
  barplot(tN, col = c('cyan', 'red'), ylim = c(0, 100), ylab = '% of observation', xlab = 'Class of Semiconductor', names.arg = c('Faulty', 'Non-Faulty'))
  title(main = t, font = 4)
}

missing_value_analysis <- function(df, framename) {
  # Function to analyse missing values in the df and visualize the missing values per feature
  message('count of total NAs: ', sum(is.na(df)))
  message('Total Number of observations: ', dim(df)[1])
  message('Total Number of observations with NA: ', sum(apply(df, 1, anyNA)))
  message('Total Number of features: ', dim(df)[2])
  message('Total Number of features with NA: ', sum(apply(df, 2, anyNA)))
  missing_per_col <- sapply(df, function(x) sum(is.na(x)) * 100 / dim(df)[1] )
  message('Summary of missing values per feature:')
  print(summary(missing_per_col))
  t <- paste("Missing Values in", framename)
  hist(missing_per_col, col = c('lightblue'), labels = TRUE,  ylim = c(0, 1000), main = t, xlab = '% of missing values per feature')
}

outliers_detection <- function(x) {
  # Function to detection outlier in each feature (x) using 3-s rule
  bv <- c(mean(x, na.rm = TRUE) - 3 * sd(x, na.rm = TRUE), mean(x, na.rm = TRUE) + 3 * sd(x, na.rm = TRUE))
  z <- (x < bv[1]) | (x > bv[2])
  sum(z, na.rm = TRUE)
}

outlier_analysis <- function(df, framename) {
  # Function to analyse outlier and visualize the outlier per feature
  outlier_per_col <- apply(df, 2, outliers_detection)
  # message('Outlier per Col', as.data.frame(outlier_per_col))
  message('Summary of outliers per feature:')
  print(summary(outlier_per_col))
  message('Total Number of observation: ', dim(df)[1])
  message('Total Number of features: ', dim(df)[2])
  message('Total Number of Outliers: ', sum(outlier_per_col))
  outlier_per_col <- outlier_per_col * 100 / dim(df)[1]
  t <- paste("Outliers in", framename)
  hist(outlier_per_col, col = c('indianred1'), labels = TRUE,  ylim = c(0, 1000), main = t, xlab = '% of outliers per feature')
}

split_data <- function(df, ratio) {
  # Function to split the data into train and test data set
  split <- sample(1:nrow(df), nrow(df) * ratio)
  train <- df[split, ]
  test <- df[-split, ]
  lst <- list('train_df' = train, 'test_df' = test)
  return(lst)
}

null_feature_removal <- function(df, null_percent) {
  # Function to remove features with more than a certain % of NAs
  df = df[, which(colMeans(!is.na(df)) >= null_percent)]
  return(df)
}

variance_feature_removal <- function(df) {
  # Function to remove features with zero variance and near zero variance. nearZeroVar removes
  # 1. one unique value i.e., zero variance features
  # 2. very few unique values relative to the number of observations 
  #		uniqueCut: the cutoff for the percentage of distinct values out of the number of total samples (Default is 10)
  # 3. Ratio of frequency of most common value to frequency of the second most common value is large
  #		freqCut: the cutoff for the ratio of the most common value to the second most common value (Default is 19)
  df <- df[-nearZeroVar(df)]
  return(df)
}

impute_outlier_3s <- function(x) {
  # Function to replace outliers with 3s boundary values
  bv <- c(mean(x, na.rm = TRUE) - 3 * sd(x, na.rm = TRUE), mean(x, na.rm = TRUE) + 3 * sd(x, na.rm = TRUE))
  x[x < bv[1]] <- bv[1]
  x[x > bv[2]] <- bv[2]
  return(x)
}

impute_outlier_NA <- function(x) {
  # Function to replace outliers with NAs or missing values
  bv <- c(mean(x, na.rm = TRUE) - 3 * sd(x, na.rm = TRUE), mean(x, na.rm = TRUE) + 3 * sd(x, na.rm = TRUE))
  x[x < bv[1]] <- NA
  x[x > bv[2]] <- NA
  return(x)
}

mean_imputation <- function(df) {
  # Function to impute missing values(NAs) by mean of the features
  for(i in 1:ncol(df)) {
    df[ , i][is.na(df[ , i])] <- mean(df[ , i], na.rm = TRUE)
  }
  return(df)
}

KNN_imputation <- function(df, train_df=NULL) {
  # Function to impute missing values(NAs) using K-Nearest Neighbour (k = 5)
  impute_KNN <- knnImputation(df, k = 5, scale = T, meth = "weighAvg", distData = train_df)
  return(impute_KNN)
}

selection_boruta <- function(df, class){
  # Function to perform feature selection using BORUTA
  boruta_df <- cbind(class, df)
  train_boruta_features <- Boruta(class ~ ., data = boruta_df, doTrace = 2, ntree = 500, maxRuns = 100)
  boruta_mF <- TentativeRoughFix(train_boruta_features)
  print("Confirmed formula after Boruta: ")
  message(getConfirmedFormula(train_boruta_features))
  plot(train_boruta_features, las = 2, cex.axis = 0.5)
  plotImpHistory(train_boruta_features)
  print(attStats(train_boruta_features))
  boruta_mF_confirmed <- names(boruta_mF$finalDecision[boruta_mF$finalDecision %in% c('Confirmed')])
  message("Number of Confirmed Features after Boruta: ", length(boruta_mF_confirmed))
  message(" Confirmed Features after Boruta: ", boruta_mF_confirmed)
  boruta_mF_tentative <- names(boruta_mF$finalDecision[boruta_mF$finalDecision %in% c('Tentative')])
  boruta_mF_rejected <- names(boruta_mF$finalDecision[boruta_mF$finalDecision %in% c('Rejected')])
  train_boruta <- df[,(names(df) %in% boruta_mF_confirmed)]
  return(train_boruta)
}


# ============================================================================= #
# READING SOURCE DATA & DATA ANALYSIS
# ============================================================================= #

semcon_original_data <- read_sav(SRC)
semcon_original_data$class[semcon_original_data$class == 1] <- "F"
semcon_original_data$class[semcon_original_data$class == 0] <- "NF"
analyse_data(semcon_original_data)
semcon_original_data <- semcon_original_data[, -c(1, 3)]
class_distribution(semcon_original_data, 'Given Sample')
missing_value_analysis(semcon_original_data, 'Given Sample')
outlier_analysis(semcon_original_data[, -c(1)], 'Given Sample')


# ============================================================================= #
# TRAIN AND TEST SPLIT
# ============================================================================= #

set.seed(666)
semcon_split_data <- split_data(semcon_original_data, TRAIN_TEST_RATIO)
semcon_train_data <- semcon_split_data[["train_df"]]
semcon_test_data <- semcon_split_data[["test_df"]]


# ============================================================================= #
# TRAIN DATA ANALYSIS
# ============================================================================= #

class_distribution(semcon_train_data, 'Train Data')
missing_value_analysis(semcon_train_data, 'train Data')
outlier_analysis(semcon_train_data[, -c(1)], 'Train Data')
# write_csv(semcon_train_data, "semcon_train_data.csv")


# ============================================================================= #
# TEST DATA ANALYSIS
# ============================================================================= #

class_distribution(semcon_test_data, 'Test Data')
missing_value_analysis(semcon_test_data, 'Test Data')
outlier_analysis(semcon_test_data[, -c(1)], 'Test Data')
rm(semcon_split_data)
# write_csv(semcon_test_data, "semcon_test_data.csv")


# ============================================================================= #
# FEATURE REDUCTION
# 01. Remove Features with more than 60% NA based on missing_value_analysis
# ============================================================================= #

class <- semcon_train_data$class
train_null_removal <- null_feature_removal(semcon_train_data[, -c(1)], NULL_THRESHOLD)
missing_value_analysis(train_null_removal, paste('Train Data after removing features with more than', NULL_THRESHOLD * 100, '% NA'))
outlier_analysis(train_null_removal, paste('train set after removing features with more than', NULL_THRESHOLD * 100, '% NA'))
# write_csv(cbind(class, train_null_removal), "train_60NA_Removal.csv")


# ============================================================================= #
# FEATURE REDUCTION
# 02. Near Zero Variance Removal
# ============================================================================= #

train_variance_removal <- variance_feature_removal(train_null_removal)
missing_value_analysis(train_variance_removal, 'train set after removing features with Near Zero Variance')
outlier_analysis(train_variance_removal, 'train set after removing features with Near Zero Variance')
# write_csv(cbind(class, train_variance_removal), "train_NZV_removal.csv")


# ============================================================================= #
# OUTLIER HANDLING
# 01. Replace by 3s boundary
# ============================================================================= #

train_outlier_3s <- apply(train_variance_removal, 2, impute_outlier_3s)
train_outlier_3s <- as.data.frame(train_outlier_3s)
missing_value_analysis(train_outlier_3s, 'train set after imputating outlier with 3s boundary')
outlier_analysis(train_outlier_3s, 'train set after imputating outlier with 3s boundary')
# write_csv(cbind(class, train_outlier_3s), "train_OTH_3s.csv")


# ============================================================================= #
# OUTLIER HANDLING
# 01. Replace by NAs
# ============================================================================= #

train_outlier_NA <- apply(train_variance_removal, 2, impute_outlier_NA)
train_outlier_NA <- as.data.frame(train_outlier_NA)
missing_value_analysis(train_outlier_NA, 'train set after imputating outlier with NA')
outlier_analysis(train_outlier_NA, 'train set after imputating outlier with NA')
# write_csv(cbind(class, train_outlier_NA), "train_OTH_NA.csv")


# ============================================================================= #
# NA HANDLING
# 01. mean imputation
# ============================================================================= #

# train_mean_imputation <- mean_imputation(train_outlier_NA)
# missing_value_analysis(train_mean_imputation, 'train set after mean imputation')
# outlier_analysis(train_mean_imputation, 'train set after mean imputation')
# write_csv(cbind(class, train_mean_imputation), "train_NAH_mean.csv")


# ============================================================================= #
# NA HANDLING
# 02. KNN imputation
# ============================================================================= #

train_knn_imputation <- KNN_imputation(train_outlier_NA)
missing_value_analysis(train_knn_imputation, 'train set after KNN imputation')
# outlier_analysis(train_knn_imputation, 'train set after KNN imputation')
# write_csv(cbind(class, train_knn_imputation), "train_NAH_knn.csv")


# ============================================================================= #
# FEATURE SELECTION AND REDUCTION
# 01. BORUTA on KNN imputed train set
# ============================================================================= #

set.seed(111)
train_FR_boruta_knn <- selection_boruta(train_knn_imputation, class)
print("Summary of selected features")
summary(train_FR_boruta_knn)
message("Variance corresponding to selected features:")
sapply(train_FR_boruta_knn, var)
missing_value_analysis(train_FR_boruta_knn, 'train set after BORUTA')
# outlier_analysis(train_FR_boruta_knn, 'train set after BORUTA')
# write_csv(cbind(train_FR_boruta_knn), "train_FR_boruta.csv")

# ============================================================================= #
# FEATURE SELECTION AND REDUCTION
# 02. Recursive Feature Elimination (RFE)
# ============================================================================= #

# set.seed(199)
# options(warn = -1)
# subsets <- c(15:25)
# ctrl <- rfeControl(functions = rfFuncs, method = 'repeatedcv', repeats = 5, verbose = TRUE)
# norm <- preProcess(train_knn_imputation)
# x <- predict(norm, train_knn_imputation)
# lmProfile <- rfe(x = x, y = as.factor(class), sizes = subsets, rfeControl = ctrl)
# lmProfile

# train <- cbind(class, train_knn_imputation[ , which(names(train_knn_imputation) %in% c(lmProfile[["optVariables"]][1:20]))])
# test <- semcon_test_data[ , which(names(semcon_test_data) %in% c(names(train)))]


# ============================================================================= #
# FINAL TRAIN DATASET
# 01. Imbalanced
# 02. Under-sampled
# 03. Over-sampled
# 04. Mixed-sampled
# 05. ROSE
# 06. SMOTE
# ============================================================================= #

imba_train_data <- cbind(class, train_FR_boruta_knn)
print(table(imba_train_data$class))
# F: 71, NF: 1182
under_train_data <- ovun.sample(class ~., data = imba_train_data, method = 'under', seed = 222, N = 142)$data
print(table(under_train_data$class))
# F: 71, NF: 71
over_train_data <- ovun.sample(class ~., data = imba_train_data, method = 'over', seed = 222, N = 2364)$data
print(table(over_train_data$class))
# F: 1182, NF: 1182
ovun_train_data <- ovun.sample(class ~., data = imba_train_data, method = 'both', p = 0.5, seed = 222, N = 1253)$data
print(table(ovun_train_data$class))
# F: 623, NF: 630
rose_train_data <- ROSE(class ~ ., data = imba_train_data, N = dim(imba_train_data)[1] * 3, p = 0.5, seed = 1)$data
print(table(rose_train_data$class))
# F: 1935, NF: 1824
smote_train_data <- SMOTE(class ~., data = imba_train_data, perc.over = 2000, perc.under = 100)
print(table(smote_train_data$class))
# F: 1491, NF: 1420


# ============================================================================= #
# PREPROCESS: TEST DATASET
# 01. Same number of features as Train after variance removal (Mirroring the features)
# 02. Outlier detection and imputation
# 03. KNN Imputation using train
# ============================================================================= #

class <- semcon_test_data$class
test_variance_removal <- semcon_test_data[ , which(names(semcon_test_data) %in% c(names(train_variance_removal)))]
test_outlier_NA <- apply(test_variance_removal, 2, impute_outlier_NA)
test_outlier_NA <- as.data.frame(test_outlier_NA)
test_knn_imputation <- KNN_imputation(test_outlier_NA, train_knn_imputation)
test <- cbind(class, test_knn_imputation)
class_distribution(test, 'Test')


# ============================================================================= #
# SAMPLING TECHNIQUES
# 01. CV
# 02. Bootstrap
# 03. K-Fold CV
# ============================================================================= #

ctrl_noTune <- trainControl(method = 'none', classProbs = TRUE, summaryFunction = twoClassSummary)
ctrl_tune <- trainControl(method = 'boot', number = 20, classProbs = TRUE, summaryFunction = twoClassSummary)
# ctrl <- trainControl(method = "cv", number = 10, savePredictions = 'final',  summaryFunction = twoClassSummary)
# ctrl <- trainControl(method = "boot632", number = 1000, savePredictions = TRUE, savePredictions = 'final', classProbs = T, summaryFunction = twoClassSummary)
# ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE, savePredictions = 'final', classProbs = T, summaryFunction = twoClassSummary)
# ctrl_under <- trainControl(method = 'repeatedcv', number = 10, verboseIter = FALSE, sampling = 'down')
# ctrl_over <- trainControl(method = 'repeatedcv', number = 10, verboseIter = FALSE, sampling = 'up')
# ctrl_rose <- trainControl(method = 'repeatedcv', number = 10, verboseIter = FALSE, sampling = 'rose')
# ctrl_smote <- trainControl(method = 'repeatedcv', number = 10, verboseIter = FALSE, sampling = 'smote')


# # ============================================================================= #
# # CLASSIFICATION MODEL WITHOUT TUNE
# # 01. Generalized Linear Model
# #	A. Imbalanced Dataset
# #	B. Undersampling
# #	C. Oversampling
# #	D. Mixedsampling
# #	E. ROSE
# #	F. SMOTE
# # ============================================================================= #

# set.seed(642)
# # Model Creation
# model_glm_imba <- train(form = class ~ ., data = imba_train_data, family = binomial(link = 'logit'), trControl = ctrl_noTune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# #exp(coef(model_glm_imba$finalModel))
# #varImp(model_glm_imba)
# # Train Evaluation
# # prob_glm_class_imba = predict(model_glm_imba, newdata = imba_train_data, type = 'prob')
# # pred_glm_class_imba = predict(model_glm_imba, newdata = imba_train_data)
# # confusionMatrix(as.factor(pred_glm_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_glm_imba <- roc.curve(imba_train_data[['class']], pred_glm_class_imba, plotit = T, main = 'ROC Curve using GLM')
# # print(roc_glm_imba$auc)
# # Test Evaluation
# prob_glm_class_imba = predict(model_glm_imba, newdata = test, type = 'prob')
# pred_glm_class_imba = predict(model_glm_imba, newdata = test)
# confusionMatrix(as.factor(pred_glm_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_glm_imba <- roc.curve(test[['class']], pred_glm_class_imba, plotit = T, main = 'ROC Curve using GLM')
# print(roc_glm_imba$auc)

# set.seed(642)
# # Model Creation
# model_glm_under <- train(form = class ~ ., data = under_train_data, family = binomial(link = 'logit'), trControl = ctrl_noTune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# # Train Evaluation
# # prob_glm_class_under = predict(model_glm_under, newdata = under_train_data, type = 'prob')
# # pred_glm_class_under = predict(model_glm_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_glm_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_glm_under <- roc.curve(under_train_data[['class']], pred_glm_class_under, plotit = T, main = 'ROC Curve using GLM & Undersmaple')
# # print(roc_glm_under$auc)
# # Test Evaluation
# prob_glm_class_under = predict(model_glm_under, newdata = test, type = 'prob')
# pred_glm_class_under = predict(model_glm_under, newdata = test)
# confusionMatrix(as.factor(pred_glm_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_glm_under <- roc.curve(test[['class']], pred_glm_class_under, plotit = T, main = 'ROC Curve using GLM & Undersmaple')
# print(roc_glm_under$auc)

# set.seed(642)
# # Model Creation
# model_glm_over <- train(form = class ~ ., data = over_train_data, family = binomial(link = 'logit'), trControl = ctrl_noTune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# # Train Evaluation
# # prob_glm_class_over = predict(model_glm_over, newdata = over_train_data, type = 'prob')
# # pred_glm_class_over = predict(model_glm_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_glm_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_glm_over <- roc.curve(over_train_data[['class']], pred_glm_class_over, plotit = T, main = 'ROC Curve using GLM & Oversmaple')
# # print(roc_glm_over$auc)
# # Test Evaluation
# prob_glm_class_over = predict(model_glm_over, newdata = test, type = 'prob')
# pred_glm_class_over = predict(model_glm_over, newdata = test)
# confusionMatrix(as.factor(pred_glm_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_glm_over <- roc.curve(test[['class']], pred_glm_class_over, plotit = T, main = 'ROC Curve using GLM & Oversmaple')
# print(roc_glm_over$auc)

# set.seed(642)
# # Model Creation
# model_glm_ovun <- train(form = class ~ ., data = ovun_train_data, family = binomial(link = 'logit'), trControl = ctrl_noTune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# # Train Evaluation
# # prob_glm_class_ovun = predict(model_glm_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_glm_class_ovun = predict(model_glm_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_glm_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_glm_under <- roc.curve(ovun_train_data[['class']], pred_glm_class_ovun, plotit = T, main = 'ROC Curve using GLM & Miexed-smaple')
# # print(roc_glm_under$auc)
# # Test Evaluation
# prob_glm_class_ovun = predict(model_glm_ovun, newdata = test, type = 'prob')
# pred_glm_class_ovun = predict(model_glm_ovun, newdata = test)
# confusionMatrix(as.factor(pred_glm_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_glm_ovun <- roc.curve(test[['class']], pred_glm_class_ovun, plotit = T, main = 'ROC Curve using GLM & Miexed-smaple')
# print(roc_glm_ovun$auc)

# set.seed(642)
# # Model Creation
# model_glm_rose <- train(form = class ~ ., data = rose_train_data, family = binomial(link = 'logit'), trControl = ctrl_noTune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# # Train Evaluation
# # prob_glm_class_rose = predict(model_glm_rose, newdata = rose_train_data, type = 'prob')
# # pred_glm_class_rose = predict(model_glm_rose, newdata = rose_train_data)
# # confusionMatrix(as.factor(pred_glm_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_glm_rose <- roc.curve(rose_train_data[['class']], pred_glm_class_rose, plotit = T, main = 'ROC Curve using GLM & ROSE')
# # print(roc_glm_rose$auc)
# # Test Evaluation
# prob_glm_class_rose = predict(model_glm_rose, newdata = test, type = 'prob')
# pred_glm_class_rose = predict(model_glm_rose, newdata = test)
# confusionMatrix(as.factor(pred_glm_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_glm_rose <- roc.curve(test[['class']], pred_glm_class_rose, plotit = T, main = 'ROC Curve using GLM & ROSE')
# print(roc_glm_rose$auc)

# set.seed(642)
# # Model Creation
# model_glm_smote <- train(form = class ~ ., data = smote_train_data, family = binomial(link = 'logit'), trControl = ctrl_noTune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# # Train Evaluation
# # prob_glm_class_smote = predict(model_glm_smote, newdata = smote_train_data, type = 'prob')
# # pred_glm_class_smote = predict(model_glm_smote, newdata = smote_train_data)
# # confusionMatrix(as.factor(pred_glm_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_glm_smote <- roc.curve(smote_train_data[['class']], pred_glm_class_smote, plotit = T, main = 'ROC Curve using GLM & SMOTE')
# # print(roc_glm_smote$auc)
# # Test Evaluation
# prob_glm_class_smote = predict(model_glm_smote, newdata = test, type = 'prob')
# pred_glm_class_smote = predict(model_glm_smote, newdata = test)
# confusionMatrix(as.factor(pred_glm_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_glm_smote <- roc.curve(test[['class']], pred_glm_class_smote, plotit = T, main = 'ROC Curve using GLM & SMOTE')
# print(roc_glm_smote$auc)


# # ============================================================================= #
# # CLASSIFICATION MODEL
# # 02. Decision Tree - ctree, chaid, C5.0, xgbTree
# #	A. Imbalanced Dataset
# #	B. Undersampling
# #	C. Oversampling
# #	D. Mixedsampling
# #	E. ROSE
# #	F. SMOTE
# # ============================================================================= #

# set.seed(642)
# # Model Creation
# model_dt_imba <- train(form = class ~ ., data = imba_train_data, method = 'ctree', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_dt_imba)
# # Train Evaluation
# # prob_dt_class_imba <- predict(model_dt_imba, newdata = imba_train_data, type = 'prob')
# # pred_dt_class_imba = predict(model_dt_imba, newdata = imba_train_data)
# # confusionMatrix(as.factor(pred_dt_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# # # Error in ROC Calculation
# # roc_dt_imba <- roc.curve(imba_train_data[['class']], pred_dt_class_imba, plotit = T, main = 'ROC Curve using DT')
# # print(roc_dt_imba$auc)
# # Test Evaluation
# prob_dt_class_imba <- predict(model_dt_imba, newdata = test, type = 'prob')
# pred_dt_class_imba = predict(model_dt_imba, newdata = test)
# confusionMatrix(as.factor(pred_dt_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_dt_imba <- roc.curve(test[['class']], pred_dt_class_imba, plotit = T, main = 'ROC Curve using DT')
# print(roc_dt_imba$auc)

# set.seed(642)
# # Model Creation
# model_dt_under <- train(form = class ~ ., data = under_train_data, method = 'ctree', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_dt_under)
# # Train Evaluation
# # prob_dt_class_under <- predict(model_dt_under, newdata = under_train_data, type = 'prob')
# # pred_dt_class_under = predict(model_dt_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_dt_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_dt_under <- roc.curve(under_train_data[['class']], pred_dt_class_under, plotit = T, main = 'ROC Curve using DT & undersample')
# # print(roc_dt_under$auc)
# # Test Evaluation
# prob_dt_class_under <- predict(model_dt_under, newdata = test, type = 'prob')
# pred_dt_class_under = predict(model_dt_under, newdata = test)
# confusionMatrix(as.factor(pred_dt_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_dt_under <- roc.curve(test[['class']], pred_dt_class_under, plotit = T, main = 'ROC Curve using DT & undersample')
# print(roc_dt_under$auc)

# set.seed(642)
# # Model Creation
# model_dt_over <- train(form = class ~ ., data = over_train_data, method = 'ctree', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_dt_over)
# # Train Evaluation
# # prob_dt_class_over <- predict(model_dt_over, newdata = over_train_data, type = 'prob')
# # pred_dt_class_over = predict(model_dt_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_dt_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_dt_over <- roc.curve(over_train_data[['class']], pred_dt_class_over, plotit = T, main = 'ROC Curve using DT & oversample')
# # print(roc_dt_over$auc)
# # Test Evaluation
# prob_dt_class_over <- predict(model_dt_over, newdata = test, type = 'prob')
# pred_dt_class_over = predict(model_dt_over, newdata = test)
# confusionMatrix(as.factor(pred_dt_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_dt_over <- roc.curve(test[['class']], pred_dt_class_over, plotit = T, main = 'ROC Curve using DT & oversampling')
# print(roc_dt_over$auc)

# set.seed(642)
# # Model Creation
# model_dt_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'ctree', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_dt_ovun)
# # Train Evaluation
# # prob_dt_class_ovun <- predict(model_dt_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_dt_class_ovun = predict(model_dt_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_dt_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_dt_ovun <- roc.curve(ovun_train_data[['class']], pred_dt_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_dt_ovun$auc)
# # Test Evaluation
# prob_dt_class_ovun <- predict(model_dt_ovun, newdata = test, type = 'prob')
# pred_dt_class_ovun = predict(model_dt_ovun, newdata = test)
# confusionMatrix(as.factor(pred_dt_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_dt_ovun <- roc.curve(test[['class']], pred_dt_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_dt_ovun$auc)

# set.seed(642)
# # Model Creation
# model_dt_rose <- train(form = class ~ ., data = rose_train_data, method = 'ctree', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_dt_rose)
# # Train Evaluation
# # prob_dt_class_rose <- predict(model_dt_rose, newdata = rose_train_data, type = 'prob')
# # pred_dt_class_rose = predict(model_dt_rose, newdata = rose_train_data)
# # confusionMatrix(as.factor(pred_dt_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_dt_rose <- roc.curve(rose_train_data[['class']], pred_dt_class_rose, plotit = T, main = 'ROC Curve using DT & ROSE')
# # print(roc_dt_rose$auc)
# # Test Evaluation
# prob_dt_class_rose <- predict(model_dt_rose, newdata = test, type = 'prob')
# pred_dt_class_rose = predict(model_dt_rose, newdata = test)
# confusionMatrix(as.factor(pred_dt_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_dt_rose <- roc.curve(test[['class']], pred_dt_class_rose, plotit = T, main = 'ROC Curve using DT & ROSE')
# print(roc_dt_rose$auc)

# set.seed(642)
# # Model Creation
# model_dt_smote <- train(form = class ~ ., data = smote_train_data, method = 'ctree', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_dt_smote)
# # Train Evaluation
# # prob_dt_class_smote <- predict(model_dt_smote, newdata = smote_train_data, type = 'prob')
# # pred_dt_class_smote = predict(model_dt_smote, newdata = smote_train_data)
# # confusionMatrix(as.factor(pred_dt_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_dt_smote <- roc.curve(smote_train_data[['class']], pred_dt_class_smote, plotit = T, main = 'ROC Curve using DT & SMOTE')
# # print(roc_dt_smote$auc)
# # Test Evaluation
# prob_dt_class_smote <- predict(model_dt_smote, newdata = test, type = 'prob')
# pred_dt_class_smote = predict(model_dt_smote, newdata = test)
# confusionMatrix(as.factor(pred_dt_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_dt_smote <- roc.curve(test[['class']], pred_dt_class_smote, plotit = T, main = 'ROC Curve using DT & SMOTE')
# print(roc_dt_smote$auc)


# # ============================================================================= #
# # CLASSIFICATION MODEL
# # 03. Random Forest - ranger, rf
# #	A. Imbalanced Dataset
# #	B. Undersampling
# #	C. Oversampling
# #	D. Mixedsampling
# #	E. ROSE
# #	F. SMOTE
# # ============================================================================= #

# set.seed(642)
# # Model Creation
# model_rf_imba <- train(form = class ~ ., data = imba_train_data, method = 'rf', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_rf_imba)
# # Train Evaluation
# # # prob_rf_class_imba <- predict(model_rf_imba, newdata = imba_train_data, type = 'prob')
# # pred_rf_class_imba = predict(model_rf_imba, newdata = imba_train_data)
# # confusionMatrix(as.factor(pred_rf_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_rf_imba <- roc.curve(imba_train_data[['class']], pred_rf_class_imba, plotit = T, main = 'ROC Curve using RF')
# # print(roc_rf_imba$auc)
# # Test Evaluation
# pred_rf_class_imba = predict(model_rf_imba, newdata = test)
# confusionMatrix(as.factor(pred_rf_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_rf_imba <- roc.curve(test[['class']], pred_rf_class_imba, plotit = T, main = 'ROC Curve using RF')
# print(roc_rf_imba$auc)

# set.seed(642)
# # Model Creation
# model_rf_under <- train(form = class ~ ., data = under_train_data, method = 'rf', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_rf_under)
# # Train Evaluation
# # pred_rf_class_under = predict(model_rf_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_rf_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_rf_under <- roc.curve(under_train_data[['class']], pred_rf_class_under, plotit = T, main = 'ROC Curve using RF & undersample')
# # print(roc_rf_under$auc)
# # Test Evaluation
# pred_rf_class_under = predict(model_rf_under, newdata = test)
# confusionMatrix(as.factor(pred_rf_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_rf_under <- roc.curve(test[['class']], pred_rf_class_under, plotit = T, main = 'ROC Curve using RF & undersample')
# print(roc_rf_under$auc)

# set.seed(642)
# # Model Creation
# model_rf_over <- train(form = class ~ ., data = over_train_data, method = 'rf', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_rf_over)
# # Train Evaluation
# # pred_rf_class_over = predict(model_rf_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_rf_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_rf_over <- roc.curve(over_train_data[['class']], pred_rf_class_over, plotit = T, main = 'ROC Curve using RF & oversample')
# # print(roc_rf_over$auc)
# # Test Evaluation
# pred_rf_class_over = predict(model_rf_over, newdata = test)
# confusionMatrix(as.factor(pred_rf_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_rf_over <- roc.curve(test[['class']], pred_rf_class_over, plotit = T, main = 'ROC Curve using RF & oversample')
# print(roc_rf_over$auc)

# set.seed(642)
# # Model Creation
# model_rf_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'rf', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_rf_ovun)
# # Train Evaluation
# # pred_rf_class_ovun = predict(model_rf_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_rf_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_rf_ovun <- roc.curve(ovun_train_data[['class']], pred_rf_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_rf_ovun$auc)
# # Test Evaluation
# pred_rf_class_ovun = predict(model_rf_ovun, newdata = test)
# confusionMatrix(as.factor(pred_rf_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_rf_ovun <- roc.curve(test[['class']], pred_rf_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_rf_ovun$auc)

# set.seed(642)
# # Model Creation
# model_rf_rose <- train(form = class ~ ., data = rose_train_data, method = 'rf', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_rf_rose)
# # Train Evaluation
# # pred_rf_class_rose = predict(model_rf_rose, newdata = rose_train_data)
# # confusionMatrix(as.factor(pred_rf_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_rf_rose <- roc.curve(rose_train_data[['class']], pred_rf_class_rose, plotit = T, main = 'ROC Curve using RF & ROSE')
# # print(roc_rf_rose$auc)
# # Test Evaluation
# pred_rf_class_rose = predict(model_rf_rose, newdata = test)
# confusionMatrix(as.factor(pred_rf_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_rf_rose <- roc.curve(test[['class']], pred_rf_class_rose, plotit = T, main = 'ROC Curve using RF & ROSE')
# print(roc_rf_rose$auc)

# set.seed(642)
# # Model Creation
# model_rf_smote <- train(form = class ~ ., data = smote_train_data, method = 'rf', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_rf_smote)
# # # Train Evaluation
# # pred_rf_class_smote = predict(model_rf_smote, newdata = smote_train_data)
# # confusionMatrix(as.factor(pred_rf_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_rf_smote <- roc.curve(smote_train_data[['class']], pred_rf_class_smote, plotit = T, main = 'ROC Curve using RF & SMOTE')
# # print(roc_rf_smote$auc)
# # Test Evaluation
# pred_rf_class_smote = predict(model_rf_smote, newdata = test)
# confusionMatrix(as.factor(pred_rf_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_rf_smote <- roc.curve(test[['class']], pred_rf_class_smote, plotit = T, main = 'ROC Curve using RF & SMOTE')
# print(roc_rf_smote$auc)


# # ============================================================================= #
# # CLASSIFICATION MODEL
# # 04. KNN
# #	A. Imbalanced Dataset
# #	B. Undersampling
# #	C. Oversampling
# #	D. Mixedsampling
# #	E. ROSE
# #	F. SMOTE
# # ============================================================================= #

# set.seed(642)
# # Model Creation
# model_knn_imba <- train(form = class ~ ., data = imba_train_data, method = 'knn', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_knn_imba)
# # Train Evaluation
# # prob_knn_class_imba <- predict(model_knn_imba, newdata = imba_train_data, type = 'prob')
# # pred_knn_class_imba = predict(model_knn_imba, newdata = imba_train_data)
# # confusionMatrix(as.factor(pred_knn_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_knn_imba <- roc.curve(imba_train_data[['class']], pred_knn_class_imba, plotit = T, main = 'ROC Curve using KNN')
# # print(roc_knn_imba$auc)
# # Test Evaluation
# prob_knn_class_imba <- predict(model_knn_imba, newdata = test, type = 'prob')
# pred_knn_class_imba = predict(model_knn_imba, newdata = test)
# confusionMatrix(as.factor(pred_knn_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_knn_imba <- roc.curve(test[['class']], pred_knn_class_imba, plotit = T, main = 'ROC Curve using KNN')
# print(roc_knn_imba$auc)

# set.seed(642)
# # Model Creation
# model_knn_under <- train(form = class ~ ., data = under_train_data, method = 'knn', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_knn_under)
# # Train Evaluation
# # prob_knn_class_under <- predict(model_knn_under, newdata = under_train_data, type = 'prob')
# # pred_knn_class_under = predict(model_knn_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_knn_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_knn_under <- roc.curve(under_train_data[['class']], pred_knn_class_under, plotit = T, main = 'ROC Curve using KNN and undersample')
# # print(roc_knn_under$auc)
# # Test Evaluation
# prob_knn_class_under <- predict(model_knn_under, newdata = test, type = 'prob')
# pred_knn_class_under = predict(model_knn_under, newdata = test)
# confusionMatrix(as.factor(pred_knn_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_knn_under <- roc.curve(test[['class']], pred_knn_class_under, plotit = T, main = 'ROC Curve using KNN and undersample')
# print(roc_knn_under$auc)

# set.seed(642)
# # Model Creation
# model_knn_over <- train(form = class ~ ., data = over_train_data, method = 'knn', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_knn_over)
# # Train Evaluation
# # prob_knn_class_over <- predict(model_knn_over, newdata = over_train_data, type = 'prob')
# # pred_knn_class_over = predict(model_knn_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_knn_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_knn_over <- roc.curve(pover_train_data[['class']], red_knn_class_over, plotit = T, main = 'ROC Curve using KNN and oversample')
# # print(roc_knn_over$auc)
# # Test Evaluation
# prob_knn_class_over <- predict(model_knn_over, newdata = test, type = 'prob')
# pred_knn_class_over = predict(model_knn_over, newdata = test)
# confusionMatrix(as.factor(pred_knn_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_knn_over <- roc.curve(test[['class']], pred_knn_class_over, plotit = T, main = 'ROC Curve using KNN and oversample')
# print(roc_knn_over$auc)

# set.seed(642)
# # Model Creation
# model_knn_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'knn', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_knn_ovun)
# # Train Evaluation
# # prob_knn_class_ovun <- predict(model_knn_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_knn_class_ovun = predict(model_knn_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_knn_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_knn_ovun <- roc.curve(ovun_train_data[['class']], pred_knn_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_knn_ovun$auc)
# # Test Evaluation
# prob_knn_class_ovun <- predict(model_knn_ovun, newdata = test, type = 'prob')
# pred_knn_class_ovun = predict(model_knn_ovun, newdata = test)
# confusionMatrix(as.factor(pred_knn_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_knn_ovun <- roc.curve(test[['class']], pred_knn_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_knn_ovun$auc)

# set.seed(642)
# # Model Creation
# model_knn_rose <- train(form = class ~ ., data = rose_train_data, method = 'knn', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_knn_rose)
# # Train Evaluation
# # prob_knn_class_rose <- predict(model_knn_rose, newdata = rose_train_data, type = 'prob')
# # pred_knn_class_rose = predict(model_knn_rose, newdata = rose_train_data)
# # confusionMatrix(as.factor(pred_knn_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_knn_rose <- roc.curve(rose_train_data[['class']], pred_knn_class_rose, plotit = T, main = 'ROC Curve using KNN and ROSE')
# # print(roc_knn_rose$auc)
# # Test Evaluation
# prob_knn_class_rose <- predict(model_knn_rose, newdata = test, type = 'prob')
# pred_knn_class_rose = predict(model_knn_rose, newdata = test)
# confusionMatrix(as.factor(pred_knn_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_knn_rose <- roc.curve(test[['class']], pred_knn_class_rose, plotit = T, main = 'ROC Curve using KNN and ROSE')
# print(roc_knn_rose$auc)

# set.seed(642)
# # Model Creation
# model_knn_smote <- train(form = class ~ ., data = smote_train_data, method = 'knn', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_knn_smote)
# # Train Evaluation
# # prob_knn_class_smote <- predict(model_knn_smote, newdata = smote_train_data, type = 'prob')
# # pred_knn_class_smote = predict(model_knn_smote, newdata = smote_train_data)
# # confusionMatrix(as.factor(pred_knn_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_knn_smote <- roc.curve(smote_train_data[['class']], pred_knn_class_smote, plotit = T, main = 'ROC Curve using KNN and SMOTE')
# # print(roc_knn_smote$auc)
# # Test Evaluation
# prob_knn_class_smote <- predict(model_knn_smote, newdata = test, type = 'prob')
# pred_knn_class_smote = predict(model_knn_smote, newdata = test)
# confusionMatrix(as.factor(pred_knn_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_knn_smote <- roc.curve(ptest[['class']], red_knn_class_smote, plotit = T, main = 'ROC Curve using KNN and SMOTE')
# print(roc_knn_smote$auc)


# # ============================================================================= #
# # CLASSIFICATION MODEL
# # 05. Naive Bayes
# #	A. Imbalanced Dataset
# #	B. Undersampling
# #	C. Oversampling
# #	D. Mixedsampling
# #	E. ROSE
# #	F. SMOTE
# # ============================================================================= #

# set.seed(642)
# # Model Creation
# model_nb_imba <- train(form = class ~ ., data = imba_train_data, method = 'naive_bayes', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_nb_imba)
# # Train Evaluation
# # prob_nb_class_imba <- predict(model_nb_imba, newdata = imba_train_data, type = 'prob')
# # pred_nb_class_imba = predict(model_nb_imba, newdata = imba_train_data)
# # confusionMatrix(as.factor(pred_nb_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# # roc_nb_imba <- roc.curve(imba_train_data[['class']], pred_nb_class_imba, plotit = T, main = 'ROC Curve using NB')
# # print(roc_nb_imba$auc)
# # Test Evaluation
# prob_nb_class_imba <- predict(model_nb_imba, newdata = test, type = 'prob')
# pred_nb_class_imba = predict(model_nb_imba, newdata = test)
# confusionMatrix(as.factor(pred_nb_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_imba <- roc.curve(test[['class']], pred_nb_class_imba, plotit = T, main = 'ROC Curve using NB')
# print(roc_nb_imba$auc)

# set.seed(642)
# # Model Creation
# model_nb_under <- train(form = class ~ ., data = under_train_data, method = 'naive_bayes', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_nb_under)
# # Train Evaluation
# # prob_nb_class_under <- predict(model_nb_under, newdata = under_train_data, type = 'prob')
# # pred_nb_class_under = predict(model_nb_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_nb_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_nb_under <- roc.curve(under_train_data[['class']], pred_nb_class_under, plotit = T, main = 'ROC Curve using NB & undersample')
# # print(roc_nb_under$auc)
# # Test Evaluation
# prob_nb_class_under <- predict(model_nb_under, newdata = test, type = 'prob')
# pred_nb_class_under = predict(model_nb_under, newdata = test)
# confusionMatrix(as.factor(pred_nb_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_under <- roc.curve(test[['class']], pred_nb_class_under, plotit = T, main = 'ROC Curve using NB & undersample')
# print(roc_nb_under)$auc

# set.seed(642)
# # Model Creation
# model_nb_over <- train(form = class ~ ., data = over_train_data, method = 'naive_bayes', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_nb_over)
# # Train Evaluation
# # prob_nb_class_over <- predict(model_nb_over, newdata = over_train_data, type = 'prob')
# # pred_nb_class_over = predict(model_nb_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_nb_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_nb_over <- roc.curve(over_train_data[['class']], pred_nb_class_over, plotit = T, main = 'ROC Curve using NB & oversample')
# # print(roc_nb_over$auc)
# # Test Evaluation
# prob_nb_class_over <- predict(model_nb_over, newdata = test, type = 'prob')
# pred_nb_class_over = predict(model_nb_over, newdata = test)
# confusionMatrix(as.factor(pred_nb_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_over <- roc.curve(test[['class']], pred_nb_class_over, plotit = T, main = 'ROC Curve using NB & oversample')
# print(roc_nb_over$auc)

# set.seed(642)
# # Model Creation
# model_nb_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'naive_bayes', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_nb_ovun)
# # Train Evaluation
# # prob_nb_class_ovun <- predict(model_nb_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_nb_class_ovun = predict(model_nb_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_nb_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_nb_ovun <- roc.curve(ovun_train_data[['class']], pred_nb_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_nb_ovun$auc)
# # Test Evaluation
# prob_nb_class_ovun <- predict(model_nb_ovun, newdata = test, type = 'prob')
# pred_nb_class_ovun = predict(model_nb_ovun, newdata = test)
# confusionMatrix(as.factor(pred_nb_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_ovun <- roc.curve(test[['class']], pred_nb_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_nb_ovun$auc)

# set.seed(642)
# # Model Creation
# model_nb_rose <- train(form = class ~ ., data = rose_train_data, method = 'naive_bayes', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_nb_rose)
# # Train Evaluation
# # prob_nb_class_rose <- predict(model_nb_rose, newdata = rose_train_data, type = 'prob')
# # pred_nb_class_rose = predict(model_nb_rose, newdata = rose_train_data)
# # confusionMatrix(as.factor(pred_nb_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_nb_rose <- roc.curve(rose_train_data[['class']], pred_nb_class_rose, plotit = T, main = 'ROC Curve using NB & ROSE')
# # print(roc_nb_rose$auc)
# # Test Evaluation
# prob_nb_class_rose <- predict(model_nb_rose, newdata = test, type = 'prob')
# pred_nb_class_rose = predict(model_nb_rose, newdata = test)
# confusionMatrix(as.factor(pred_nb_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_rose <- roc.curve(test[['class']], pred_nb_class_rose, plotit = T, main = 'ROC Curve using NB & ROSE')
# print(roc_nb_rose$auc)

# set.seed(642)
# # Model Creation
# model_nb_smote <- train(form = class ~ ., data = smote_train_data, method = 'naive_bayes', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_nb_smote)
# # Train Evaluation
# # prob_nb_class_smote <- predict(model_nb_smote, newdata = smote_train_data, type = 'prob')
# # pred_nb_class_smote = predict(model_nb_smote, newdata = smote_train_data)
# # confusionMatrix(as.factor(pred_nb_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_nb_smote <- roc.curve(smote_train_data[['class']], pred_nb_class_smote, plotit = T, main = 'ROC Curve using NB & SMOTE')
# # print(roc_nb_smote$auc)
# # Test Evaluation
# prob_nb_class_smote <- predict(model_nb_smote, newdata = test, type = 'prob')
# pred_nb_class_smote = predict(model_nb_smote, newdata = test)
# confusionMatrix(as.factor(pred_nb_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_smote <- roc.curve(test[['class']], pred_nb_class_smote, plotit = T, main = 'ROC Curve using NB & SMOTE')
# print(roc_nb_smote$auc)


# # ============================================================================= #
# # CLASSIFICATION MODEL
# # 06. SVM
# #	A. Imbalanced Dataset
# #	B. Undersampling
# #	C. Oversampling
# #	D. Mixedsampling
# #	E. ROSE
# #	F. SMOTE
# # ============================================================================= #

# set.seed(642)
# # Model Creation
# model_svm_imba <- train(form = class ~ ., data = imba_train_data, method = 'svmRadial', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_svm_imba)
# # Train Evaluation
# # prob_svm_class_imba <- predict(model_svm_imba, newdata = imba_train_data, type = 'prob')
# # pred_svm_class_imba = predict(model_svm_imba, newdata = imba_train_data)
# # confusionMatrix(as.factor(pred_svm_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_svm_imba <- roc.curve(imba_train_data[['class']], pred_svm_class_imba, plotit = T, main = 'ROC Curve using SVM')
# # print(roc_svm_imba$auc)
# # Test Evaluation
# prob_svm_class_imba <- predict(model_svm_imba, newdata = test, type = 'prob')
# pred_svm_class_imba = predict(model_svm_imba, newdata = test)
# confusionMatrix(as.factor(pred_svm_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_svm_imba <- roc.curve(test[['class']], pred_svm_class_imba, plotit = T, main = 'ROC Curve using SVM')
# print(roc_svm_imba$auc)

# set.seed(642)
# # Model Creation
# model_svm_under <- train(form = class ~ ., data = under_train_data, method = 'svmRadial', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_svm_under)
# # Train Evaluation
# # prob_svm_class_under <- predict(model_svm_under, newdata = under_train_data, type = 'prob')
# # pred_svm_class_under = predict(model_svm_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_svm_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_svm_under <- roc.curve(under_train_data[['class']], pred_svm_class_under, plotit = T, main = 'ROC Curve using SVM & undersample')
# # print(roc_svm_under$auc)
# # Test Evaluation
# prob_svm_class_under <- predict(model_svm_under, newdata = test, type = 'prob')
# pred_svm_class_under = predict(model_svm_under, newdata = test)
# confusionMatrix(as.factor(pred_svm_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_svm_under <- roc.curve(test[['class']], pred_svm_class_under, plotit = T, main = 'ROC Curve using SVM & undersample')
# print(roc_svm_under$auc)

# set.seed(642)
# # Model Creation
# model_svm_over <- train(form = class ~ ., data = over_train_data, method = 'svmRadial', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_svm_over)
# # Train Evaluation
# # prob_svm_class_over <- predict(model_svm_over, newdata = over_train_data, type = 'prob')
# # pred_svm_class_over = predict(model_svm_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_svm_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_svm_over <- roc.curve(over_train_data[['class']], pred_svm_class_over, plotit = T, main = 'ROC Curve using SVM & oversample')
# # print(roc_svm_over$auc)
# # Test Evaluation
# prob_svm_class_over <- predict(model_svm_over, newdata = test, type = 'prob')
# pred_svm_class_over = predict(model_svm_over, newdata = test)
# confusionMatrix(as.factor(pred_svm_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_svm_over <- roc.curve(test[['class']], pred_svm_class_over, plotit = T, main = 'ROC Curve using SVM & oversample')
# print(roc_svm_over$auc)

# set.seed(642)
# # Model Creation
# model_svm_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'svmRadial', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_svm_ovun)
# # Train Evaluation
# # prob_svm_class_ovun <- predict(model_svm_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_svm_class_ovun = predict(model_svm_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_svm_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_svm_ovun <- roc.curve(ovun_train_data[['class']], pred_svm_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_svm_ovun$auc)
# # Test Evaluation
# prob_svm_class_ovun <- predict(model_svm_ovun, newdata = test, type = 'prob')
# pred_svm_class_ovun = predict(model_svm_ovun, newdata = test)
# confusionMatrix(as.factor(pred_svm_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_svm_ovun <- roc.curve(test[['class']], pred_svm_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_svm_ovun$auc)

# set.seed(642)
# # Model Creation
# model_svm_rose <- train(form = class ~ ., data = rose_train_data, method = 'svmRadial', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_svm_rose)
# # Train Evaluation
# # prob_svm_class_rose <- predict(model_svm_rose, newdata = rose_train_data, type = 'prob')
# # pred_svm_class_rose = predict(model_svm_rose, newdata = rose_train_data)
# # confusionMatrix(as.factor(pred_svm_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_svm_rose <- roc.curve(rose_train_data[['class']], pred_svm_class_rose, plotit = T, main = 'ROC Curve using SVM & ROSE')
# # print(roc_svm_rose$auc)
# # Test Evaluation
# prob_svm_class_rose <- predict(model_svm_rose, newdata = test, type = 'prob')
# pred_svm_class_rose = predict(model_svm_rose, newdata = test)
# confusionMatrix(as.factor(pred_svm_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_svm_rose <- roc.curve(test[['class']], pred_svm_class_rose, plotit = T, main = 'ROC Curve using SVM & ROSE')
# print(roc_svm_rose$auc)

# set.seed(642)
# # Model Creation
# model_svm_smote <- train(form = class ~ ., data = smote_train_data, method = 'svmRadial', trControl = ctrl_noTune, preProcess = c("center", "scale"), metric = 'ROC')
# # plot(model_svm_smote)
# # Train Evaluation
# # prob_svm_class_smote <- predict(model_svm_smote, newdata = smote_train_data, type = 'prob')
# # pred_svm_class_smote = predict(model_svm_smote, newdata = smote_train_data)
# # confusionMatrix(as.factor(pred_svm_class_smote), as.factor(smote_train_data), positive = 'F', mode = 'everything')
# # roc_svm_smote <- roc.curve(smote_train_data[['class']], pred_svm_class_smote, plotit = T, main = 'ROC Curve using SVM & SMOTE')
# # print(roc_svm_smote$auc)
# # Test Evaluation
# prob_svm_class_smote <- predict(model_svm_smote, newdata = test, type = 'prob')
# pred_svm_class_smote = predict(model_svm_smote, newdata = test)
# confusionMatrix(as.factor(pred_svm_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_svm_smote <- roc.curve(test[['class']], pred_svm_class_smote, plotit = T, main = 'ROC Curve using SVM & SMOTE')
# print(roc_svm_smote$auc)




# ============================================================================= #
# CLASSIFICATION MODEL WITH TUNE
# 01. Generalized Linear Model
#	A. Imbalanced Dataset
#	B. Undersampling
#	C. Oversampling
#	D. Mixedsampling
#	E. ROSE
#	F. SMOTE
# ============================================================================= #


# Model Creation
model_glm_imba <- train(form = class ~ ., data = imba_train_data, family = binomial(link = 'logit'), trControl = ctrl_tune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
#exp(coef(model_glm_imba$finalModel))
#varImp(model_glm_imba)
# Train Evaluation
# prob_glm_class_imba = predict(model_glm_imba, newdata = imba_train_data, type = 'prob')
# pred_glm_class_imba = predict(model_glm_imba, newdata = imba_train_data)
# confusionMatrix(as.factor(pred_glm_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_glm_imba <- roc.curve(imba_train_data[['class']], pred_glm_class_imba, plotit = T, main = 'ROC Curve using GLM')
# print(roc_glm_imba$auc)
# Test Evaluation
prob_glm_class_imba = predict(model_glm_imba, newdata = test, type = 'prob')
pred_glm_class_imba = predict(model_glm_imba, newdata = test)
confusionMatrix(as.factor(pred_glm_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_glm_imba <- roc.curve(test[['class']], pred_glm_class_imba, plotit = T, main = 'ROC Curve using GLM')
print(roc_glm_imba$auc)

# 
# # Model Creation
# model_glm_under <- train(form = class ~ ., data = under_train_data, family = binomial(link = 'logit'), trControl = ctrl_tune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# # Train Evaluation
# # prob_glm_class_under = predict(model_glm_under, newdata = under_train_data, type = 'prob')
# # pred_glm_class_under = predict(model_glm_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_glm_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_glm_under <- roc.curve(under_train_data[['class']], pred_glm_class_under, plotit = T, main = 'ROC Curve using GLM & Undersmaple')
# # print(roc_glm_under$auc)
# # Test Evaluation
# prob_glm_class_under = predict(model_glm_under, newdata = test, type = 'prob')
# pred_glm_class_under = predict(model_glm_under, newdata = test)
# confusionMatrix(as.factor(pred_glm_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_glm_under <- roc.curve(test[['class']], pred_glm_class_under, plotit = T, main = 'ROC Curve using GLM & Undersmaple')
# print(roc_glm_under$auc)
# 
# 
# # Model Creation
# model_glm_over <- train(form = class ~ ., data = over_train_data, family = binomial(link = 'logit'), trControl = ctrl_tune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# # Train Evaluation
# # prob_glm_class_over = predict(model_glm_over, newdata = over_train_data, type = 'prob')
# # pred_glm_class_over = predict(model_glm_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_glm_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_glm_over <- roc.curve(over_train_data[['class']], pred_glm_class_over, plotit = T, main = 'ROC Curve using GLM & Oversmaple')
# # print(roc_glm_over$auc)
# # Test Evaluation
# prob_glm_class_over = predict(model_glm_over, newdata = test, type = 'prob')
# pred_glm_class_over = predict(model_glm_over, newdata = test)
# confusionMatrix(as.factor(pred_glm_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_glm_over <- roc.curve(test[['class']], pred_glm_class_over, plotit = T, main = 'ROC Curve using GLM & Oversmaple')
# print(roc_glm_over$auc)
# 
# 
# # Model Creation
# model_glm_ovun <- train(form = class ~ ., data = ovun_train_data, family = binomial(link = 'logit'), trControl = ctrl_tune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# # Train Evaluation
# # prob_glm_class_ovun = predict(model_glm_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_glm_class_ovun = predict(model_glm_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_glm_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_glm_under <- roc.curve(ovun_train_data[['class']], pred_glm_class_ovun, plotit = T, main = 'ROC Curve using GLM & Miexed-smaple')
# # print(roc_glm_under$auc)
# # Test Evaluation
# prob_glm_class_ovun = predict(model_glm_ovun, newdata = test, type = 'prob')
# pred_glm_class_ovun = predict(model_glm_ovun, newdata = test)
# confusionMatrix(as.factor(pred_glm_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_glm_ovun <- roc.curve(test[['class']], pred_glm_class_ovun, plotit = T, main = 'ROC Curve using GLM & Miexed-smaple')
# print(roc_glm_ovun$auc)


# Model Creation
model_glm_rose <- train(form = class ~ ., data = rose_train_data, family = binomial(link = 'logit'), trControl = ctrl_tune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# Train Evaluation
# prob_glm_class_rose = predict(model_glm_rose, newdata = rose_train_data, type = 'prob')
# pred_glm_class_rose = predict(model_glm_rose, newdata = rose_train_data)
# confusionMatrix(as.factor(pred_glm_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_glm_rose <- roc.curve(rose_train_data[['class']], pred_glm_class_rose, plotit = T, main = 'ROC Curve using GLM & ROSE')
# print(roc_glm_rose$auc)
# Test Evaluation
prob_glm_class_rose = predict(model_glm_rose, newdata = test, type = 'prob')
pred_glm_class_rose = predict(model_glm_rose, newdata = test)
confusionMatrix(as.factor(pred_glm_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_glm_rose <- roc.curve(test[['class']], pred_glm_class_rose, plotit = T, main = 'ROC Curve using GLM & ROSE')
print(roc_glm_rose$auc)


# Model Creation
model_glm_smote <- train(form = class ~ ., data = smote_train_data, family = binomial(link = 'logit'), trControl = ctrl_tune, method = 'glm', preProcess = c("center", "scale"), metric = 'ROC')
# Train Evaluation
# prob_glm_class_smote = predict(model_glm_smote, newdata = smote_train_data, type = 'prob')
# pred_glm_class_smote = predict(model_glm_smote, newdata = smote_train_data)
# confusionMatrix(as.factor(pred_glm_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_glm_smote <- roc.curve(smote_train_data[['class']], pred_glm_class_smote, plotit = T, main = 'ROC Curve using GLM & SMOTE')
# print(roc_glm_smote$auc)
# Test Evaluation
prob_glm_class_smote = predict(model_glm_smote, newdata = test, type = 'prob')
pred_glm_class_smote = predict(model_glm_smote, newdata = test)
confusionMatrix(as.factor(pred_glm_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_glm_smote <- roc.curve(test[['class']], pred_glm_class_smote, plotit = T, main = 'ROC Curve using GLM & SMOTE')
print(roc_glm_smote$auc)


# ============================================================================= #
# CLASSIFICATION MODEL
# 02. Decision Tree - ctree, chaid, C5.0, xgbTree
#	A. Imbalanced Dataset
#	B. Undersampling
#	C. Oversampling
#	D. Mixedsampling
#	E. ROSE
#	F. SMOTE
# ============================================================================= #


# Model Creation
model_dt_imba <- train(form = class ~ ., data = imba_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_dt_imba)
# Train Evaluation
# prob_dt_class_imba <- predict(model_dt_imba, newdata = imba_train_data, type = 'prob')
# pred_dt_class_imba = predict(model_dt_imba, newdata = imba_train_data)
# confusionMatrix(as.factor(pred_dt_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# # Error in ROC Calculation
# roc_dt_imba <- roc.curve(imba_train_data[['class']], pred_dt_class_imba, plotit = T, main = 'ROC Curve using DT')
# print(roc_dt_imba$auc)
# Test Evaluation
prob_dt_class_imba <- predict(model_dt_imba, newdata = test, type = 'prob')
pred_dt_class_imba = predict(model_dt_imba, newdata = test)
confusionMatrix(as.factor(pred_dt_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_dt_imba <- roc.curve(test[['class']], pred_dt_class_imba, plotit = T, main = 'ROC Curve using DT')
print(roc_dt_imba$auc)


# # Model Creation
# model_dt_under <- train(form = class ~ ., data = under_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_dt_under)
# # Train Evaluation
# # prob_dt_class_under <- predict(model_dt_under, newdata = under_train_data, type = 'prob')
# # pred_dt_class_under = predict(model_dt_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_dt_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_dt_under <- roc.curve(under_train_data[['class']], pred_dt_class_under, plotit = T, main = 'ROC Curve using DT & undersample')
# # print(roc_dt_under$auc)
# # Test Evaluation
# prob_dt_class_under <- predict(model_dt_under, newdata = test, type = 'prob')
# pred_dt_class_under = predict(model_dt_under, newdata = test)
# confusionMatrix(as.factor(pred_dt_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_dt_under <- roc.curve(test[['class']], pred_dt_class_under, plotit = T, main = 'ROC Curve using DT & undersample')
# print(roc_dt_under$auc)
# 
# 
# # Model Creation
# model_dt_over <- train(form = class ~ ., data = over_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_dt_over)
# # Train Evaluation
# # prob_dt_class_over <- predict(model_dt_over, newdata = over_train_data, type = 'prob')
# # pred_dt_class_over = predict(model_dt_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_dt_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_dt_over <- roc.curve(over_train_data[['class']], pred_dt_class_over, plotit = T, main = 'ROC Curve using DT & oversample')
# # print(roc_dt_over$auc)
# # Test Evaluation
# prob_dt_class_over <- predict(model_dt_over, newdata = test, type = 'prob')
# pred_dt_class_over = predict(model_dt_over, newdata = test)
# confusionMatrix(as.factor(pred_dt_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_dt_over <- roc.curve(test[['class']], pred_dt_class_over, plotit = T, main = 'ROC Curve using DT & oversampling')
# print(roc_dt_over$auc)
# 
# 
# # Model Creation
# model_dt_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_dt_ovun)
# # Train Evaluation
# # prob_dt_class_ovun <- predict(model_dt_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_dt_class_ovun = predict(model_dt_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_dt_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_dt_ovun <- roc.curve(ovun_train_data[['class']], pred_dt_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_dt_ovun$auc)
# # Test Evaluation
# prob_dt_class_ovun <- predict(model_dt_ovun, newdata = test, type = 'prob')
# pred_dt_class_ovun = predict(model_dt_ovun, newdata = test)
# confusionMatrix(as.factor(pred_dt_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_dt_ovun <- roc.curve(test[['class']], pred_dt_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_dt_ovun$auc)


# Model Creation
model_dt_rose <- train(form = class ~ ., data = rose_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_dt_rose)
# Train Evaluation
# prob_dt_class_rose <- predict(model_dt_rose, newdata = rose_train_data, type = 'prob')
# pred_dt_class_rose = predict(model_dt_rose, newdata = rose_train_data)
# confusionMatrix(as.factor(pred_dt_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_dt_rose <- roc.curve(rose_train_data[['class']], pred_dt_class_rose, plotit = T, main = 'ROC Curve using DT & ROSE')
# print(roc_dt_rose$auc)
# Test Evaluation
prob_dt_class_rose <- predict(model_dt_rose, newdata = test, type = 'prob')
pred_dt_class_rose = predict(model_dt_rose, newdata = test)
confusionMatrix(as.factor(pred_dt_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_dt_rose <- roc.curve(test[['class']], pred_dt_class_rose, plotit = T, main = 'ROC Curve using DT & ROSE')
print(roc_dt_rose$auc)


# Model Creation
model_dt_smote <- train(form = class ~ ., data = smote_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_dt_smote)
# Train Evaluation
# prob_dt_class_smote <- predict(model_dt_smote, newdata = smote_train_data, type = 'prob')
# pred_dt_class_smote = predict(model_dt_smote, newdata = smote_train_data)
# confusionMatrix(as.factor(pred_dt_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_dt_smote <- roc.curve(smote_train_data[['class']], pred_dt_class_smote, plotit = T, main = 'ROC Curve using DT & SMOTE')
# print(roc_dt_smote$auc)
# Test Evaluation
prob_dt_class_smote <- predict(model_dt_smote, newdata = test, type = 'prob')
pred_dt_class_smote = predict(model_dt_smote, newdata = test)
confusionMatrix(as.factor(pred_dt_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_dt_smote <- roc.curve(test[['class']], pred_dt_class_smote, plotit = T, main = 'ROC Curve using DT & SMOTE')
print(roc_dt_smote$auc)


# ============================================================================= #
# CLASSIFICATION MODEL
# 03. Random Forest - ranger, rf
#	A. Imbalanced Dataset
#	B. Undersampling
#	C. Oversampling
#	D. Mixedsampling
#	E. ROSE
#	F. SMOTE
# ============================================================================= #


# Model Creation
model_rf_imba <- train(form = class ~ ., data = imba_train_data, method = 'rf', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_rf_imba)
# Train Evaluation
# # prob_rf_class_imba <- predict(model_rf_imba, newdata = imba_train_data, type = 'prob')
# pred_rf_class_imba = predict(model_rf_imba, newdata = imba_train_data)
# confusionMatrix(as.factor(pred_rf_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_rf_imba <- roc.curve(pred_rf_class_imba, imba_train_data[['class']], plotit = T, main = 'ROC Curve using RF')
# print(roc_rf_imba$auc)
# Test Evaluation
pred_rf_class_imba = predict(model_rf_imba, newdata = test)
confusionMatrix(as.factor(pred_rf_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_rf_imba <- roc.curve(test[['class']], pred_rf_class_imba, plotit = T, main = 'ROC Curve using RF')
print(roc_rf_imba$auc)

# 
# # Model Creation
# model_rf_under <- train(form = class ~ ., data = under_train_data, method = 'rf', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_rf_under)
# # Train Evaluation
# # pred_rf_class_under = predict(model_rf_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_rf_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_rf_under <- roc.curve(under_train_data[['class']], pred_rf_class_under, plotit = T, main = 'ROC Curve using RF & undersample')
# # print(roc_rf_under$auc)
# # Test Evaluation
# pred_rf_class_under = predict(model_rf_under, newdata = test)
# confusionMatrix(as.factor(pred_rf_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_rf_under <- roc.curve(test[['class']], pred_rf_class_under, plotit = T, main = 'ROC Curve using RF & undersample')
# print(roc_rf_under$auc)
# 
# 
# # Model Creation
# model_rf_over <- train(form = class ~ ., data = over_train_data, method = 'rf', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_rf_over)
# # Train Evaluation
# # pred_rf_class_over = predict(model_rf_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_rf_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_rf_over <- roc.curve(over_train_data[['class']], pred_rf_class_over, plotit = T, main = 'ROC Curve using RF & oversample')
# # print(roc_rf_over$auc)
# # Test Evaluation
# pred_rf_class_over = predict(model_rf_over, newdata = test)
# confusionMatrix(as.factor(pred_rf_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_rf_over <- roc.curve(test[['class']], pred_rf_class_over, plotit = T, main = 'ROC Curve using RF & oversample')
# print(roc_rf_over$auc)
# 
# 
# # Model Creation
# model_rf_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'rf', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_rf_ovun)
# # Train Evaluation
# # pred_rf_class_ovun = predict(model_rf_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_rf_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_rf_ovun <- roc.curve(ovun_train_data[['class']], pred_rf_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_rf_ovun$auc)
# # Test Evaluation
# pred_rf_class_ovun = predict(model_rf_ovun, newdata = test)
# confusionMatrix(as.factor(pred_rf_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_rf_ovun <- roc.curve(test[['class']], pred_rf_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_rf_ovun$auc)


# Model Creation
model_rf_rose <- train(form = class ~ ., data = rose_train_data, method = 'rf', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_rf_rose)
# Train Evaluation
# pred_rf_class_rose = predict(model_rf_rose, newdata = rose_train_data)
# confusionMatrix(as.factor(pred_rf_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_rf_rose <- roc.curve(rose_train_data[['class']], pred_rf_class_rose, plotit = T, main = 'ROC Curve using RF & ROSE')
# print(roc_rf_rose$auc)
# Test Evaluation
pred_rf_class_rose = predict(model_rf_rose, newdata = test)
confusionMatrix(as.factor(pred_rf_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_rf_rose <- roc.curve(test[['class']], pred_rf_class_rose, plotit = T, main = 'ROC Curve using RF & ROSE')
print(roc_rf_rose$auc)


# Model Creation
model_rf_smote <- train(form = class ~ ., data = smote_train_data, method = 'rf', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_rf_smote)
# # Train Evaluation
# pred_rf_class_smote = predict(model_rf_smote, newdata = smote_train_data)
# confusionMatrix(as.factor(pred_rf_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_rf_smote <- roc.curve(smote_train_data[['class']], pred_rf_class_smote, plotit = T, main = 'ROC Curve using RF & SMOTE')
# print(roc_rf_smote$auc)
# Test Evaluation
pred_rf_class_smote = predict(model_rf_smote, newdata = test)
confusionMatrix(as.factor(pred_rf_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_rf_smote <- roc.curve(test[['class']], pred_rf_class_smote, plotit = T, main = 'ROC Curve using RF & SMOTE')
print(roc_rf_smote$auc)


# ============================================================================= #
# CLASSIFICATION MODEL
# 04. KNN
#	A. Imbalanced Dataset
#	B. Undersampling
#	C. Oversampling
#	D. Mixedsampling
#	E. ROSE
#	F. SMOTE
# ============================================================================= #


# Model Creation
model_knn_imba <- train(form = class ~ ., data = imba_train_data, method = 'knn', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_knn_imba)
# Train Evaluation
# prob_knn_class_imba <- predict(model_knn_imba, newdata = imba_train_data, type = 'prob')
# pred_knn_class_imba = predict(model_knn_imba, newdata = imba_train_data)
# confusionMatrix(as.factor(pred_knn_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_knn_imba <- roc.curve(imba_train_data[['class']], pred_knn_class_imba, plotit = T, main = 'ROC Curve using KNN')
# print(roc_knn_imba$auc)
# Test Evaluation
prob_knn_class_imba <- predict(model_knn_imba, newdata = test, type = 'prob')
pred_knn_class_imba = predict(model_knn_imba, newdata = test)
confusionMatrix(as.factor(pred_knn_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_knn_imba <- roc.curve(test[['class']], pred_knn_class_imba, plotit = T, main = 'ROC Curve using KNN')
print(roc_knn_imba$auc)

# 
# # Model Creation
# model_knn_under <- train(form = class ~ ., data = under_train_data, method = 'knn', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_knn_under)
# # Train Evaluation
# # prob_knn_class_under <- predict(model_knn_under, newdata = under_train_data, type = 'prob')
# # pred_knn_class_under = predict(model_knn_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_knn_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_knn_under <- roc.curve(under_train_data[['class']], pred_knn_class_under, plotit = T, main = 'ROC Curve using KNN and undersample')
# # print(roc_knn_under$auc)
# # Test Evaluation
# prob_knn_class_under <- predict(model_knn_under, newdata = test, type = 'prob')
# pred_knn_class_under = predict(model_knn_under, newdata = test)
# confusionMatrix(as.factor(pred_knn_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_knn_under <- roc.curve(test[['class']], pred_knn_class_under, plotit = T, main = 'ROC Curve using KNN and undersample')
# print(roc_knn_under$auc)
# 
# 
# # Model Creation
# model_knn_over <- train(form = class ~ ., data = over_train_data, method = 'knn', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_knn_over)
# # Train Evaluation
# # prob_knn_class_over <- predict(model_knn_over, newdata = over_train_data, type = 'prob')
# # pred_knn_class_over = predict(model_knn_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_knn_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_knn_over <- roc.curve(over_train_data[['class']], pred_knn_class_over, plotit = T, main = 'ROC Curve using KNN and oversample')
# # print(roc_knn_over$auc)
# # Test Evaluation
# prob_knn_class_over <- predict(model_knn_over, newdata = test, type = 'prob')
# pred_knn_class_over = predict(model_knn_over, newdata = test)
# confusionMatrix(as.factor(pred_knn_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_knn_over <- roc.curve(test[['class']], pred_knn_class_over, plotit = T, main = 'ROC Curve using KNN and oversample')
# print(roc_knn_over$auc)
# 
# 
# # Model Creation
# model_knn_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'knn', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_knn_ovun)
# # Train Evaluation
# # prob_knn_class_ovun <- predict(model_knn_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_knn_class_ovun = predict(model_knn_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_knn_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_knn_ovun <- roc.curve(ovun_train_data[['class']], pred_knn_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_knn_ovun$auc)
# # Test Evaluation
# prob_knn_class_ovun <- predict(model_knn_ovun, newdata = test, type = 'prob')
# pred_knn_class_ovun = predict(model_knn_ovun, newdata = test)
# confusionMatrix(as.factor(pred_knn_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_knn_ovun <- roc.curve(test[['class']], pred_knn_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_knn_ovun$auc)


# Model Creation
model_knn_rose <- train(form = class ~ ., data = rose_train_data, method = 'knn', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_knn_rose)
# Train Evaluation
# prob_knn_class_rose <- predict(model_knn_rose, newdata = rose_train_data, type = 'prob')
# pred_knn_class_rose = predict(model_knn_rose, newdata = rose_train_data)
# confusionMatrix(as.factor(pred_knn_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_knn_rose <- roc.curve(rose_train_data[['class']], pred_knn_class_rose, plotit = T, main = 'ROC Curve using KNN and ROSE')
# print(roc_knn_rose$auc)
# Test Evaluation
prob_knn_class_rose <- predict(model_knn_rose, newdata = test, type = 'prob')
pred_knn_class_rose = predict(model_knn_rose, newdata = test)
confusionMatrix(as.factor(pred_knn_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_knn_rose <- roc.curve(test[['class']], pred_knn_class_rose, plotit = T, main = 'ROC Curve using KNN and ROSE')
print(roc_knn_rose$auc)


# Model Creation
model_knn_smote <- train(form = class ~ ., data = smote_train_data, method = 'knn', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_knn_smote)
# Train Evaluation
# prob_knn_class_smote <- predict(model_knn_smote, newdata = smote_train_data, type = 'prob')
# pred_knn_class_smote = predict(model_knn_smote, newdata = smote_train_data)
# confusionMatrix(as.factor(pred_knn_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_knn_smote <- roc.curve(smote_train_data[['class']], pred_knn_class_smote, plotit = T, main = 'ROC Curve using KNN and SMOTE')
# print(roc_knn_smote$auc)
# Test Evaluation
prob_knn_class_smote <- predict(model_knn_smote, newdata = test, type = 'prob')
pred_knn_class_smote = predict(model_knn_smote, newdata = test)
confusionMatrix(as.factor(pred_knn_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_knn_smote <- roc.curve(test[['class']], pred_knn_class_smote, plotit = T, main = 'ROC Curve using KNN and SMOTE')
print(roc_knn_smote$auc)


# ============================================================================= #
# CLASSIFICATION MODEL
# 05. Naive Bayes
#	A. Imbalanced Dataset
#	B. Undersampling
#	C. Oversampling
#	D. Mixedsampling
#	E. ROSE
#	F. SMOTE
# ============================================================================= #


# Model Creation
model_nb_imba <- train(form = class ~ ., data = imba_train_data, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_nb_imba)
# Train Evaluation
# prob_nb_class_imba <- predict(model_nb_imba, newdata = imba_train_data, type = 'prob')
# pred_nb_class_imba = predict(model_nb_imba, newdata = imba_train_data)
# confusionMatrix(as.factor(pred_nb_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_imba <- roc.curve(imba_train_data[['class']], pred_nb_class_imba, plotit = T, main = 'ROC Curve using NB')
# print(roc_nb_imba$auc)
# Test Evaluation
prob_nb_class_imba <- predict(model_nb_imba, newdata = test, type = 'prob')
pred_nb_class_imba = predict(model_nb_imba, newdata = test)
confusionMatrix(as.factor(pred_nb_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_nb_imba <- roc.curve(test[['class']], pred_nb_class_imba, plotit = T, main = 'ROC Curve using NB')
print(roc_nb_imba$auc)

# 
# # Model Creation
# model_nb_under <- train(form = class ~ ., data = under_train_data, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_nb_under)
# # Train Evaluation
# # prob_nb_class_under <- predict(model_nb_under, newdata = under_train_data, type = 'prob')
# # pred_nb_class_under = predict(model_nb_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_nb_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_nb_under <- roc.curve(under_train_data[['class']], pred_nb_class_under, plotit = T, main = 'ROC Curve using NB & undersample')
# # print(roc_nb_under$auc)
# # Test Evaluation
# prob_nb_class_under <- predict(model_nb_under, newdata = test, type = 'prob')
# pred_nb_class_under = predict(model_nb_under, newdata = test)
# confusionMatrix(as.factor(pred_nb_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_under <- roc.curve(test[['class']], pred_nb_class_under, plotit = T, main = 'ROC Curve using NB & undersample')
# print(roc_nb_under$auc)
# 
# 
# # Model Creation
# model_nb_over <- train(form = class ~ ., data = over_train_data, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_nb_over)
# # Train Evaluation
# # prob_nb_class_over <- predict(model_nb_over, newdata = over_train_data, type = 'prob')
# # pred_nb_class_over = predict(model_nb_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_nb_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_nb_over <- roc.curve(over_train_data[['class']], pred_nb_class_over, plotit = T, main = 'ROC Curve using NB & oversample')
# # print(roc_nb_over$auc)
# # Test Evaluation
# prob_nb_class_over <- predict(model_nb_over, newdata = test, type = 'prob')
# pred_nb_class_over = predict(model_nb_over, newdata = test)
# confusionMatrix(as.factor(pred_nb_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_over <- roc.curve(test[['class']], pred_nb_class_over, plotit = T, main = 'ROC Curve using NB & oversample')
# print(roc_nb_over$auc)
# 
# 
# # Model Creation
# model_nb_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_nb_ovun)
# # Train Evaluation
# # prob_nb_class_ovun <- predict(model_nb_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_nb_class_ovun = predict(model_nb_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_nb_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_nb_ovun <- roc.curve(ovun_train_data[['class']], pred_nb_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_nb_ovun$auc)
# # Test Evaluation
# prob_nb_class_ovun <- predict(model_nb_ovun, newdata = test, type = 'prob')
# pred_nb_class_ovun = predict(model_nb_ovun, newdata = test)
# confusionMatrix(as.factor(pred_nb_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_nb_ovun <- roc.curve(test[['class']], pred_nb_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_nb_ovun$auc)


# Model Creation
model_nb_rose <- train(form = class ~ ., data = rose_train_data, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_nb_rose)
# Train Evaluation
# prob_nb_class_rose <- predict(model_nb_rose, newdata = rose_train_data, type = 'prob')
# pred_nb_class_rose = predict(model_nb_rose, newdata = rose_train_data)
# confusionMatrix(as.factor(pred_nb_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_nb_rose <- roc.curve(rose_train_data[['class']], pred_nb_class_rose, plotit = T, main = 'ROC Curve using NB & ROSE')
# print(roc_nb_rose$auc)
# Test Evaluation
prob_nb_class_rose <- predict(model_nb_rose, newdata = test, type = 'prob')
pred_nb_class_rose = predict(model_nb_rose, newdata = test)
confusionMatrix(as.factor(pred_nb_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_nb_rose <- roc.curve(test[['class']], pred_nb_class_rose, plotit = T, main = 'ROC Curve using NB & ROSE')
print(roc_nb_rose$auc)


# Model Creation
model_nb_smote <- train(form = class ~ ., data = smote_train_data, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_nb_smote)
# Train Evaluation
# prob_nb_class_smote <- predict(model_nb_smote, newdata = smote_train_data, type = 'prob')
# pred_nb_class_smote = predict(model_nb_smote, newdata = smote_train_data)
# confusionMatrix(as.factor(pred_nb_class_smote), as.factor(smote_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_nb_smote <- roc.curve(smote_train_data[['class']], pred_nb_class_smote, plotit = T, main = 'ROC Curve using NB & SMOTE')
# print(roc_nb_smote$auc)
# Test Evaluation
prob_nb_class_smote <- predict(model_nb_smote, newdata = test, type = 'prob')
pred_nb_class_smote = predict(model_nb_smote, newdata = test)
confusionMatrix(as.factor(pred_nb_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_nb_smote <- roc.curve(test[['class']], pred_nb_class_smote, plotit = T, main = 'ROC Curve using NB & SMOTE')
print(roc_nb_smote$auc)


# ============================================================================= #
# CLASSIFICATION MODEL
# 06. SVM
#	A. Imbalanced Dataset
#	B. Undersampling
#	C. Oversampling
#	D. Mixedsampling
#	E. ROSE
#	F. SMOTE
# ============================================================================= #


# Model Creation
model_svm_imba <- train(form = class ~ ., data = imba_train_data, method = 'svmRadial', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_svm_imba)
# Train Evaluation
# prob_svm_class_imba <- predict(model_svm_imba, newdata = imba_train_data, type = 'prob')
# pred_svm_class_imba = predict(model_svm_imba, newdata = imba_train_data)
# confusionMatrix(as.factor(pred_svm_class_imba), as.factor(imba_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_svm_imba <- roc.curve(imba_train_data[['class']], pred_svm_class_imba, plotit = T, main = 'ROC Curve using SVM')
# print(roc_svm_imba$auc)
# Test Evaluation
prob_svm_class_imba <- predict(model_svm_imba, newdata = test, type = 'prob')
pred_svm_class_imba = predict(model_svm_imba, newdata = test)
confusionMatrix(as.factor(pred_svm_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_svm_imba <- roc.curve(test[['class']], pred_svm_class_imba, plotit = T, main = 'ROC Curve using SVM')
print(roc_svm_imba$auc)

# 
# # Model Creation
# model_svm_under <- train(form = class ~ ., data = under_train_data, method = 'svmRadial', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_svm_under)
# # Train Evaluation
# # prob_svm_class_under <- predict(model_svm_under, newdata = under_train_data, type = 'prob')
# # pred_svm_class_under = predict(model_svm_under, newdata = under_train_data)
# # confusionMatrix(as.factor(pred_svm_class_under), as.factor(under_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_svm_under <- roc.curve(under_train_data[['class']], pred_svm_class_under, plotit = T, main = 'ROC Curve using SVM & undersample')
# # print(roc_svm_under$auc)
# # Test Evaluation
# prob_svm_class_under <- predict(model_svm_under, newdata = test, type = 'prob')
# pred_svm_class_under = predict(model_svm_under, newdata = test)
# confusionMatrix(as.factor(pred_svm_class_under), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_svm_under <- roc.curve(test[['class']], pred_svm_class_under, plotit = T, main = 'ROC Curve using SVM & undersample')
# print(roc_svm_under$auc)
# 
# 
# # Model Creation
# model_svm_over <- train(form = class ~ ., data = over_train_data, method = 'svmRadial', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_svm_over)
# # Train Evaluation
# # prob_svm_class_over <- predict(model_svm_over, newdata = over_train_data, type = 'prob')
# # pred_svm_class_over = predict(model_svm_over, newdata = over_train_data)
# # confusionMatrix(as.factor(pred_svm_class_over), as.factor(over_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_svm_over <- roc.curve(over_train_data[['class']], pred_svm_class_over, plotit = T, main = 'ROC Curve using SVM & oversample')
# # print(roc_svm_over$auc)
# # Test Evaluation
# prob_svm_class_over <- predict(model_svm_over, newdata = test, type = 'prob')
# pred_svm_class_over = predict(model_svm_over, newdata = test)
# confusionMatrix(as.factor(pred_svm_class_over), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_svm_over <- roc.curve(test[['class']], pred_svm_class_over, plotit = T, main = 'ROC Curve using SVM & oversample')
# print(roc_svm_over$auc)
# 
# 
# # Model Creation
# model_svm_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'svmRadial', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
# plot(model_svm_ovun)
# # Train Evaluation
# # prob_svm_class_ovun <- predict(model_svm_ovun, newdata = ovun_train_data, type = 'prob')
# # pred_svm_class_ovun = predict(model_svm_ovun, newdata = ovun_train_data)
# # confusionMatrix(as.factor(pred_svm_class_ovun), as.factor(ovun_train_data[['class']]), positive = 'F', mode = 'everything')
# # roc_svm_ovun <- roc.curve(ovun_train_data[['class']], pred_svm_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# # print(roc_svm_ovun$auc)
# # Test Evaluation
# prob_svm_class_ovun <- predict(model_svm_ovun, newdata = test, type = 'prob')
# pred_svm_class_ovun = predict(model_svm_ovun, newdata = test)
# confusionMatrix(as.factor(pred_svm_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
# roc_svm_ovun <- roc.curve(test[['class']], pred_svm_class_ovun, plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_svm_ovun$auc)


# Model Creation
model_svm_rose <- train(form = class ~ ., data = rose_train_data, method = 'svmRadial', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_svm_rose)
# Train Evaluation
# prob_svm_class_rose <- predict(model_svm_rose, newdata = rose_train_data, type = 'prob')
# pred_svm_class_rose = predict(model_svm_rose, newdata = rose_train_data)
# confusionMatrix(as.factor(pred_svm_class_rose), as.factor(rose_train_data[['class']]), positive = 'F', mode = 'everything')
# roc_svm_rose <- roc.curve(rose_train_data[['class']], pred_svm_class_rose, plotit = T, main = 'ROC Curve using SVM & ROSE')
# print(roc_svm_rose$auc)
# Test Evaluation
prob_svm_class_rose <- predict(model_svm_rose, newdata = test, type = 'prob')
pred_svm_class_rose = predict(model_svm_rose, newdata = test)
confusionMatrix(as.factor(pred_svm_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_svm_rose <- roc.curve(test[['class']], pred_svm_class_rose, plotit = T, main = 'ROC Curve using SVM & ROSE')
print(roc_svm_rose$auc)


# Model Creation
model_svm_smote <- train(form = class ~ ., data = smote_train_data, method = 'svmRadial', tuneLength = 5, trControl = ctrl_tune, preProcess = c("center", "scale"), metric = 'ROC')
plot(model_svm_smote)
# Train Evaluation
# prob_svm_class_smote <- predict(model_svm_smote, newdata = smote_train_data, type = 'prob')
# pred_svm_class_smote = predict(model_svm_smote, newdata = smote_train_data)
# confusionMatrix(as.factor(pred_svm_class_smote), as.factor(smote_train_data), positive = 'F', mode = 'everything')
# roc_svm_smote <- roc.curve(smote_train_data[['class']], pred_svm_class_smote, plotit = T, main = 'ROC Curve using SVM & SMOTE')
# print(roc_svm_smote$auc)
# Test Evaluation
prob_svm_class_smote <- predict(model_svm_smote, newdata = test, type = 'prob')
pred_svm_class_smote = predict(model_svm_smote, newdata = test)
confusionMatrix(as.factor(pred_svm_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_svm_smote <- roc.curve(test[['class']], pred_svm_class_smote, plotit = T, main = 'ROC Curve using SVM & SMOTE')
print(roc_svm_smote$auc)









# ============================================================================= #
# MODEL COMPARISON FINAL
# ============================================================================= #

# # Compare model performances using resample()
# models_compare <- resamples(list(GLM = model_glm, DT = model_dt, RF = model_rf, KNN = model_knn, NB = model_nb, MARS = model_mars, NN = model_nn, SVM = model_svm))
# # Summary of the models performances
# summary(models_compare)
# # Draw box plots to compare models
# scales <- list(x=list(relation="free"), y=list(relation="free"))
# bwplot(models_compare, scales=scales)



