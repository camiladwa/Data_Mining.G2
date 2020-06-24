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
# 04. MAXITR: Maximum number of iterations to perform for MICE
# 05. EXTRACTITR: Extract imputated data set after performing MICE 
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

KNN_imputation <- function(df) {
  # Function to impute missing values(NAs) using K-Nearest Neighbour (k = 5)
  #impute_KNN <- knnImputation(df, k = 5, scale = T, meth = "weighAvg", distData = NULL)
  KNNImpute <- preProcess(df, method="knnImpute")
  impute_KNN <- predict(KNNImpute, newdata = df)
  return(impute_KNN)
}

BAG_imputation <- function(df) {
  # Function to impute missing values(NAs) using K-Nearest Neighbour (k = 5)
  BAGImpute <- preProcess(df, method="bagImpute")
  impute_BAG <- predict(BAGImpute, newdata = df)
  return(impute_BAG)
}

# mice_imputation <- function(df, itr = 5, extract = 5){
#   # Function to impute missing values(NAs) using MICE, default 5 iterations and extract data after completion of 5th iteration
#   impute_mice <- mice(df, m = itr, method = 'pmm', seed = 500)
#   summary(impute_mice)
#   mice_imputate <- complete(impute_mice, extract)
#   return(mice_imputate)
# }

# PCA_Analysis <- function(df) {
#   # Function to check the PCs using Scree Plot and Kaiser-Guttman Rule.
#   pca <- prcomp(df, center = TRUE, scale. = TRUE)
#   pca_var <- pca$sdev ^ 2
#   pca_var_per <- round(pca_var * 100 / sum(pca_var) , 2)
#   plot(1:length(pca_var), pca_var, type="b", col='blue', ylab="Eigenvalue", xlab="Component Number", main = 'Scree Plot') 
#   abline(h = 1,lty = 2,col = "red")
#   message('Principal components as per Kaiser-Guttman rule: ', length(pca_var[pca_var >= 1]))
#   plot(1:length(pca_var_per), pca_var_per, type="b", col='red', ylab="Proportion of Variance Explained", xlab="Component Number", main = 'Variance Explained by Components')
#   message('Total Variance explained by 118 PCs: ', sum(pca_var_per[1:length(pca_var[pca_var >= 1])]))
# }

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
# FEATURE REDUCTION (WIP)
# 03. Highly correlated features Removal (corr < 0.99)
# ============================================================================= #

# corr <- cor(train_variance_removal, use = "pairwise.complete.obs")
# highCor <- c(findCorrelation(corr, names = TRUE, cutoff = 0.99))
# train_highcorr_removal <- train_variance_removal[ , -which(names(train_variance_removal) %in% c(highCor))]
# missing_value_analysis(train_highcorr_removal, 'train set after removing Highly Correlated features')
# outlier_analysis(train_highcorr_removal, 'train set after removing highly correlated features')
# write_csv(cbind(class, train_highcorr_removal), "train_HC_removal.csv")


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
# NA HANDLING
# 03. MICE - DO NOT RUN, 
# 	Time consuming (5-7 hours for 5 iterations, 15 hours for 15 iterations)
#	NA reduced from 11183 to 121 after 5 iterations	
# ============================================================================= #

# train_mice_imputation <- mice_imputation(train_outlier_NA, MAXITR, EXTRACTITR)
# missing_value_analysis(train_mice_imputation, 'train set after MICE imputation')
# outlier_analysis(train_mice_imputation, 'train set after MICE imputation')
# write_csv(cbind(class, train_mice_imputation), "train_NAH_mice.csv")


# ============================================================================= #
# NA HANDLING
# 04. PCA
#	estimate the number of components from incomplete data
#	iterativePCA algorithm
# ============================================================================= #

# nPCs <- estim_ncpPCA(train_outlier_NA, method.cv = "Kfold", verbose = FALSE)
# res_comp <- imputePCA(train_outlier_NA, ncp = nPCs$ncp)
# train_NAH_pca <- as.data.frame(res_comp$completeObs)


# ============================================================================= #
# NA HANDLING
# 05. Bagged Tree
# ============================================================================= #

# train_bag_Imputation <- BAG_imputation(train_outlier_NA)
# missing_value_analysis(train_bag_Imputation, 'train set after Bagged Tree imputation')
# outlier_analysis(train_bag_Imputation, 'train set after Bagged Tree imputation')
# write_csv(cbind(class, train_bag_Imputation), "train_NAH_bt.csv")


# ============================================================================= #
# FEATURE SELECTION AND REDUCTION
# 01. PCA
# ============================================================================= #

# PCA_Analysis(train_knn_imputation)
# PCA_Analysis(train_NAH_pca)


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
# 02. BORUTA on PCA imputed train set
# ============================================================================= #

# train_FR_boruta_pca <- selection_boruta(train_NAH_pca, class)
# print("Summary of selected features")
# summary(train_FR_boruta_pca)
# message("Variance corresponding to selected features:")
# sapply(train_FR_boruta_pca, var)
# missing_value_analysis(train_FR_boruta_pca, 'train set after BORUTA')
# outlier_analysis(train_FR_boruta_pca, 'train set after BORUTA')
# write_csv(cbind(train_FR_boruta_pca), "train_FR_boruta.csv")


# ============================================================================= #
# FEATURE SELECTION AND REDUCTION
# 03. BORUTA on Bag Tree imputed train set
# ============================================================================= #

# train_FR_boruta_bt <- selection_boruta(train_bag_Imputation, class)
# print("Summary of selected features")
# summary(train_FR_boruta_bt)
# message("Variance corresponding to selected features:")
# sapply(train_FR_boruta_bt, var)
# missing_value_analysis(train_FR_boruta_bt, 'train set after BORUTA')
# outlier_analysis(train_FR_boruta_bt, 'train set after BORUTA')
# write_csv(cbind(train_FR_boruta_bt), "train_FR_boruta.csv")


# ============================================================================= #
# FEATURE SELECTION AND REDUCTION
# 04. Recursive Feature Elimination (RFE)
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
# ============================================================================= #

train <- cbind(class, train_FR_boruta_knn)


# ============================================================================= #
# PREPROCESS: TEST DATASET
# 01. KNN Imputation
# ============================================================================= #

test <- semcon_test_data[ , which(names(semcon_test_data) %in% c(names(train)))]
class_distribution(test, 'Test Sample')
missing_value_analysis(test, 'Test Sample')
outlier_analysis(test, 'Test Sample')
test <- KNN_imputation(test)
missing_value_analysis(test, 'Test set after KNN imputation')
outlier_analysis(test, 'Test Sample after KNN imputation')


# ============================================================================= #
# PREPROCESS: TEST DATASET
# 02. Bag Tree Imputation
# ============================================================================= #

# bagImpute <- preProcess(test, method="bagImpute")
# test <- predict(bagImpute, newdata = test)
# missing_value_analysis(test, 'Test set after KNN imputation')

# ============================================================================= #
# SAMPLING TECHNIQUES
# 01. CV
# 02. Bootstrap
# 03. K-Fold CV
# 04. Undersampling
# 05. Oversampling
# 06. ROSE
# 07. SMOTE
# ============================================================================= #

ctrl <- trainControl(method = 'repeatedcv', number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
# ctrl <- trainControl(method = "cv", number = 10, savePredictions = 'final',  summaryFunction = twoClassSummary)
# ctrl <- trainControl(method = "boot632", number = 1000, savePredictions = TRUE, savePredictions = 'final', classProbs = T, summaryFunction = twoClassSummary)
# ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE, savePredictions = 'final', classProbs = T, summaryFunction = twoClassSummary)
ctrl_under <- trainControl(method = 'repeatedcv', number = 10, verboseIter = FALSE, sampling = 'down')
ctrl_over <- trainControl(method = 'repeatedcv', number = 10, verboseIter = FALSE, sampling = 'up')
ctrl_rose <- trainControl(method = 'repeatedcv', number = 10, verboseIter = FALSE, sampling = 'rose')
ctrl_smote <- trainControl(method = 'repeatedcv', number = 10, verboseIter = FALSE, sampling = 'smote')

ovun_train_data <- ovun.sample(class ~., data = train, method = 'both', p = 0.5, seed = 222, N = 1253)$data

# ============================================================================= #
# CLASSIFICATION MODEL
# 01. Generalized Linear Model
#	A. Imbalanced Dataset
#	B. Undersampling
#	C. Oversampling
#	D. Mixedsampling
#	E. ROSE
#	F. SMOTE
# ============================================================================= #

set.seed(642)
# Train Evaluation
# model_glm <- train(form = class ~ ., data = train, family = binomial(link = 'logit'), trControl = ctrl, method = 'glm', preProcess = c("center", "scale"))
# #exp(coef(model_glm$finalModel))
# #varImp(model_glm)
# prob_glm_class = predict(model_glm, newdata = train, type = 'prob')
# pred_glm_class = predict(model_glm, newdata = train)
# confusionMatrix(table(pred_glm_class, train$class), positive = 'F', mode = 'everything')
# roc_glm <- roc.curve(pred_glm_class, train[['class']], plotit = T, main = 'ROC Curve using GLM')
# print(roc_glm)
# Test Evaluation
model_glm <- train(form = class ~ ., data = test, family = binomial(link = 'logit'), trControl = ctrl, method = 'glm', preProcess = c("center", "scale"))
prob_glm_class = predict(model_glm, newdata = test, type = 'prob')
pred_glm_class = predict(model_glm, newdata = test)
confusionMatrix(table(pred_glm_class, test$class), positive = 'F', mode = 'everything')
roc_glm <- roc.curve(pred_glm_class, test[['class']], plotit = T, main = 'ROC Curve using GLM')
print(roc_glm)

set.seed(642)
# Train Evaluation
# model_glm_under <- train(form = class ~ ., data = train, family = binomial(link = 'logit'), trControl = ctrl_under, method = 'glm', preProcess = c("center", "scale"))
# prob_glm_class_under = predict(model_glm_under, newdata = train, type = 'prob')
# pred_glm_class_under = predict(model_glm_under, newdata = train)
# confusionMatrix(table(pred_glm_class_under, train[['class']]), positive = 'F', mode = 'everything')
# roc_glm_under <- roc.curve(pred_glm_class_under, train[['class']], plotit = T, main = 'ROC Curve using GLM & Undersmaple')
# print(roc_glm_under)
# Test Evaluation
model_glm_under <- train(form = class ~ ., data = test, family = binomial(link = 'logit'), trControl = ctrl_under, method = 'glm', preProcess = c("center", "scale"))
prob_glm_class_under = predict(model_glm_under, newdata = test, type = 'prob')
pred_glm_class_under = predict(model_glm_under, newdata = test)
confusionMatrix(table(pred_glm_class_under, test[['class']]), positive = 'F', mode = 'everything')
roc_glm_under <- roc.curve(pred_glm_class_under, test[['class']], plotit = T, main = 'ROC Curve using GLM & Undersmaple')
print(roc_glm_under)

set.seed(642)
# Train Evaluation
# model_glm_over <- train(form = class ~ ., data = train, family = binomial(link = 'logit'), trControl = ctrl_over, method = 'glm', preProcess = c("center", "scale"))
# prob_glm_class_over = predict(model_glm_over, newdata = train, type = 'prob')
# pred_glm_class_over = predict(model_glm_over, newdata = train)
# confusionMatrix(table(pred_glm_class_over, train[['class']]), positive = 'F', mode = 'everything')
# roc_glm_over <- roc.curve(pred_glm_class_over, train[['class']], plotit = T, main = 'ROC Curve using GLM & Oversmaple')
# print(roc_glm_over)
# Test Evaluation
model_glm_over <- train(form = class ~ ., data = test, family = binomial(link = 'logit'), trControl = ctrl_over, method = 'glm', preProcess = c("center", "scale"))
prob_glm_class_over = predict(model_glm_over, newdata = test, type = 'prob')
pred_glm_class_over = predict(model_glm_over, newdata = test)
confusionMatrix(table(pred_glm_class_over, test[['class']]), positive = 'F', mode = 'everything')
roc_glm_over <- roc.curve(pred_glm_class_over, test[['class']], plotit = T, main = 'ROC Curve using GLM & Oversmaple')
print(roc_glm_over)

set.seed(642)
# Train Evaluation
# model_glm_ovun <- train(form = class ~ ., data = ovun_train_data, family = binomial(link = 'logit'), trControl = ctrl, method = 'glm', preProcess = c("center", "scale"))
# prob_glm_class_ovun = predict(model_glm_ovun, newdata = train, type = 'prob')
# pred_glm_class_ovun = predict(model_glm_ovun, newdata = train)
# confusionMatrix(as.factor(pred_glm_class_ovun), as.factor(train[['class']]), positive = 'F', mode = 'everything')
# roc_glm_under <- roc.curve(pred_glm_class_ovun, train[['class']], plotit = T, main = 'ROC Curve using GLM & Miexed-smaple')
# print(roc_glm_under)
# Test Evaluation
model_glm_ovun <- train(form = class ~ ., data = test, family = binomial(link = 'logit'), trControl = ctrl, method = 'glm', preProcess = c("center", "scale"))
prob_glm_class_ovun = predict(model_glm_ovun, newdata = test, type = 'prob')
pred_glm_class_ovun = predict(model_glm_ovun, newdata = test)
confusionMatrix(as.factor(pred_glm_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_glm_ovun <- roc.curve(pred_glm_class_ovun, test[['class']], plotit = T, main = 'ROC Curve using GLM & Miexed-smaple')
print(roc_glm_ovun)

set.seed(642)
# Train Evaluation
# model_glm_rose <- train(form = class ~ ., data = train, family = binomial(link = 'logit'), trControl = ctrl_rose, method = 'glm', preProcess = c("center", "scale"))
# prob_glm_class_rose = predict(model_glm_rose, newdata = train, type = 'prob')
# pred_glm_class_rose = predict(model_glm_rose, newdata = train)
# confusionMatrix(table(pred_glm_class_rose, train[['class']]), positive = 'F', mode = 'everything')
# roc_glm_rose <- roc.curve(pred_glm_class_rose, train[['class']], plotit = T, main = 'ROC Curve using GLM & ROSE')
# print(roc_glm_rose)
# Test Evaluation
model_glm_rose <- train(form = class ~ ., data = test, family = binomial(link = 'logit'), trControl = ctrl_rose, method = 'glm', preProcess = c("center", "scale"))
prob_glm_class_rose = predict(model_glm_rose, newdata = test, type = 'prob')
pred_glm_class_rose = predict(model_glm_rose, newdata = test)
confusionMatrix(table(pred_glm_class_rose, test[['class']]), positive = 'F', mode = 'everything')
roc_glm_rose <- roc.curve(pred_glm_class_rose, test[['class']], plotit = T, main = 'ROC Curve using GLM & ROSE')
print(roc_glm_rose)

set.seed(642)
# Train Evaluation
# model_glm_smote <- train(form = class ~ ., data = train, family = binomial(link = 'logit'), trControl = ctrl_smote, method = 'glm', preProcess = c("center", "scale"))
# prob_glm_class_smote = predict(model_glm_smote, newdata = train, type = 'prob')
# pred_glm_class_smote = predict(model_glm_smote, newdata = train)
# confusionMatrix(table(pred_glm_class_smote, train[['class']]), positive = 'F', mode = 'everything')
# roc_glm_smote <- roc.curve(pred_glm_class_smote, train[['class']], plotit = T, main = 'ROC Curve using GLM & SMOTE')
# print(roc_glm_smote)
# Test Evaluation
model_glm_smote <- train(form = class ~ ., data = test, family = binomial(link = 'logit'), trControl = ctrl_smote, method = 'glm', preProcess = c("center", "scale"))
prob_glm_class_smote = predict(model_glm_smote, newdata = test, type = 'prob')
pred_glm_class_smote = predict(model_glm_smote, newdata = test)
confusionMatrix(table(pred_glm_class_smote, test[['class']]), positive = 'F', mode = 'everything')
roc_glm_smote <- roc.curve(pred_glm_class_smote, test[['class']], plotit = T, main = 'ROC Curve using GLM & SMOTE')
print(roc_glm_smote)


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

set.seed(642)
# Train Evaluation
# model_dt <- train(form = class ~ ., data = train, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_dt)
# prob_dt_class <- predict(model_dt, newdata = train, type = 'prob')
# pred_dt_class = predict(model_dt, newdata = train)
# confusionMatrix(table(pred_dt_class, train[['class']]), positive = 'F', mode = 'everything')
# roc_dt <- roc.curve(pred_dt_class, train[['class']], plotit = T, main = 'ROC Curve using DT')
# print(roc_dt)
# Test Evaluation
model_dt <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_dt)
prob_dt_class <- predict(model_dt, newdata = test, type = 'prob')
pred_dt_class = predict(model_dt, newdata = test)
confusionMatrix(table(pred_dt_class, test[['class']]), positive = 'F', mode = 'everything')
roc_dt <- roc.curve(pred_dt_class, test[['class']], plotit = T, main = 'ROC Curve using DT')
print(roc_dt)

set.seed(642)
# Train Evaluation
# model_dt_under <- train(form = class ~ ., data = train, method = 'ctree', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
# plot(model_dt_under)
# prob_dt_class_under <- predict(model_dt_under, newdata = train, type = 'prob')
# pred_dt_class_under = predict(model_dt_under, newdata = train)
# confusionMatrix(table(pred_dt_class_under, train[['class']]), positive = 'F', mode = 'everything')
# roc_dt_under <- roc.curve(pred_dt_class_under, train[['class']], plotit = T, main = 'ROC Curve using DT & undersample')
# print(roc_dt_under)
# Test Evaluation
model_dt_under <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
plot(model_dt_under)
prob_dt_class_under <- predict(model_dt_under, newdata = test, type = 'prob')
pred_dt_class_under = predict(model_dt_under, newdata = test)
confusionMatrix(table(pred_dt_class_under, test[['class']]), positive = 'F', mode = 'everything')
roc_dt_under <- roc.curve(pred_dt_class_under, test[['class']], plotit = T, main = 'ROC Curve using DT & undersample')
print(roc_dt_under)

set.seed(642)
# Train Evaluation
# model_dt_over <- train(form = class ~ ., data = train, method = 'ctree', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
# plot(model_dt_over)
# prob_dt_class_over <- predict(model_dt_over, newdata = train, type = 'prob')
# pred_dt_class_over = predict(model_dt_over, newdata = train)
# confusionMatrix(table(pred_dt_class_over, train[['class']]), positive = 'F', mode = 'everything')
# roc_dt_under <- roc.curve(pred_dt_class_under, train[['class']], plotit = T, main = 'ROC Curve using DT & oversample')
# print(roc_dt_under)
# Test Evaluation
model_dt_over <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
plot(model_dt_over)
prob_dt_class_over <- predict(model_dt_over, newdata = test, type = 'prob')
pred_dt_class_over = predict(model_dt_over, newdata = test)
confusionMatrix(table(pred_dt_class_over, test[['class']]), positive = 'F', mode = 'everything')
roc_dt_over <- roc.curve(pred_dt_class_over, test[['class']], plotit = T, main = 'ROC Curve using DT & oversampling')
print(roc_dt_over)

set.seed(642)
# Train Evaluation
# model_dt_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_dt_ovun)
# prob_dt_class_ovun <- predict(model_dt_ovun, newdata = train, type = 'prob')
# pred_dt_class_ovun = predict(model_dt_ovun, newdata = train)
# confusionMatrix(as.factor(pred_dt_class_ovun), as.factor(train[['class']]), positive = 'F', mode = 'everything')
# roc_dt_ovun <- roc.curve(pred_dt_class_ovun, train[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_dt_ovun)
# Test Evaluation
model_dt_ovun <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_dt_ovun)
prob_dt_class_ovun <- predict(model_dt_ovun, newdata = test, type = 'prob')
pred_dt_class_ovun = predict(model_dt_ovun, newdata = test)
confusionMatrix(as.factor(pred_dt_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_dt_ovun <- roc.curve(pred_dt_class_ovun, test[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
print(roc_dt_ovun)

set.seed(642)
# Train Evaluation
# model_dt_rose <- train(form = class ~ ., data = train, method = 'ctree', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
# plot(model_dt_rose)
# prob_dt_class_rose <- predict(model_dt_rose, newdata = train, type = 'prob')
# pred_dt_class_rose = predict(model_dt_rose, newdata = train)
# confusionMatrix(table(pred_dt_class_rose, train[['class']]), positive = 'F', mode = 'everything')
# roc_dt_rose <- roc.curve(pred_dt_class_rose, train[['class']], plotit = T, main = 'ROC Curve using DT & ROSE')
# print(roc_dt_rose)
# Test Evaluation
model_dt_rose <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
plot(model_dt_rose)
prob_dt_class_rose <- predict(model_dt_rose, newdata = test, type = 'prob')
pred_dt_class_rose = predict(model_dt_rose, newdata = test)
confusionMatrix(table(pred_dt_class_rose, test[['class']]), positive = 'F', mode = 'everything')
roc_dt_rose <- roc.curve(pred_dt_class_rose, test[['class']], plotit = T, main = 'ROC Curve using DT & ROSE')
print(roc_dt_rose)

set.seed(642)
# Train Evaluation
# model_dt_smote <- train(form = class ~ ., data = train, method = 'ctree', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
# plot(model_dt_smote)
# prob_dt_class_smote <- predict(model_dt_smote, newdata = train, type = 'prob')
# pred_dt_class_smote = predict(model_dt_smote, newdata = train)
# confusionMatrix(table(pred_dt_class_smote, train[['class']]), positive = 'F', mode = 'everything')
# roc_dt_smote <- roc.curve(pred_dt_class_smote, train[['class']], plotit = T, main = 'ROC Curve using DT & SMOTE')
# print(roc_dt_smote)
# Test Evaluation
model_dt_smote <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
plot(model_dt_smote)
prob_dt_class_smote <- predict(model_dt_smote, newdata = test, type = 'prob')
pred_dt_class_smote = predict(model_dt_smote, newdata = test)
confusionMatrix(table(pred_dt_class_smote, test[['class']]), positive = 'F', mode = 'everything')
roc_dt_smote <- roc.curve(pred_dt_class_smote, test[['class']], plotit = T, main = 'ROC Curve using DT & SMOTE')
print(roc_dt_smote)


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

set.seed(642)
# Train Evaluation
# model_rf <- train(form = class ~ ., data = train, method = 'ranger', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_rf)
# #prob_rf_class <- predict(model_rf, newdata = train, type = 'prob')
# pred_rf_class = predict(model_rf, newdata = train)
# confusionMatrix(table(pred_rf_class, train[['class']]), positive = 'F', mode = 'everything')
# roc_rf <- roc.curve(pred_rf_class, train[['class']], plotit = T, main = 'ROC Curve using RF')
# print(roc_rf)
# Test Evaluation
model_rf <- train(form = class ~ ., data = test, method = 'ranger', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_rf)
pred_rf_class = predict(model_rf, newdata = test)
confusionMatrix(table(pred_rf_class, test[['class']]), positive = 'F', mode = 'everything')
roc_rf <- roc.curve(pred_rf_class, test[['class']], plotit = T, main = 'ROC Curve using RF')
print(roc_rf)

set.seed(642)
# Train Evaluation
# model_rf_under <- train(form = class ~ ., data = train, method = 'ranger', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
# plot(model_rf_under)
# pred_rf_class_under = predict(model_rf_under, newdata = train)
# confusionMatrix(table(pred_rf_class_under, train[['class']]), positive = 'F', mode = 'everything')
# roc_rf_under <- roc.curve(pred_rf_class_under, train[['class']], plotit = T, main = 'ROC Curve using RF & undersample')
# print(roc_rf_under)
# Test Evaluation
model_rf_under <- train(form = class ~ ., data = test, method = 'ranger', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
plot(model_rf_under)
pred_rf_class_under = predict(model_rf_under, newdata = test)
confusionMatrix(table(pred_rf_class_under, test[['class']]), positive = 'F', mode = 'everything')
roc_rf_under <- roc.curve(pred_rf_class_under, test[['class']], plotit = T, main = 'ROC Curve using RF & undersample')
print(roc_rf_under)

set.seed(642)
# Train Evaluation
# model_rf_over <- train(form = class ~ ., data = train, method = 'ranger', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
# plot(model_rf_over)
# pred_rf_class_over = predict(model_rf_over, newdata = train)
# confusionMatrix(table(pred_rf_class_over, train[['class']]), positive = 'F', mode = 'everything')
# roc_rf_over <- roc.curve(pred_rf_class_over, train[['class']], plotit = T, main = 'ROC Curve using RF & oversample')
# print(roc_rf_over)
# Test Evaluation
model_rf_over <- train(form = class ~ ., data = test, method = 'ranger', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
plot(model_rf_over)
pred_rf_class_over = predict(model_rf_over, newdata = test)
confusionMatrix(table(pred_rf_class_over, test[['class']]), positive = 'F', mode = 'everything')
roc_rf_over <- roc.curve(pred_rf_class_over, test[['class']], plotit = T, main = 'ROC Curve using RF & oversample')
print(roc_rf_over)

set.seed(642)
# Train Evaluation
# model_rf_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_rf_ovun)
# pred_rf_class_ovun = predict(model_rf_ovun, newdata = train)
# confusionMatrix(as.factor(pred_rf_class_ovun), as.factor(train[['class']]), positive = 'F', mode = 'everything')
# roc_rf_ovun <- roc.curve(pred_rf_class_ovun, train[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_rf_ovun)
# Test Evaluation
model_rf_ovun <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_rf_ovun)
pred_rf_class_ovun = predict(model_rf_ovun, newdata = test)
confusionMatrix(as.factor(pred_rf_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_rf_ovun <- roc.curve(pred_rf_class_ovun, test[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
print(roc_rf_ovun)

set.seed(642)
# Train Evaluation
# model_rf_rose <- train(form = class ~ ., data = train, method = 'ranger', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
# plot(model_rf_rose)
# pred_rf_class_rose = predict(model_rf_rose, newdata = train)
# confusionMatrix(table(pred_rf_class_rose, train[['class']]), positive = 'F', mode = 'everything')
# roc_rf_rose <- roc.curve(pred_rf_class_rose, train[['class']], plotit = T, main = 'ROC Curve using RF & ROSE')
# print(roc_rf_rose)
# Test Evaluation
model_rf_rose <- train(form = class ~ ., data = test, method = 'ranger', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
plot(model_rf_rose)
pred_rf_class_rose = predict(model_rf_rose, newdata = test)
confusionMatrix(table(pred_rf_class_rose, test[['class']]), positive = 'F', mode = 'everything')
roc_rf_rose <- roc.curve(pred_rf_class_rose, test[['class']], plotit = T, main = 'ROC Curve using RF & ROSE')
print(roc_rf_rose)

set.seed(642)
# Train Evaluation
# model_rf_smote <- train(form = class ~ ., data = train, method = 'ranger', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
# plot(model_rf_smote)
# pred_rf_class_smote = predict(model_rf_smote, newdata = train)
# confusionMatrix(table(pred_rf_class_smote, train[['class']]), positive = 'F', mode = 'everything')
# roc_rf_smote <- roc.curve(pred_rf_class_smote, train[['class']], plotit = T, main = 'ROC Curve using RF & SMOTE')
# print(roc_rf_smote)
# Test Evaluation
model_rf_smote <- train(form = class ~ ., data = test, method = 'ranger', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
plot(model_rf_smote)
pred_rf_class_smote = predict(model_rf_smote, newdata = test)
confusionMatrix(table(pred_rf_class_smote, test[['class']]), positive = 'F', mode = 'everything')
roc_rf_smote <- roc.curve(pred_rf_class_smote, test[['class']], plotit = T, main = 'ROC Curve using RF & SMOTE')
print(roc_rf_smote)


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

set.seed(642)
# Train Evaluation
# model_knn <- train(form = class ~ ., data = train, method = 'knn', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_knn)
# prob_knn_class <- predict(model_knn, newdata = train, type = 'prob')
# pred_knn_class = predict(model_knn, newdata = train)
# confusionMatrix(table(pred_knn_class, train[['class']]), positive = 'F', mode = 'everything')
# roc_knn <- roc.curve(pred_knn_class, train[['class']], plotit = T, main = 'ROC Curve using KNN')
# print(roc_knn)
# Test Evaluation
model_knn <- train(form = class ~ ., data = test, method = 'knn', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_knn)
prob_knn_class <- predict(model_knn, newdata = test, type = 'prob')
pred_knn_class = predict(model_knn, newdata = test)
confusionMatrix(table(pred_knn_class, test[['class']]), positive = 'F', mode = 'everything')
roc_knn <- roc.curve(pred_knn_class, test[['class']], plotit = T, main = 'ROC Curve using KNN')
print(roc_knn)

set.seed(642)
# Train Evaluation
# model_knn_under <- train(form = class ~ ., data = train, method = 'knn', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
# plot(model_knn_under)
# prob_knn_class_under <- predict(model_knn_under, newdata = train, type = 'prob')
# pred_knn_class_under = predict(model_knn_under, newdata = train)
# confusionMatrix(table(pred_knn_class_under, train[['class']]), positive = 'F', mode = 'everything')
# roc_knn_under <- roc.curve(pred_knn_class_under, train[['class']], plotit = T, main = 'ROC Curve using KNN and undersample')
# print(roc_knn_under)
# Test Evaluation
model_knn_under <- train(form = class ~ ., data = test, method = 'knn', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
plot(model_knn_under)
prob_knn_class_under <- predict(model_knn_under, newdata = test, type = 'prob')
pred_knn_class_under = predict(model_knn_under, newdata = test)
confusionMatrix(table(pred_knn_class_under, test[['class']]), positive = 'F', mode = 'everything')
roc_knn_under <- roc.curve(pred_knn_class_under, test[['class']], plotit = T, main = 'ROC Curve using KNN and undersample')
print(roc_knn_under)

set.seed(642)
# Train Evaluation
# model_knn_over <- train(form = class ~ ., data = train, method = 'knn', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
# plot(model_knn_over)
# prob_knn_class_over <- predict(model_knn_over, newdata = train, type = 'prob')
# pred_knn_class_over = predict(model_knn_over, newdata = train)
# confusionMatrix(table(pred_knn_class_over, train[['class']]), positive = 'F', mode = 'everything')
# roc_knn_over <- roc.curve(pred_knn_class_over, train[['class']], plotit = T, main = 'ROC Curve using KNN and oversample')
# print(roc_knn_over)
# Test Evaluation
model_knn_over <- train(form = class ~ ., data = test, method = 'knn', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
plot(model_knn_over)
prob_knn_class_over <- predict(model_knn_over, newdata = test, type = 'prob')
pred_knn_class_over = predict(model_knn_over, newdata = test)
confusionMatrix(table(pred_knn_class_over, test[['class']]), positive = 'F', mode = 'everything')
roc_knn_over <- roc.curve(pred_knn_class_over, test[['class']], plotit = T, main = 'ROC Curve using KNN and oversample')
print(roc_knn_over)

set.seed(642)
# Train Evaluation
# model_knn_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_knn_ovun)
# prob_knn_class_ovun <- predict(model_knn_ovun, newdata = train, type = 'prob')
# pred_knn_class_ovun = predict(model_knn_ovun, newdata = train)
# confusionMatrix(as.factor(pred_knn_class_ovun), as.factor(train[['class']]), positive = 'F', mode = 'everything')
# roc_knn_ovun <- roc.curve(pred_knn_class_ovun, train[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_knn_ovun)
# Test Evaluation
model_knn_ovun <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_knn_ovun)
prob_knn_class_ovun <- predict(model_knn_ovun, newdata = test, type = 'prob')
pred_knn_class_ovun = predict(model_knn_ovun, newdata = test)
confusionMatrix(as.factor(pred_knn_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_knn_ovun <- roc.curve(pred_knn_class_ovun, test[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
print(roc_knn_ovun)

set.seed(642)
# Train Evaluation
# model_knn_rose <- train(form = class ~ ., data = train, method = 'knn', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
# plot(model_knn_rose)
# prob_knn_class_rose <- predict(model_knn_rose, newdata = train, type = 'prob')
# pred_knn_class_rose = predict(model_knn_rose, newdata = train)
# confusionMatrix(table(pred_knn_class_rose, train[['class']]), positive = 'F', mode = 'everything')
# roc_knn_rose <- roc.curve(pred_knn_class_rose, train[['class']], plotit = T, main = 'ROC Curve using KNN and ROSE')
# print(roc_knn_rose)
# Test Evaluation
model_knn_rose <- train(form = class ~ ., data = test, method = 'knn', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
plot(model_knn_rose)
prob_knn_class_rose <- predict(model_knn_rose, newdata = test, type = 'prob')
pred_knn_class_rose = predict(model_knn_rose, newdata = test)
confusionMatrix(table(pred_knn_class_rose, test[['class']]), positive = 'F', mode = 'everything')
roc_knn_rose <- roc.curve(pred_knn_class_rose, test[['class']], plotit = T, main = 'ROC Curve using KNN and ROSE')
print(roc_knn_rose)

set.seed(642)
# Train Evaluation
# model_knn_smote <- train(form = class ~ ., data = train, method = 'knn', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
# plot(model_knn_smote)
# prob_knn_class_smote <- predict(model_knn_smote, newdata = train, type = 'prob')
# pred_knn_class_smote = predict(model_knn_smote, newdata = train)
# confusionMatrix(table(pred_knn_class_smote, train[['class']]), positive = 'F', mode = 'everything')
# roc_knn_smote <- roc.curve(pred_knn_class_smote, train[['class']], plotit = T, main = 'ROC Curve using KNN and SMOTE')
# print(roc_knn_smote)
# Test Evaluation
model_knn_smote <- train(form = class ~ ., data = test, method = 'knn', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
plot(model_knn_smote)
prob_knn_class_smote <- predict(model_knn_smote, newdata = test, type = 'prob')
pred_knn_class_smote = predict(model_knn_smote, newdata = test)
confusionMatrix(table(pred_knn_class_smote, test[['class']]), positive = 'F', mode = 'everything')
roc_knn_smote <- roc.curve(pred_knn_class_smote, test[['class']], plotit = T, main = 'ROC Curve using KNN and SMOTE')
print(roc_knn_smote)


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

set.seed(642)
# Train Evaluation
# model_nb <- train(form = class ~ ., data = train, method = 'naive_bayes', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_nb)
# prob_nb_class <- predict(model_nb, newdata = train, type = 'prob')
# pred_nb_class = predict(model_nb, newdata = train)
# confusionMatrix(table(pred_nb_class, train[['class']]), positive = 'F', mode = 'everything')
# roc_nb <- roc.curve(pred_nb_class, train[['class']], plotit = T, main = 'ROC Curve using NB')
# print(roc_nb)
# Test Evaluation
model_nb <- train(form = class ~ ., data = test, method = 'naive_bayes', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_nb)
prob_nb_class <- predict(model_nb, newdata = test, type = 'prob')
pred_nb_class = predict(model_nb, newdata = test)
confusionMatrix(table(pred_nb_class, test[['class']]), positive = 'F', mode = 'everything')
roc_nb <- roc.curve(pred_nb_class, test[['class']], plotit = T, main = 'ROC Curve using NB')
print(roc_nb)

set.seed(642)
# Train Evaluation
# model_nb_under <- train(form = class ~ ., data = train, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
# plot(model_nb_under)
# prob_nb_class_under <- predict(model_nb_under, newdata = train, type = 'prob')
# pred_nb_class_under = predict(model_nb_under, newdata = train)
# confusionMatrix(table(pred_nb_class_under, train[['class']]), positive = 'F', mode = 'everything')
# roc_nb_under <- roc.curve(pred_nb_class_under, train[['class']], plotit = T, main = 'ROC Curve using NB & undersample')
# print(roc_nb_under)
# Test Evaluation
model_nb_under <- train(form = class ~ ., data = test, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
plot(model_nb_under)
prob_nb_class_under <- predict(model_nb_under, newdata = test, type = 'prob')
pred_nb_class_under = predict(model_nb_under, newdata = test)
confusionMatrix(table(pred_nb_class_under, test[['class']]), positive = 'F', mode = 'everything')
roc_nb_under <- roc.curve(pred_nb_class_under, test[['class']], plotit = T, main = 'ROC Curve using NB & undersample')
print(roc_nb_under)

set.seed(642)
# Train Evaluation
# model_nb_over <- train(form = class ~ ., data = train, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
# plot(model_nb_over)
# prob_nb_class_over <- predict(model_nb_over, newdata = train, type = 'prob')
# pred_nb_class_over = predict(model_nb_over, newdata = train)
# confusionMatrix(table(pred_nb_class_over, train[['class']]), positive = 'F', mode = 'everything')
# roc_nb_over <- roc.curve(pred_nb_class_over, train[['class']], plotit = T, main = 'ROC Curve using NB & oversample')
# print(roc_nb_over)
# Test Evaluation
model_nb_over <- train(form = class ~ ., data = test, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
plot(model_nb_over)
prob_nb_class_over <- predict(model_nb_over, newdata = test, type = 'prob')
pred_nb_class_over = predict(model_nb_over, newdata = test)
confusionMatrix(table(pred_nb_class_over, test[['class']]), positive = 'F', mode = 'everything')
roc_nb_over <- roc.curve(pred_nb_class_over, test[['class']], plotit = T, main = 'ROC Curve using NB & oversample')
print(roc_nb_over)

set.seed(642)
# Train Evaluation
# model_nb_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_nb_ovun)
# prob_nb_class_ovun <- predict(model_nb_ovun, newdata = train, type = 'prob')
# pred_nb_class_ovun = predict(model_nb_ovun, newdata = train)
# confusionMatrix(as.factor(pred_nb_class_ovun), as.factor(train[['class']]), positive = 'F', mode = 'everything')
# roc_nb_ovun <- roc.curve(pred_nb_class_ovun, train[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_nb_ovun)
# Test Evaluation
model_nb_ovun <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_nb_ovun)
prob_nb_class_ovun <- predict(model_nb_ovun, newdata = test, type = 'prob')
pred_nb_class_ovun = predict(model_nb_ovun, newdata = test)
confusionMatrix(as.factor(pred_nb_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_nb_ovun <- roc.curve(pred_nb_class_ovun, test[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
print(roc_nb_ovun)

set.seed(642)
# Train Evaluation
# model_nb_rose <- train(form = class ~ ., data = train, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
# plot(model_nb_rose)
# prob_nb_class_rose <- predict(model_nb_rose, newdata = train, type = 'prob')
# pred_nb_class_rose = predict(model_nb_rose, newdata = train)
# confusionMatrix(table(pred_nb_class_rose, train[['class']]), positive = 'F', mode = 'everything')
# roc_nb_rose <- roc.curve(pred_nb_class_rose, train[['class']], plotit = T, main = 'ROC Curve using NB & ROSE')
# print(roc_nb_rose)
# Test Evaluation
model_nb_rose <- train(form = class ~ ., data = test, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
plot(model_nb_rose)
prob_nb_class_rose <- predict(model_nb_rose, newdata = test, type = 'prob')
pred_nb_class_rose = predict(model_nb_rose, newdata = test)
confusionMatrix(table(pred_nb_class_rose, test[['class']]), positive = 'F', mode = 'everything')
roc_nb_rose <- roc.curve(pred_nb_class_rose, test[['class']], plotit = T, main = 'ROC Curve using NB & ROSE')
print(roc_nb_rose)

set.seed(642)
# Train Evaluation
# model_nb_smote <- train(form = class ~ ., data = train, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
# plot(model_nb_smote)
# prob_nb_class_smote <- predict(model_nb_smote, newdata = train, type = 'prob')
# pred_nb_class_smote = predict(model_nb_smote, newdata = train)
# confusionMatrix(table(pred_nb_class_smote, train[['class']]), positive = 'F', mode = 'everything')
# roc_nb_smote <- roc.curve(pred_nb_class_smote, train[['class']], plotit = T, main = 'ROC Curve using NB & SMOTE')
# print(roc_nb_smote)
# Test Evaluation
model_nb_smote <- train(form = class ~ ., data = test, method = 'naive_bayes', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
plot(model_nb_smote)
prob_nb_class_smote <- predict(model_nb_smote, newdata = test, type = 'prob')
pred_nb_class_smote = predict(model_nb_smote, newdata = test)
confusionMatrix(table(pred_nb_class_smote, test[['class']]), positive = 'F', mode = 'everything')
roc_nb_smote <- roc.curve(pred_nb_class_smote, test[['class']], plotit = T, main = 'ROC Curve using NB & SMOTE')
print(roc_nb_smote)


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

set.seed(642)
# Train Evaluation
# model_svm <- train(form = class ~ ., data = train, method = 'svmRadial', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_svm)
# prob_svm_class <- predict(model_svm, newdata = train, type = 'prob')
# pred_svm_class = predict(model_svm, newdata = train)
# confusionMatrix(table(pred_svm_class, train[['class']]), positive = 'F', mode = 'everything')
# roc_svm <- roc.curve(pred_svm_class, train[['class']], plotit = T, main = 'ROC Curve using SVM')
# print(roc_svm)
# Test Evaluation
model_svm <- train(form = class ~ ., data = test, method = 'svmRadial', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_svm)
prob_svm_class <- predict(model_svm, newdata = test, type = 'prob')
pred_svm_class = predict(model_svm, newdata = test)
confusionMatrix(table(pred_svm_class, test[['class']]), positive = 'F', mode = 'everything')
roc_svm <- roc.curve(pred_svm_class, test[['class']], plotit = T, main = 'ROC Curve using SVM')
print(roc_svm)

set.seed(642)
# Train Evaluation
# model_svm_under <- train(form = class ~ ., data = train, method = 'svmRadial', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
# plot(model_svm_under)
# prob_svm_class_under <- predict(model_svm_under, newdata = train, type = 'prob')
# pred_svm_class_under = predict(model_svm_under, newdata = train)
# confusionMatrix(table(pred_svm_class_under, train[['class']]), positive = 'F', mode = 'everything')
# roc_svm_under <- roc.curve(pred_svm_class_under, train[['class']], plotit = T, main = 'ROC Curve using SVM & undersample')
# print(roc_svm_under)
# Test Evaluation
model_svm_under <- train(form = class ~ ., data = test, method = 'svmRadial', tuneLength = 5, trControl = ctrl_under, preProcess = c("center", "scale"))
plot(model_svm_under)
prob_svm_class_under <- predict(model_svm_under, newdata = test, type = 'prob')
pred_svm_class_under = predict(model_svm_under, newdata = test)
confusionMatrix(table(pred_svm_class_under, test[['class']]), positive = 'F', mode = 'everything')
roc_svm_under <- roc.curve(pred_svm_class_under, test[['class']], plotit = T, main = 'ROC Curve using SVM & undersample')
print(roc_svm_under)

set.seed(642)
# Train Evaluation
# model_svm_over <- train(form = class ~ ., data = train, method = 'svmRadial', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
# plot(model_svm_over)
# prob_svm_class_over <- predict(model_svm_over, newdata = train, type = 'prob')
# pred_svm_class_over = predict(model_svm_over, newdata = train)
# confusionMatrix(table(pred_svm_class_over, train[['class']]), positive = 'F', mode = 'everything')
# roc_svm_over <- roc.curve(pred_svm_class_over, train[['class']], plotit = T, main = 'ROC Curve using SVM & oversample')
# print(roc_svm_over)
# Test Evaluation
model_svm_over <- train(form = class ~ ., data = test, method = 'svmRadial', tuneLength = 5, trControl = ctrl_over, preProcess = c("center", "scale"))
plot(model_svm_over)
prob_svm_class_over <- predict(model_svm_over, newdata = test, type = 'prob')
pred_svm_class_over = predict(model_svm_over, newdata = test)
confusionMatrix(table(pred_svm_class_over, test[['class']]), positive = 'F', mode = 'everything')
roc_svm_over <- roc.curve(pred_svm_class_over, test[['class']], plotit = T, main = 'ROC Curve using SVM & oversample')
print(roc_svm_over)

set.seed(642)
# Train Evaluation
# model_svm_ovun <- train(form = class ~ ., data = ovun_train_data, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
# plot(model_svm_ovun)
# prob_svm_class_ovun <- predict(model_svm_ovun, newdata = train, type = 'prob')
# pred_svm_class_ovun = predict(model_svm_ovun, newdata = train)
# confusionMatrix(as.factor(pred_svm_class_ovun), as.factor(train[['class']]), positive = 'F', mode = 'everything')
# roc_svm_ovun <- roc.curve(pred_svm_class_ovun, train[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
# print(roc_svm_ovun)
# Test Evaluation
model_svm_ovun <- train(form = class ~ ., data = test, method = 'ctree', tuneLength = 5, trControl = ctrl, preProcess = c("center", "scale"))
plot(model_svm_ovun)
prob_svm_class_ovun <- predict(model_svm_ovun, newdata = test, type = 'prob')
pred_svm_class_ovun = predict(model_svm_ovun, newdata = test)
confusionMatrix(as.factor(pred_svm_class_ovun), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_svm_under <- roc.curve(pred_svm_class_ovun, test[['class']], plotit = T, main = 'ROC Curve using DT & mixed-sample')
print(roc_svm_ovun)

set.seed(642)
# Train Evaluation
# model_svm_rose <- train(form = class ~ ., data = train, method = 'svmRadial', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
# plot(model_svm_rose)
# prob_svm_class_rose <- predict(model_svm_rose, newdata = train, type = 'prob')
# pred_svm_class_rose = predict(model_svm_rose, newdata = train)
# confusionMatrix(table(pred_svm_class_rose, train[['class']]), positive = 'F', mode = 'everything')
# roc_svm_rose <- roc.curve(pred_svm_class_rose, train[['class']], plotit = T, main = 'ROC Curve using SVM & ROSE')
# print(roc_svm_rose)
# Test Evaluation
model_svm_rose <- train(form = class ~ ., data = test, method = 'svmRadial', tuneLength = 5, trControl = ctrl_rose, preProcess = c("center", "scale"))
plot(model_svm_rose)
prob_svm_class_rose <- predict(model_svm_rose, newdata = test, type = 'prob')
pred_svm_class_rose = predict(model_svm_rose, newdata = test)
confusionMatrix(table(pred_svm_class_rose, test[['class']]), positive = 'F', mode = 'everything')
roc_svm_rose <- roc.curve(pred_svm_class_rose, test[['class']], plotit = T, main = 'ROC Curve using SVM & ROSE')
print(roc_svm_rose)

set.seed(642)
# Train Evaluation
# model_svm_smote <- train(form = class ~ ., data = train, method = 'svmRadial', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
# plot(model_svm_smote)
# prob_svm_class_smote <- predict(model_svm_smote, newdata = train, type = 'prob')
# pred_svm_class_smote = predict(model_svm_smote, newdata = train)
# confusionMatrix(table(pred_svm_class_smote, train[['class']]), positive = 'F', mode = 'everything')
# roc_svm_smote <- roc.curve(pred_svm_class_smote, train[['class']], plotit = T, main = 'ROC Curve using SVM & SMOTE')
# print(roc_svm_smote)
# Test Evaluation
model_svm_smote <- train(form = class ~ ., data = test, method = 'svmRadial', tuneLength = 5, trControl = ctrl_smote, preProcess = c("center", "scale"))
plot(model_svm_smote)
prob_svm_class_smote <- predict(model_svm_smote, newdata = test, type = 'prob')
pred_svm_class_smote = predict(model_svm_smote, newdata = test)
confusionMatrix(table(pred_svm_class_smote, test[['class']]), positive = 'F', mode = 'everything')
roc_svm_smote <- roc.curve(pred_svm_class_smote, test[['class']], plotit = T, main = 'ROC Curve using SVM & SMOTE')
print(roc_svm_smote)


# ============================================================================= #
# MODEL COMPARISON FINAL
# ============================================================================= #

# Compare model performances using resample()
models_compare <- resamples(list(GLM = model_glm, DT = model_dt, RF = model_rf, KNN = model_knn, NB = model_nb, MARS = model_mars, NN = model_nn, SVM = model_svm))
# Summary of the models performances
summary(models_compare)
# Draw box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales)
