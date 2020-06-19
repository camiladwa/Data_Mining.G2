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
# ============================================================================= #

load_lib <- c('grDevices', 'foreign', 'caret', 'DMwR', 'mice', 'Boruta', 'randomForest', 'ggplot2', 'pacman', 'missMDA', 'caTools', 'party', 'earth', 'ROSE', 'kernlab')
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
  barplot(tN, col = c('cyan', 'red'), ylim = c(0, 100), ylab = '% of observation', xlab = 'Class of Semiconductor', names.arg = c('Non-Faulty', 'Faulty'))
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
  set.seed(666)
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
  impute_KNN <- knnImputation(df, k = 5, scale = T, meth = "weighAvg", distData = NULL)
  return(impute_KNN)
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
  set.seed(777)
  train_boruta_features <- Boruta(class ~ ., data = boruta_df, doTrace = 2, ntree = 500, maxRuns = 150)
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
analyse_data(semcon_original_data)
semcon_original_data <- semcon_original_data[, -c(1, 3)]
class_distribution(semcon_original_data, 'Given Sample')
missing_value_analysis(semcon_original_data, 'Given Sample')
outlier_analysis(semcon_original_data[, -c(1)], 'Given Sample')


# ============================================================================= #
# TRAIN AND TEST SPLIT INTO 80:20
# ============================================================================= #

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
# write_csv(semcon_test_data, "semcon_test_data.csv")


# ============================================================================= #
# FEATURE REDUCTION
# 01. Remove Features with more than 60% NA based on missing_value_analysis
# ============================================================================= #

class <- as.factor(semcon_train_data$class)
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

# bagImpute <- preProcess(train_outlier_NA, method="bagImpute", verbose = TRUE)
# train_bag_Impute <- predict(bagImpute, newdata = train_outlier_NA)
# missing_value_analysis(train_bag_Impute, 'train set after Bagged Tree imputation')
# outlier_analysis(train_bag_Impute, 'train set after Bagged Tree imputation')
# write_csv(cbind(class, train_knn_imputation), "train_NAH_knn.csv")


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

# train_FR_boruta_bt <- selection_boruta(train_bag_Impute, class)
# print("Summary of selected features")
# summary(train_FR_boruta_bt)
# message("Variance corresponding to selected features:")
# sapply(train_FR_boruta_bt, var)
# missing_value_analysis(train_FR_boruta_bt, 'train set after BORUTA')
# outlier_analysis(train_FR_boruta_knn, 'train set after BORUTA')
# write_csv(cbind(train_FR_boruta_knn), "train_FR_boruta.csv")


# ============================================================================= #
# FEATURE SELECTION AND REDUCTION
# 04. Recursive Feature Elimination (RFE)
# ============================================================================= #

# set.seed(199)
# options(warn = -1)
# subsets <- c(15:25)
# ctrl <- rfeControl(functions = rfFuncs, method = 'repeatedcv', repeats = 5, verbose = TRUE)
# lmProfile <- rfe(x = train_knn_imputation, y = class, sizes = subsets, rfeControl = ctrl)
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

ctrl <- trainControl(method = 'repeatedcv', number = 10, repeats = 10)
# ctrl <- trainControl(method = "cv", number = 5, savePredictions = 'final',  summaryFunction = twoClassSummary)
# ctrl <- trainControl(method = "boot632", number = 1000, savePredictions = TRUE, savePredictions = 'final', classProbs = T, summaryFunction = twoClassSummary)
# ctrl <- trainControl(method = "repeatedcv", number = 5, savePredictions = TRUE, savePredictions = 'final', classProbs = T, summaryFunction = twoClassSummary)
ctrl_under <- trainControl(method = 'repeatedcv', number = 10, repeats = 10, verboseIter = FALSE, sampling = 'down')
ctrl_over <- trainControl(method = 'repeatedcv', number = 10, repeats = 10, verboseIter = FALSE, sampling = 'up')
ctrl_rose <- trainControl(method = 'repeatedcv', number = 10, repeats = 10, verboseIter = FALSE, sampling = 'rose')
ctrl_smote <- trainControl(method = 'repeatedcv', number = 10, repeats = 10, verboseIter = FALSE, sampling = 'smote')


# ============================================================================= #
# CLASSIFICATION MODEL
# 01. Generalized Linear Model
# ============================================================================= #

set.seed(892)
model_glm <- train(form = class ~ ., data = train, family = binomial(link = 'logit'), trControl = ctrl, method = 'glm')
#exp(coef(model_glm$finalModel))
#varImp(model_glm)
prob_glm_class = predict(model_glm, newdata = test, type = 'prob')
pred_glm_class = predict(model_glm, newdata = test)
confusionMatrix(table(pred_glm_class, test[['class']]), positive = '1', mode = 'everything')
colAUC(prob_glm_class, test[['class']], plotROC = TRUE)


# ============================================================================= #
# CLASSIFICATION MODEL
# 02. Decision Tree - ctree, chaid, C5.0, xgbTree
# ============================================================================= #

set.seed(892)
model_dt <- train(form = class ~ ., data = train, method = 'ctree', tuneLength = 5, trControl = ctrl)
plot(model_dt)
prob_dt_class <- predict(model_dt, newdata = test, type = 'prob')
pred_dt_class = predict(model_dt, newdata = test)
confusionMatrix(table(pred_dt_class, test[['class']]), positive = '1', mode = 'everything')
colAUC(prob_dt_class, test[['class']], plotROC = TRUE)


# ============================================================================= #
# CLASSIFICATION MODEL
# 03. Random Forest - ranger, rf
# ============================================================================= #

set.seed(892)
model_rf <- train(form = class ~ ., data = train, method = 'ranger', tuneLength = 5, trControl = ctrl)
plot(model_rf)
#prob_rf_class <- predict(model_rf, newdata = test, type = 'prob')
pred_rf_class = predict(model_rf, newdata = test)
confusionMatrix(table(pred_rf_class, test[['class']]), positive = '1', mode = 'everything')
# colAUC(prob_rf_class, train[['class']], plotROC = TRUE)


# ============================================================================= #
# CLASSIFICATION MODEL
# 04. KNN
# ============================================================================= #

model_knn <- train(form = class ~ ., data = train, method = 'knn', tuneLength = 5, trControl = ctrl)
plot(model_knn)
prob_knn_class <- predict(model_knn, newdata = test, type = 'prob')
pred_knn_class = predict(model_knn, newdata = test)
confusionMatrix(table(pred_knn_class, test[['class']]), positive = '1', mode = 'everything')
colAUC(prob_knn_class, test[['class']], plotROC = TRUE)


# ============================================================================= #
# CLASSIFICATION MODEL
# 05. Naive Bayes
# ============================================================================= #

model_nb <- train(form = class ~ ., data = train, method = 'naive_bayes', tuneLength = 5, trControl = ctrl)
plot(model_nb)
prob_nb_class <- predict(model_nb, newdata = test, type = 'prob')
pred_nb_class = predict(model_nb, newdata = test)
confusionMatrix(table(pred_nb_class, test[['class']]), positive = '1', mode = 'everything')
colAUC(prob_nb_class, test[['class']], plotROC = TRUE)


# ============================================================================= #
# CLASSIFICATION MODEL
# 06. Multivariate Adaptive Regression Splines (MARS)
# ============================================================================= #

set.seed(40)
model_mars <- train(form = class ~ ., data = train, method = 'earth', tuneLength = 5, trControl = ctrl)
plot(model_mars)
prob_mars_class <- predict(model_mars, newdata = test, type = 'prob')
pred_mars_class = predict(model_mars, newdata = test)
confusionMatrix(table(pred_mars_class, test[['class']]), positive = '1', mode = 'everything')
colAUC(prob_mars_class, test[['class']], plotROC = TRUE)


# ============================================================================= #
# CLASSIFICATION MODEL
# 07. Neural Network: Penalized Multinomial Regression
# ============================================================================= #

model_nn <- train(form = class ~ ., data = train, method = 'multinom', tuneLength = 5, trControl = ctrl)
plot(model_nn)
prob_nn_class <- predict(model_nn, newdata = test, type = 'prob')
pred_nn_class = predict(model_nn, newdata = test)
confusionMatrix(table(pred_nn_class, test[['class']]), positive = '1', mode = 'everything')
colAUC(prob_nn_class, test[['class']], plotROC = TRUE)


# ============================================================================= #
# CLASSIFICATION MODEL
# 08. SVM
# ============================================================================= #

model_svm <- train(form = class ~ ., data = train, method = 'svmRadial', tuneLength = 5, trControl = ctrl)
plot(model_svm)
prob_svm_class <- predict(model_svm, newdata = test, type = 'prob')
pred_svm_class = predict(model_svm, newdata = test)
confusionMatrix(table(pred_svm_class, test[['class']]), positive = '1', mode = 'everything')
colAUC(prob_svm_class, test[['class']], plotROC = TRUE)


# ============================================================================= #
# CLASSIFICATION MODEL
# 09. Adaboost
# ============================================================================= #

model_adb <- train(form = class ~ ., data = train, method = 'adaboost', tuneLength = 5, trControl = ctrl)
plot(model_adb)
prob_adb_class <- predict(model_adb, newdata = test, type = 'prob')
pred_adb_class = predict(model_adb, newdata = test)
confusionMatrix(table(pred_adb_class, test[['class']]), positive = '1', mode = 'everything')
colAUC(prob_adb_class, test[['class']], plotROC = TRUE)


# ============================================================================= #
# CLASSIFICATION MODEL
# 10. xgBoost Dart
# ============================================================================= #

model_xgb <- train(form = class ~ ., data = train, method = 'xgbDART', tuneLength = 5, trControl = ctrl)
plot(model_xgb)
prob_xgb_class <- predict(model_xgb, newdata = test, type = 'prob')
pred_xgb_class = predict(model_xgb, newdata = test)
confusionMatrix(table(pred_xgb_class, test[['class']]), positive = '1', mode = 'everything')
colAUC(prob_xgb_class, test[['class']], plotROC = TRUE)


# ============================================================================= #
# MODEL COMPARISON
# ============================================================================= #

# Compare model performances using resample()
models_compare <- resamples(list(GLM = model_glm, DT = model_dt, RF = model_rf, KNN = model_knn, NB = model_nb, MARS = model_mars, NN = model_nn, SVM = model_svm, ADABOOST = model_adb))
# Summary of the models performances
summary(models_compare)
# Draw box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales)
