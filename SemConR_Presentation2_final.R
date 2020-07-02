# ============================================================================= #
# INSTALL & LOAD REQUIRED LIBRARIES
# 01. foreign: Reading SPSS Source
# 02. caret: nearZeroVar
# 03. DMwR: knnImputation
# 04. Boruta
# 05. randomForest
# 06. ggplot2
# 07. plotly
# 08. ROSE
# 09. summarytools
# ============================================================================= #

load_lib <- c('foreign', 'caret', 'DMwR', 'Boruta', 'randomForest', 'ggplot2', 'plotly', 'ROSE', 'summarytools')
install_lib <- load_lib[!load_lib %in% installed.packages()]
for(lib in install_lib) install.packages(lib, dependencies = TRUE)
sapply(load_lib, require, character=TRUE)

# ============================================================================= #
# FILE SETUP
# - SRC: Source file location with file name
# ============================================================================= #

SRC <- 'C:/Users/YDS/Desktop/R/SemCom/SemCon/Source/secom_mod.SAV'


# ============================================================================= #
# FUNCTION DEFINITIONS (Repeatedly used block of code lines are converted into functions)
# 01. class_distribution
# 02. missing_value_analysis
# 03. outliers_detection
# 04. outlier_analysis
# 05. impute_outlier_3s
# 06. impute_outlier_NA
# 07. KNN_imputation
# ============================================================================= #

class_distribution <- function(df, framename) {
  # Function to create barplot and visualize the distribution of class in the df
  tN <- table(df$class) * 100 / nrow(df)
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

KNN_imputation <- function(df, train_df=NULL) {
  # Function to impute missing values(NAs) using K-Nearest Neighbour (k = 5)
  impute_KNN <- knnImputation(df, k = 5, scale = T, meth = "weighAvg", distData = train_df)
  return(impute_KNN)
}


# ============================================================================= #
# DATA LOAD
# - Load data from Source file into the dataframe
# - EDA on original dataset
# ============================================================================= #

semcon_original_data <- read.spss(SRC, to.data.frame = TRUE)
semcon_original_data <- semcon_original_data[, -c(1, 3)] 				# Remove ID and timestamp columns
semcon_original_data$class[semcon_original_data$class == 1] <- 'F'		# Define 1 as F (Faulty)
semcon_original_data$class[semcon_original_data$class == 0] <- 'NF' 	# Define 0 as NF (Non-faulty)
dfSummary(semcon_original_data$class)
class_distribution(semcon_original_data, 'Given Sample')
missing_value_analysis(semcon_original_data, 'Given Sample')
outlier_analysis(semcon_original_data[, -c(1)], 'Given Sample')			# checking outliers after removing class variable


# ============================================================================= #
# TRAIN AND TEST SPLIT
# ============================================================================= #

TRAIN_TEST_RATIO <- 0.8
set.seed(666)
split <- sample(1:nrow(semcon_original_data), nrow(semcon_original_data) * TRAIN_TEST_RATIO)
semcon_train_data <- semcon_original_data[split, ]
semcon_test_data <- semcon_original_data[-split, ]


# ============================================================================= #
# TRAIN DATA ANALYSIS
# ============================================================================= #

dfSummary(semcon_train_data$class)
class_distribution(semcon_train_data, 'Given Sample')
missing_value_analysis(semcon_train_data, 'Given Sample')
outlier_analysis(semcon_train_data[, -c(1)], 'Given Sample')			# checking outliers after removing class variable


# ============================================================================= #
# TEST DATA ANALYSIS
# ============================================================================= #

dfSummary(semcon_test_data$class)
class_distribution(semcon_test_data, 'Test Data')
missing_value_analysis(semcon_test_data, 'Test Data')
outlier_analysis(semcon_test_data[, -c(1)], 'Test Data')				# checking outliers after removing class variable


# ============================================================================= #
# FEATURE REDUCTION
# - Remove Features with more than 55% NA based on missing_value_analysis
# ============================================================================= #

NULL_THRESHOLD <- 0.55
class <- semcon_train_data$class
semcon_train_data <- semcon_train_data[, -c(1)]
train_null_removal = semcon_train_data[, which(colMeans(!is.na(semcon_train_data)) >= NULL_THRESHOLD)]
cat("Number of feature removed after NA removal with threshold ",NULL_THRESHOLD*100,"% is :", length(semcon_train_data)-length(train_null_removal)-1)


# ============================================================================= #
# FEATURE REDUCTION
# - Near Zero Variance Removal
# ============================================================================= #

train_variance_removal <- train_null_removal[-nearZeroVar(train_null_removal)]
cat("Number of feature removed after near zero variance is :", length(train_null_removal)-length(train_variance_removal))


# ============================================================================= #
# OUTLIER HANDLING
# - Replace by 3s boundary
# 	- EDA after feature reduction and outlier handling
# ============================================================================= #

train_outlier_3s <- apply(train_variance_removal, 2, impute_outlier_3s)
train_outlier_3s <- as.data.frame(train_outlier_3s)
missing_value_analysis(train_outlier_3s, 'train set after imputating outlier with 3s boundary')
outlier_analysis(train_outlier_3s, 'train set after imputating outlier with 3s boundary')


# ============================================================================= #
# OUTLIER HANDLING
# - Replace by NAs
# 	- EDA after feature reduction and outlier handling
# ============================================================================= #

train_outlier_NA <- apply(train_variance_removal, 2, impute_outlier_NA)
train_outlier_NA <- as.data.frame(train_outlier_NA)
missing_value_analysis(train_outlier_NA, 'train set after imputating outlier with NA')
outlier_analysis(train_outlier_NA, 'train set after imputating outlier with NA')


# ============================================================================= #
# NA HANDLING
# - KNN imputation
# ============================================================================= #

train_knn_imputation <- KNN_imputation(train_outlier_NA)
missing_value_analysis(train_knn_imputation, 'train set after KNN imputation')


# ============================================================================= #
# FEATURE SELECTION AND REDUCTION
# - BORUTA on KNN imputed train set
# ============================================================================= #

set.seed(111)
# Perform feature selection using BORUTA
boruta_df <- cbind(class, train_knn_imputation)
train_boruta_features <- Boruta(class ~ ., data = boruta_df, doTrace = 2, ntree = 500, maxRuns = 100)
boruta_mF <- TentativeRoughFix(train_boruta_features)
print("Confirmed formula after Boruta: ")
message(getConfirmedFormula(train_boruta_features))
plot(train_boruta_features, las = 2, cex.axis = 0.5)
plotImpHistory(train_boruta_features)
print(attStats(train_boruta_features))
boruta_mF_confirmed <- names(boruta_mF$finalDecision[boruta_mF$finalDecision %in% c('Confirmed')])
cat("Number of Confirmed Features after Boruta: ", length(boruta_mF_confirmed))
cat(" Confirmed Features after Boruta: ", boruta_mF_confirmed)
train_FR_boruta_knn <- train_knn_imputation[,(names(train_knn_imputation) %in% boruta_mF_confirmed)]
train_FR_boruta_knn <- cbind(class, train_FR_boruta_knn)
summarytools::view(dfSummary(train_FR_boruta_knn))


# ============================================================================= #
# BALANCING METHODS
# - Imbalanced
# - ROSE
# - SMOTE
# ============================================================================= #

imba_train_data <- train_FR_boruta_knn
print(table(imba_train_data$class))
dfSummary(imba_train_data$class)
# F: 71, NF: 1182
rose_train_data <- ROSE(class ~ ., data = imba_train_data, N = dim(imba_train_data)[1] * 3, p = 0.5, seed = 1)$data
print(table(rose_train_data$class))
dfSummary(rose_train_data$class)
# F: 1935, NF: 1824
smote_train_data <- SMOTE(class ~., data = imba_train_data, perc.over = 2000, perc.under = 100)
print(table(smote_train_data$class))
dfSummary(smote_train_data$class)
# F: 1491, NF: 1420


# ============================================================================= #
# PARAMETER OPTIMZATION
# - Bootstrap
# - grid search
# ============================================================================= #

ctrl_tune <- trainControl(method = 'boot', number = 20, classProbs = TRUE, summaryFunction = twoClassSummary, search = 'grid')


# ============================================================================= #
# PREPROCESS: TEST DATASET
# - Same number of features as Train after variance removal (Mirroring the features)
# - Outlier detection and handling
# - KNN Imputation on test using train
# ============================================================================= #

class <- semcon_test_data$class
test_variance_removal <- semcon_test_data[ , which(names(semcon_test_data) %in% c(names(train_variance_removal)))]
test_outlier_NA <- apply(test_variance_removal, 2, impute_outlier_NA)
test_outlier_NA <- as.data.frame(test_outlier_NA)
test_knn_imputation <- KNN_imputation(test_outlier_NA, train_knn_imputation)
test <- cbind(class, test_knn_imputation)
dfSummary(test$class)
class_distribution(test, 'Test')


# ============================================================================= #
# CLASSIFICATION MODEL
# - Random Forest using rf
#	- Imbalanced Dataset
#	- ROSE
#	- SMOTE
# ============================================================================= #

# Model Creation
model_rf_imba <- train(form = class ~ ., data = imba_train_data, method = 'rf', tuneLength = 5, trControl = ctrl_tune, metric = 'ROC')
plot(model_rf_imba)
# Test Evaluation
pred_rf_class_imba = predict(model_rf_imba, newdata = test)
confusionMatrix(as.factor(pred_rf_class_imba), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_rf_imba <- roc.curve(test[['class']], pred_rf_class_imba, plotit = T, main = 'ROC Curve using RF')
print(roc_rf_imba$auc)

# Model Creation
model_rf_rose <- train(form = class ~ ., data = rose_train_data, method = 'rf', tuneLength = 5, trControl = ctrl_tune, metric = 'ROC')
plot(model_rf_rose)
# Test Evaluation
pred_rf_class_rose = predict(model_rf_rose, newdata = test)
confusionMatrix(as.factor(pred_rf_class_rose), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_rf_rose <- roc.curve(test[['class']], pred_rf_class_rose, plotit = T, main = 'ROC Curve using RF & ROSE')
print(roc_rf_rose$auc)


# Model Creation
model_rf_smote <- train(form = class ~ ., data = smote_train_data, method = 'rf', tuneLength = 5, trControl = ctrl_tune, metric = 'ROC')
plot(model_rf_smote)
# Test Evaluation
pred_rf_class_smote = predict(model_rf_smote, newdata = test)
confusionMatrix(as.factor(pred_rf_class_smote), as.factor(test[['class']]), positive = 'F', mode = 'everything')
roc_rf_smote <- roc.curve(test[['class']], pred_rf_class_smote, plotit = T, main = 'ROC Curve using RF & SMOTE')
print(roc_rf_smote$auc)
