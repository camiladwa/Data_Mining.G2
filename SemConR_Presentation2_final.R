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
# ============================================================================= #

load_lib <- c('grDevices', 'foreign', 'caret', 'DMwR', 'mice', 'Boruta', 'randomForest', 'ggplot2', 'pacman', 'missMDA')
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

SRC <- 'C:/Users/YDS/Desktop/R/SemCom/SemCon/Source/secom_mod.SAV'
# SRC <- 'C:/Users/unnat/Desktop/YDS/YS_SEMCON/Source/secom_mod.SAV'
TRAIN_TEST_RATIO <- 0.80
NULL_THRESHOLD <- 0.6
# MAXITR <- 15
# EXTRACTITR <- 15


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
  dir <- "C:/Users/YDS/Desktop/R/SemCom/SemCon/IntermediateFiles"
  # dir <- "C:/Users/unnat/Desktop/YDS/YS_SEMCON"
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
  df <- train_null_removal[-nearZeroVar(df)] 
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

mice_imputation <- function(df, itr = 5, extract = 5){
  # Function to impute missing values(NAs) using MICE, default 5 iterations and extract data after completion of 5th iteration
  impute_mice <- mice(df, m = itr, method = 'pmm', seed = 500)
  summary(impute_mice)
  mice_imputate <- complete(impute_mice, extract)
  return(mice_imputate)
}

PCA_Analysis <- function(df) {
  # Function to check the PCs using Scree Plot and Kaiser-Guttman Rule.
  pca <- prcomp(df, center = TRUE, scale. = TRUE)
  pca_var <- pca$sdev ^ 2
  pca_var_per <- round(pca_var * 100 / sum(pca_var) , 2)
  plot(1:length(pca_var), pca_var, type="b", col='blue', ylab="Eigenvalue", xlab="Component Number", main = 'Scree Plot') 
  abline(h = 1,lty = 2,col = "red")
  message('Principal components as per Kaiser-Guttman rule: ', length(pca_var[pca_var >= 1]))
  plot(1:length(pca_var_per), pca_var_per, type="b", col='red', ylab="Proportion of Variance Explained", xlab="Component Number", main = 'Variance Explained by Components')
  message('Total Variance explained by 118 PCs: ', sum(pca_var_per[1:length(pca_var[pca_var >= 1])]))
}

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
# READING SOURCE DATA
# ============================================================================= #

semcon_original_data <- read_sav(SRC)


# ============================================================================= #
# SOURCE DATA ANALYSIS
# ============================================================================= #

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

class_distribution(semcon_train_data, 'Test Data')
missing_value_analysis(semcon_test_data, 'Test Data')
outlier_analysis(semcon_test_data[, -c(1)], 'Test Data')
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

train_mean_imputation <- mean_imputation(train_outlier_NA)
missing_value_analysis(train_mean_imputation, 'train set after mean imputation')
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
#	
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

nPCs <- estim_ncpPCA(train_outlier_NA, method.cv = "Kfold", verbose = FALSE)
res_comp <- imputePCA(train_outlier_NA, ncp = nPCs$ncp)
train_NAH_pca <- as.data.frame(res_comp$completeObs)

# ============================================================================= #
# FEATURE SELECTION AND REDUCTION
# 01. PCA
# ============================================================================= #

PCA_Analysis(train_knn_imputation)
PCA_Analysis(train_NAH_pca)


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

train_FR_boruta_pca <- selection_boruta(train_NAH_pca, class)
print("Summary of selected features")
summary(train_FR_boruta_pca)
message("Variance corresponding to selected features:")
sapply(train_FR_boruta_pca, var)
missing_value_analysis(train_FR_boruta_pca, 'train set after BORUTA')
# outlier_analysis(train_FR_boruta_pca, 'train set after BORUTA')
# write_csv(cbind(train_FR_boruta_pca), "train_FR_boruta.csv")


# ============================================================================= #
# CLASSIFICATION MODEL - WIP (DO NOT RUN)
# 01. Random Forest
# ============================================================================= #

set.seed(333)
#rf_all <- randomForest(class ~ ., data = )
rf_model <- randomForest(class ~ ., data = train_FR_boruta_knn)
rf_prob <- predict(rf_model, train_FR_boruta_knn)
rf_pred = ifelse(rf_prob > 0.5, 1, 0)
confusionMatrix(table(rf_pred, semcon_train_data$class))
