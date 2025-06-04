source("TreatmentClassificationFunctions.R")

load("../Rout/Rdata/StandardizedByDATWidened.Rdata")

# Note on data shape: The data used in this script is "widened", i.e. each
# sample has a row and each trait at each time point has a column. This is the
# typical format for caret-trained models.

################################################################
#-------------- FULL DATA CLASSIFICATION ----------------------#
################################################################
# Using predictors from all time points (most comprehensive model)

# Set up repeated 3-fold cross-validation
control = trainControl(method = "repeatedcv",  # Repeated cross-validation
                        number = 3,             # Number of folds
                        repeats = 5,            # Number of repeats
                        verboseIter = TRUE)     # To see progress

# Train the Random Forest model
set.seed(123)  # Set a seed for reproducibility
rf_model = train(Treatment ~ .,                # Formula to predict Treatment using all other variables
                  data = AllPredictorVars,                  # The dataset
                  method = "rf",                # Random Forest method
                  trControl = control,          # Cross-validation control
                  importance = TRUE)            # To calculate variable importance


if (!dir.exists("../Rout/mlModels")) {
  dir.create("../Rout/mlModels", recursive = TRUE)
}
save(rf_model,file = "../Rout/mlModels/TreatmentPrediction_RF.Rdata")

#optional: load
#load("../Rout/mlModels/TreatmentPrediction_RF.Rdata")
# Output the results
print(rf_model)

# View the importance of features
varImportance = varImp(rf_model)$importance
varImportance = varImportance[order(varImp(rf_model)$importance$Drought, decreasing = TRUE),]

if (!dir.exists("../Rout/treatmentClassification")) {
  dir.create("../Rout/treatmentClassification", recursive = TRUE)
}
pdf("../Rout/treatmentClassification/treatmentClassificationVarImp.pdf", width = 10, height = 7, onefile = T)
fullClusters = plot_rf_importance(rf_model, main_title = "")
dev.off()

################################################################
#--------------- BY WEEK CLASSIFICATION -----------------------#
################################################################
# Using data from each week separately
folds = 3
repeats = 1
nCVs = folds*repeats
nTopPredictors = 10

# Set up repeated 3-fold cross-validation
control = trainControl(method = "repeatedcv",  # Repeated cross-validation
                        number = folds,             # Number of folds
                        repeats = repeats,            # Number of repeats
                        verboseIter = TRUE)     # To see progress
# Train the Random Forest model
set.seed(123)  # Set a seed for reproducibility

# time ranges (weeks by DAT)
timeRanges = list(24:29,30:36,37:43,44:49,50:57,58:64,65:71,72:78,79:85,86:92)

modelAccuracies = data.frame()
featureImportances = data.frame()

modelList = list()

for (timeRange in timeRanges) {
  print(timeRange)
  subsetData = subsetPredictors(AllPredictorVars,byTime = timeRange)
  
  rf_model = train(Treatment ~ .,                # Formula to predict Treatment using all other variables
                    data = subsetData,                  # The dataset
                    method = "rf",                # Random Forest method
                    trControl = control,          # Cross-validation control
                    importance = TRUE)            # To calculate variable importance
  
  accuracies = rf_model$resample$Accuracy  # 'rf_model' is the Random Forest model object created with caret
  varImportance = varImp(rf_model)$importance
  varImportance = varImportance[order(varImp(rf_model)$importance$Drought, decreasing = TRUE),]
  
  thisRow = c(ncol(subsetData)-1,accuracies)
  
  modelAccuracies = rbind(modelAccuracies,thisRow)
  featureImportances = rbind(featureImportances,rownames(varImportance)[1:nTopPredictors])
  
  modelList[[paste0("model", min(timeRange), "to", max(timeRange))]] = rf_model
  
}

rownames(modelAccuracies) = as.character(timeRanges)
colnames(modelAccuracies)[1] = "nPredictors"
colnames(modelAccuracies)[2:(nCVs+1)] = paste("accuracy_",1:nCVs,sep = "")
colnames(featureImportances) = paste("importantFeature_",1:nTopPredictors,sep = "")

# PLOTTING ACCURACIES
# add full model accuracies
# load("../Rout/mlModels/TreatmentPrediction_RF.Rdata")
modelList[["model24to93"]] = rf_model
accuracies = rf_model$resample$Accuracy
thisRow = c(ncol(rf_model$trainingData)-1,accuracies)

modelAccuracies = rbind(modelAccuracies,thisRow)

# rename the new row
rownames(modelAccuracies)[nrow(modelAccuracies)] = "24:93"


pdf("../Rout/treatmentClassification/weekWiseAccuracy.pdf", width = 10, height = 6)
boxplot(t(as.matrix(modelAccuracies[,-1])),xlab = "Time Range of Predictors (DAT)", ylab = "Classification Accuracy", col = "lightblue")

# or using weekly labels
rownames(modelAccuracies) = c(as.character(0:9),"all weeks")
boxplot(t(as.matrix(modelAccuracies[,-1])),xlab = "Time Range of Predictors (week)", ylab = "Classification Accuracy", col = "lightblue")
dev.off()

rowMeans(modelAccuracies[,-1])

