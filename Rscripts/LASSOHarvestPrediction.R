source("LassoHarvestPredictionFunctions.R")

if (!dir.exists("../Rout/")) {
  dir.create("../Rout/", recursive = TRUE)
}
# For outputs:
if (!dir.exists("../Rout/LASSOHarvestPrediction")) {
  dir.create("../Rout/LASSOHarvestPrediction", recursive = TRUE)
}
# For outputs:
if (!dir.exists("../Rout/mlModels/")) {
  dir.create("../Rout/mlModels/", recursive = TRUE)
}

# Note on data shape: The data used in this script is "widened", i.e. each
# sample has a row and each trait at each time point has a column. This is the
# typical format for caret-trained models.

### LOADING DATA ###
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# For cross-validation: we also load the non-reduces variable data
load("../Rout/Rdata/standardized_normalized_widened.Rdata")
load("../Rout/Rdata/StandardizedByDATWidened.Rdata")
AllPredictorVars$Treatment = NULL

# To use data aggregated by week uncomment lines below
# We also rename data frame to work with script
# load("../Rout/Rdata/StandardizedByDATWidenedAggregated.Rdata")
# load("../Rout/Rdata/standardized_normalized_widened_aggregated.Rdata")
# 
# AllPredictorVars = AllPredictorVarsAggregated
# AllPredictorVars$Treatment = NULL
# AllControlVars = AllControlVarsAggregated
# AllDroughtVars = AllDroughtVarsAggregated


### SCRIPT SETTINGS ###
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
yieldVars = c("Total.Biomass.DW","Total.Spike.weight")
# Uncomment to predict all yield variables
#yieldVars = setdiff(colnames(harvest),c("Plant.ID","Treatment","Plant.Name"))

num_folds = 3
# we used num_repeats 15, I set it to 1 here for testing speed
num_repeats = 1

set.seed(42)  # Set a seed for reproducibility

# Specify is you want to save the models. This will save each model in a separate file, but we can also save them together later
saveModels = F

# Specify if you want plots to be output to a pdf
savePlots = T

# PREDICTING X-RAY DATA 
# If you want to predict x-ray data, the easiest way to do this is just replace the harvest data with harvest_xray (uncomment below)
#harvest = harvest_xray


# CONTROL DATA

controlOut = trainAndEvaluateModels(AllControlVars, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                        saveModels = saveModels, savePlots = savePlots, modelListName = "control",
                                        mainTitle = "R^2 of Lasso Models, Control Data",
                                        Algorithm = "LASSO", DataGroup = "control", Subsetting = "None"
                                        )




# DROUGHT DATA
droughtOut = trainAndEvaluateModels(AllDroughtVars, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                        saveModels = saveModels, savePlots = savePlots, modelListName = "drought",
                                        mainTitle = "R^2 of Lasso Models, Drought Data",
                                        Algorithm = "LASSO", DataGroup = "drought", Subsetting = "None"
                                        )


# POOLED DATA
# Add a treatment column
AllPredictorVars$Treatment = as.factor(ifelse(grepl("L[0-9]+_C_[0-9]+", rownames(AllPredictorVars)),
                                               "Control",
                                               ifelse(grepl("L[0-9]+_D_[0-9]+", rownames(AllPredictorVars)),
                                                      "Drought",
                                                      NA)))  # Use NA for any unmatched rows (optional)

pooledOut = trainAndEvaluateStratifiedModels(AllPredictorVars, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                     stratify_by = "Treatment", saveModels = saveModels, savePlots = savePlots, modelListName = "pooled",
                                     mainTitle = "R^2 of Lasso Models, Pooled Data",
                                     Algorithm = "LASSO", DataGroup = "pooled", Subsetting = "None")



#############################
#### NON RGB PREDICTORS #####
#############################

### SUBSETTING PREDICTORS ###
AllDroughtVarsNoRGB = subsetPredictors(AllDroughtVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = T)
AllControlVarsNoRGB = subsetPredictors(AllControlVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = T)
AllPredictorVarsNoRGB = subsetPredictors(AllPredictorVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = T)

# CONTROL DATA

controlOutNoRGB = trainAndEvaluateModels(AllControlVarsNoRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                     saveModels = saveModels, savePlots = savePlots, modelListName = "controlNoRGB",
                                     mainTitle = "R^2 of Lasso Models, Control Data, NoRGB",
                                     Algorithm = "LASSO", DataGroup = "control", Subsetting = "NoRGB")




# DROUGHT DATA

droughtOutNoRGB = trainAndEvaluateModels(AllDroughtVarsNoRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                     saveModels = saveModels, savePlots = savePlots, modelListName = "droughtNoRGB", 
                                     mainTitle = "R^2 of Lasso Models, Drought Data, NoRGB",
                                     Algorithm = "LASSO", DataGroup = "drought", Subsetting = "NoRGB")


# POOLED DATA
# Add a treatment column
AllPredictorVarsNoRGB$Treatment = as.factor(ifelse(grepl("L[0-9]+_C_[0-9]+", rownames(AllPredictorVarsNoRGB)),
                                               "Control",
                                               ifelse(grepl("L[0-9]+_D_[0-9]+", rownames(AllPredictorVarsNoRGB)),
                                                      "Drought",
                                                      NA)))  # Use NA for any unmatched rows (optional)

pooledOutNoRGB = trainAndEvaluateStratifiedModels(AllPredictorVarsNoRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                              stratify_by = "Treatment", saveModels = saveModels, savePlots = savePlots, modelListName = "pooledNoRGB", 
                                              mainTitle = "R^2 of Lasso Models, Pooled Data, NoRGB",
                                              Algorithm = "LASSO", DataGroup = "pooled", Subsetting = "NoRGB")





##############################
#### ONLY RGB PREDICTORS #####
##############################

### SUBSETTING PREDICTORS ###
AllDroughtVarsOnlyRGB = subsetPredictors(AllDroughtVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = F)
AllControlVarsOnlyRGB = subsetPredictors(AllControlVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = F)
AllPredictorVarsOnlyRGB = subsetPredictors(AllPredictorVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = F)

# CONTROL DATA

controlOutOnlyRGB = trainAndEvaluateModels(AllControlVarsOnlyRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                          saveModels = saveModels, savePlots = savePlots, modelListName = "controlOnlyRGB",
                                          mainTitle = "R^2 of Lasso Models, Control Data, OnlyRGB",
                                          Algorithm = "LASSO", DataGroup = "control", Subsetting = "OnlyRGB")




# DROUGHT DATA

droughtOutOnlyRGB = trainAndEvaluateModels(AllDroughtVarsOnlyRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                          saveModels = saveModels, savePlots = savePlots, modelListName = "droughtOnlyRGB", 
                                          mainTitle = "R^2 of Lasso Models, Drought Data, OnlyRGB",
                                          Algorithm = "LASSO", DataGroup = "drought", Subsetting = "OnlyRGB")


# POOLED DATA
# Add a treatment column
AllPredictorVarsOnlyRGB$Treatment = as.factor(ifelse(grepl("L[0-9]+_C_[0-9]+", rownames(AllPredictorVarsOnlyRGB)),
                                                    "Control",
                                                    ifelse(grepl("L[0-9]+_D_[0-9]+", rownames(AllPredictorVarsOnlyRGB)),
                                                           "Drought",
                                                           NA)))  # Use NA for any unmatched rows (optional)

pooledOutOnlyRGB = trainAndEvaluateStratifiedModels(AllPredictorVarsOnlyRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                                   stratify_by = "Treatment", saveModels = saveModels, savePlots = savePlots, modelListName = "pooledOnlyRGB", 
                                                   mainTitle = "R^2 of Lasso Models, Pooled Data, OnlyRGB",
                                                   Algorithm = "LASSO", DataGroup = "pooled", Subsetting = "OnlyRGB")







##############################
### PRE DAT 51 PREDICTORS ####
##############################

### SUBSETTING PREDICTORS ###
AllDroughtVars1to51 = subsetPredictors(AllDroughtVars, byTime = 1:51)
AllControlVars1to51 = subsetPredictors(AllControlVars, byTime = 1:51)
AllPredictorVars1to51 = subsetPredictors(AllPredictorVars, byTime = 1:51)

# CONTROL DATA

controlOut1to51 = trainAndEvaluateModels(AllControlVars1to51, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                            saveModels = saveModels, savePlots = savePlots, modelListName = "control1to51",
                                            mainTitle = "R^2 of Lasso Models, Control Data, 1to51",
                                            Algorithm = "LASSO", DataGroup = "control", Subsetting = "1to51")




# DROUGHT DATA

droughtOut1to51 = trainAndEvaluateModels(AllDroughtVars1to51, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                            saveModels = saveModels, savePlots = savePlots, modelListName = "drought1to51", 
                                            mainTitle = "R^2 of Lasso Models, Drought Data, 1to51",
                                            Algorithm = "LASSO", DataGroup = "drought", Subsetting = "1to51")


# POOLED DATA
# Add a treatment column
AllPredictorVars1to51$Treatment = as.factor(ifelse(grepl("L[0-9]+_C_[0-9]+", rownames(AllPredictorVars1to51)),
                                                      "Control",
                                                      ifelse(grepl("L[0-9]+_D_[0-9]+", rownames(AllPredictorVars1to51)),
                                                             "Drought",
                                                             NA)))  # Use NA for any unmatched rows (optional)

pooledOut1to51 = trainAndEvaluateStratifiedModels(AllPredictorVars1to51, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                                     stratify_by = "Treatment", saveModels = saveModels, savePlots = savePlots, modelListName = "pooled1to51", 
                                                     mainTitle = "R^2 of Lasso Models, Pooled Data, 1to51",
                                                     Algorithm = "LASSO", DataGroup = "pooled", Subsetting = "1to51")


#######################################
### PRE DAT 51 PREDICTORS, NON RGB ####
#######################################

### SUBSETTING PREDICTORS ###
AllDroughtVars1to51NoRGB = subsetPredictors(AllDroughtVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = T)
AllDroughtVars1to51NoRGB = subsetPredictors(AllDroughtVars1to51NoRGB, byTime = 1:51)

AllControlVars1to51NoRGB = subsetPredictors(AllControlVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = T)
AllControlVars1to51NoRGB = subsetPredictors(AllControlVars1to51NoRGB, byTime = 1:51)

AllPredictorVars1to51NoRGB = subsetPredictors(AllPredictorVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = T)
AllPredictorVars1to51NoRGB = subsetPredictors(AllPredictorVars1to51NoRGB, byTime = 1:51)

# CONTROL DATA

controlOut1to51NoRGB = trainAndEvaluateModels(AllControlVars1to51NoRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                            saveModels = saveModels, savePlots = savePlots, modelListName = "control1to51NoRGB",
                                            mainTitle = "R^2 of Lasso Models, Control Data, 1to51NoRGB",
                                            Algorithm = "LASSO", DataGroup = "control", Subsetting = "1to51NoRGB")




# DROUGHT DATA

droughtOut1to51NoRGB = trainAndEvaluateModels(AllDroughtVars1to51NoRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                            saveModels = saveModels, savePlots = savePlots, modelListName = "drought1to51NoRGB", 
                                            mainTitle = "R^2 of Lasso Models, Drought Data, 1to51NoRGB",
                                            Algorithm = "LASSO", DataGroup = "drought", Subsetting = "1to51NoRGB")


# POOLED DATA
# Add a treatment column
AllPredictorVars1to51NoRGB$Treatment = as.factor(ifelse(grepl("L[0-9]+_C_[0-9]+", rownames(AllPredictorVars1to51NoRGB)),
                                                      "Control",
                                                      ifelse(grepl("L[0-9]+_D_[0-9]+", rownames(AllPredictorVars1to51NoRGB)),
                                                             "Drought",
                                                             NA)))  # Use NA for any unmatched rows (optional)

pooledOut1to51NoRGB = trainAndEvaluateStratifiedModels(AllPredictorVars1to51NoRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                                     stratify_by = "Treatment", saveModels = saveModels, savePlots = savePlots, modelListName = "pooled1to51NoRGB", 
                                                     mainTitle = "R^2 of Lasso Models, Pooled Data, 1to51NoRGB",
                                                     Algorithm = "LASSO", DataGroup = "pooled", Subsetting = "1to51NoRGB")

#######################################
### PRE DAT 51 PREDICTORS, ONLY RGB ###
#######################################

### SUBSETTING PREDICTORS ###
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
AllDroughtVars1to51OnlyRGB = subsetPredictors(AllDroughtVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = F)
AllDroughtVars1to51OnlyRGB = subsetPredictors(AllDroughtVars1to51OnlyRGB, byTime = 1:51)

AllControlVars1to51OnlyRGB = subsetPredictors(AllControlVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = F)
AllControlVars1to51OnlyRGB = subsetPredictors(AllControlVars1to51OnlyRGB, byTime = 1:51)

AllPredictorVars1to51OnlyRGB = subsetPredictors(AllPredictorVars, bySensor = c("RGB1","RGB2","RGB1.RGB2"), exclude = F)
AllPredictorVars1to51OnlyRGB = subsetPredictors(AllPredictorVars1to51OnlyRGB, byTime = 1:51)

# CONTROL DATA

controlOut1to51OnlyRGB = trainAndEvaluateModels(AllControlVars1to51OnlyRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                            saveModels = saveModels, savePlots = savePlots, modelListName = "control1to51OnlyRGB",
                                            mainTitle = "R^2 of Lasso Models, Control Data, 1to51OnlyRGB",
                                            Algorithm = "LASSO", DataGroup = "control", Subsetting = "1to51OnlyRGB")




# DROUGHT DATA

droughtOut1to51OnlyRGB = trainAndEvaluateModels(AllDroughtVars1to51OnlyRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                            saveModels = saveModels, savePlots = savePlots, modelListName = "drought1to51OnlyRGB", 
                                            mainTitle = "R^2 of Lasso Models, Drought Data, 1to51OnlyRGB",
                                            Algorithm = "LASSO", DataGroup = "drought", Subsetting = "1to51OnlyRGB")


# POOLED DATA
# Add a treatment column
AllPredictorVars1to51OnlyRGB$Treatment = as.factor(ifelse(grepl("L[0-9]+_C_[0-9]+", rownames(AllPredictorVars1to51OnlyRGB)),
                                                      "Control",
                                                      ifelse(grepl("L[0-9]+_D_[0-9]+", rownames(AllPredictorVars1to51OnlyRGB)),
                                                             "Drought",
                                                             NA)))  # Use NA for any unmatched rows (optional)

pooledOut1to51OnlyRGB = trainAndEvaluateStratifiedModels(AllPredictorVars1to51OnlyRGB, outcomeDF = harvest, yieldVars = yieldVars, num_folds = num_folds, num_repeats = num_repeats,
                                                     stratify_by = "Treatment", saveModels = saveModels, savePlots = savePlots, modelListName = "pooled1to51OnlyRGB", 
                                                     mainTitle = "R^2 of Lasso Models, Pooled Data, 1to51OnlyRGB",
                                                     Algorithm = "LASSO", DataGroup = "pooled", Subsetting = "1to51OnlyRGB")




##########################################
############ JOIN ACCURACY DFs ###########
##########################################
modelsList = c("controlOut", "droughtOut", "pooledOut", 
               "controlOut1to51", "droughtOut1to51", "pooledOut1to51", 
               "controlOut1to51NoRGB", "droughtOut1to51NoRGB", "pooledOut1to51NoRGB", 
               "controlOut1to51OnlyRGB", "droughtOut1to51OnlyRGB", "pooledOut1to51OnlyRGB", 
               "controlOutOnlyRGB", "droughtOutOnlyRGB", "pooledOutOnlyRGB", 
               "controlOutNoRGB", "droughtOutNoRGB", "pooledOutNoRGB")

combinedLASSOAccuracies = combine_results(modelsList)

save(combinedLASSOAccuracies, file = "../Rout/LASSOHarvestPrediction/combinedLASSOAccuracies.Rdata")

### OPTIONAL: SAVE ALL MODELS TOGETHER ###
save(list = modelsList, file = "../Rout/mlModels/allLASSOHarvestPredictions.Rdata")
