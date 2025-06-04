# List and check required packages
required_packages = c("caret", "randomForest")
missing_packages = required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# Prompt to install if needed
if (length(missing_packages) > 0) {
  message("The following packages are missing: ", paste(missing_packages, collapse = ", "))
  install = readline(prompt = "Do you want to install them now? [y/n]: ")
  if (tolower(install) == "y") {
    install.packages(missing_packages, dependencies = TRUE)
    rm(install)
  } else {
    stop("Required packages not installed. Script cannot continue.")
  }
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Remove temporary variables
rm(required_packages, missing_packages)

#' Subset Predictor Columns by Time or Sensor
#'
#' Filters a data frame to include (or exclude) predictors based on sensor prefix and/or time suffix.
#'
#' @param df A data frame containing predictor variables.
#' @param byTime Optional numeric vector; keeps only columns with time suffixes matching these values.
#' @param bySensor Optional character vector; matches column names starting with these sensor names.
#' @param keep Character vector of column names to always keep (e.g., "Treatment").
#' @param exclude Logical; if TRUE, excludes columns matching the `bySensor` prefix instead of including them.
#'
#' @return A filtered data frame containing only the selected columns.
#'
#' @examples
#' subsetPredictors(mydata, byTime = c(35:49), bySensor = c("VNIR", "RGB"))
subsetPredictors = function(df, byTime = NULL, bySensor = NULL, keep = c("Treatment"), exclude = FALSE) {
  # Initialize the columns to keep with those specified in 'keep'
  columns_to_keep = keep[keep %in% colnames(df)]  # Ensure they exist in the data frame
  
  # Check if byTime is provided
  if (!is.null(byTime)) {
    # Extract the last numeric part from each column name
    time_values = as.numeric(sub(".*_(\\d+)$", "\\1", colnames(df)))
    
    # Find the columns that have the last numeric part within the byTime range
    time_columns = colnames(df)[time_values %in% byTime]
    
    # Combine 'time_columns' with 'keep' columns, ensuring no duplicates
    columns_to_keep = unique(c(columns_to_keep, time_columns))
  }
  
  # Check if bySensor is provided
  if (!is.null(bySensor)) {
    # Create a regex pattern for matching the prefixes
    sensor_pattern = paste0("^(", paste(bySensor, collapse = "|"), ")")
    
    # Find columns that match any sensor prefix in the bySensor vector
    sensor_columns = colnames(df)[grepl(sensor_pattern, colnames(df))]
    
    if (exclude) {
      # Exclude the sensor columns and keep the others
      columns_to_keep = setdiff(columns_to_keep, sensor_columns)
      all_other_columns = setdiff(colnames(df), sensor_columns)
      columns_to_keep = unique(c(columns_to_keep, all_other_columns))
    } else {
      # Include only the specified sensor columns
      columns_to_keep = unique(c(columns_to_keep, sensor_columns))
    }
  }
  
  # Subset the data frame to keep only the selected columns
  df = df[, columns_to_keep, drop = FALSE]
  
  # Return the filtered data frame
  return(df)
}


#' Plot Model Accuracy Distributions
#'
#' Generates boxplots comparing the accuracy of models across different outcome variables.
#'
#' @param models_list A named list of trained caret models.
#' @param accuracy_measure A character string indicating which accuracy metric to plot ("RMSE", "Rsquared", or "MAE").
#' @param models_suffix A suffix string to remove from model names for cleaner x-axis labels.
#' @param main A character string specifying the plot title.
#'
#' @return A boxplot is printed to the graphics device.
#'
#' @examples
#' plot_model_accuracies(models_list, accuracy_measure = "Rsquared", models_suffix = "_model", main = "RF Performance")
plot_model_accuracies = function(models_list, accuracy_measure = "Rsquared", models_suffix, main) {
  # Check if the accuracy measure is valid
  if (!accuracy_measure %in% c("RMSE", "Rsquared", "MAE")) {
    stop("Invalid accuracy measure. Choose from 'RMSE', 'Rsquared', or 'MAE'.")
  }
  
  # List to store accuracies for each model
  accuracies_list = list()
  
  # Loop over each model in the list
  for (model_name in names(models_list)) {
    # Get the trained model
    model = models_list[[model_name]]
    
    # Extract the resamples for the optimal lambda
    resamples = model$resample
    
    # Extract the desired accuracy measure
    accuracies = resamples[[accuracy_measure]]
    
    # Calculate mean and standard deviation
    mean_accuracy = mean(accuracies, na.rm = TRUE)
    sd_accuracy = sd(accuracies, na.rm = TRUE)
    
    # Print mean and standard deviation for each model
    cat(sprintf("Model: %s - Mean %s: %.4f, SD %s: %.4f\n", model_name, accuracy_measure, mean_accuracy, accuracy_measure, sd_accuracy))
    
    # Store accuracies in the list for later use
    accuracies_list[[sub(models_suffix, "", model_name)]] = accuracies
  }
  
  # Combine all accuracies into a single data frame for plotting
  combined_accuracies = stack(accuracies_list)
  
  # Adjust the margins and plot size for the boxplot
  par(mar = c(10, 5, 4, 2) + 0.1)  # Adjust inner margins: bottom, left, top, right
  par(oma = c(2, 2, 2, 2))  # Adjust outer margins: bottom, left, top, right
  
  # Create boxplots of all accuracies for each model
  boxplot(values ~ ind, data = combined_accuracies,
          main = main,
          xlab = "",
          ylab = paste(accuracy_measure, "Accuracy"),
          col = "lightblue",
          las = 2)  # Rotate x-axis labels for better readability
}


#' Train and Evaluate Random Forest Models for Harvest Trait Prediction
#'
#' Trains random forest models using cross-validation for each trait in `harvestVars`, evaluates performance, and optionally saves outputs.
#'
#' @param predictorDF A data frame of predictor variables (row names should match Plant.IDs).
#' @param outcomeDF A data frame containing outcome/yield variables and a "Plant.ID" column.
#' @param harvestVars A character vector of outcome variable names to model.
#' @param num_folds Number of folds for cross-validation.
#' @param num_repeats Number of repeats for cross-validation.
#' @param saveModels Logical; whether to save the trained models.
#' @param savePlots Logical; whether to save a summary plot of model accuracy.
#' @param modelListName Base name for saved plot/model files.
#' @param mainTitle Title for the accuracy plot.
#' @param Algorithm A label to include in the result metadata.
#' @param DataGroup A label to include in the result metadata.
#' @param Subsetting A label to include in the result metadata.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{models}{A list of trained model objects.}
#'   \item{results}{A list of resample accuracy data frames for each model.}
#' }
#'
#' @examples
#' trainAndEvaluateModels(predictorDF, outcomeDF, harvestVars = c("Yield"), 5, 3, TRUE, TRUE, "ModelList", "Accuracy Plot", "RF", "Group1", "SubsetA")
trainAndEvaluateModels = function(predictorDF, outcomeDF, harvestVars, num_folds, num_repeats, 
                                   saveModels, savePlots, modelListName, mainTitle, 
                                   Algorithm, DataGroup, Subsetting) {
  models_list = list()
  results_list = list()
  
  missing_vars = setdiff(harvestVars, names(outcomeDF))
  if (length(missing_vars) > 0) {
    stop("The following yield variables are missing in outcomeDF: ", paste(missing_vars, collapse=", "))
  }
  
  for (outcomeVar in harvestVars) {
    model_data = merge(predictorDF, outcomeDF[c("Plant.ID", outcomeVar)], by.x = "row.names", by.y = "Plant.ID")
    rownames(model_data) = model_data$Row.names
    model_data$Row.names = NULL
    
    control = trainControl(method = "repeatedcv", number = num_folds, repeats = num_repeats, 
                            search = "random", savePredictions = "final", classProbs = FALSE, verboseIter = FALSE)
    
    rf_model = train(as.formula(paste(outcomeVar, "~ .")), data = model_data, 
                      method = "rf", trControl = control, ntree = 500)
    
    model_name = paste(outcomeVar, "model", sep = "_")
    models_list[[model_name]] = rf_model
    res = rf_model$resample
    res$Algorithm = Algorithm
    res$DataGroup = DataGroup
    res$Subsetting = Subsetting
    results_list[[model_name]] = res
    
    if (saveModels) {
      save(rf_model, file = paste0("../Rout/mlModels/yieldPrediction_RF_", model_name, ".Rdata"))
    }
  }
  
  if (savePlots) {
    pdf(paste0("../Rout/RFHarvestPrediction/", modelListName, ".pdf"), width = 10, height = 8)
    plot_model_accuracies(models_list, models_suffix = "_model", main = mainTitle)
    dev.off()
  }
  
  return(list(models = models_list, results = results_list))
}


#' Train Random Forest Models with Stratified Cross-Validation
#'
#' Trains random forest models for each trait, using stratified repeated cross-validation on a user-specified grouping variable.
#'
#' @param predictorDF A data frame of predictors with row names as Plant.ID.
#' @param outcomeDF A data frame with outcomes and a "Plant.ID" column.
#' @param harvestVars A vector of outcome variable names to model.
#' @param num_folds Number of folds for cross-validation.
#' @param num_repeats Number of repeats for cross-validation.
#' @param stratify_by Column name in outcomeDF used to create stratified folds.
#' @param saveModels Logical; if TRUE, saves trained models to files.
#' @param savePlots Logical; if TRUE, saves summary accuracy plots to PDF.
#' @param modelListName Base name for output files.
#' @param mainTitle Title for the accuracy plot.
#' @param Algorithm Label to annotate results.
#' @param DataGroup Label to annotate results.
#' @param Subsetting Label to annotate results.
#'
#' @return A list with:
#' \describe{
#'   \item{models}{A list of trained models.}
#'   \item{results}{A list of resampling results per model.}
#' }
#'
#' @examples
#' trainAndEvaluateStratifiedModels(predictorDF, outcomeDF, c("Yield"), 5, 3, "Genotype", TRUE, TRUE, "StratifiedModels", "Stratified Accuracy", "RF", "DataGroup", "SubsetX")
trainAndEvaluateStratifiedModels = function(predictorDF, outcomeDF, harvestVars, num_folds, num_repeats, 
                                             stratify_by, saveModels, savePlots, modelListName, mainTitle, 
                                             Algorithm, DataGroup, Subsetting) {
  # Initialize the list to store models and results
  models_list = list()
  results_list = list()
  
  # Check if all specified yield variables exist in the outcome data frame
  missing_vars = setdiff(harvestVars, names(outcomeDF))
  if (length(missing_vars) > 0) {
    stop("The following yield variables are missing in outcomeDF: ", paste(missing_vars, collapse=", "))
  }
  
  for (outcomeVar in harvestVars) {
    # Merge predictor data and outcome data based on Plant.ID
    model_data = merge(predictorDF, outcomeDF[c("Plant.ID", outcomeVar)], by.x = "row.names", by.y = "Plant.ID")
    rownames(model_data) = model_data$Row.names
    model_data$Row.names = NULL
    
    # Setup repeated stratified cross-validation folds
    folds = createMultiFolds(model_data[[stratify_by]], k = num_folds, times = num_repeats)
    control = trainControl(method = "repeatedcv", number = num_folds, repeats = num_repeats,
                            savePredictions = "final", classProbs = FALSE, verboseIter = FALSE, index = folds)
    
    # Remove the stratification column from the data frame
    model_data[[stratify_by]] = NULL
    
    # Train the Random Forest model
    rf_model = train(as.formula(paste(outcomeVar, "~ .")), data = model_data, 
                      method = "rf", trControl = control, ntree = 500)
    
    # Save the trained model in the list with its associated name
    model_name = paste(outcomeVar, "model", sep = "_")
    models_list[[model_name]] = rf_model
    res = rf_model$resample
    res$Algorithm = Algorithm
    res$DataGroup = DataGroup
    res$Subsetting = Subsetting
    results_list[[model_name]] = res
    
    # Optionally save the models to files
    if (saveModels) {
      save(rf_model, file = paste0("../Rout/mlModels/yieldPrediction_RF_", model_name, ".Rdata"))
    }
  }
  
  # Optionally create plots of model accuracies
  if (savePlots) {
    pdf(paste0("../Rout/RFHarvestPrediction/", modelListName, ".pdf"), width = 10, height = 8)
    plot_model_accuracies(models_list, models_suffix = "_model", main = mainTitle)
    dev.off()
  }
  
  return(list(models = models_list, results = results_list))
}


#' Combine Resampling Results from Multiple Model Runs
#'
#' Aggregates resampling results from several modeling runs (stored as named lists in the global environment).
#'
#' @param result_names A character vector of object names containing results from `trainAndEvaluateModels` or similar.
#'
#' @return A named list of combined data frames, one per outcome variable.
#'
#' @examples
#' combined <- combine_results(c("modelRun1", "modelRun2", "modelRun3"))
combine_results = function(result_names) {
  combined_results = list()
  
  # Iterate through each outcome variable
  for (outcomeVar in harvestVars) {
    all_dfs = list()
    
    # Gather all data frames for the current outcome variable
    for (list_name in result_names) {
      # Check if the list and the specific model result exist in the global environment
      if (exists(list_name) && !is.null(get(list_name)$results[[paste0(outcomeVar, "_model")]])) {
        all_dfs[[list_name]] = get(list_name)$results[[paste0(outcomeVar, "_model")]]
      }
    }
    
    # Combine all gathered data frames into one, if not empty
    if (length(all_dfs) > 0) {
      combined_results[[outcomeVar]] = do.call(rbind, all_dfs)
    }
  }
  
  return(combined_results)
}


