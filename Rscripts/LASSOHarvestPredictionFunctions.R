# List and check required packages
required_packages = c("caret", "glmnet","ggplot2")
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
subsetPredictors <- function(df, byTime = NULL, bySensor = NULL, keep = c("Treatment"), exclude = FALSE) {
  # Initialize the columns to keep with those specified in 'keep'
  columns_to_keep <- keep[keep %in% colnames(df)]  # Ensure they exist in the data frame
  
  # Check if byTime is provided
  if (!is.null(byTime)) {
    # Extract the last numeric part from each column name
    time_values <- as.numeric(sub(".*_(\\d+)$", "\\1", colnames(df)))
    
    # Find the columns that have the last numeric part within the byTime range
    time_columns <- colnames(df)[time_values %in% byTime]
    
    # Combine 'time_columns' with 'keep' columns, ensuring no duplicates
    columns_to_keep <- unique(c(columns_to_keep, time_columns))
  }
  
  # Check if bySensor is provided
  if (!is.null(bySensor)) {
    # Create a regex pattern for matching the prefixes
    sensor_pattern <- paste0("^(", paste(bySensor, collapse = "|"), ")")
    
    # Find columns that match any sensor prefix in the bySensor vector
    sensor_columns <- colnames(df)[grepl(sensor_pattern, colnames(df))]
    
    if (exclude) {
      # Exclude the sensor columns and keep the others
      columns_to_keep <- setdiff(columns_to_keep, sensor_columns)
      all_other_columns <- setdiff(colnames(df), sensor_columns)
      columns_to_keep <- unique(c(columns_to_keep, all_other_columns))
    } else {
      # Include only the specified sensor columns
      columns_to_keep <- unique(c(columns_to_keep, sensor_columns))
    }
  }
  
  # Subset the data frame to keep only the selected columns
  df <- df[, columns_to_keep, drop = FALSE]
  
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
plot_model_accuracies <- function(models_list, accuracy_measure = "Rsquared", models_suffix, main = "Accuracies for Each Model", ...) {
  # Check if the accuracy measure is valid
  if (!accuracy_measure %in% c("RMSE", "Rsquared", "MAE")) {
    stop("Invalid accuracy measure. Choose from 'RMSE', 'Rsquared', or 'MAE'.")
  }
  
  # List to store accuracies for each model
  accuracies_list <- list()
  
  # Loop over each model in the list
  for (model_name in names(models_list)) {
    # Get the trained model
    model <- models_list[[model_name]]
    
    # Extract the resamples for the optimal lambda
    resamples <- model$resample
    
    # Extract the desired accuracy measure
    accuracies <- resamples[[accuracy_measure]]
    
    # Calculate mean and standard deviation
    mean_accuracy <- mean(accuracies, na.rm = TRUE)
    sd_accuracy <- sd(accuracies, na.rm = TRUE)
    
    # Print mean and standard deviation for each model
    cat(sprintf("Model: %s - Mean %s: %.4f, SD %s: %.4f\n", model_name, accuracy_measure, mean_accuracy, accuracy_measure, sd_accuracy))
    
    # Store accuracies in the list for later use
    accuracies_list[[sub(models_suffix, "", model_name)]] <- accuracies
  }
  
  # Combine all accuracies into a single data frame for plotting
  combined_accuracies <- stack(accuracies_list)
  
  # Adjust the margins and plot size for the boxplot
  par(mar = c(10, 5, 4, 2) + 0.1)  # Adjust inner margins: bottom, left, top, right
  par(oma = c(6, 2, 2, 2))  # Adjust outer margins: bottom, left, top, right
  
  # Create boxplots of all accuracies for each model
  boxplot(values ~ ind, data = combined_accuracies,
          main = main,
          xlab = "",
          ylab = paste(accuracy_measure, "Accuracy"),
          col = "lightblue",
          las = 2,
          ...)  # Rotate x-axis labels for better readability
}



#' Plot Non-Zero Lasso Coefficients
#'
#' Visualizes the non-zero coefficients from Lasso regression models stored in a list.
#' Each model's coefficients are extracted at the optimal lambda and plotted using `ggplot2`.
#'
#' @param models_list A named list of `caret::train` model objects trained with method `"glmnet"` (Lasso).
#'
#' @return A set of coefficient bar plots (printed to output).
#'
#' @examples
#' plot_non_zero_coefficients(my_lasso_models)
plot_non_zero_coefficients <- function(models_list) {
  # Loop over each model in the list
  for (model_name in names(models_list)) {
    # Get the trained model
    lasso_model <- models_list[[model_name]]
    
    # Extract coefficients for the optimal lambda
    optimal_lambda <- lasso_model$bestTune$lambda
    coef_optimal <- coef(lasso_model$finalModel, s = optimal_lambda)  # Extract coefficients
    
    # Convert coefficients to a data frame
    coef_df <- as.data.frame(as.matrix(coef_optimal))
    coef_df$Predictor <- rownames(coef_df)
    names(coef_df)[1] <- "Coefficient"
    
    # Remove the intercept term
    coef_df <- coef_df[coef_df$Predictor != "(Intercept)", ]
    
    # Filter out zero coefficients
    coef_df <- coef_df[coef_df$Coefficient != 0, ]
    
    # Sort by absolute coefficient value for better visualization
    coef_df <- coef_df[order(abs(coef_df$Coefficient), decreasing = TRUE), ]
    
    # Check if there are any non-zero coefficients to plot
    if (nrow(coef_df) > 0) {
      # Plot the coefficients using ggplot2
      p <- ggplot(coef_df, aes(x = reorder(Predictor, Coefficient), y = Coefficient)) +
        geom_bar(stat = "identity", fill = "lightblue") +
        coord_flip() +  # Flip the axes for better readability
        theme_minimal() +
        labs(title = paste("Lasso Coefficients for", model_name),
             x = "Predictor",
             y = "Coefficient") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      # Print the plot to ensure it shows in the output
      print(p)
    } else {
      cat(sprintf("No non-zero coefficients to plot for model: %s\n", model_name))
    }
  }
}

#' Train and Evaluate Lasso Models
#'
#' Trains Lasso regression models for one or more yield traits using repeated cross-validation.
#' Results and plots (accuracy + non-zero coefficients) are returned and optionally saved.
#'
#' @param predictorDF Data frame containing predictor variables (rows = plants).
#' @param outcomeDF Data frame containing outcome/yield variables with column `Plant.ID`.
#' @param yieldVars Character vector of outcome variable names to model.
#' @param num_folds Number of folds for cross-validation.
#' @param num_repeats Number of repeats for cross-validation.
#' @param saveModels Logical; whether to save the trained models to disk.
#' @param savePlots Logical; whether to generate and save performance plots.
#' @param modelListName String used in saved file names.
#' @param mainTitle Title for the model accuracy plot.
#' @param Algorithm Label for model type, saved with results (e.g., "Lasso").
#' @param DataGroup Optional label (e.g., sensor, stage).
#' @param Subsetting Description of subsetting/filtering applied.
#'
#' @return A list with elements:
#' \describe{
#'   \item{models}{Named list of trained `glmnet` models.}
#'   \item{results}{List of data frames with cross-validated resample performance.}
#' }
#'
#' @examples
#' fit <- trainAndEvaluateModels(predictors, outcomes, c("Yield"), 5, 3, TRUE, TRUE, "lassoSet1", "Lasso Accuracies", "Lasso", "Stage1", "NDVI_only")
trainAndEvaluateModels <- function(predictorDF, outcomeDF, yieldVars, num_folds, num_repeats, 
                                   saveModels, savePlots, modelListName, mainTitle, 
                                   Algorithm, DataGroup, Subsetting) {
  # Initialize the list to store models and the list to store results data frames
  models_list <- list()
  results_list <- list()
  
  # Verify that all specified yieldVars exist in the outcomeDF
  missing_vars <- setdiff(yieldVars, names(outcomeDF))
  if (length(missing_vars) > 0) {
    stop("The following yield variables are missing in outcomeDF: ", paste(missing_vars, collapse=", "))
  }
  
  # Loop over each outcome variable in yieldVars
  for (outcomeVar in yieldVars) {
    # Ensure outcomeVar is part of the columns
    if (!outcomeVar %in% names(outcomeDF)) {
      stop("Outcome variable '", outcomeVar, "' not found in outcomeDF")
    }
    
    # Prepare the data by merging predictor and outcome data frames based on the Plant.ID
    model_data <- merge(predictorDF, outcomeDF[c("Plant.ID", outcomeVar)], by.x = "row.names", by.y = "Plant.ID")
    rownames(model_data) <- model_data$Row.names
    model_data$Row.names <- NULL
    
    # Prepare predictor matrix and outcome vector
    x <- model.matrix(~ . - 1, data = model_data[, !(names(model_data) %in% c("Plant.ID", outcomeVar))])
    y <- model_data[[outcomeVar]]
    
    # Set up repeated cross-validation with stratification
    control <- trainControl(method = "repeatedcv", number = num_folds, repeats = num_repeats,
                            summaryFunction = defaultSummary, savePredictions = "final", verboseIter = FALSE)
    
    # Train the Lasso model
    lasso_model <- train(x, y, method = "glmnet", trControl = control,
                         tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.1, length = 100)),
                         preProcess = c("center", "scale"))
    
    # Save the trained model in the list
    model_name = paste(outcomeVar, "model", sep = "_")
    models_list[[model_name]] = lasso_model
    
    # Collect resample results
    res = lasso_model$resample
    res$Algorithm = Algorithm
    res$DataGroup = DataGroup
    res$Subsetting = Subsetting
    results_list[[model_name]] = res
    
    # Optionally save models
    if (saveModels) {
      save(lasso_model, file = paste0("../Rout/mlModels/yieldPrediction_", model_name, ".Rdata"))
    }
  }
  
  
  # Optionally create plots
  if (savePlots) {
    pdf(paste0("../Rout/LASSOHarvestPrediction/", modelListName, ".pdf"), width = 10, height = 8)
    plot_model_accuracies(models_list, models_suffix = "model", main = mainTitle)
    plot_non_zero_coefficients(models_list)  # Added this line to plot non-zero coefficients after accuracies
    dev.off()
  }
  
  return(list(models = models_list, results = results_list))
}


#' Train and Evaluate Stratified Lasso Models
#'
#' Similar to `trainAndEvaluateModels`, but uses stratified resampling based on a grouping variable.
#' This ensures that each fold maintains class distribution across repeats.
#'
#' @param predictorDF Data frame with predictors.
#' @param outcomeDF Data frame with yield outcomes and `Plant.ID`.
#' @param yieldVars Character vector of outcome variable names.
#' @param num_folds Number of cross-validation folds.
#' @param num_repeats Number of repeats for cross-validation.
#' @param stratify_by Column name in `model_data` used for stratification.
#' @param saveModels Logical; whether to save models to disk.
#' @param savePlots Logical; whether to save performance plots.
#' @param modelListName String used in file names.
#' @param mainTitle Title for plots.
#' @param Algorithm Label for model algorithm (e.g., "Lasso").
#' @param DataGroup Optional string to tag data source/group.
#' @param Subsetting Description of data filtering.
#'
#' @return A list with:
#' \describe{
#'   \item{models}{List of trained `glmnet` models.}
#'   \item{results}{List of data frames with performance metrics.}
#' }
#'
#' @examples
#' fit <- trainAndEvaluateStratifiedModels(predictors, outcomes, c("Yield"), 5, 3, "Treatment", TRUE, TRUE, "lassoStrat", "Lasso Accuracy", "Lasso", "Group1", "NDVI")
trainAndEvaluateStratifiedModels <- function(predictorDF, outcomeDF, yieldVars, num_folds, num_repeats, 
                                             stratify_by, saveModels, savePlots, modelListName, mainTitle, 
                                             Algorithm, DataGroup, Subsetting) {
  # Initialize the list to store models and the list to store results data frames
  models_list <- list()
  results_list <- list()
  
  # Verify that all specified yieldVars exist in the outcomeDF
  missing_vars <- setdiff(yieldVars, names(outcomeDF))
  if (length(missing_vars) > 0) {
    stop("The following yield variables are missing in outcomeDF: ", paste(missing_vars, collapse=", "))
  }
  
  # Loop over each outcome variable in yieldVars
  for (outcomeVar in yieldVars) {
    # Ensure outcomeVar is part of the columns
    if (!outcomeVar %in% names(outcomeDF)) {
      stop("Outcome variable '", outcomeVar, "' not found in outcomeDF")
    }
    
    # Prepare the data by merging predictor and outcome data frames based on the Plant.ID
    model_data <- merge(predictorDF, outcomeDF[c("Plant.ID", outcomeVar)], by.x = "row.names", by.y = "Plant.ID")
    rownames(model_data) <- model_data$Row.names
    model_data$Row.names <- NULL
    
    # Set up repeated stratified cross-validation
    folds <- createMultiFolds(model_data[[stratify_by]], k = num_folds, times = num_repeats)
    
    # Delete the stratification column from model_data before modeling
    model_data[[stratify_by]] <- NULL
    
    # Prepare predictor matrix and outcome vector
    x <- model.matrix(~ . - 1, data = model_data[, !(names(model_data) %in% c("Plant.ID", outcomeVar))])
    y <- model_data[[outcomeVar]]
    
    control <- trainControl(method = "repeatedcv", number = num_folds, repeats = num_repeats,
                            summaryFunction = defaultSummary, savePredictions = "final", verboseIter = FALSE,
                            index = folds)
    
    # Train the model using the specified method and parameters
    model <- train(x, y, method = "glmnet", trControl = control,
                   tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.1, length = 100)),
                   preProcess = c("center", "scale"))
    
    # Save the trained model in the list
    model_name <- paste(outcomeVar, "model", sep = "_")
    models_list[[model_name]] <- model
    
    # Collect resample results
    res <- model$resample
    res$Algorithm <- Algorithm
    res$DataGroup <- DataGroup
    res$Subsetting <- Subsetting
    results_list[[model_name]] <- res
    
    # Optionally save models
    if (saveModels) {
      save(model, file = paste0("../Rout/mlModels/yieldPrediction_", model_name, ".Rdata"))
    }
  }
  
  # Optionally create plots
  if (savePlots) {
    pdf(paste0("../Rout/LASSOHarvestPrediction/", modelListName,".pdf"), width = 10, height = 8)
    plot_model_accuracies(models_list, models_suffix = "model", main = mainTitle)
    plot_non_zero_coefficients(models_list)
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
combine_results <- function(result_names) {
  combined_results <- list()
  
  # Iterate through each outcome variable
  for (outcomeVar in yieldVars) {
    all_dfs <- list()
    
    # Gather all data frames for the current outcome variable
    for (list_name in result_names) {
      # Check if the list and the specific model result exist in the global environment
      if (exists(list_name) && !is.null(get(list_name)$results[[paste0(outcomeVar, "_model")]])) {
        all_dfs[[list_name]] <- get(list_name)$results[[paste0(outcomeVar, "_model")]]
      }
    }
    
    # Combine all gathered data frames into one, if not empty
    if (length(all_dfs) > 0) {
      combined_results[[outcomeVar]] <- do.call(rbind, all_dfs)
    }
  }
  
  return(combined_results)
}






