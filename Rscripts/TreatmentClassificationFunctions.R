# List and check required packages
required_packages = c("xtable", "caret", "randomForest", "ggplot2")
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

# For nicer predictor names
clean_label = function(x) {
  # Remove "MM" (case-sensitive) for unit in RGB
  x = gsub("_MM_", "_", x)
  # Remove the word "plant_Avg" (case-insensitive)
  x = gsub("_(?i)plant_Avg_", " ", x, perl = TRUE)
  # Replace underscore-number endings with ", DAT xx"
  x = gsub("_(\\d+)$", ", \\1 DAT", x, perl = TRUE)
  # Remove the word "plant" (case-insensitive)
  x = gsub("_(?i)plant_", " ", x, perl = TRUE)
  # Replace underscores with spaces
  x = gsub("_", " ", x)
  # Replace literal dots with spaces
  x = gsub("\\.", " ", x)
  # Trim extra spaces
  x = trimws(x)
  # Remove trailing "MM" (case-sensitive) for unit in RGB
  x = gsub("MM$", "", x)
  # Capitalize only the first letter
  x = sub("^(\\w)", "\\U\\1", x, perl = TRUE)
  # Reintroduce dots between rgb color codes
  x = gsub("138 139 97", "138.139.97", x)
  x = gsub("68 85 63", "68.85.63", x)
  x = gsub("(?i)deltaT", "CTD", x, perl = TRUE)
  # Take care of caps-locked RGB traits
  x = gsub("AREA", "area", x)
  x = gsub("HEIGHT", "height", x)
  x = gsub("WIDTH", "width", x)
  x = gsub("COMPACTNESS", "compactness", x)
  x = gsub("PERIMETER", "perimeter", x)
  
  return(x)
}

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
subsetPredictors = function(df, byTime = NULL, bySensor = NULL, keep = c("Treatment")) {
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
    
    # Combine 'sensor_columns' with 'columns_to_keep', ensuring no duplicates
    columns_to_keep = unique(c(columns_to_keep, sensor_columns))
  }
  
  # Subset the data frame to keep only the selected columns
  df = df[, columns_to_keep, drop = FALSE]
  
  # Return the filtered data frame
  return(df)
}


#' Pairwise Wilcoxon Rank-Sum Tests Between Rows
#'
#' Performs pairwise Wilcoxon rank-sum tests between all unique row pairs in a numeric data frame. Useful for testing whether rows differ significantly in their distributions.
#'
#' @param df A data frame where each row represents a sample or group and each column is a numeric feature.
#' @param p.adjust.method Method used to adjust p-values for multiple testing (default: "BH"). See \code{\link{p.adjust.methods}}.
#'
#' @return A data frame with each row representing a pairwise comparison. Includes raw and adjusted p-values.
#'
#' @examples
#' results <- pairwise_wilcox(my_data_matrix)
pairwise_wilcox = function(df, p.adjust.method = "BH") {
  n = nrow(df)
  results = data.frame(
    Row1 = character(),
    Row2 = character(),
    P.Value = numeric(),
    Adjusted.P.Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Iterate through all row pairs
  pairs = combn(1:n, 2)
  
  for (i in 1:ncol(pairs)) {
    r1 = pairs[1, i]
    r2 = pairs[2, i]
    
    # Perform Wilcoxon test
    test = wilcox.test(
      as.numeric(df[r1, ]), 
      as.numeric(df[r2, ])
    )
    
    results = rbind(
      results,
      data.frame(
        Row1 = paste0("Row_", r1),
        Row2 = paste0("Row_", r2),
        P.Value = test$p.value,
        Adjusted.P.Value = NA
      )
    )
  }
  
  # Adjust p-values for multiple testing
  results$Adjusted.P.Value = p.adjust(results$P.Value, method = p.adjust.method)
  
  return(results)
}

#' Plot Random Forest Variable Importance
#'
#' Creates a horizontal bar plot of the top `k` most important variables from a random forest model trained via the `caret` package.
#'
#' @param rf_model A trained random forest model object from \code{\link[caret]{train}}.
#' @param top_k Integer; number of top variables to display (default: 20).
#' @param main_title Character string; title for the plot (default: "Top K Important Variables").
#' @param show_mean_accuracy Logical; whether to display the mean accuracy ± SD on the plot (default: TRUE).
#'
#' @return A base R barplot is displayed.
#'
#' @note Uses the `clean_label` function to simplify variable names for display. You may remove or replace this step if not applicable.
#'
#' @examples
#' plot_rf_importance(rf_model, top_k = 15)
plot_rf_importance = function(rf_model, 
                                               top_k = 20, 
                                               main_title = "Top K Important Variables", 
                                               show_mean_accuracy = TRUE) {
  # 1) Extract variable importance
  var_importance = varImp(rf_model)$importance
  
  if ("Overall" %in% colnames(var_importance)) {
    importance_df = data.frame(
      Variable = rownames(var_importance), 
      Importance = var_importance$Overall, 
      stringsAsFactors = FALSE
    )
  } else {
    importance_df = data.frame(
      Variable = rownames(var_importance), 
      Importance = rowMeans(var_importance), 
      stringsAsFactors = FALSE
    )
  }
  
  # 2) Sort and select top K
  top_variables = head(importance_df[order(-importance_df$Importance), ], top_k)
  
  # 3) Apply label cleaning (optional — remove this line if you don't want it)
  top_variables$CleanedVariable = sapply(top_variables$Variable, clean_label)
  
  # 4) Set factor levels for plotting
  top_variables$CleanedVariable = factor(top_variables$CleanedVariable, 
                                          levels = top_variables$CleanedVariable[order(top_variables$Importance)])
  
  # 5) Plot settings
  par(mar = c(5, 12, 4, 4), oma = c(1, 3, 1, 3))
  
  # 6) Barplot
  x_max = max(top_variables$Importance)
  barplot(top_variables$Importance, 
          names.arg = top_variables$CleanedVariable, 
          las = 2, horiz = TRUE,
          main = main_title, 
          cex.main = 1,
          xlab = "Variable Importance", 
          col = "lightblue", 
          border = "black", 
          xlim = c(0, x_max * 1.19))  # Extend x-axis
  
  # 7) Display mean accuracy (optional)
  if (show_mean_accuracy && "Accuracy" %in% colnames(rf_model$resample)) {
    mean_accuracy = mean(rf_model$resample$Accuracy, na.rm = TRUE)
    sd_accuracy = sd(rf_model$resample$Accuracy, na.rm = TRUE)
    accuracy_label = sprintf("Mean Accuracy: %.2f ± %.2f", mean_accuracy, sd_accuracy)
    mtext(side = 3, line = -0.5, accuracy_label, adj = 0, cex = 1)
  }
}
