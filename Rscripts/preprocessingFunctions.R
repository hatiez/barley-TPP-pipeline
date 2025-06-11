# List of required packages
required_packages = c(
  "ggplot2", "dplyr", "tidyr", "gridExtra", "grid", "cowplot", "lme4",
  "reshape2", "stringr", "missForest", "MASS", "car", "readr", "openxlsx","purrr"
)

# Identify missing packages
missing_packages = required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# Prompt user to install missing packages
if (length(missing_packages) > 0) {
  message("The following packages are missing: ", paste(missing_packages, collapse = ", "))
  install = readline(prompt = "Do you want to install them now? [y/n]: ")
  if (tolower(install) == "y") {
    install.packages(missing_packages, dependencies = TRUE)
  } else {
    stop("Required packages not installed. Script cannot continue.")
  }
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))



#' Merge and Filter Data Frames
#'
#' Function to remove unnecessary columns and merge data frames from multiple files.
#' @param mergeVars Variables to merge by (usually DAT and Plant.ID).
#' @param removeVars Variables to remove (usually Date, Plant.Name, Treatment, PID).
#' @param ... Data frames to merge.
#' @return A merged data frame.
#' @examples
#' merge_and_filter_dfs(c("DAT","Plant.ID"), c("Date", "PID"), df1, df2, df3)
merge_and_filter_dfs = function(mergeVars=c("DAT","Plant.ID","Plant.Name", "Treatment"), removeVars=c("Date", "PID"), ...) {
  dfs = list(...)
  df_names = as.character(substitute(list(...)))[-1L]
  all_col_names = unlist(lapply(dfs, names))
  duplicate_col_names = all_col_names[duplicated(all_col_names)]
  duplicate_col_names = setdiff(duplicate_col_names, mergeVars)
  
  remove_columns = function(df, removeVars) {
    df[, !names(df) %in% removeVars]
  }
  
  rename_columns = function(df, df_name, duplicate_col_names) {
    for (col in duplicate_col_names) {
      if (col %in% names(df)) {
        names(df)[names(df) == col] = paste0(col, "_", df_name)
      }
    }
    return(df)
  }
  
  reduced_dfs = lapply(dfs, remove_columns, removeVars = removeVars)
  
  for (i in seq_along(reduced_dfs)) {
    reduced_dfs[[i]] = rename_columns(reduced_dfs[[i]], df_names[i], duplicate_col_names)
  }
  
  merge_dfs = function(x, y, mergeVars) {
    merge(x, y, by = mergeVars, all = TRUE)
  }
  
  merged_df = Reduce(function(x, y) merge_dfs(x, y, mergeVars), reduced_dfs)
  
  return(merged_df)
}

#' Group-Wise Imputation Using missForest
#'
#' Function to apply missForest imputation on missing values within groups of data and return imputed data and out-of-bag (OOB) error.
#' @param df The input data frame.
#' @param groupVars Variables to group by.
#' @param ignoreVars Variables to ignore during imputation.
#' @param ... Additional arguments passed to missForest.
#' @return A list containing the imputed data frame and a data frame with OOB error for each group.
#' @examples
#' groupWiseMissForest(df, groupVars = c("DAT", "Treatment", "Plant.Name"), ignoreVars = c("Plant.ID"))
groupWiseMissForest = function(df, groupVars = c("DAT", "Treatment", "Plant.Name"), ignoreVars = c("Plant.ID"), ...) {
  if (!all(groupVars %in% names(df))) {
    stop("All groupVars must be columns in the data frame")
  }
  
  result_df = df
  oob_error_df = data.frame(Group = character(), OOBError = numeric(), stringsAsFactors = FALSE)
  group_combinations = unique(df[groupVars])
  
  for (i in 1:nrow(group_combinations)) {
    subset_condition = rep(TRUE, nrow(df))
    
    for (var in groupVars) {
      df_var = as.character(df[[var]])
      group_var = as.character(group_combinations[i, var])
      subset_condition = subset_condition & (df_var == group_var)
    }
    
    group_data = df[subset_condition, ]
    ignore_data = group_data[, ignoreVars, drop = FALSE]
    group_data = group_data[, !(names(group_data) %in% ignoreVars)]
    
    if (!any(is.na(group_data))) {
      next
    }
    
    for (var in groupVars) {
      group_data[[var]] = as.factor(group_data[[var]])
    }
    
    missForest_result = suppressWarnings(missForest(as.data.frame(group_data), ...))
    imputed_data = missForest_result$ximp
    
    oob_error = missForest_result$OOBerror[["NRMSE"]]
    group_label = paste(group_combinations[i, ], collapse = ", ")
    oob_error_df = rbind(oob_error_df, data.frame(Group = group_label, OOBError = oob_error, stringsAsFactors = FALSE))
    
    result_df[subset_condition, setdiff(names(imputed_data),groupVars)] = imputed_data[,!(colnames(imputed_data) %in% groupVars)]
  }
  
  return(list(ImputedData = result_df, OOBError = oob_error_df))
}



#' Z-Score Normalization
#'
#' Apply z-score normalization across the entire column or within groups using dplyr.
#' @param df The input data frame.
#' @param groupVars Variables to group by. If NULL, standardization is applied across the whole column.
#' @param ignoreVars Variables to ignore during scaling.
#' @return A data frame with z-score normalized values.
#' @examples
#' zScore(df, groupVars = c("DAT", "Treatment"), ignoreVars = c("Plant.Name","Plant.ID"))
zScore = function(df, groupVars = c("DAT", "Treatment"), ignoreVars = c("Plant.Name","Plant.ID")) {
  
  if (!is.null(groupVars) && !all(groupVars %in% names(df))) {
    stop("groupVars must be columns in the data frame")
  }
  
  scaleVars = setdiff(names(df), c(groupVars, ignoreVars))
  
  if (!is.null(groupVars)) {
    df = df %>%
      group_by(across(all_of(groupVars))) %>%
      mutate(across(all_of(scaleVars), ~ scale(.)[, 1]))
    df = ungroup(df)
  } else {
    df = df %>%
      mutate(across(all_of(scaleVars), ~ scale(.)[, 1]))
  }
  
  return(df)
}




#' Min-Max Normalization
#'
#' Apply Min-Max normalization across the entire column or within groups using dplyr.
#' @param df The input data frame.
#' @param groupVars Variables to group by. If NULL, normalization is applied across the whole column.
#' @param ignoreVars Variables to ignore during scaling.
#' @return A data frame with Min-Max normalized values.
#' @examples
#' minMaxNormalize(df, groupVars = c("DAT", "Treatment"), ignoreVars = c("Plant.Name","Plant.ID"))
minMaxNormalize = function(df, groupVars = c("DAT", "Treatment"), ignoreVars = c("Plant.Name", "Plant.ID")) {
  
  if (!is.null(groupVars) && !all(groupVars %in% names(df))) {
    stop("groupVars must be columns in the data frame")
  }
  
  # Identify the columns to be normalized (excluding groupVars and ignoreVars)
  scaleVars = setdiff(names(df), c(groupVars, ignoreVars))
  
  if (!is.null(groupVars)) {
    # Perform Min-Max normalization within groups
    df = df %>%
      group_by(across(all_of(groupVars))) %>%
      mutate(across(all_of(scaleVars), ~ (. - min(.)) / (max(.) - min(.)))) %>%
      ungroup()
  } else {
    # Perform Min-Max normalization across the entire column
    df = df %>%
      mutate(across(all_of(scaleVars), ~ (. - min(.)) / (max(.) - min(.))))
  }
  
  return(df)
}


#' Calculate NA Rate for Numeric Columns
#'
#' Function to calculate the rate of NA values in numeric columns of a data frame.
#' @param df The input data frame.
#' @param ignoreVars Variables to ignore during NA calculation.
#' @return A list with the number of NAs, total entries, and NA rate.
#' @examples
#' numericNArate(df, ignoreVars = c("DAT"))
numericNArate = function(df, ignoreVars = c("DAT")){
  numeric_df = df[sapply(df, is.numeric)]
  numeric_df = numeric_df[, !(names(numeric_df) %in% ignoreVars)]
  
  NAs = sum(colSums(is.na(numeric_df)))
  entries = nrow(numeric_df) * ncol(numeric_df)
  
  return(list(nNAs = NAs, nEntries = entries, NArate = (NAs / entries)))
}

#' Calculate NA Rate Per Variable
#'
#' Function to calculate the NA rate per variable in a data frame.
#' @param df The input data frame.
#' @param ignoreVars Variables to ignore during NA rate calculation.
#' @return A data frame with NA counts and rates for each variable.
#' @examples
#' numericNAratePerVar(df, ignoreVars = c("DAT"))
numericNAratePerVar = function(df, ignoreVars = c("DAT")) {
  numeric_df = df[sapply(df, is.numeric)]
  numeric_df = numeric_df[, !(names(numeric_df) %in% ignoreVars)]
  
  na_rate_df = data.frame(Column = character(),
                           nNAs = integer(),
                           nEntries = integer(),
                           NArate = numeric(),
                           stringsAsFactors = FALSE)
  
  for (col in names(numeric_df)) {
    nNAs = sum(is.na(numeric_df[[col]]))
    nEntries = nrow(numeric_df)
    NArate = nNAs / nEntries
    
    na_rate_df = rbind(na_rate_df, data.frame(Column = col, nNAs = nNAs, nEntries = nEntries, NArate = NArate))
  }
  
  return(na_rate_df)
}

#' Calculate NA Rate Per Replicate
#'
#' Function to calculate NA counts and rates per replicate in a data frame.
#' @param df The input data frame.
#' @param ignoreVars Variables to ignore during NA rate calculation.
#' @param groupingColumns Columns to group by (e.g., Plant.Name, Treatment).
#' @param IDcolumn Column identifying replicates.
#' @return A list containing data frames for NA rate and NA count by variable.
#' @examples
#' numericNAratePerReplicate(df, ignoreVars = c("DAT"), groupingColumns = c("Plant.Name", "Treatment"))
numericNAratePerReplicate = function(df, ignoreVars = c("DAT"), groupingColumns = c("Plant.Name", "Treatment"), IDcolumn = "Plant.ID") {
  
  add_replicate_column = function(data, column_name = IDcolumn) {
    data$Replicate = as.character(sub(".*_(\\d+)$", "\\1", data[[column_name]]))
    return(data)
  }
  
  df = add_replicate_column(df, column_name = IDcolumn)
  
  df[[IDcolumn]] = as.character(df[[IDcolumn]])
  df[groupingColumns] = lapply(df[groupingColumns], as.character)
  
  numeric_df = df[sapply(df, is.numeric)]
  numeric_df = numeric_df[, !(names(numeric_df) %in% ignoreVars)]
  
  na_rate_df = data.frame(
    Genotype = character(),
    Treatment = character(),
    Replicate = character(),
    nNAs = integer(),
    nEntries = integer(),
    NArate = numeric(),
    stringsAsFactors = FALSE
  )
  
  na_count_by_variable_df = data.frame(
    Genotype = character(),
    Treatment = character(),
    Replicate = character(),
    Variable = character(),
    nNAs = integer(),
    stringsAsFactors = FALSE
  )
  
  groups = unique(df[, c(groupingColumns, "Replicate")])
  
  for (i in 1:nrow(groups)) {
    genotype = groups[i, "Plant.Name"]
    treatment = groups[i, "Treatment"]
    replicate = groups[i, "Replicate"]
    subset_df = numeric_df[df[["Plant.Name"]] == genotype & df[["Treatment"]] == treatment & df$Replicate == replicate, ]
    
    nNAs = sum(is.na(subset_df))
    nEntries = nrow(subset_df) * ncol(subset_df)
    NArate = ifelse(nEntries == 0, NA, nNAs / nEntries)
    
    na_rate_df = rbind(na_rate_df, data.frame(
      Genotype = genotype,
      Treatment = treatment,
      Replicate = replicate,
      nNAs = nNAs,
      nEntries = nEntries,
      NArate = NArate
    ))
    
    for (col in names(subset_df)) {
      nNAs_var = sum(is.na(subset_df[[col]]))
      na_count_by_variable_df = rbind(na_count_by_variable_df, data.frame(
        Genotype = genotype,
        Treatment = treatment,
        Replicate = replicate,
        Variable = col,
        nNAs = nNAs_var
      ))
    }
  }
  
  return(list(na_rate_df = na_rate_df, na_count_by_variable_df = na_count_by_variable_df))
}

#' Summed Outlier Counts Per Plant
#'
#' Function to sum outlier counts across multiple data frames and return a merged data frame with total outliers per plant.
#' @param outlier_list A list of data frames containing outlier counts.
#' @param merge_keys Columns to merge by (e.g., Genotype, Treatment, Replicate).
#' @param VarNameCol Name of the column containing variable names.
#' @param outlierCol Name of the column containing outlier counts.
#' @return A data frame with summed outlier counts per plant.
#' @examples
#' plantwise_summed_outlier_counts(outlier_list, merge_keys = c("Genotype", "Treatment", "Replicate"))
plantwise_summed_outlier_counts = function(outlier_list, merge_keys = c("Genotype", "Treatment", "Replicate"), VarNameCol = "Variable", outlierCol = "nNAs") {
  results_list = list()
  all_keys = unique(do.call(rbind, lapply(outlier_list, function(df) df[, merge_keys])))
  
  merged_df = all_keys
  rownames(merged_df) = NULL
  
  for (df_name in names(outlier_list)) {
    df = outlier_list[[df_name]]
    merged_df[[df_name]] = 0
    
    for (i in 1:nrow(merged_df)) {
      matching_rows = df
      for (key in merge_keys) {
        matching_rows = matching_rows[matching_rows[[key]] == merged_df[i, key], ]
      }
      
      if (nrow(matching_rows) > 0) {
        merged_df[i, df_name] = sum(matching_rows[[outlierCol]], na.rm = TRUE)
      }
    }
  }
  
  outlier_columns = setdiff(names(merged_df), merge_keys)
  merged_df$TotalOutliers = rowSums(merged_df[, outlier_columns], na.rm = TRUE)
  
  return(merged_df)
}

#' Impute Missing Values by Group
#'
#' Function to impute missing values based on the mean value of replicates within groups.
#' @param df The input data frame.
#' @param identifiers Variables to group by.
#' @return A data frame with imputed values.
#' @examples
#' imputeMissingValuesByGroup(df, identifiers = c("DAT", "Plant.ID"))
imputeMissingValuesByGroup = function(df, identifiers) {
  if (!all(identifiers %in% names(df))) {
    stop("Identifiers must be columns in the data frame")
  }
  
  df[identifiers] = lapply(df[identifiers], as.factor)
  
  df = df %>%
    group_by(across(all_of(identifiers))) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
    ungroup()
  
  return(df)
}

#' Remove Z-Score Outliers
#'
#' Function to remove outliers based on z-score threshold.
#' @param df The input data frame.
#' @param threshold The z-score threshold to identify outliers.
#' @param groupVars Variables to group by.
#' @param ignoreVars Variables to ignore during outlier detection.
#' @return A data frame with outliers removed.
#' @examples
#' removeZscoreOutliers(df, threshold = 3, groupVars = c("DAT", "Treatment"))
removeZscoreOutliers = function(df, threshold = 3, groupVars = c("DAT", "Treatment"), ignoreVars = c("Date", "Plant.ID", "PID", "Plant.Name","DAT")) {
  ZscoreDF = zScore(df, groupVars = groupVars, ignoreVars = ignoreVars)
  for (col in names(ZscoreDF)) {
    if (is.numeric(ZscoreDF[[col]]) && !col %in% ignoreVars) {
      df[which(ZscoreDF[[col]] >= threshold), col] = NA
      df[which(ZscoreDF[[col]] <= threshold * (-1)), col] = NA
    }
  }
  return(df)
}

#' Remove IQR Outliers
#'
#' Function to remove outliers based on IQR method.
#' @param df The input data frame.
#' @param groupVars Variables to group by.
#' @param ignoreVars Variables to ignore during outlier detection.
#' @param multiplier The IQR multiplier to identify outliers.
#' @return A data frame with outliers removed.
#' @examples
#' removeIQROutliers(df, groupVars = c("DAT", "Treatment"), multiplier = 1.5)
removeIQROutliers = function(df, groupVars = c("DAT", "Treatment"), ignoreVars = c("Plant.Name", "Plant.ID"), multiplier = 1.5) {
  if (!all(groupVars %in% names(df))) {
    stop("groupVars must be columns in the data frame")
  }
  
  outlierVars = setdiff(names(df), c(groupVars, ignoreVars))
  
  df = df %>%
    group_by(across(all_of(groupVars))) %>%
    mutate(across(all_of(outlierVars), ~ {
      Q1 = quantile(., 0.25, na.rm = TRUE)
      Q3 = quantile(., 0.75, na.rm = TRUE)
      IQR = Q3 - Q1
      lower_bound = Q1 - multiplier * IQR
      upper_bound = Q3 + multiplier * IQR
      replace(., . < lower_bound | . > upper_bound, NA)
    })) %>%
    ungroup()
  
  return(df)
}

#' Calculate NA Rate Per Variable Per Group
#'
#' Function to calculate NA rates per variable and group by a specified column.
#' @param df The input data frame.
#' @param groupVar The column to group by for NA rate calculation.
#' @return A data frame with NA counts per variable per group.
#' @examples
#' numericNAratePerVarPerGroup(df, groupVar = "DAT")
numericNAratePerVarPerGroup = function(df, groupVar) {
  na_count_df = df %>%
    group_by(across(all_of(groupVar))) %>%
    summarise(across(where(is.numeric), ~ sum(is.na(.)), .names = "nNAs_{col}")) %>%
    pivot_longer(cols = starts_with("nNAs_"), names_to = "Variable", values_to = "nNAs") %>%
    mutate(Variable = gsub("nNAs_", "", Variable)) %>%
    ungroup()
  return(na_count_df)
}

#' Group-wise Shapiro-Wilk Test
#'
#' Function to perform Shapiro-Wilk test for normality within groups or across entire columns if groupVars is NULL.
#' @param df The input data frame.
#' @param groupVars Variables to group by. If NULL, the test is applied across the entire column.
#' @param ignoreVars Variables to ignore during the test.
#' @param p.adjust.method The method to adjust p-values.
#' @return A list with unadjusted and adjusted p-value tables.
#' @examples
#' groupwiseShapiroTest(df, groupVars = c("DAT", "Treatment", "Plant.Name"), p.adjust.method = "holm")
groupwiseShapiroTest = function(df, groupVars = c("DAT", "Treatment", "Plant.Name"), ignoreVars = c("DAT", "Treatment", "Plant.ID", "Plant.Name"), p.adjust.method = "holm") {
  
  if (!is.null(groupVars) && !all(groupVars %in% names(df))) {
    stop("groupVars must be columns in the data frame")
  }
  
  testVars = setdiff(names(df), c(groupVars, ignoreVars))
  
  safe_shapiro_test = function(x) {
    result = tryCatch(shapiro.test(x)$p.value, error = function(e) NA)
    return(result)
  }
  
  if (!is.null(groupVars)) {
    df %>%
      group_by(across(all_of(groupVars))) %>%
      summarise(across(all_of(testVars), list(p.value = ~ safe_shapiro_test(na.omit(.))), .names = "{.col}_p.value")) %>%
      ungroup() %>%
      pivot_longer(cols = -all_of(groupVars), names_to = "variable", values_to = "p.value") %>%
      mutate(variable = gsub("_p.value", "", variable)) %>%
      group_by(variable) %>%
      mutate(p.value.adjusted = p.adjust(p.value, method = p.adjust.method)) %>%
      pivot_wider(names_from = variable, values_from = c(p.value, p.value.adjusted)) %>%
      as.data.frame() -> results_table
  } else {
    df %>%
      summarise(across(all_of(testVars), list(p.value = ~ safe_shapiro_test(na.omit(.))), .names = "{.col}_p.value")) %>%
      pivot_longer(everything(), names_to = "variable", values_to = "p.value") %>%
      mutate(variable = gsub("_p.value", "", variable)) %>%
      mutate(p.value.adjusted = p.adjust(p.value, method = p.adjust.method)) %>%
      as.data.frame() -> results_table
  }
  
  unadjusted_table = results_table %>%
    select(all_of(groupVars), starts_with("p.value_")) %>%
    rename_with(~ gsub("p.value_", "", .), starts_with("p.value_"))
  
  adjusted_table = results_table %>%
    select(all_of(groupVars), starts_with("p.value.adjusted_")) %>%
    rename_with(~ gsub("p.value.adjusted_", "", .), starts_with("p.value.adjusted_"))
  
  return(list(unadjusted_table = unadjusted_table, adjusted_table = adjusted_table))
}

#' Group-wise Normality Test
#'
#' Function to perform normality test (Shapiro-Wilk or Kolmogorov-Smirnov) within groups.
#' @param df The input data frame.
#' @param groupVars Variables to group by.
#' @param ignoreVars Variables to ignore during the test.
#' @param test The type of normality test ("shapiro" or "ks").
#' @param p.adjust.method The method to adjust p-values.
#' @return A list with unadjusted and adjusted p-value tables.
#' @examples
#' groupWiseNormalityTest(df, groupVars = c("DAT", "Treatment", "Plant.Name"), test = "shapiro", p.adjust.method = "holm")
groupWiseNormalityTest = function(df, groupVars = c("DAT", "Treatment", "Plant.Name"), ignoreVars = c("DAT", "Treatment", "Plant.ID", "Plant.Name"), test = "shapiro", p.adjust.method = "holm") {
  if (!all(groupVars %in% names(df))) {
    stop("groupVars must be columns in the data frame")
  }
  
  testVars = setdiff(names(df), c(groupVars, ignoreVars))
  
  safe_shapiro_test = function(x) {
    result = tryCatch(shapiro.test(x)$p.value, error = function(e) NA)
    return(result)
  }
  
  safe_ks_test = function(x) {
    result = tryCatch(ks.test(x, y = "pnorm")$p.value, error = function(e) NA)
    return(result)
  }
  
  test_function = ifelse(test == "shapiro", safe_shapiro_test, safe_ks_test)
  
  df %>%
    group_by(across(all_of(groupVars))) %>%
    summarise(across(all_of(testVars), list(p.value = ~ test_function(na.omit(.))), .names = "{.col}_p.value")) %>%
    ungroup() %>%
    pivot_longer(cols = -all_of(groupVars), names_to = "variable", values_to = "p.value") %>%
    mutate(variable = gsub("_p.value", "", variable)) %>%
    group_by(variable) %>%
    mutate(p.value.adjusted = p.adjust(p.value, method = p.adjust.method)) %>%
    pivot_wider(names_from = variable, values_from = c(p.value, p.value.adjusted)) %>%
    as.data.frame() -> results_table
  
  unadjusted_table = results_table %>%
    dplyr::select(all_of(groupVars), starts_with("p.value_")) %>%
    dplyr::rename_with(~ gsub("p.value_", "", .))
  
  adjusted_table = results_table %>%
    dplyr::select(all_of(groupVars), starts_with("p.value.adjusted_")) %>%
    dplyr::rename_with(~ gsub("p.value.adjusted_", "", .))
  
  return(list(unadjusted_table = unadjusted_table, adjusted_table = adjusted_table))
}

#' Box-Cox Transformation and Redistribution
#'
#' Function to apply Box-Cox transformation on specified groups within a data frame.
#' @param df The input data frame.
#' @param nonNormalGroups Data frame specifying the groups that are non-normal.
#' @param GroupVars Variables to group by.
#' @param numericVar The variable to transform.
#' @return A data frame with transformed values.
#' @examples
#' boxcoxRedistribute(df, nonNormalGroups, GroupVars = c("Treatment", "DAT"))
boxcoxRedistribute = function(df, nonNormalGroups, GroupVars = c("Treatment", "DAT"), numericVar = "variable") {
  for (i in 1:nrow(nonNormalGroups)) {
    current_group = nonNormalGroups[i, ]
    
    subset_condition = Reduce(`&`, lapply(GroupVars, function(col) df[[col]] == current_group[[col]]))
    variable_name = current_group[[numericVar]]
    numeric_vector = df[subset_condition, variable_name]
    
    if (all(is.na(numeric_vector))) {
      warning(paste("All values are NA for", variable_name, "in group:", paste(current_group, collapse = ",")))
      next
    }
    
    if (any(numeric_vector <= 0, na.rm = TRUE)) {
      shift_amount = -min(numeric_vector, na.rm = TRUE) + 1
      numeric_vector = numeric_vector + shift_amount
      #print(shift_amount)
    }
    
    if (all(numeric_vector > 0, na.rm = TRUE)) {
      boxcox_result = boxcox(numeric_vector ~ 1, plotit = FALSE)
      lambda = boxcox_result$x[which.max(boxcox_result$y)]
      #print(paste(lambda))
      transformed_vector = if (lambda == 0) log(numeric_vector) else (numeric_vector^lambda - 1) / lambda
      df[subset_condition, variable_name] = transformed_vector
    } else {
      warning(paste("Non-positive values found in", variable_name, "for group:", paste(current_group, collapse = ",")))
    }
  }
  
  return(df)
}

#' Box-Cox Transformation without Grouping
#'
#' Function to apply Box-Cox transformation across the entire data frame without grouping.
#' @param df The input data frame.
#' @param nonNormalGroups Data frame specifying the non-normal variables.
#' @param numericVar The variable to transform.
#' @return A data frame with transformed values.
#' @examples
#' boxcoxRedistributeNongrouping(df, nonNormalGroups)
boxcoxRedistributeNongrouping = function(df, nonNormalGroups, numericVar = "variable") {
  
  for (i in 1:nrow(nonNormalGroups)) {
    variable_name = nonNormalGroups[i, numericVar, drop = TRUE]
    numeric_vector = df[[variable_name]]
    
    if (all(is.na(numeric_vector))) {
      warning(paste("All values are NA for", variable_name))
      next
    }
    
    if (any(numeric_vector <= 0, na.rm = TRUE)) {
      shift_amount = -min(numeric_vector, na.rm = TRUE) + 1
      numeric_vector = numeric_vector + shift_amount
    }
    
    if (all(numeric_vector > 0, na.rm = TRUE)) {
      boxcox_result = boxcox(numeric_vector ~ 1, plotit = FALSE)
      lambda = boxcox_result$x[which.max(boxcox_result$y)]
      transformed_vector = if (lambda == 0) log(numeric_vector) else (numeric_vector^lambda - 1) / lambda
      df[[variable_name]] = transformed_vector
    } else {
      warning(paste("Non-positive values found in", variable_name))
    }
  }
  
  return(df)
}

#' Global Box-Cox Transformation
#'
#' Function to apply Box-Cox transformation globally across groups in a data frame.
#' @param df The input data frame.
#' @param groupVars Variables to group by.
#' @param ignoreVars Variables to ignore during transformation.
#' @return A list containing the transformed data frame and transformation details.
#' @examples
#' globalBoxCox(df, groupVars = c("DAT", "Treatment"))
globalBoxCox = function(df, groupVars = c("DAT", "Treatment"), ignoreVars = c("DAT", "Treatment", "Plant.ID", "Plant.Name")) {
  transformations = data.frame(
    Group = character(),
    Variable = character(),
    Lambda = numeric(),
    stringsAsFactors = FALSE
  )
  
  numeric_df = df[sapply(df, is.numeric)]
  numeric_df = numeric_df[, !(names(numeric_df) %in% ignoreVars)]
  
  grouped_df = split(numeric_df, df[groupVars], drop = TRUE)
  
  for (group_name in names(grouped_df)) {
    group_data = grouped_df[[group_name]]
    
    for (variable in names(group_data)) {
      numeric_vector = group_data[[variable]]
      
      if (all(numeric_vector > 0, na.rm = TRUE)) {
        boxcox_result = boxcox(numeric_vector ~ 1, plotit = FALSE)
        lambda = boxcox_result$x[which.max(boxcox_result$y)]
        
        if (abs(lambda - 1) > 0.01) {
          group_details = unlist(strsplit(group_name, ".", fixed = TRUE))
          transformations = rbind(transformations, data.frame(
            Group = paste(groupVars, group_details, sep = "=", collapse = ", "),
            Variable = variable,
            Lambda = lambda,
            stringsAsFactors = FALSE
          ))
        }
        
        transformed_vector = if (lambda == 0) log(numeric_vector) else (numeric_vector^lambda - 1) / lambda
        df[match(row.names(group_data), row.names(df)), variable] = transformed_vector
      } else {
        warning(paste("Non-positive values found in", variable, "for group:", group_name))
      }
    }
  }
  
  return(list(transformed_df = df, transformations = transformations))
}

#' Plant ID Correction
#'
#' Function to correct treatment column errors based on Plant.ID information.
#' @param df The input data frame.
#' @param silent If TRUE, suppresses output messages.
#' @return A data frame with corrected treatment columns.
#' @examples
#' plantIDCorrection(df)
plantIDCorrection = function(df, silent = FALSE) {
  current_df = df
  current_df = cbind(current_df, TreatmentCode = str_extract(current_df$Plant.ID, "_[D,C]_"))
  current_df$TreatmentCode = dplyr::recode(current_df$TreatmentCode, "_C_" = "Control", "_D_" = "Drought")
  
  errorCount = 0
  for (i in 1:nrow(current_df)) {
    if (current_df$Treatment[i] != current_df$TreatmentCode[i]) {
      errorCount = errorCount + 1
      df$Treatment[current_df$DAT[i] == df$DAT & current_df$Plant.ID[i] == df$Plant.ID & current_df$Treatment == df$Treatment] = current_df$TreatmentCode[i]
    }
  }
  
  if (!silent) {
    print(paste("Corrected", errorCount, "treatment column errors."))
  }
  
  return(df)
}

#' Pivot Data Frame for Analysis
#'
#' Function to pivot a data frame to a wider format, one row per replicate and one column per trait at each time point.
#' @param df The input data frame.
#' @param replicate_col Column identifying replicates.
#' @param ignoreVars Variables to ignore during pivoting.
#' @return A pivoted data frame.
#' @examples
#' pivot_data_frame(df, replicate_col = "Plant.ID")
pivot_data_frame = function(df, replicate_col = "Plant.ID", ignoreVars = c("Plant.ID", "DAT", "Plant.Name", "Treatment")) {
  # Identify trait columns to pivot
  trait_cols = setdiff(colnames(df), ignoreVars)
  
  # Check if there is only one trait column
  if (length(trait_cols) == 1) {
    # Handle case with only one variable column
    trait_col_name = trait_cols[1]
    
    # Pivot wider using the trait name and DAT for new column names
    df_wide = df %>%
      pivot_wider(
        names_from = "DAT",
        values_from = all_of(trait_col_name),
        names_prefix = paste0(trait_col_name, "_")  # Prefix with the variable name
      )
  } else {
    # Default case: multiple trait columns, pivot as before
    df_wide = df %>%
      pivot_wider(
        names_from = c("DAT"),
        values_from = all_of(trait_cols),
        names_sep = "_"
      )
  }
  
  return(df_wide)
}



#' Widen and Merge Data Frames
#'
#' Function to widen multiple data frames and merge them on a common column.
#' @param DFlist A list of data frames to widen and merge.
#' @param mergeby Columns to merge by.
#' @param removeVars Variables to remove during merging.
#' @return A merged data frame.
#' @examples
#' widen_and_merge(list(df1, df2, df3), mergeby = c("Plant.ID"))
widen_and_merge = function(DFlist,mergeby = c("Plant.ID"),removeVars = c("DAT","Treatment","Plant.Name")){
  output = pivot_data_frame(DFlist[[1]], replicate_col = mergeby)
  
  output = output[ , -which(names(output) %in% removeVars)]
  
  for (i in 2:length(DFlist)) {
    df = DFlist[[i]]
    df = pivot_data_frame(df)
    df = df[ , -which(names(df) %in% removeVars)]
    output = merge(output,df, by = mergeby)
  }
  return(as.data.frame(output))
}

#' Aggregate Weekly Trait Data by Plant
#'
#' This function aggregates temporal trait data into weekly summaries for each plant.
#' For each plant and each week (as defined by `weekDATs`), it calculates min, mean, and max
#' values for each trait. If the aggregated values are globally constant (i.e., min = mean = max
#' for all traits across the entire dataset), only the mean values are retained.
#'
#' @param df A data frame containing temporal phenotypic data.
#' @param timeVar The name of the column representing time (default: "DAT").
#' @param idVar The name of the column identifying unique plants (default: "Plant.ID").
#' @param weekDATs A list of numeric vectors, each representing DAT values for one experimental week.
#' @param ignoreVars A character vector of column names to exclude from aggregation (e.g., metadata).
#'
#' @return A data frame containing aggregated weekly trait data with a column `valType`
#'         ("min", "mean", "max") and time column renamed to `timeVar`.
#'         If all trait values are constant across valTypes, only "mean" rows are kept.
#'
#' @examples
#' weekly_aggregate(my_data)
weekly_aggregate <- function(df,
                             timeVar    = "DAT",
                             idVar      = "Plant.ID",
                             weekDATs   = list(24:29, 30:36, 37:43, 44:49, 
                                               50:57, 58:64, 65:71, 72:78,
                                               79:85, 86:92, 93:98),
                             ignoreVars = c("Plant.Name", "Treatment")) {
  # 1) Tag each row with its corresponding week number
  #    This maps each DAT to the index of the list it belongs to in weekDATs
  df2 <- df %>%
    mutate(Week = purrr::map_int(.data[[timeVar]],
                                 ~ which(sapply(weekDATs, function(w) .x %in% w))[1]
    )) %>%
    filter(!is.na(Week))  # Drop rows that don't match any week
  
  # 2) Identify which columns are trait measurements (everything except meta-columns)
  traitCols <- setdiff(names(df2), c(timeVar, idVar, ignoreVars, "Week"))
  
  # 3) For each valType ("min", "mean", "max"), compute aggregated values by plant and week
  summaries <- lapply(c("min","mean","max"), function(valType) {
    fn <- match.fun(valType)
    df2 %>%
      group_by(across(all_of(c(idVar, ignoreVars, "Week")))) %>%
      summarise(across(all_of(traitCols), ~ fn(.x, na.rm = TRUE)),
                .groups = "drop") %>%
      mutate(valType = valType)
  })
  
  # Combine all summaries into one long-format data frame
  combined <- bind_rows(summaries)
  
  # 4) Pivot to wide format so that min, mean, max for each trait are in separate columns
  wide <- combined %>%
    pivot_wider(
      id_cols    = all_of(c(idVar, ignoreVars, "Week")),
      names_from = valType,
      values_from= all_of(traitCols),
      names_glue = "{.value}_{valType}"
    )
  
  # 5) Check globally (across the entire dataset) whether min == mean == max for all traits
  all_equal_global <- all(
    unlist(lapply(traitCols, function(tr) {
      min_col  <- paste0(tr, "_min")
      mean_col <- paste0(tr, "_mean")
      max_col  <- paste0(tr, "_max")
      all(wide[[min_col]] == wide[[mean_col]], na.rm = TRUE) &&
        all(wide[[max_col]] == wide[[mean_col]], na.rm = TRUE)
    }))
  )
  
  # 6) If all traits are constant across valTypes, keep only "mean" rows; else keep all
  cleaned <- if (all_equal_global) {
    combined %>% filter(valType == "mean")
  } else {
    combined
  }
  
  # 7) Rename "Week" column back to the original timeVar name (usually "DAT")
  cleaned %>%
    rename(!!timeVar := Week)
}


#' Widen and Merge a List of Weekly-Aggregated Data Frames
#'
#' Each df in DFlist must have:
#'  - an ID column (e.g. "Plant.ID")
#'  - any ignoreVars (e.g. "Plant.Name","Treatment") you want to carry over
#'  - a "DAT" column containing the week number (1,2,3…)
#'  - a "valType" column ("min","mean","max")
#'  - one column per trait
#'
#' @param DFlist     List of those data frames
#' @param mergeby    Character vector of columns to merge by (e.g. "Plant.ID")
#' @param ignoreVars Character vector of other columns to carry along (e.g. "Plant.Name","Treatment")
#' @return           A single data.frame with one row per ID, and columns named
#'                   trait_valType_W{DAT} for every combination found.
#' @examples
#' # Suppose agg1, agg2 are two weekly-aggregated tibbles:
#' widen_and_merge(list(agg1, agg2), mergeby="Plant.ID", ignoreVars=c("Plant.Name","Treatment"))
widen_and_merge_weekly <- function(DFlist,
                            mergeby    = "Plant.ID",
                            ignoreVars = c("Plant.Name","Treatment")) {
  # A helper that widens one df
  widen_one <- function(df) {
    # identify trait columns
    traitCols <- setdiff(names(df), c(mergeby, ignoreVars, "DAT", "valType"))
    
    df %>%
      pivot_wider(
        id_cols      = c(all_of(mergeby), all_of(ignoreVars)),
        names_from   = c("valType","DAT"),
        values_from  = traitCols,
        names_glue   = "{.value}_{valType}_{DAT}W"
      )
  }
  
  # 1) Widen each
  wide_list <- purrr::map(DFlist, widen_one)
  
  # 2) Merge all together by 'mergeby' + ignoreVars
  output <- reduce(wide_list, full_join, by = c(mergeby, ignoreVars))
  
  as.data.frame(output)
}


compute_stats = function(df, group_vars, trait_prefix = "", treatment_label = NULL, include_dat = TRUE) {
  # If we're NOT grouping by Treatment, drop it from the data before computing anything
  if (!"Treatment" %in% group_vars && "Treatment" %in% names(df)) {
    df = df %>% dplyr::select(-Treatment)
  }
  
  # Recalculate predictor columns after potentially dropping Treatment
  predictor_cols = setdiff(names(df), c(group_vars, ignoreVars))
  
  # Summarise: compute mean and sd
  df_grouped = df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(across(all_of(predictor_cols), list(mean = mean, sd = sd),
                     .names = "{.col}_{.fn}"),
              .groups = "drop")
  
  # Pivot to long format
  df_long = df_grouped %>%
    pivot_longer(
      cols = -all_of(group_vars),
      names_to = c("Trait", ".value"),
      names_pattern = "(.*)_(mean|sd)$"
    )
  
  # WideName and Treatment tagging
  if ("DAT" %in% group_vars && include_dat) {
    df_long = df_long %>%
      mutate(DAT = as.character(DAT),
             WideName = paste0(Trait, "_", DAT))
  } else {
    df_long = df_long %>%
      mutate(DAT = NA_character_,
             WideName = Trait)
  }
  
  if (!is.null(treatment_label)) {
    df_long$Treatment = treatment_label
  }
  
  df_long %>%
    dplyr::select(WideName, Treatment, DAT, Trait, Mean = mean, SD = sd)
}

