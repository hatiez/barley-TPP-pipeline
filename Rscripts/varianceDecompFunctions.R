# Define required packages
.required_pkgs = c("reshape2", "lme4", "ggplot2", "lmerTest", "patchwork", "purrr", "tibble")

# Identify missing ones
.missing_pkgs = .required_pkgs[!(.required_pkgs %in% installed.packages()[, "Package"])]

# Prompt to install if necessary
if (length(.missing_pkgs) > 0) {
  message("The following packages are missing: ", paste(.missing_pkgs, collapse = ", "))
  .user_input = readline(prompt = "Do you want to install them now? [y/n]: ")
  if (tolower(.user_input) == "y") {
    install.packages(.missing_pkgs, dependencies = TRUE)
    remove( .user_input)
  } else {
    stop("Missing required packages. Script execution halted.")
  }
}

# Load all packages
invisible(lapply(.required_pkgs, library, character.only = TRUE))

# Clean up temporary variables
rm(.required_pkgs, .missing_pkgs)

# Function for making predictor traits more readable
clean_label_TemporalTraits = function(x) {
  # Remove the word "plant_Avg" (case-insensitive)
  x = gsub("_(?i)plant_Avg_", " ", x, perl = TRUE)
  # remove MM unit
  x = gsub("_MM_", "_", x)
  # Replace underscore-number endings with ", DAT xx"
  x = gsub("_(\\d+)$", ", \\1 DAT", x, perl = TRUE)
  # Remove the word "plant" (case-insensitive)
  x = gsub("_(?i)plant_", " ", x, perl = TRUE)
  # Replace underscores with spaces
  x = gsub("_", " ", x)
  # Replace literal dots with spaces
  x = gsub("\\.", " ", x)
  # Remove trailing "MM" (case-sensitive) for unit in RGB
  x = gsub("MM$", "", x)
  # Trim extra spaces
  x = trimws(x)
  # Capitalize only the first letter
  x = sub("^(\\w)", "\\U\\1", x, perl = TRUE)
  # Reintroduce dots between rgb color codes
  x = gsub("138 139 97", "138.139.97", x)
  x = gsub("68 85 63", "68.85.63", x)
  # Take care of caps-locked RGB traits
  x = gsub("AREA", "area", x)
  x = gsub("HEIGHT", "height", x)
  x = gsub("WIDTH", "width", x)
  x = gsub("COMPACTNESS", "compactness", x)
  x = gsub("PERIMETER", "perimeter", x)
  x = gsub("(?i)deltaT", "CTD", x, perl = TRUE)
  
  return(x)
}

# Function for making harvest trait names more readable
clean_label_HarvestTraits= function(x) {
  # Define replacement mapping as a named vector
  replacements = c(
    "Tiller.number.per.pot" = "Tiller No",
    "Spike.number.per.pot" = "Spike No",
    "Total.Biomass.DW" = "Biomass DW",
    "Total.Spike.weight" = "Total Spike W",
    "N5.spike.weight" = "Spike W.5Spk",
    "Grain.number.per.5.spikes" = "Grain No.5Spk",
    "Grain.weight.per.5.spikes" = "Grain W.5Spk",
    "Overall.spike.length" = "Spike L overall",
    "spikelet.density.to.overall.spike.length" = "Spkl D overall",
    "Spike.length.to.base.of.top.seed" = "Spike L Top",
    "spikelet.density.to.base.of.top.spike" = "Spkl D Top",
    "Infertal.tip.length" = "Infert Tip L",
    "N20.seed.weight" = "N20 Grain W",
    "N.spikelets" = "Spkl No",
    "N.Seeds.per.spike" = "Grain No.1Spk",
    "Seeds.per.spike" = "Grain Spk ratio"
  )
  
  # Replace labels using the mapping
  x = ifelse(x %in% names(replacements), replacements[x], x)
  
  return(x)
}


#' Compute Temporal Repeatability
#'
#' This function computes the temporal repeatability (TR) of a given variable using mixed-effects models.
#' Temporal repeatability is calculated based on variance components extracted from the model.
#'
#' @param data The input data frame containing the measurements.
#' @param variable The variable name (as a string) for which temporal repeatability is to be computed.
#' @param n_replicates_range A numeric vector of length 2 specifying the range of the number of replicates to consider for the TR calculation.
#' @param return_model Logical; if TRUE, the function returns the mixed-effects model along with the TR values.
#' @return A list containing TR_min, TR_max, and optionally the fitted model.
#' @examples
#' compute_temporal_repeatability(data = df, variable = "Height")
compute_temporal_repeatability = function(data, variable, n_replicates_range = c(8,20), return_model = FALSE) {
  formula = as.formula(paste(variable, "~ 1 + (1 | DAT) + (1 | DAT:Treatment) + (1 | DAT:Plant.Name)"))
  
  # Fit the mixed model
  model = lmer(formula, data = data)
  
  # Extract variance components
  var_components = VarCorr(model)
  
  # Extract the specific variances
  genotypic_variance = as.numeric(var_components$`DAT:Plant.Name`[1])
  residual_variance = attr(var_components, "sc")^2
  
  # Calculate the number of replications
  replications = length(unique(data$DAT))
  
  # Compute temporal repeatability (TR)
  TR_min = genotypic_variance / (genotypic_variance + (residual_variance/min(n_replicates_range)))
  TR_max = genotypic_variance / (genotypic_variance + (residual_variance/max(n_replicates_range)))
  
  # Return the result
  if (return_model) {
    return(list(TR_min = TR_min, TR_max = TR_max, model = model))
  } else {
    return(list(TR_min = TR_min, TR_max = TR_max))
  }
}

#' Temporal Repeatability Table
#'
#' This function generates a table summarizing the variance components and temporal repeatability for multiple variables.
#' The table includes the variance components for DAT, Plant.Name, Treatment, Residual, and the computed temporal repeatability.
#'
#' @param data The input data frame containing the measurements.
#' @param excludeVars A vector of variable names to exclude from the analysis (e.g., "DAT", "Plant.Name", "Treatment").
#' @param varCompByPercentage Logical; if TRUE, the variance components are expressed as a percentage of the total variance.
#' @return A data frame summarizing the variance components and temporal repeatability for each variable.
#' @examples
#' temporal_repeatability_table(data = df)
temporal_repeatability_table = function(data, excludeVars = c("DAT", "Plant.Name", "Plant.ID", "Treatment", "Replicate"), varCompByPercentage = TRUE) {
  # Get list of variables
  measurement_vars = setdiff(names(data), excludeVars)
  n_variables = length(measurement_vars)
  
  # Create results df
  results = data.frame(
    Variable = measurement_vars,
    Plant.Name.var = rep(0, n_variables),
    DAT.var = rep(0, n_variables),
    Treatment.var = rep(0, n_variables),
    Residual.var = rep(0, n_variables),
    TR_min = rep(0, n_variables),
    TR_max = rep(0, n_variables)
  )
  
  for (i in 1:n_variables) {
    MM = compute_temporal_repeatability(data = data, variable = measurement_vars[i], n_replicates_range = c(8, 20), return_model = TRUE)
    
    # Get and insert variance components
    variance_comps = VarCorr(MM$model)
    results$DAT.var[i] = variance_comps$`DAT`[1]
    results$Plant.Name.var[i] = variance_comps$`DAT:Plant.Name`[1]
    results$Treatment.var[i] = variance_comps$`DAT:Treatment`[1]
    results$Residual.var[i] = attr(variance_comps, "sc")^2
    
    # If requested, give them as percentage
    if (varCompByPercentage) {
      results[i, 2:5] = 100 * (results[i, 2:5] / sum(results[i, 2:5]))
    }
    
    # Insert temporal repeatability
    results$TR_min[i] = MM$TR_min
    results$TR_max[i] = MM$TR_max
  }
  
  return(results)
}


#' Plot Temporal Repeatability
#'
#' This function generates a stacked bar plot showing the explained variance components 
#' for multiple temporal traits. It optionally overlays error bars indicating the range 
#' of temporal repeatability (TR_min to TR_max) and supports displaying a secondary axis.
#'
#' @param data A data frame containing variance components and temporal repeatability values.
#' @param varCol The name of the column in `data` representing the trait/variable names. Default is "Variable".
#' @param TRminCol The column name representing the minimum temporal repeatability. Default is "TR_min".
#' @param TRmaxCol The column name representing the maximum temporal repeatability. Default is "TR_max".
#' @param legendLabels Optional vector of labels for the legend. If NULL, default component names are used.
#' @param main The main title for the plot. If NULL, an empty title is used.
#' @param cex.axis Numeric scaling factor for axis text size. Default is 0.9.
#' @param include_error_bars Logical; if TRUE, error bars are drawn from TR_min to TR_max. Default is FALSE.
#' @param include_secondary_axis Logical; if TRUE, a secondary y-axis showing repeatability is added. Default is FALSE.
#'
#' @return A ggplot object displaying the variance decomposition and repeatability.
#'
#' @examples
#' plot_temporal_repeatability(data = results_df)
plot_temporal_repeatability = function(data, 
                                        varCol = "Variable", 
                                        TRminCol = "TR_min", 
                                        TRmaxCol = "TR_max", 
                                        legendLabels = NULL, 
                                        main = NULL, 
                                        cex.axis = 0.9,
                                        include_error_bars = FALSE,
                                        include_secondary_axis = FALSE) {
    

  # Convert the data to long format for ggplot
  data_long = melt(data, id.vars = c(varCol, TRminCol, TRmaxCol), 
                    variable.name = "Var.Component", 
                    value.name = "Percent.Variation")
  
  # Set factor levels based on input order
  data_long[[varCol]] = factor(data_long[[varCol]], levels = unique(data[[varCol]]))
  
  if (is.null(main)) {
    main = ""
  }
  
  # Default legend labels if none provided
  if (is.null(legendLabels)) {
    legendLabels = unique(data_long$Var.Component)
  }
  
  # Define bar colors
  bar_colors = c("lightgreen", "lightblue", "yellow", "orange")
  
  # Base plot
  p = ggplot(data_long, aes_string(x = varCol, y = "Percent.Variation", fill = "Var.Component")) +
    geom_bar(stat = "identity", position = "stack", width = 0.7)
  
  # Add error bars if enabled
  if (include_error_bars) {
    p = p + geom_errorbar(aes_string(x = varCol, ymin = paste0(TRminCol, "*100"), ymax = paste0(TRmaxCol, "*100")), 
                           color = "black", width = 0.2, size = 0.6, show.legend = FALSE)
  }
  
  # Add secondary y-axis if enabled
  if (include_secondary_axis) {
    p = p + scale_y_continuous(sec.axis = sec_axis(~ . / 100, name = "Temporal repeatability"), 
                                name = "Variance component (%)")
  } else {
    p = p + scale_y_continuous(name = "Variance component (%)")
  }
  
  p = p + scale_fill_manual(values = bar_colors, labels = legendLabels) +
    theme_classic() +  
    theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 10 * cex.axis, color = "black"),  
          axis.text.y = element_text(size = 10),  
          axis.title.x = element_text(size = 12, margin = ggplot2::margin(t = 15)),  
          axis.title.y = element_text(size = 12, margin = ggplot2::margin(r = 10)),  
          axis.title.y.right = element_text(size = 12, margin = ggplot2::margin(l = 10)),  
          plot.title = element_text(size = 14, hjust = 0.5),  
          legend.position = "right",  
          legend.background = element_rect(fill = "white", color = "black"),
          legend.text = element_text(size = 11)) +  
    labs(title = clean_label_TemporalTraits(main), x = "Temporal Trait", fill = "Var. Component") +
    guides(fill = guide_legend(override.aes = list(shape = NA)))
  
  p = p + scale_x_discrete(labels = Vectorize(clean_label_TemporalTraits))
  print(p)
}





#' Compute Harvest Repeatability
#'
#' This function computes the repeatability (R) of a given variable using mixed-effects models for harvest data.
#' Repeatability is calculated based on variance components extracted from the model.
#'
#' @param data The input data frame containing the harvest measurements.
#' @param variable The variable name (as a string) for which repeatability is to be computed.
#' @param n_replicates_range A numeric vector of length 2 specifying the range of the number of replicates to consider for the R calculation.
#' @param return_model Logical; if TRUE, the function returns the mixed-effects model along with the R values.
#' @return A list containing R_min, R_max, and optionally the fitted model.
#' @examples
#' harvest_repeatability(data = df, variable = "Yield")
harvest_repeatability = function(data, variable, n_replicates_range = c(8,20), return_model = FALSE) {
  formula = as.formula(paste(variable, "~ 1 + (1 | Treatment) + (1 | Plant.Name)"))
  
  # Fit the mixed model
  model = lmer(formula, data = data)
  
  # Extract variance components
  var_components = VarCorr(model)
  
  # Extract the specific variances
  genotypic_variance = as.numeric(var_components$`Plant.Name`[1])
  residual_variance = attr(var_components, "sc")^2
  
  # Compute repeatability
  R_min = genotypic_variance / (genotypic_variance + (residual_variance/min(n_replicates_range)))
  R_max = genotypic_variance / (genotypic_variance + (residual_variance/max(n_replicates_range)))
  
  # Return the result
  if (return_model) {
    return(list(R_min = R_min, R_max = R_max, model = model))
  } else {
    return(list(R_min = R_min, R_max = R_max))
  }
}


#' Harvest Repeatability Table
#'
#' This function generates a table summarizing the variance components and repeatability for multiple harvest variables.
#' The table includes the variance components for Plant.Name, Treatment, Residual, and the computed repeatability.
#'
#' @param data The input data frame containing the harvest measurements.
#' @param excludeVars A vector of variable names to exclude from the analysis (e.g., "DAT", "Plant.Name", "Treatment").
#' @param n_replicates_range A numeric vector specifying the range of the number of replicates to consider for the repeatability calculation.
#' @param varCompByPercentage Logical; if TRUE, the variance components are expressed as a percentage of the total variance.
#' @return A data frame summarizing the variance components and repeatability for each harvest variable.
#' @examples
#' harvest_repeatability_table(data = df)
harvest_repeatability_table = function(data, excludeVars = c("DAT", "Plant.Name", "Plant.ID", "Treatment", "Replicate"), n_replicates_range = c(8,20), varCompByPercentage = TRUE) {
  # Get list of variables
  measurement_vars = setdiff(names(data), excludeVars)
  n_variables = length(measurement_vars)
  
  # Create results df
  results = data.frame(
    Variable = measurement_vars,
    Plant.Name.var = rep(0, n_variables),
    Treatment.var = rep(0, n_variables),
    Residual.var = rep(0, n_variables),
    R_min = rep(0, n_variables),
    R_max = rep(0, n_variables)
  )
  
  for (i in 1:n_variables) {
    MM = harvest_repeatability(data = data, variable = measurement_vars[i], n_replicates_range = n_replicates_range, return_model = TRUE)
    
    # Get and insert variance components
    variance_comps = VarCorr(MM$model)
    results$Plant.Name.var[i] = variance_comps$`Plant.Name`[1]
    results$Treatment.var[i] = variance_comps$`Treatment`[1]
    results$Residual.var[i] = attr(variance_comps, "sc")^2
    
    # If requested, give them as percentage
    if (varCompByPercentage) {
      results[i, 2:5] = 100 * (results[i, 2:5] / sum(results[i, 2:5]))
    }
    
    # Insert repeatability
    results$R_min[i] = MM$R_min
    results$R_max[i] = MM$R_max
  }
  
  return(results)
}


#' Plot Harvest Repeatability
#'
#' This function generates a bar plot showing the explained variance components and repeatability for multiple harvest variables.
#' The bar plot includes error bars representing the range of repeatability (R_min to R_max).
#'
#' @param data The data frame containing variance components and repeatability values.
#' @param varCol The column name in `data` representing the variables (default: "Variable").
#' @param RminCol The column name in `data` representing the minimum repeatability (default: "R_min").
#' @param RmaxCol The column name in `data` representing the maximum repeatability (default: "R_max").
#' @param legendLabels A vector of labels for the legend; if NULL, default labels are used.
#' @param main The main title of the plot; if NULL, a default title is used.
#' @return A ggplot object representing the bar plot.
#' @examples
#' plot_harvest_repeatability(data = df)
plot_harvest_repeatability = function(data, 
                                       varCol = "Variable", 
                                       RminCol = "R_min", 
                                       RmaxCol = "R_max", 
                                       legendLabels = NULL, 
                                       main = NULL, 
                                       cex.axis = 0.9) {
  # Convert the data to long format for ggplot
  data_long = melt(data, id.vars = c(varCol, RminCol, RmaxCol), 
                    variable.name = "Var.Component", 
                    value.name = "Percent.Variation")
  
  # Set the factor levels of varCol according to the order in the input data
  data_long[[varCol]] = factor(data_long[[varCol]], levels = unique(data[[varCol]]))
  
  if (is.null(main)) {
    main = ""
  }
  
  # Default legend labels if none provided
  if (is.null(legendLabels)) {
    legendLabels = unique(data_long$Var.Component)
  }
  
  # Define bar colors
  bar_colors = c("lightgreen", "yellow", "orange")
  
  # Create plot
  p = ggplot(data_long, aes_string(x = varCol, y = "Percent.Variation", fill = "Var.Component")) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Adjust width to align with x-axis
    geom_errorbar(aes_string(x = varCol, ymin = paste0(RminCol, "*100"), ymax = paste0(RmaxCol, "*100")), 
                  color = "gray40", width = 0.2, size = 0.6, show.legend = FALSE) +
    scale_y_continuous(sec.axis = sec_axis(~ . / 100, name = "Repeatability"), 
                       name = "Variance component (%)") +
    scale_fill_manual(values = bar_colors, labels = legendLabels) +
    
    # Base R-like theme with proper spacing
    theme_classic() +  
    theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 10 * cex.axis),  # Fix label alignment
          axis.text.y = element_text(size = 10),  
          axis.title.x = element_text(size = 12, margin = ggplot2::margin(t = 15)),  # Space between x-axis title and labels
          axis.title.y = element_text(size = 12, margin = ggplot2::margin(r = 10)),  # Space between left y-axis title and labels
          axis.title.y.right = element_text(size = 12, margin = ggplot2::margin(l = 10)),  # Space for right y-axis title
          plot.title = element_text(size = 14, hjust = 0.5),  
          legend.position = "right",  # Restore legend
          legend.background = element_rect(fill = "white", color = "black")) +  
    
    labs(title = clean_label_HarvestTraits(main), x = "Trait", fill = "Var. Component") +
    guides(fill = guide_legend(override.aes = list(shape = NA)))  # Keep legend formatting
  
  p = p + scale_x_discrete(labels = clean_label_HarvestTraits)
  
  print(p)
}






#' Fit Random Regression Models for Longitudinal Traits
#'
#' This function fits random regression mixed-effects models for multiple traits measured over time. 
#' It compares models with and without random slopes to assess whether slope variability is statistically significant.
#' Optionally, it includes treatment effects in the model structure.
#'
#' @param data A data frame containing longitudinal trait data, including time, ID, and trait variables.
#' @param timeVar A string specifying the column name for time (default: "DAT").
#' @param idVar A string specifying the column name for the subject or plant identifier (default: "Plant.ID").
#' @param excludeVars A character vector of column names to exclude from trait modeling (default includes ID and meta-data).
#' @param pAdjustMethod A string specifying the method for p-value adjustment (default: "BH").
#' @param includeTreatment Logical; if TRUE, includes random slopes by treatment in the model (default: FALSE).
#' @param treatmentVar A string specifying the column name for the treatment variable (default: "Treatment").
#'
#' @return A list with two elements:
#' \describe{
#'   \item{results}{A data frame with raw and adjusted p-values for each trait, indicating slope significance.}
#'   \item{plots}{A named list of ggplot objects for visualizing fitted models per trait.}
#' }
#'
#' @examples
#' \dontrun{
#' rr_output <- fit_random_regression(data = growth_df, includeTreatment = TRUE)
#' rr_output$results  # View significance of slope variation
#' rr_output$plots[["TraitName"]]  # View plot for a specific trait
#' }
fit_random_regression = function(data, 
                                            timeVar = "DAT", 
                                            idVar = "Plant.ID", 
                                            excludeVars = c("DAT", "Plant.Name", "Plant.ID", "Replicate"), 
                                            pAdjustMethod = "BH", 
                                            includeTreatment = FALSE, 
                                            treatmentVar = "Treatment") {
  # Ensure Treatment is a factor if included
  if (includeTreatment) {
    data[[treatmentVar]] = as.factor(data[[treatmentVar]])
  }
  
  # Determine traits and convert to numeric
  exclude_all = c(excludeVars, timeVar, idVar, treatmentVar)
  traits = setdiff(names(data), exclude_all)
  data[traits] = lapply(data[traits], function(x) as.numeric(as.character(x)))
  # If trait should be scaled
  #data[traits] = lapply(data[traits], scale)
  
  
  # Prepare results dataframe
  results = data.frame(Trait = traits, 
                        p_value_raw = rep(NA, length(traits)), 
                        p_value_adj = rep(NA, length(traits)), 
                        stringsAsFactors = FALSE)
  
  # List for storing plots
  plot_list = list()
  
  # Define custom colors
  custom_colors = c("Drought" = "orange", "Control" = "lightblue")
  
  # Calculate raw p-values via likelihood ratio tests
  p_values = numeric(length(traits))
  for (i in seq_along(traits)) {
    trait = traits[i]
    
    # Define model formulas
    if (includeTreatment) {
      full_formula = as.formula(paste(trait, "~", timeVar, 
                                       "+ (1 +", timeVar, "|", idVar, ")",
                                       "+ (1 +", timeVar, "|", treatmentVar, ")"))
      reduced_formula = as.formula(paste(trait, "~", timeVar, 
                                          "+ (1 |", idVar, ")",
                                          "+ (1 +", timeVar, "|", treatmentVar, ")"))
    } else {
      full_formula = as.formula(paste(trait, "~", timeVar, "+ (1 +", timeVar, "|", idVar, ")"))
      reduced_formula = as.formula(paste(trait, "~", timeVar, "+ (1 |", idVar, ")"))
    }
    
    # Fit models
    model_full = lmerTest::lmer(full_formula, data = data, REML = FALSE)
    model_reduced = lmerTest::lmer(reduced_formula, data = data, REML = FALSE)
    model_comparison = anova(model_reduced, model_full)
    p_value_slope = model_comparison$`Pr(>Chisq)`[2]
    p_values[i] = p_value_slope
    
    # Store raw p-value
    results$p_value_raw[i] = p_value_slope
  }
  
  # Adjust p-values
  results$p_value_adj = p.adjust(p_values, method = pAdjustMethod)
  
  # Generate plots with adjusted p-values
  for (i in seq_along(traits)) {
    trait = traits[i]
    
    # Create a separate data frame to avoid overwriting issues
    plot_data = data[, c(timeVar, idVar, treatmentVar, trait)]
    plot_data$Predicted = predict(lmerTest::lmer(as.formula(paste(trait, "~", timeVar, "+ (1 +", timeVar, "|", idVar, ")")), 
                                                  data = data, REML = FALSE), re.form = NULL)
    
    # Format adjusted p-value for annotation
    p_adjust_label = ifelse(results$p_value_adj[i] < 0.001, "Slope p < 0.001", 
                             paste0("Slope p = ", signif(results$p_value_adj[i], 3)))
    
    # Determine the position for annotation (top-right corner)
    x_max = max(plot_data[[timeVar]], na.rm = TRUE)  # Max time value
    y_max = max(plot_data[[trait]], na.rm = TRUE)    # Max trait value
    
    # Create plot with annotation
    p = ggplot(plot_data, aes_string(x = timeVar, y = trait, group = idVar, color = treatmentVar)) +
      geom_point(alpha = 0.6, size = 1) +  # Raw data points, now colored by treatment
      geom_line(aes(y = Predicted), size = 0.8) +  # Treatment-based colored slopes
      scale_color_manual(values = custom_colors) +  # Custom colors for treatments
      labs(title = paste("Random Regression Model for", trait),
           x = "Time", y = trait,
           color = "Treatment") +  # Legend label
      annotate("text", x = x_max, y = y_max, label = p_adjust_label, hjust = 1.1, vjust = 1.2, size = 4.5) +  # Add p-value annotation
      theme_classic() +
      theme(legend.position = "right")
    
    print(p)  # Display the plot
    plot_list[[trait]] = p  # Store for later use
  }
  
  return(list(results = results, plots = plot_list))
}








