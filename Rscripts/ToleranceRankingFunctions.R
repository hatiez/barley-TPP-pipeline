# List of required packages
required_packages = c(
  "effectsize", "effsize", "ggplot2", "biotools", "vegan", "rcompanion",
  "ggpubr", "ggfortify", "scales", "vcd", "dplyr"
)

# Identify missing packages
missing_packages = required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# Prompt user to install missing packages
if (length(missing_packages) > 0) {
  message("The following packages are missing: ", paste(missing_packages, collapse = ", "))
  install = readline(prompt = "Do you want to install them now? [y/n]: ")
  if (tolower(install) == "y") {
    install.packages(missing_packages, dependencies = TRUE)
    remove(install)
  } else {
    stop("Required packages not installed. Script cannot continue.")
    remove(install)
  }
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

remove(missing_packages,required_packages)


harvest_testByTreatment = function(df, 
                                    genotypeCol = "Plant.Name", 
                                    genotypes = c("L1","L2","L3","L4","L5","L6","L7","L8","L9"), 
                                    treatmentCol = "Treatment", 
                                    ignoreCols = c("Treatment", "Plant.Name", "Plant.ID"),
                                    effect_size_method = "cohensD",
                                    test_method = "wilcox.test",
                                    p_adjust_method = "none") {
  
  # Identify the trait columns
  traitCols = setdiff(colnames(df), ignoreCols)
  
  # Initialize output data frames for effect sizes and p-values
  effect_size_table = data.frame(matrix(nrow = length(genotypes), ncol = length(traitCols)))
  test_statistic_table = data.frame(matrix(nrow = length(genotypes), ncol = length(traitCols)))
  colnames(effect_size_table) = colnames(test_statistic_table) = traitCols
  rownames(effect_size_table) = rownames(test_statistic_table) = genotypes
  
  # Store p-values for adjustment
  p_values = numeric()
  index_map = list()
  
  # Loop over each genotype and each trait to calculate the effect size and test statistic
  counter = 1
  for (genotype in genotypes) {
    for (trait in traitCols) {
      # Extract the data for the current genotype and trait
      control_values = df[df[[treatmentCol]] == "Control" & df[[genotypeCol]] == genotype, trait]
      drought_values = df[df[[treatmentCol]] == "Drought" & df[[genotypeCol]] == genotype, trait]
      
      if(length(control_values == length(drought_values))){
        # Calculate the effect size
        effect_size = switch(effect_size_method,
                             "rank_biserial" = rank_biserial(control_values, drought_values)$r_rank_biserial,
                             "cohensD" = cohen.d(control_values, drought_values)$estimate,
                             NA)
        
        # Calculate the test statistic (p-value)
        test_statistic = switch(test_method,
                                "wilcox.test" = wilcox.test(control_values, drought_values)$p.value,
                                "t.test" = t.test(control_values, drought_values)$p.value,
                                NA)
        
        # Store the p-value for correction later
        p_values = c(p_values, test_statistic)
        index_map[[counter]] = c(genotype, trait)
        counter = counter + 1
      }
      else{
        effect_size = NaN
        test_statistic = NaN
      }
      
      # Store the effect size in the output tables
      effect_size_table[genotype, trait] = effect_size
    }
  }
  
  # Adjust p-values if needed
  adjusted_p_values = p.adjust(p_values, method = p_adjust_method)
  
  # Map adjusted p-values back to the test_statistic_table
  for (i in seq_along(adjusted_p_values)) {
    genotype = index_map[[i]][1]
    trait = index_map[[i]][2]
    test_statistic_table[genotype, trait] = adjusted_p_values[i]
  }
  
  return(list(effect_sizes = effect_size_table, test_statistics = test_statistic_table))
}

significance_level = function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}


plot_effect_size_barplot = function(effect_size_df, 
                                    test_statistic_df, 
                                    trait, 
                                    include_significance = TRUE) {
  
  # Check if the specified trait is in the data frames
  if (!(trait %in% colnames(effect_size_df))) {
    stop(paste("Trait", trait, "not found in the data frame"))
  }
  
  # Extract and prepare data for the specified trait
  plot_data = data.frame(
    Genotype = rownames(effect_size_df),
    EffectSize = effect_size_df[, trait],
    PValue = test_statistic_df[, trait]
  )
  
  # Remove rows with NA values
  plot_data = na.omit(plot_data)
  
  # Order the data by increasing effect size
  plot_data = plot_data[order(plot_data$EffectSize), ]
  
  # Determine significance levels
  if (include_significance) {
    plot_data$Significance = sapply(plot_data$PValue, function(p) {
      significance_level(p)
    })
  } else {
    plot_data$Significance = ""
  }
  
  # Calculate max/min effect size for proper y-axis scaling
  max_effect_size = max(plot_data$EffectSize, na.rm = TRUE)
  min_effect_size = min(plot_data$EffectSize, na.rm = TRUE)
  if (min_effect_size > 0) min_effect_size = 0
  
  # Plotting
  p = ggplot(plot_data, aes(x = reorder(Genotype, EffectSize), y = EffectSize)) +
    geom_bar(stat = "identity", fill = "gray80", color = "black") +  # Neutral fill
    geom_text(aes(label = Significance), vjust = -0.3, size = 5, color = "black") +
    scale_y_continuous(limits = c(min_effect_size, max_effect_size * 1.2)) +  # Extend y-axis
    theme_classic() +
    labs(title = paste("Effect Sizes for", trait),
         x = "Line",
         y = "Effect Size") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}



PERMANOVA_per_Genotype = function(df, 
                                   genotypeCol = "Plant.Name", 
                                   genotypes = c("L1","L2","L3","L4","L5","L6","L7","L8","L9"), 
                                   treatmentCol = "Treatment", 
                                   ignoreCols = c("Treatment", "Plant.Name", "Plant.ID"),
                                   p_adjust_method = "none",
                                   permutations = 3000) {
  # Load the vegan package if not already loaded
  if (!require("vegan")) install.packages("vegan", dependencies = TRUE)
  library(vegan)
  
  # Identify the trait columns
  traitCols = setdiff(colnames(df), ignoreCols)
  
  # Initialize output data frames for effect sizes and p-values
  effect_size_table = data.frame(matrix(nrow = length(genotypes), ncol = 1))  # One column for PERMANOVA effect size
  test_statistic_table = data.frame(matrix(nrow = length(genotypes), ncol = 1))  # One column for PERMANOVA p-value
  colnames(effect_size_table) = colnames(test_statistic_table) = c("PERMANOVA")
  rownames(effect_size_table) = rownames(test_statistic_table) = genotypes
  
  # Store p-values for adjustment
  p_values = numeric()
  index_map = list()
  
  # Loop over each genotype to perform PERMANOVA and calculate the effect size and test statistic
  counter = 1
  for (genotype in genotypes) {
    print(genotype)
    # Extract the data for the current genotype
    genotype_data = df[df[[genotypeCol]] == genotype, ]
    
    # Ensure the data is balanced (contains both "Control" and "Drought")
    if (length(unique(genotype_data[[treatmentCol]])) < 2) {
      effect_size_table[genotype, "PERMANOVA"] = NA
      test_statistic_table[genotype, "PERMANOVA"] = NA
      next
    }
    
    # Ensure all dependent variables are numeric
    numeric_data = as.matrix(genotype_data[, traitCols])
    numeric_data = apply(numeric_data, 2, as.numeric)  # Convert to numeric
    
    # Combine numeric matrix back with the Treatment column in a data frame
    genotype_data_numeric = data.frame(numeric_data, Treatment = genotype_data[[treatmentCol]])
    
    # Perform PERMANOVA using adonis2
    perma_result = adonis2(numeric_data ~ Treatment, data = genotype_data_numeric, method = "euclidean", permutations = permutations)
    
    # Extract R² as effect size and p-value from the PERMANOVA results
    R2 = perma_result$R2[1]  # R² for the Treatment term
    p_value = perma_result$`Pr(>F)`[1]  # p-value for the Treatment term
    
    # Store effect size (R²) and p-value
    effect_size_table[genotype, "PERMANOVA"] = R2
    test_statistic_table[genotype, "PERMANOVA"] = p_value
    
    # Store the p-value for correction later
    p_values = c(p_values, p_value)
    index_map[[counter]] = genotype
    counter = counter + 1
  }
  
  # Adjust p-values if needed
  adjusted_p_values = p.adjust(p_values, method = p_adjust_method)
  
  # Map adjusted p-values back to the test_statistic_table
  for (i in seq_along(adjusted_p_values)) {
    genotype = index_map[[i]]
    test_statistic_table[genotype, "PERMANOVA"] = adjusted_p_values[i]
  }
  
  return(list(effect_sizes = effect_size_table, test_statistics = test_statistic_table))
}


plot_permanova_effect_size_barplot = function(effect_size_df, 
                                              test_statistic_df, 
                                              include_significance = TRUE,
                                              main = "Effect Sizes from MANOVA by Genotype") {
  
  # Extract and prepare data
  plot_data = data.frame(
    Genotype = rownames(effect_size_df),
    EffectSize = effect_size_df$PERMANOVA,
    PValue = test_statistic_df$PERMANOVA
  )
  
  # Remove rows with NA values
  plot_data = na.omit(plot_data)
  
  # Stop if there's no valid data
  if (nrow(plot_data) == 0) {
    stop("No valid data to plot. All effect sizes are NA or missing.")
  }
  
  # Order the data by increasing effect size
  plot_data = plot_data[order(plot_data$EffectSize), ]
  
  # Add significance labels
  if (include_significance) {
    plot_data$Significance = sapply(plot_data$PValue, function(p) {
      if (p < 0.001) return("***")
      else if (p < 0.01) return("**")
      else if (p < 0.05) return("*")
      else return("")
    })
  } else {
    plot_data$Significance = ""
  }
  
  # Axis scaling
  max_effect_size = max(plot_data$EffectSize, na.rm = TRUE)
  min_effect_size = min(plot_data$EffectSize, na.rm = TRUE)
  if (min_effect_size > 0) min_effect_size = 0
  
  # Plotting
  p = ggplot(plot_data, aes(x = reorder(Genotype, EffectSize), y = EffectSize)) +
    geom_bar(stat = "identity", fill = "lightblue", color = "black") +
    geom_text(aes(label = Significance), vjust = -0.3, size = 5, color = "black") +
    scale_y_continuous(limits = c(min_effect_size, max_effect_size * 1.2)) +
    theme_classic() +
    labs(title = main,
         x = "Line",
         y = "Effect Size (generalized Eta Squared)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}





check_normality = function(data, variable, group1, group2) {
  # Create an empty list to store the results
  results = list()
  
  # Combine the grouping factors into a single group
  data$group = interaction(data[[group1]], data[[group2]])
  
  # Get the unique combinations of the groups
  unique_groups = unique(data$group)
  
  # Loop over each group and perform the Shapiro-Wilk test
  for (g in unique_groups) {
    # Subset the data for this group
    group_data = subset(data, group == g)[[variable]]
    
    # Check if there are enough data points to perform the Shapiro-Wilk test
    if (length(group_data) >= 3) {
      # Perform the Shapiro-Wilk test
      shapiro_test = shapiro.test(group_data)
      
      # Store the p-value along with the group name
      results[[as.character(g)]] = list(
        "p_value" = shapiro_test$p.value,
        "is_gaussian" = ifelse(shapiro_test$p.value > 0.05, TRUE, FALSE)
      )
    } else {
      # If not enough data, store a message
      results[[as.character(g)]] = "Not enough data for Shapiro-Wilk test"
    }
  }
  
  return(results)
}






