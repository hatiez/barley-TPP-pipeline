# Load functions
source("HarvestPredictionEvaluationFunctions.R")

if (!dir.exists("../Rout/")) {
  dir.create("../Rout/", recursive = TRUE)
}

# Load RF and LASSO accuracy data frames
load("../Rout/RFHarvestPrediction/combinedRFAccuracies.Rdata")
RFAccs = combinedRFAccuracies
rm(combinedRFAccuracies)

load("../Rout/LassoHarvestPrediction/combinedLASSOAccuracies.Rdata")
lassoAccs = combinedLASSOAccuracies
rm(combinedLASSOAccuracies)

# ------------------------------------------------------------------------------
# Make master accuracy df
allAccs = bind_rows(
  lapply(list(lassoAccs = lassoAccs, RFAccs = RFAccs), function(list_item) bind_rows(list_item, .id = "OutcomeVar")),
  .id = "ListName"
)
# remove redundant listName column and resample column
allAccs$ListName = NULL
allAccs$Resample = NULL

# Make a table with mean accuracy values
means = aggregate(cbind(RMSE, Rsquared, MAE) ~ OutcomeVar + Algorithm + DataGroup + Subsetting, data = allAccs, FUN = mean)

# save excel
if (!dir.exists("../Rout/tables")) {
  dir.create("../Rout/tables", recursive = TRUE)
}
write.xlsx(means,"../Rout/tables/predictionAccuracyMeans.xlsx")


# ------------------------------------------------------------------------------
# Define the desired order for boxplots
desired_order = c(
  "Total.Biomass.DW", 
  "Spike.number.per.pot",
  "Total.Spike.weight", 
  "N5.spike.weight", 
  "Tiller.number.per.pot"
)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ACCURACY PLOTS

# ------------------------------------------------------------------------------
# OVERVIEW OF ALL OUTCOMES 

# Define algorithms and treatment groups
algorithms = c("LASSO", "RF")
treatmentGroups = c("control", "drought", "pooled")

# Define the output directory
output_dir = "../Rout/modelComparison/"

if (!dir.exists("../Rout/modelComparison")) {
  dir.create("../Rout/modelComparison", recursive = TRUE)
}

# Initialize an empty list to store ggplot objects
plot_list = list()

# Loop over both algorithms and treatment groups
for (algorithm in algorithms) {
  for (treatment in treatmentGroups) {
    
    # Subset data
    temp = allAccs[allAccs$Algorithm == algorithm & 
                      allAccs$Subsetting == "None" & 
                      allAccs$DataGroup == treatment, ]
    
    # Create a new factor with levels reordered
    temp$OutcomeVar = factor(temp$OutcomeVar, 
                                 levels = c(desired_order, 
                                            setdiff(unique(temp$OutcomeVar), desired_order)))
    
    # Reorder rows based on the new factor order
    temp = temp[order(temp$OutcomeVar), ]
    
    
    # Clean names
    temp = temp %>%
      mutate(OutcomeVar = clean_label_HarvestTraits(OutcomeVar))
    
    # Generate plot
    p = ggboxplot(temp, x = "OutcomeVar", y = "Rsquared",
                   fill = I("lightblue"),   # Set a constant fill color (lightblue)
                   add = "none") +
      labs(title = "",
           x = "Outcome Variable",
           y = expression(R^2)) +
      theme_classic() +
      theme(
        axis.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1),  # Rotate x-ticks
        axis.title = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = ggplot2::margin(t = 15, r = 0, b = 0, l = 0)), 
        plot.title = element_text(hjust = 0.5, size = 16),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.6)
      )
    
    # Define output filename dynamically
    output_filename = paste0(output_dir, "allResponses", algorithm, toTitleCase(treatment), ".pdf")
    
    # Save plot to PDF
    pdf(output_filename, width = 10, height = 6)
    print(p)
    dev.off()
    
    # Store the plot in a list using a structured key
    plot_key = paste0(algorithm, "_", treatment)
    plot_list[[plot_key]] = p  # Store ggplot object
    
    # Print status message
    message("Saved: ", output_filename)
  }
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# COMPARISON PLOTS
# ------------------------------------------------------------------------------
# DATA FOR COMPARISONS
# get variables of interest 
SubsetCompAccs = allAccs[allAccs$OutcomeVar %in% c("Total.Biomass.DW","Total.Spike.weight","N5.spike.weight","Spike.number.per.pot"),]

desired_order = c("Total.Biomass.DW","Total.Spike.weight","N5.spike.weight","Spike.number.per.pot")

# Create a new factor with levels reordered
SubsetCompAccs$OutcomeVar = factor(SubsetCompAccs$OutcomeVar, 
                          levels = c(desired_order, 
                                     setdiff(unique(SubsetCompAccs$OutcomeVar), desired_order)))

# Reorder rows based on the new factor order
SubsetCompAccs = SubsetCompAccs[order(SubsetCompAccs$OutcomeVar), ]


# Clean names
SubsetCompAccs = SubsetCompAccs %>%
  mutate(OutcomeVar = clean_label_HarvestTraits(OutcomeVar))


#-------------------------------------------------------------------------------
# ALGORITHM COMPARISONS 

temp = SubsetCompAccs[SubsetCompAccs$Subsetting == "None",]

pdf("../Rout/modelComparison/fullPredAlgorithmComparisonsPagePerTreatment.pdf", width = 10, height = 6)
compare_model_accuracies_byFactor(temp, accuracy_measure = "Rsquared", color_by = "Algorithm", 
                         data_group_col = "OutcomeVar", pageBy = "DataGroup",
                         compare_scope = "within",
                         annotate_comparisons = F)
dev.off()



#-------------------------------------------------------------------------------
# DATA GROUP COMPARISONS (NON-SUBSET)

temp = SubsetCompAccs[SubsetCompAccs$Subsetting == "None",]

pdf("../Rout/modelComparison/fullPredDataGroupComparisonsPagePerAlgorithm.pdf", width = 10, height = 6)
compare_model_accuracies_byFactor(temp, accuracy_measure = "Rsquared", color_by = "DataGroup", 
                                  data_group_col = "OutcomeVar", pageBy = "Algorithm",
                                  compare_scope = "within",
                                  annotate_comparisons = F)
dev.off()

#-------------------------------------------------------------------------------
# SUBSET BY TIME COMPARISON PLOTS: LASSO
temp = SubsetCompAccs[SubsetCompAccs$Subsetting %in% c("None", "1to51") & SubsetCompAccs$Algorithm == "LASSO",]
pdf("../Rout/modelComparison/TimeSubsetComparisonsPagePerTreatment.pdf", width = 10, height = 6)
timeSubsCompLasso = compare_model_accuracies_byFactor(temp,accuracy_measure = "Rsquared", data_group_col = "OutcomeVar",
                                                      pageBy = "DataGroup",
                                             annotate_comparisons = F,
                                             color_by = "Subsetting",
                                             compare_scope = "within",
                                             test_type = "wilcox",
                                             main_title = "LASSO -")



# SUBSET BY TIME COMPARISON PLOTS: RF
temp = SubsetCompAccs[SubsetCompAccs$Subsetting %in% c("None", "1to51") & SubsetCompAccs$Algorithm == "RF",]
timeSubsCompRF = compare_model_accuracies_byFactor(temp,accuracy_measure = "Rsquared", data_group_col = "OutcomeVar",
                                          pageBy = "DataGroup",
                                          annotate_comparisons = F,
                                          color_by = "Subsetting",
                                          compare_scope = "within",
                                          test_type = "wilcox",
                                          main_title = "RF -")
dev.off()


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FEATURE IMPORTANCE 


# ------------------------------------------------------------------------------
# RF

# load models
listElems = load("../Rout/mlModels/allRFHarvestPredictions.Rdata")
RFmodels = mget(listElems)
rm(list = listElems)

# get a list of outcome variables
outVars = c("Total.Biomass.DW_model","Total.Spike.weight_model","N5.spike.weight_model","Spike.number.per.pot_model")

if (!dir.exists("../Rout/modelComparison/FeatureImportance")) {
  dir.create("../Rout/modelComparison/FeatureImportance", recursive = TRUE)
}

pdf("../Rout/modelComparison/FeatureImportance/RF.pdf", width = 12, height = 10)

for(trainSet in names(RFmodels)){
  for(outVar in outVars){
    varName = sub("_model$", "", outVar)
    
    plot_rf_importance(RFmodels[[trainSet]][["models"]][[outVar]], top_k = 20,
                                       main_title = paste("RF",varName,"model feature importance,",trainSet),
                                       show_mean_r2 = T)
  }
}

dev.off()



# ------------------------------------------------------------------------------
# LASSO FEATURE IMPORTANCE

# load models
listElems = load("../Rout/mlModels/allLASSOHarvestPredictions.Rdata")
lassoModels = mget(listElems)
rm(list = listElems)

# get a list of outcome variables
outVars = c("Total.Biomass.DW_model","Total.Spike.weight_model","N5.spike.weight_model","Spike.number.per.pot_model")

pdf("../Rout/modelComparison/FeatureImportance/lasso.pdf", width = 12, height = 10)

for (trainSet in names(lassoModels)) {
  for (outVar in outVars) {
    varName = sub("_model$", "", outVar)
    plot_lasso_importance(
      lassoModels[[trainSet]][["models"]][[outVar]],
      minAbsCoef = 0.01,
      main_title = paste("LASSO", varName, "Coefficients,", trainSet)
    )
  }
}

dev.off()





# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# TESTING ACROSS TREATMENT GROUPS

# Let's test how well models do when used to predict harvest traits using data
# from other treatment groups. We neeed the original means and standard 
# deviations, so that we can undo the z-score transformation and reapply the
# fitting one from the model we are testing.

# We just do this for LASSO models here


### --- Load the necessary objects --- ###
# Load pooled predictors and harvest data
listElems = load("../Rout/Rdata/StandardizedByDATWidened.Rdata")
pooledData = mget(listElems)
rm(list = listElems)

# Load the summary statistics for the traits
load("../Rout/Rdata/traitMeans.Rdata")

# Create subsets of pooled predictors based on Treatment
pooledData$controlPreds = pooledData$AllPredictorVars[ pooledData$AllPredictorVars$Treatment == "Control", ]
pooledData$droughtPreds = pooledData$AllPredictorVars[ pooledData$AllPredictorVars$Treatment == "Drought", ]

# Initialize an empty results data frame
results_df = data.frame(ResponseTrait = character(),
                         DataGroup = character(),
                         Rsq_Control = numeric(),
                         Rsq_Drought = numeric(),
                         Rsq_Pooled = numeric(),
                         stringsAsFactors = FALSE)

# Loop over each treatment group that models were trained on:
# "pooled": trained on pooled predictors (no re-transformation needed)
# "control" and "drought": models trained on treatment-specific z-scores
for (group in c("pooled", "control", "drought")) {
  
  # Extract the model list for the current group (e.g., lassoModels$pooledOut$models)
  model_list = lassoModels[[paste0(group, "Out")]]$models
  
  # Iterate through each trait model for this group
  for (trait in names(model_list)) {
    # Remove the "_model" suffix from the trait name
    clean_trait = sub("_model$", "", trait)
    
    # Match the observed response values from harvest data.
    # The harvest data are not z-transformed and remain as originally measured.
    obs_control = pooledData$harvest[ match(rownames(pooledData$controlPreds), pooledData$harvest$Plant.ID),
                                       clean_trait, drop = TRUE ]
    obs_drought = pooledData$harvest[ match(rownames(pooledData$droughtPreds), pooledData$harvest$Plant.ID),
                                       clean_trait, drop = TRUE ]
    obs_pooled  = pooledData$harvest[ match(rownames(pooledData$AllPredictorVars), pooledData$harvest$Plant.ID),
                                       clean_trait, drop = TRUE ]
    
    # Get the current trait model
    model = model_list[[trait]]
    
    # For models trained on Control or Drought,
    # we need to transform the pooled predictors to the correct z-scale.
    if (group == "control") {
      # Re-standardize predictors so that they are on the Control scale
      control_transformed = restandardize(pooledData$controlPreds, summary_stats, target_treatment = "Control")
      pooled_transformed  = restandardize(pooledData$AllPredictorVars, summary_stats, target_treatment = "Control")
      drought_transformed = restandardize(pooledData$droughtPreds, summary_stats, target_treatment = "Control")
      
      Rsq_Control = caret::postResample(predict(model, control_transformed), obs_control)["Rsquared"]
      Rsq_Drought = caret::postResample(predict(model, drought_transformed), obs_drought)["Rsquared"]
      Rsq_Pooled  = caret::postResample(predict(model, pooled_transformed), obs_pooled)["Rsquared"]
      
    } else if (group == "drought") {
      # Re-standardize predictors on the Drought scale
      drought_transformed = restandardize(pooledData$droughtPreds, summary_stats, target_treatment = "Drought")
      pooled_transformed  = restandardize(pooledData$AllPredictorVars, summary_stats, target_treatment = "Drought")
      control_transformed = restandardize(pooledData$controlPreds, summary_stats, target_treatment = "Drought")
      
      Rsq_Control = caret::postResample(predict(model, control_transformed), obs_control)["Rsquared"]
      Rsq_Drought = caret::postResample(predict(model, drought_transformed), obs_drought)["Rsquared"]
      Rsq_Pooled  = caret::postResample(predict(model, pooled_transformed), obs_pooled)["Rsquared"]
      
    } else if (group == "pooled") {
      # For pooled models, the predictors are already on the correct (pooled) scale.
      Rsq_Control = caret::postResample(predict(model, pooledData$controlPreds), obs_control)["Rsquared"]
      Rsq_Drought = caret::postResample(predict(model, pooledData$droughtPreds), obs_drought)["Rsquared"]
      Rsq_Pooled  = caret::postResample(predict(model, pooledData$AllPredictorVars), obs_pooled)["Rsquared"]
    }
    
    # Append the results for this trait and group
    results_df = rbind(results_df, data.frame(
      ResponseTrait = clean_trait,
      DataGroup = group,
      Rsq_Control = Rsq_Control,
      Rsq_Drought = Rsq_Drought,
      Rsq_Pooled = Rsq_Pooled,
      stringsAsFactors = FALSE
    ))
  }
}


# Clean names
results_df$ResponseTrait = sapply(results_df$ResponseTrait, clean_label_HarvestTraits)

# MAKE LATEX TABLE
# Define traits of interest and row order
selected_traits = c("Spike No", "Biomass DW", "Total Spike W", "Spike W.5Spk")
group_order = c("control", "drought", "pooled")

# Prepare data
latex_df = results_df %>%
  filter(ResponseTrait %in% selected_traits) %>%
  mutate(ResponseTrait = factor(ResponseTrait, levels = selected_traits),
         DataGroup = factor(DataGroup, levels = group_order)) %>%
  arrange(ResponseTrait, DataGroup) %>%
  mutate(across(starts_with("Rsq_"), ~ round(., 3)))

# Initialize empty vector to store LaTeX row strings
rows = c()

# Loop through traits and build LaTeX-formatted rows
for (trait in selected_traits) {
  sub_df = latex_df %>% filter(ResponseTrait == trait)
  
  rows = c(
    rows,
    sprintf("\\multirow{3}{*}{%s} & %s & %.2f & %.2f & %.2f \\\\",
            trait,
            as.character(sub_df$DataGroup[1]),
            sub_df$Rsq_Control[1], sub_df$Rsq_Drought[1], sub_df$Rsq_Pooled[1]),
    sprintf(" & %s & %.2f & %.2f & %.2f \\\\",
            as.character(sub_df$DataGroup[2]),
            sub_df$Rsq_Control[2], sub_df$Rsq_Drought[2], sub_df$Rsq_Pooled[2]),
    sprintf(" & %s & %.2f & %.2f & %.2f \\\\",
            as.character(sub_df$DataGroup[3]),
            sub_df$Rsq_Control[3], sub_df$Rsq_Drought[3], sub_df$Rsq_Pooled[3]),
    "\\arrayrulecolor{gray}\\hline\\arrayrulecolor{black}"
  )
}

# Paste the full table output
cat("\\begin{table}[ht]\n",
    "\\centering\n",
    "\\begin{tabular}{llrrr}\n",
    "  \\hline\n",
    "ResponseTrait & DataGroup & $R^2_{\\text{Control}}$ & $R^2_{\\text{Drought}}$ & $R^2_{\\text{Pooled}}$ \\\\\n",
    "  \\hline\n",
    paste(rows, collapse = "\n"),
    "\n\\end{tabular}\n",
    "\\caption{Cross-validation $R^2$ values for LASSO models trained on control, drought, or pooled data.}\n",
    "\\label{tab:lassoRsqClean}\n",
    "\\end{table}\n",
    sep = "")

# save excel
excel_df = latex_df %>%
  rename(
    Trait = ResponseTrait,
    ModelGroup = DataGroup,
    R2_Control = Rsq_Control,
    R2_Drought = Rsq_Drought,
    R2_Pooled = Rsq_Pooled
  )

# Save to Excel
write.xlsx(excel_df, "../Rout/tables/LASSOCrossValidation.xlsx", rowNames = FALSE)

