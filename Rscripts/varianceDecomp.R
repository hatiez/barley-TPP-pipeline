source("varianceDecompFunctions.R")

# load Data
load("../Rout/Rdata/reimputed_outliers.Rdata")
# make list of data frames
DFs = names(sapply(ls(), function(x) is.data.frame(get(x)))[sapply(ls(), function(x) is.data.frame(get(x))) == TRUE])

if (!dir.exists("../Rout/")) {
  dir.create("../Rout/", recursive = TRUE)
}

if (!dir.exists("../Rout/VarianceComponents")) {
  dir.create("../Rout/VarianceComponents", recursive = TRUE)
}

#-------------------------------------------------------------------------------
# REMOVE MIXED GENOTYPES
plants_to_remove =  c("L7", "L8", "L9")

# Loop through each data frame in DFs
for (df_name in DFs) {
  # Check if the data frame exists in the global environment
  if (exists(df_name, envir = .GlobalEnv)) {
    # Get the data frame
    df =  get(df_name, envir = .GlobalEnv)
    
    # Check if 'Plant.Name' column exists in the data frame
    if ("Plant.Name" %in% names(df)) {
      # Remove rows where Plant.Name matches any in the list
      df =  df[!df$Plant.Name %in% plants_to_remove, ]
      
      # Rename "G8_494" to "G7_493"
      df$Plant.Name =  ifelse(df$Plant.Name == "G8_494", "G7_493", df$Plant.Name)
      
      # Assign the cleaned data frame back to the global environment
      assign(df_name, df, envir = .GlobalEnv)
    }
  }
}


#-------------------------------------------------------------------------------
## Variance Decomposition: Harvest

# Initialize a list to store the count of GR scores and the mean GR
countRepScores = list(
  below_0_5 = 0,
  between_0_5_and_0_75 = 0,
  between_0_75_and_0_9 = 0,
  above_0_9 = 0,
  total_sum_for_mean = 0,
  mean_GR = 0
)
minRep = c()
maxRep = c()

# Generate plot
tempRepTable = harvest_repeatability_table(harvest)
for (i in 1:nrow(tempRepTable)){
  if (tempRepTable$R_min[i]<0.5){
    countRepScores$below_0_5 = countRepScores$below_0_5 + 1
  }
  else if (tempRepTable$R_min[i]<0.75){
    countRepScores$between_0_5_and_0_75 = countRepScores$between_0_5_and_0_75 + 1
  }
  else if (tempRepTable$R_min[i]<0.90){
    countRepScores$between_0_75_and_0_9 = countRepScores$between_0_75_and_0_9 + 1
  }
  else{
    countRepScores$above_0_9 = countRepScores$above_0_9 + 1
  }
}
# Sum them for the mean
countRepScores$total_sum_for_mean = countRepScores$total_sum_for_mean + sum(as.matrix(tempRepTable[,5:6]))
minRep = c(minRep,min(tempRepTable$R_min))
maxRep = c(maxRep,max(tempRepTable$R_max))
# Divide mean by total
# Calculate the mean TR
countRepScores$mean_GR = (countRepScores$total_sum_for_mean / sum(countRepScores$below_0_5, countRepScores$between_0_5_and_0_75, countRepScores$between_0_75_and_0_9, countRepScores$above_0_9)) / 2
# Open PDF device
pdf("../Rout/VarianceComponents/harvest_repeatability_plots.pdf", width = 10, height = 7)
plot_harvest_repeatability(tempRepTable,
                           legendLabels = c("Genetic Line","Treatment","Residual"),
                           main = "Harvest Variance Components and Repeatability",
                           cex.axis = 1.2)

# Close PDF device
dev.off()

print(countRepScores)

# ------------------------------------------------------------------------------
## Temporal Repeatability: Testing Assumptions

# Preparing data

# Remove harvest data from DFs list
DFs = DFs[!(DFs %in% c("harvest","harvest_xray"))]

#### Delete nonsensical columns ####
# this needs to be done in weight and IR
weight = weight[ , -which(names(weight) %in% c("Weight","Weight.After.Watering","CumulativeWater"))]
IR = IR[ , -which(names(IR) %in% c("IR1_Plant_Temp","Probes_T4.unit"))]



# Testing assumptions: We want to use temporal repeatability to quantify how 
# stable temporal traits are within genotypes across time points. This requires
# that traits are procedurally homogenous, which we test for here.
# for more info check https://doi.org/10.1016/j.anbehav.2015.04.008

# Testing assumptions: Genotype-wise Slope

# Testing assumptions
pdf("../Rout/VarianceComponents/traitGenotypeSlopesWithTreatmentEffectUnscaled.pdf", width = 11, height = 8.5)

# Initialize an empty dataframe to store slopes
GenotypeSlopesTreatmentEff =  data.frame()

for (sensor in DFs) {
  if(!(sensor %in% c("harvest","harvest_xray"))){
    df =  get(sensor)
    df$GenotypeTreatment = paste(df$Plant.Name,df$Treatment,sep = "_")
    # remove data from first week
    df = df[df$DAT > (24+7),]
    new_slopes =  fit_random_regression(df, includeTreatment = T, idVar = "GenotypeTreatment")  # Compute slopes

    # Append new_slopes to slopes
    GenotypeSlopesTreatmentEff =  rbind(GenotypeSlopesTreatmentEff, new_slopes$results)
  }
}

dev.off()

write.csv(GenotypeSlopesTreatmentEff[which(GenotypeSlopesTreatmentEff$p_value_adj > 0.05),],file = "../Rout/tables/tempTraitsProcedurallyHomogenousByGenotypeUnscaled.csv")

# Conclusion: CHange in trait across time has different slopes between genotypes
# in many traits.
# We should not use temporal repeatability, and thus leave it out from our plots


# ------------------------------------------------------------------------------
## Temporal Repeatability (TR): Variance Decomposition

# Initialize a list to store the count of TR scores and the mean TR.
# Note that TR is invalid and was thus left out from the publication, see above.
countRepScores = list(
  below_0_5 = 0,
  between_0_5_and_0_75 = 0,
  between_0_75_and_0_9 = 0,
  above_0_9 = 0,
  total_sum_for_mean = 0,
  mean_TR = 0
)
minRep = c()
maxRep = c()

# To make variance decomp plots combining a custom selection of traits, we make
# a few lists of traits to plot together

# For a custom variance decomp plot showing a wide selection of traits
selected_vars =  c("IR1.Probes_Plant_deltaT",
                   "FC1_Plant_HL.ind_QY_L2",
                   "FC1_Plant_HL.LL_QY_Lss1",
                   "FC1_Plant_HL.LL_QY_Lss2",
                   "FC1_Plant_HL.LL_QY_Lss2.Lss1",
                   "RGB1.RGB2_Plant_Volume",
                   "RGB2_Plant_AREA_MM",
                   "RGB1_Plant_Avg_HEIGHT_MM")

# Initialize an empty data frame to store selected rows
selectedTempReps =  data.frame()

# For a plot showing top 20 traits for TPC
treatmentPred_vars =  c("IR1.Probes_Plant_deltaT",
                   "RGB2_Plant_AREA_MM",
                   "FC1_Plant_HL.LL_QY_Lss1",
                   "FC1_Plant_HL.LL_QY_Lss2.Lss1",
                   "RGB1.RGB2_Plant_Volume",
                   "RGB1_Plant_Avg_HEIGHT_MM",
                   "RGB2_Plant_PERIMETER_MM"
                   )

# Initialize an empty data frame to store selected rows
treatmentPredReps =  data.frame()

# Traits used in biomass prediction (drought treatment)
BiomassDroughtLassoPred_vars =  c("RGB2_Plant_AREA_MM",
                                  "RGB2_Plant_PERIMETER_MM",
                                  "RGB1.RGB2_Plant_Volume",
                                  "RGB2_Plant_COMPACTNESS",
                                  "SWIR_Plant_LWVI2",
                                  "FC1_Plant_HL.ind_FvFm_L1",
                                  "RGB1_Plant_RGB.68.85.63",
                                  "FC1_Plant_CL.ind_QY_max",
                                  "FC1_Plant_ChlCont_SFR_R",
                                  "FC1_Plant_ChlCont_SFR_G",
                                  "FC1_Plant_HL.ind_qL_Lss",
                                  "FC1_Plant_CL.ind_qL_L1")

# Initialize an empty data frame to store selected rows
BiomassDroughtLassoPredReps =  data.frame()

# Traits used in biomass prediction (drought treatment), early prediction
BiomassDroughtLassoPred1to51_vars =  c("RGB2_Plant_AREA_MM",
                                  "VNIR_Plant_MCARI1",
                                  "FC1_Plant_CL.ind_QY_max",
                                  "RGB1_Plant_Avg_WIDTH_MM",
                                  "FC1_Plant_HL.LL_QY_Lss2.Lss1",
                                  "RGB1_Plant_Avg_HEIGHT_MM",
                                  "RGB2_Plant_COMPACTNESS",
                                  "FC1_Plant_CL.ind_qL_L1",
                                  "FC1_Plant_HL.ind_qL_Lss",
                                  "SWIR_Plant_WATER")


# Initialize an empty data frame to store selected rows
BiomassDroughtLassoPredReps1to51 =  data.frame()

# POOLED BIOMASS PRED
# Traits used in biomass prediction (pooled treatment)
BiomassPooledLassoPred_vars =  c("RGB1_Plant_Avg_AREA_MM", "RGB2_Plant_AREA_MM", "FC1_Plant_HL.ind_qP_Lss", "RGB1_Plant_Avg_HEIGHT_MM", "RGB2_Plant_PERIMETER_MM", "SWIR_Plant_LWVI2", "VNIR_Plant_NDVI", "FC1_Plant_HL.LL_QY_Lss2.Lss1", "RGB2_Plant_COMPACTNESS", "FC1_Plant_HL.ind_QY_max", "FC1_Plant_HL.ind_qN_Lss", "IR1.Probes_Plant_deltaT", "FC1_Plant_CL.ind_qL_Lss")

# Initialize an empty data frame to store selected rows
BiomassPooledLassoPredReps =  data.frame()

# Traits used in biomass prediction (pooled treatment), early prediction
BiomassPooledLassoPred1to51_vars =  c("RGB1_Plant_Avg_AREA_MM", "RGB2_Plant_AREA_MM", "VNIR_Plant_MCARI1", "SWIR_Plant_LWVI2", "VNIR_Plant_NDVI", "RGB1_Plant_Avg_WIDTH_MM", "FC1_Plant_HL.LL_QY_Lss1", "VNIR_Plant_PRI", "FC1_Plant_CL.ind_qL_L1", "IR1.Probes_Plant_deltaT", "FC1_Plant_HL.ind_NPQ_Lss", "RGB2_Plant_COMPACTNESS", "SWIR_Plant_WATER", "FC1_Plant_HL.LL_QY_Lss2.Lss1", "RGB1.RGB2_Plant_Volume")

# Initialize an empty data frame to store selected rows
BiomassPooledLassoPredReps1to51 =  data.frame()

# Open PDF device
pdf("../Rout/VarianceComponents/temporal_varianceComponents.pdf", width = 10, height = 7)

for(df in DFs){
  print(paste("processing", df))
  current_df =  get(df)
  tempRepTable = temporal_repeatability_table(current_df)
  
  for (i in 1:nrow(tempRepTable)){
    if (tempRepTable$TR_min[i] < 0.5){
      countRepScores$below_0_5 = countRepScores$below_0_5 + 1
    } else if (tempRepTable$TR_min[i] < 0.75){
      countRepScores$between_0_5_and_0_75 = countRepScores$between_0_5_and_0_75 + 1
    } else if (tempRepTable$TR_min[i] < 0.90){
      countRepScores$between_0_75_and_0_9 = countRepScores$between_0_75_and_0_9 + 1
    } else {
      countRepScores$above_0_9 = countRepScores$above_0_9 + 1
    }
  }
  countRepScores$total_sum_for_mean = countRepScores$total_sum_for_mean + sum(as.matrix(tempRepTable[,6:7]))
  
  selectedTempReps =  rbind(selectedTempReps, tempRepTable[tempRepTable$Variable %in% selected_vars, ])
  treatmentPredReps =  rbind(treatmentPredReps, tempRepTable[tempRepTable$Variable %in% treatmentPred_vars, ])
  BiomassDroughtLassoPredReps =  rbind(BiomassDroughtLassoPredReps, tempRepTable[tempRepTable$Variable %in% BiomassDroughtLassoPred_vars, ])
  BiomassDroughtLassoPredReps1to51 =  rbind(BiomassDroughtLassoPredReps1to51, tempRepTable[tempRepTable$Variable %in% BiomassDroughtLassoPred1to51_vars, ])
  BiomassPooledLassoPredReps =  rbind(BiomassPooledLassoPredReps, tempRepTable[tempRepTable$Variable %in% BiomassPooledLassoPred_vars, ])
  BiomassPooledLassoPredReps1to51 =  rbind(BiomassPooledLassoPredReps1to51, tempRepTable[tempRepTable$Variable %in% BiomassPooledLassoPred1to51_vars, ])
  
  plot_temporal_repeatability(tempRepTable,
                              legendLabels = c("Genetic Line(DAT)","DAT","Treatment(DAT)","Residual"),
                              main = paste(df, "Variance Components and Repeatability"),
                              cex.axis = 0.9)
  minRep = c(minRep, min(tempRepTable$TR_min))
  maxRep = c(maxRep, max(tempRepTable$TR_max))
}
dev.off()

# Calculate the mean TR
countRepScores$mean_TR = (countRepScores$total_sum_for_mean / sum(countRepScores$below_0_5, countRepScores$between_0_5_and_0_75, countRepScores$between_0_75_and_0_9, countRepScores$above_0_9)) / 2

# Print the result
print(countRepScores)


# Print or store selectedTempReps for plotting later
plot_temporal_repeatability(selectedTempReps,
                            legendLabels = c("Genetic Line(DAT)","DAT","Treatment(DAT)","Residual"),
                            main = paste("Selected Variables Variance Components"),
                            cex.axis = 1.2)


# Open PDF device
pdf("../Rout/VarianceComponents/VarianceDecompTreatmentPred.pdf", width = 8, height = 7)

plot_temporal_repeatability(treatmentPredReps,
                            legendLabels = c("Genetic Line(DAT)","DAT","Treatment(DAT)","Residual"),
                            #main = paste("Treatment Prediction Variance Components"),
                            cex.axis = 1.2)

# Close PDF device
dev.off()

# Open PDF device
pdf("../Rout/VarianceComponents/VarianceDecompBiomassDroughtLassoPred.pdf", width = 8, height = 7)

plot_temporal_repeatability(BiomassDroughtLassoPredReps,
                            legendLabels = c("Genetic Line(DAT)","DAT","Treatment(DAT)","Residual"),
                            #main = paste("Treatment Prediction Variance Components"),
                            cex.axis = 1.2)

# Close PDF device
dev.off()

# Open PDF device
pdf("../Rout/VarianceComponents/VarianceDecompBiomassDroughtLassoPred1to51.pdf", width = 8, height = 7)

plot_temporal_repeatability(BiomassDroughtLassoPredReps1to51,
                            legendLabels = c("Genetic Line(DAT)","DAT","Treatment(DAT)","Residual"),
                            #main = paste("Treatment Prediction Variance Components"),
                            cex.axis = 1.2)

# Close PDF device
dev.off()

pdf("../Rout/VarianceComponents/VarianceDecompBiomassPooledLassoPred.pdf", width = 8, height = 7)

plot_temporal_repeatability(BiomassPooledLassoPredReps,
                            legendLabels = c("Genetic Line(DAT)", "DAT", "Treatment(DAT)", "Residual"),
                            cex.axis = 1.2)

dev.off()

pdf("../Rout/VarianceComponents/VarianceDecompBiomassPooledLasso1to51Pred.pdf", width = 8, height = 6)

plot_temporal_repeatability(BiomassPooledLassoPredReps1to51,
                            legendLabels = c("Genetic Line(DAT)", "DAT", "Treatment(DAT)", "Residual"),
                            cex.axis = 1.2)

dev.off()


#------------------------------------------------------------------------------
# COEFFICIENT OF VARIATION
# Here we get the coefficient of variation for all traits. This was not included
# in the publication.

# Identify numeric columns (excluding the treatment column)
numeric_vars =  harvest %>% dplyr::select(where(is.numeric)) %>% names()

#Compute the coefficient of variation for each numeric column
cv_table =  map_dfr(numeric_vars, function(var) { control_cv =  harvest %>% filter(Treatment == "Control") %>% summarise(cv = sd(.data[[var]], na.rm = TRUE) / mean(.data[[var]], na.rm = TRUE)) %>% pull(cv)

drought_cv =  harvest %>% filter(Treatment == "Drought") %>% summarise(cv = sd(.data[[var]], na.rm = TRUE) / mean(.data[[var]], na.rm = TRUE)) %>% pull(cv)

pooled_cv =  harvest %>% summarise(cv = sd(.data[[var]], na.rm = TRUE) / mean(.data[[var]], na.rm = TRUE)) %>% pull(cv)

tibble(variable = var, control = control_cv, drought = drought_cv, pooled = pooled_cv) })

# improve row name
cv_table$variable = sapply(cv_table$variable, clean_label_TemporalTraits)
