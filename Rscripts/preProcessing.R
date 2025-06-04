source("preprocessingFunctions.R")

if (!dir.exists("../Rout/")) {
  dir.create("../Rout/", recursive = TRUE)
}

################################################################
#------------------------- LOAD DATA --------------------------#
################################################################

### RGB ###

RGB1 = read.csv("../PublData/RGB1.csv")
RGB2 = read.csv("../PublData/RGB2.csv")
RGB1_colorSegmentation = read.csv("../PublData/RGB1_colorSegmentation.csv")
RGB1RGB2 = read.csv("../PublData/RGB1RGB2.csv")

### HYPERSPECTRAL ###
SWIR = read.csv("../PublData/SWIR.csv")
VNIR = read.csv("../PublData/VNIR.csv")

### THERMAL IR ###
IR = read.csv("../PublData/IR.csv")


### CHL FLUORESCENCE ###
FcHL_LL = read.csv("../PublData/FcHL_LL.csv")
FcHLind = read.csv("../PublData/FcHLind.csv")
FcCLind = read.csv("../PublData/FcCLind.csv")
FcChlCont = read.csv("../PublData/FcChlCont.csv")


### WATERING AND WEIGHING ###
weight = read.csv("../PublData/weight.csv")

### HARVEST DATA ###
harvest = read.csv("../PublData/harvest.csv")
harvest_xray = read.csv("../PublData/harvest_xray.csv")



###################################################################
#--------------------------- IMPUTING ----------------------------#
###################################################################

# Get the names of all currently open data frames
DFs = names(sapply(ls(), function(x) is.data.frame(get(x)))[sapply(ls(), function(x) is.data.frame(get(x))) == TRUE])

# Checking for missingness
for(df_name in DFs){
  df = get(df_name)
  print(paste("Missing rows in ",df_name,": ",sum(apply(df, 1, function(row) any(is.na(row))))," out of ",nrow(df),sep = ""))
}
remove(df)


# Adressing missing values in FcHL_LL
rows_with_na = which(apply(FcHL_LL, 1, function(row) any(is.na(row))))

# All missing values in FcHL_LL are from the same day:
print(paste("Missing values in FcHL_LL present in DAT:",unique(FcHL_LL[rows_with_na,"DAT"])))

# Therefore, it makes most sense to just entirely remove this day from the data set
FcHL_LL = FcHL_LL[FcHL_LL$DAT != "54",]


# HARVEST DATA
# One sample was removed from the experiment but still has an entry in the 
# harvest data frame. Let's remove it.
print(harvest[harvest$Plant.ID == "L6_C_17",])
harvest = harvest[harvest$Plant.ID != "L6_C_17",]

# There are a few missing values in the Total.Spike.weight column
# We can impute these using missForest. They will be imputed groupwise (per genotype)
imputation = groupWiseMissForest(harvest,groupVars = c("Plant.Name","Plant.Info"), ignoreVars = c("Plant.ID"), ntree = 20)
harvest = imputation$ImputedData

# Print the out-of-bag error of imputations
print("Out-of-bag error of imputations for Total.Spike.Weight in harvest data:")
print(imputation$OOBError)

# Ensure the 'tables' directory exists
if (!dir.exists("../Rout/tables")) {
  dir.create("../Rout/tables", recursive = TRUE)
}

# Save the CSV file
write.csv(imputation$OOBError, "../Rout/tables/missingValOOB.csv")


###################################################################
#------------------------- REORGANIZING --------------------------#
###################################################################

### RGB ###
RGB = merge_and_filter_dfs(RGB1=RGB1,RGB2=RGB2,RGB1RGB2=RGB1RGB2)
remove(RGB1,RGB2,RGB1RGB2)

RGB1_colorSegmentation = merge_and_filter_dfs(RGB1_colorSegmentation=RGB1_colorSegmentation)


### HYPERSPECTRAL ###
SWIR = merge_and_filter_dfs(SWIR = SWIR)
VNIR = merge_and_filter_dfs(VNIR = VNIR)


### THERMAL IR ###
IR = merge_and_filter_dfs(IR=IR)

### CHL FLUORESCENCE ###
#simplifying datasets
FcHL_LL = merge_and_filter_dfs(FcHL_LL=FcHL_LL)
# set last column to numeric, as it seems to have been interpreted as non-numeric
FcHL_LL$FC1_Plant_HL.LL_QY_Lss2.Lss1 = as.numeric(FcHL_LL$FC1_Plant_HL.LL_QY_Lss2.Lss1)

FcHLind = merge_and_filter_dfs(FcHLind=FcHLind)
FcCLind = merge_and_filter_dfs(FcCLind=FcCLind)
FcChlCont = merge_and_filter_dfs(FcChlCont=FcChlCont)


### WATERING AND WEIGHING ###
weight = merge_and_filter_dfs(weight = weight)

### HARVEST DATA ###
colnames(harvest)[colnames(harvest)=="X5.spike.weight"] = "N5.spike.weight"
colnames(harvest)[colnames(harvest)=="X20.seed.weight"] = "N20.seed.weight"
colnames(harvest)[colnames(harvest)=="Plant.Info"] = "Treatment"

colnames(harvest_xray)[colnames(harvest_xray)=="X..spikelets"] = "N.Spikelets"
colnames(harvest_xray)[colnames(harvest_xray)=="X..Seeds.per.spike"] = "N.Seeds.per.spike"
colnames(harvest_xray)[colnames(harvest_xray)=="Seeds.spikes"] = "Seeds.per.spike"
colnames(harvest_xray)[colnames(harvest_xray)=="Plant.Info"] = "Treatment"


# Get the names of all currently open data frames
DFs = names(sapply(ls(), function(x) is.data.frame(get(x)))[sapply(ls(), function(x) is.data.frame(get(x))) == TRUE])

# Correcting mislabelled treatment columns
for (df in DFs) {
  current_df = get(df)
  print(paste(df,":",sep = ""))
  assign(df, plantIDCorrection(current_df))
  print("")
}


###################################################################
#--------------------------- OUTLIERS ----------------------------#
###################################################################
#Detection method is either IQR or Zscore 
chosenMethod = "IQR"

# Depending on the method chosen, you need to pick a threshold for what is considered an outlier (higher threshold -> less values considered outliers)
zScoreThreshold = 3
IQRmultipler = 3

# If removeOutliers is false, the next chunk will detect and count outliers, but not replace them with NAs
removeOutliers = TRUE

# Groupvars decide at what level we group the data before checking for outliers
groupVarsTemporal = c("DAT", "Plant.Name", "Treatment")
groupVarsHarvest = c("Plant.Name", "Treatment")


# Store outlier rates in a results table
outlier_rate_table = data.frame(
  DataFrame = character(),
  OutlierRate = numeric(),
  stringsAsFactors = FALSE
)

# Initialize list to store outlier count tables for each plant
plant_outlier_count = list()

# Initialize data frame to store outlier count for traits
trait_outlier_count = data.frame(
  DataFrame = character(),
  Trait = character(),
  OutlierCount = integer(),
  stringsAsFactors = FALSE
)

# Initialize data frame to store outlier count per trait per DAT
trait_outlier_count_perDAT = data.frame(
  DataFrame = character(),
  DAT = character(),
  Trait = character(),
  OutlierCount = integer(),
  stringsAsFactors = FALSE
)

# Initialize data frame to store outlier count per trait per Plant.ID
trait_outlier_count_perReplicate = data.frame(
  DataFrame = character(),
  Plant.ID = character(),
  Trait = character(),
  OutlierCount = integer(),
  stringsAsFactors = FALSE
)

# Initialize data frame to store outlier count per replicate per day
outlier_count_perReplicate_perDay = data.frame(
  DataFrame = character(),
  DAT = character(),
  Plant.ID = character(),
  OutlierCount = integer(),
  stringsAsFactors = FALSE
)

for (df_name in DFs) {
  df = get(df_name)
  
  if (df_name %in% c("harvest","harvest_xray")) {
    groupVars = groupVarsHarvest
  } else {
    groupVars = groupVarsTemporal
  }
  
  # Remove outliers based on the chosen method
  if (chosenMethod == "Zscore") {
    df_outliersRemoved = removeZscoreOutliers(df, threshold = zScoreThreshold, groupVars = groupVars)
  } else if (chosenMethod == "IQR") {
    df_outliersRemoved = removeIQROutliers(df, multiplier = IQRmultipler, groupVars = groupVars)
  }
  
  # Calculate the NA (outlier) rates
  removedValueRate = numericNArate(df_outliersRemoved)
  
  # Store the rates in the table
  outlier_rate_table = rbind(outlier_rate_table, 
                              data.frame(DataFrame = df_name, 
                                         OutlierRate = removedValueRate$NArate))
  
  # Store the cases of outliers for each plant
  outliersPerPlant = numericNAratePerReplicate(df_outliersRemoved)$na_count_by_variable_df
  outliersPerPlant = outliersPerPlant[outliersPerPlant$nNAs != 0,]
  
  # Append the data frame to the list with the data frame name as the key
  plant_outlier_count[[df_name]] = outliersPerPlant
  
  # Store the cases of outliers for each variable
  outliersPerVariable = numericNAratePerVar(df_outliersRemoved)
  outliersPerVariable = data.frame(DataFrame = df_name, outliersPerVariable)
  
  # Stack the data frames for all variables
  trait_outlier_count = rbind(trait_outlier_count, outliersPerVariable)
  
  # Calculate and store the cases of outliers for each variable per DAT
  if ("DAT" %in% names(df_outliersRemoved)) {
    outliersPerVarPerDAT = numericNAratePerVarPerGroup(df_outliersRemoved, groupVar = "DAT")
    outliersPerVarPerDAT = data.frame(DataFrame = df_name, outliersPerVarPerDAT)
    trait_outlier_count_perDAT = rbind(trait_outlier_count_perDAT, outliersPerVarPerDAT)
  }
  
  # Calculate and store the cases of outliers for each variable per Plant.ID
  if ("Plant.ID" %in% names(df_outliersRemoved)) {
    outliersPerVarPerReplicate = numericNAratePerVarPerGroup(df_outliersRemoved, groupVar = "Plant.ID")
    outliersPerVarPerReplicate = data.frame(DataFrame = df_name, outliersPerVarPerReplicate)
    trait_outlier_count_perReplicate = rbind(trait_outlier_count_perReplicate, outliersPerVarPerReplicate)
  }
  
  # Calculate and store the cases of outliers per replicate per day (DAT)
  if (all(c("DAT", "Plant.ID") %in% names(df_outliersRemoved))) {
    outliersPerReplicatePerDay = numericNAratePerVarPerGroup(df_outliersRemoved, groupVar = c("DAT", "Plant.ID"))
    outliersPerReplicatePerDay = data.frame(DataFrame = df_name, outliersPerReplicatePerDay)
    outlier_count_perReplicate_perDay = rbind(outlier_count_perReplicate_perDay, outliersPerReplicatePerDay)
  }
  
  # Optionally, assign the modified data frame back to the environment
  if (removeOutliers) {
    # Except for weight data, this one has too instable outlier detection
    if(df_name != "weight"){
      assign(df_name, df_outliersRemoved)
    }
  }
}


# Remove 'weight' from plant_outlier_count if it exists
plant_outlier_count$weight = NULL

# Combine all outlier counts into a single data frame
summed_outlier_df = do.call(rbind, plant_outlier_count)

# Ensure the 'outlier_counts' directory exists
if (!dir.exists("../Rout/tables/outlier_counts")) {
  dir.create("../Rout/tables/outlier_counts", recursive = TRUE)
}

# Save the CSV file
write.csv(summed_outlier_df, "../Rout/tables/outlier_counts/outliersPerPlant.csv")

# Rename columns
colnames(trait_outlier_count)[colnames(trait_outlier_count) == "nNAs"] = "nOutliers"
colnames(trait_outlier_count)[colnames(trait_outlier_count) == "NArate"] = "outlierRate"

# Write the trait outlier counts to CSV
write.csv(trait_outlier_count, "../Rout/tables/outlier_counts/outliersPerTrait.csv")

# Remove cases where no outlier was found
trait_outlier_count_perDAT = trait_outlier_count_perDAT[trait_outlier_count_perDAT$nNAs > 0,]

# Rename columns for trait_outlier_count_perDAT
colnames(trait_outlier_count_perDAT)[colnames(trait_outlier_count_perDAT) == "nNAs"] = "nOutliers"

# Write the trait outlier counts per DAT to CSV
write.csv(trait_outlier_count_perDAT, "../Rout/tables/outlier_counts/outliersPerTraitPerDAT.csv")

# Per trait outlier tables
# remove 0 values
trait_outlier_count_perReplicate = trait_outlier_count_perReplicate[trait_outlier_count_perReplicate$nNAs != 0,]
outlier_count_perReplicate_perDay = outlier_count_perReplicate_perDay[outlier_count_perReplicate_perDay$nNAs != 0,]

# Save or print the new tables
write.csv(trait_outlier_count_perReplicate, "../Rout/tables/outlier_counts/outliersPerReplicatePerTrait.csv")
write.csv(outlier_count_perReplicate_perDay, "../Rout/tables/outlier_counts/outliersPerReplicatePerDAT.csv")


# REIMPUTING VALUES
# Hyperparameter: number of trees used by missForest. 
# Set to 1 for test run
nTree = 1
# Set to higher for actual run
# nTree = 20

reimputationOOB = list()
for (df_name in DFs) {
  df = as.data.frame(get(df_name))
  if (df_name %in% c("harvest","harvest_xray")){
    imputation = groupWiseMissForest(df,groupVars = c("Plant.Name","Treatment"), ignoreVars = c("Plant.ID"),ntree = nTree)
    df = imputation$ImputedData
  }
  # OTHER CASES
  else {
    imputation = groupWiseMissForest(df,groupVars = c("Plant.Name","Treatment","DAT"), ignoreVars = c("Plant.ID"),ntree = nTree)
    df = imputation$ImputedData
  }
  # print mean OOB
  print(paste(df_name,"oob:",mean(imputation$OOBError$OOBError,na.rm = TRUE),"+/-",sd(imputation$OOBError$OOBError,na.rm = TRUE)))
  reimputationOOB[[df_name]] = imputation$OOBError
  # Assign the modified data frame back to the original name
  assign(df_name, df)
}

# Storing the reimputation OOBs
# Disable scientific notation
options(scipen = 999)

# Ensure the 'reimputationOOBs' directory exists
if (!dir.exists("../Rout/tables/reimputationOOBs")) {
  dir.create("../Rout/tables/reimputationOOBs", recursive = TRUE)
}

# Loop through each data frame in the reimputationOOB list
for (df_name in names(reimputationOOB)) {
  # Get the data frame
  df = reimputationOOB[[df_name]]
  
  # Construct the file path
  file_path = file.path("../Rout/tables/reimputationOOBs/", paste0(df_name, ".csv"))
  
  # Write the data frame to a CSV file
  write.csv(df, file = file_path, row.names = FALSE)
}

# Making a less detailed overbiew of the reimputation OOBs
## 1.  summarise the OOB errors 
oob_summary = lapply(names(reimputationOOB), function(nm) {
  data.frame(
    DataFrame = nm,
    OOB_mean  = mean(reimputationOOB[[nm]]$OOBError, na.rm = TRUE),
    OOB_sd    =  sd  (reimputationOOB[[nm]]$OOBError, na.rm = TRUE)
  )
}) |>
  bind_rows()

## 2.  combine with the previously-built outlier-rate table 
df_quality_overview = outlier_rate_table |>
  left_join(oob_summary, by = "DataFrame") |>
  arrange(DataFrame)          # optional, keeps rows in a predictable order

write.csv(df_quality_overview, "../Rout/tables/outliersOverview.csv")

# STORING: store data after reimputation
# Ensure the 'Rdata' directory exists 
if (!dir.exists("../Rout/Rdata")) {
  dir.create("../Rout/Rdata", recursive = TRUE)
}
save(list = DFs, file = "../Rout/Rdata/reimputed_outliers.Rdata")


###################################################################
#--------------------- STORING MEANS AND SDs ---------------------#
###################################################################
# We need these values to reverse z-transformations later to test prediction models

# Define harvest and time-series predictor data frames separately
harvest_DFs = c("harvest", "harvest_xray")
predictor_DFs = setdiff(DFs, harvest_DFs)

ignoreVars = c("Plant.Name", "Plant.ID")

# Final output
summary_stats = data.frame(WideName = character(),
                            Treatment = character(),
                            DAT = character(),
                            Trait = character(),
                            Mean = numeric(),
                            SD = numeric(),
                            stringsAsFactors = FALSE)

### --- Process predictor data frames (with DAT) --- ###
for (df_name in predictor_DFs) {
  df = get(df_name)
  
  # 1. Control/Drought treatment-specific stats
  treatment_stats = compute_stats(df, group_vars = c("Treatment", "DAT"))
  
  # 2. Pooled stats: ignore Treatment, group only by DAT
  pooled_stats = compute_stats(df, group_vars = "DAT", treatment_label = "Pooled")
  
  # Bind both
  summary_stats = bind_rows(summary_stats, treatment_stats, pooled_stats)
}

### --- Process harvest data frames --- ###
for (df_name in harvest_DFs) {
  df = get(df_name)
  
  # 1. Control/Drought stats (group by Treatment only)
  treatment_stats = compute_stats(df, group_vars = "Treatment", include_dat = FALSE)
  
  # 2. Pooled stats: no grouping
  pooled_stats = compute_stats(df, group_vars = NULL, treatment_label = "Pooled", include_dat = FALSE)
  
  # Bind both
  summary_stats = bind_rows(summary_stats, treatment_stats, pooled_stats)
}

save(summary_stats, file = "../Rout/Rdata/traitMeans.Rdata")



###################################################################
#----------------------- BOXCOX TRANSFORM ------------------------#
###################################################################
# First, we detect which groups have a non-normal distribution
# make a data frame to count non-normality per data frame processed.
non_normality_counts = data.frame(DataFrame = character(), Count = numeric(), stringsAsFactors = FALSE)

# Initialize the list to store data frames with non-normality information (non-normality per group)
non_normality_list = list()

p.adjustment = "bonferroni"

for (df_name in DFs) {
  df = get(df_name)
  if (df_name %in% c("harvest","harvest_xray")) {
    groupVars = c("Treatment")
  } else {
    groupVars = c("Treatment","DAT")
  }
  ignoreVars = c("DAT", "Treatment", "Plant.ID", "Plant.Name","Weight","Weight.After.Watering","Plant_Temp","Env_Temp")
  
  shapiroTest = groupWiseNormalityTest(df, groupVars = groupVars, p.adjust.method = p.adjustment, ignoreVars = ignoreVars,test = "shapiro")
  
  # Check if any numerical column in adjusted_table has a value below 0.05
  if (any(shapiroTest$adjusted_table %>% dplyr::select(-all_of(groupVars)) < 0.05, na.rm = TRUE)) {
    rows_below_0_05 = shapiroTest$adjusted_table %>%
      pivot_longer(cols = -all_of(groupVars), names_to = "variable", values_to = "p.value") %>%
      filter(p.value < 0.05) %>%
      dplyr::select(all_of(groupVars), variable, p.value)
    
    # Store the result in the list
    non_normality_list[[df_name]] = rows_below_0_05
  }
  else{
     # State if all groups were normal
    print(paste(df_name,"all gaussian!"))
    # add empty df to non_normality_list
    non_normality_list[[df_name]] = data.frame()
  }
}

# Print the non-normality counts table
non_normality_counts = data.frame(
  DataFrame = names(non_normality_list),
  Count = sapply(non_normality_list, nrow)
)

print(non_normality_counts)

# TRANSFORM: Use Box-Cox to transform non-normal groups
for (df_name in DFs) {
  df = get(df_name)
  if (df_name %in% c("harvest","harvest_xray")){
    GroupVars = c("Treatment")
  }
  else if (df_name == "weight"){
    # skip this round
    next
  }
  else{
    GroupVars = c("Treatment","DAT")
  }
  # check if the data frame has any non-gaussian variables
  if (nrow(non_normality_list[[df_name]]) > 0){
    redistributed_df = boxcoxRedistribute(df, non_normality_list[[df_name]], GroupVars = GroupVars)
    assign(df_name, redistributed_df)
  }
}


# Check non-normality for DFs post transformation
redistributed_non_normality_counts = data.frame(DataFrame = character(), Count = numeric(), stringsAsFactors = FALSE)

# Initialize the list to store data frames with non-normality information
redistributed_non_normality_list = list()

p.adjustment = "bonferroni"

for (df_name in DFs) {
  if (df_name %in% c("harvest","harvest_xray")) {
    groupVars = c("Treatment")
  } else {
    groupVars = c("Treatment","DAT")
  }
  df = get(df_name)
  ignoreVars = c("DAT", "Treatment", "Plant.ID", "Plant.Name","Weight","Weight.After.Watering","Plant_Temp","Env_Temp")
  
  shapiroTest = groupWiseNormalityTest(df, groupVars = groupVars, p.adjust.method = p.adjustment, ignoreVars = ignoreVars,test = "shapiro")
  
  # Check if any numerical column in adjusted_table has a value below 0.05
  if (any(shapiroTest$adjusted_table %>% dplyr::select(-all_of(groupVars)) < 0.05, na.rm = TRUE)) {
    rows_below_0_05 = shapiroTest$adjusted_table %>%
      pivot_longer(cols = -all_of(groupVars), names_to = "variable", values_to = "p.value") %>%
      filter(p.value < 0.05) %>%
      dplyr::select(all_of(groupVars), variable, p.value)
    
    # Store the result in the list
    redistributed_non_normality_list[[paste(df_name,"_redistributed",sep = "")]] = rows_below_0_05
  }
  else{
    # State if all groups were normal
    print(paste(df_name,"all gaussian!"))
    # add empty df to non_normality_list
    redistributed_non_normality_list[[paste(df_name,"_redistributed",sep = "")]] = data.frame()
  }
}

# Print the non-normality counts table
redistributed_non_normality_counts = data.frame(
  DataFrame = names(redistributed_non_normality_list),
  Count = sapply(redistributed_non_normality_list, nrow)
)

print(redistributed_non_normality_counts)


# Ensure the 'non_normality_counts' directory exists
if (!dir.exists("../Rout/tables/non_normality_counts")) {
  dir.create("../Rout/tables/non_normality_counts", recursive = TRUE)
}
write.csv(redistributed_non_normality_counts, "../Rout/tables/non_normality_counts/DF_wise_redistributed_NNcounts.csv")
write.csv(non_normality_counts, "../Rout/tables/non_normality_counts/DF_wise_NNcounts.csv")


# Ensure the 'groupWise' directory exists
if (!dir.exists("../Rout/tables/non_normality_counts/groupWise")) {
  dir.create("../Rout/tables/non_normality_counts/groupWise", recursive = TRUE)
}
for (df_name in names(non_normality_list)) {
  if(nrow(non_normality_list[[df_name]]) > 0){
    write.csv(non_normality_list[[df_name]], paste("../Rout/tables/non_normality_counts/groupWise/",df_name,".csv",sep = ""))
  }
}

for (df_name in names(redistributed_non_normality_list)) {
  if(nrow(redistributed_non_normality_list[[df_name]]) > 0){
    write.csv(redistributed_non_normality_list[[df_name]], paste("../Rout/tables/non_normality_counts/groupWise/",df_name,".csv",sep = ""))
  }
}


###################################################################
#----------------------- STANDARDISATION -------------------------#
###################################################################
# BoxCox causes a shift in values, so standardisation is required
# Also, it is needed for many of the methods used. We use zscore.
transformMethod = "zscore"

for (df_name in DFs){
  df = get(df_name)
  if (df_name %in% c("harvest","harvest_xray")) {
    groupVars = c("Treatment")
  } else {
    groupVars = c("Treatment","DAT")
  }
  if(transformMethod == "zscore"){
    assign(df_name,zScore(df,groupVars = groupVars,ignoreVars = c("Plant.Name","Plant.ID","DAT")))
  }
  if(transformMethod == "minmax"){
    assign(df_name,minMaxNormalize(df,groupVars = groupVars))
  }
}

# Weight data: On the first day, some of the columns contain only 0s. This causes issues with transformation, introduction NaNs. Let's put back 0s.
weight[weight$DAT == 25,c("WaterConsumption","CumulativeWater")] = 0

if(transformMethod == "zscore"){
  save(list = DFs, file = "../Rout/Rdata/normalized_standardized.Rdata")
}
if(transformMethod == "minmax"){
  save(list = DFs, file = "../Rout/Rdata/normalized_minMaxTrasnformed.Rdata")
}

###################################################################
#--------------------- REMOVING TRAITS ---------------------------#
###################################################################
# Some traits were used only to calculate others and are not needed
# E.g. environmental temperature for canopy temperature depression
traits_to_remove = read_csv("../PublData/Tier3_Traitstoremove.csv", col_names = FALSE)$X1

# Initialize a counter for the total number of removed variables
total_removed_count = 0

# Loop over each data frame and remove the specified columns
for (df_name in DFs) {
  df = get(df_name)
  
  # Identify columns to be removed
  cols_to_remove = names(df) %in% traits_to_remove
  
  # Count the number of columns to be removed
  removed_count = sum(cols_to_remove)
  
  # Update the total removed count
  total_removed_count = total_removed_count + removed_count
  
  # Remove the columns
  df = df[, !cols_to_remove]
  
  # Assign the modified data frame back to its original name
  assign(df_name, df)
  
  # Print the number of removed columns for the current data frame
  print(paste("Removed", removed_count, "variables from", df_name))
}

# Print the total number of removed variables
print(paste("Total number of removed variables:", total_removed_count))


if(transformMethod == "zscore"){
  save(list = DFs, file = "../Rout/Rdata/normalized_standardized_reducedVar.Rdata")
}
if(transformMethod == "minmax"){
  save(list = DFs, file = "../Rout/Rdata/normalized_minmaxTransformed_reducedVar.Rdata")
}

###################################################################
#------------------- WIDENING AND MERGING ------------------------#
###################################################################
# To get a data shape appropriate for machine learning methods

# Get a list of all the predictor data frames
PredictorDFs = setdiff(DFs,c("harvest","weight","harvest_xray"))

# We'll split the data into two lists
controlDFs = list()
droughtDFs = list()

for (df_name in PredictorDFs){
  df = get(df_name)
  controlDFs[[df_name]] = df[df$Treatment == "Control",]
  droughtDFs[[df_name]] = df[df$Treatment == "Drought",]
}


# the function widen and merge applies the reshaping operation to all data frames in the list and then fuses them into one
AllControlVars = widen_and_merge(controlDFs)
AllDroughtVars = widen_and_merge(droughtDFs)

# We'll remove the plant.ID column and use these as row names instead, just leaving the predictor variables 
rownames(AllControlVars) = AllControlVars$Plant.ID
AllControlVars$Plant.ID = NULL

rownames(AllDroughtVars) = AllDroughtVars$Plant.ID
AllDroughtVars$Plant.ID = NULL

# This operation introduced some NAs, due to missing values present as missing rows 
# This is likely due accidentally skipping plants
# CONTROL: Find indices of NA values
nans = which(is.na(AllControlVars), arr.ind = TRUE)

# Combine the results into a data frame
na_table = data.frame(Row = rownames(AllControlVars)[nans[, 1]], Column = colnames(AllControlVars)[nans[, 2]])

# Ensure the 'skipped_measurements' directory exists
if (!dir.exists("../Rout/tables/skipped_measurements")) {
  dir.create("../Rout/tables/skipped_measurements", recursive = TRUE)
}

# Store the table
write.csv(na_table, "../Rout/tables/skipped_measurements/skippedMeasurementsControl.csv")

# View the result
print(na_table)

# REPEAT WITH DROUGH: Find indices of NA values
nans = which(is.na(AllDroughtVars), arr.ind = TRUE)

# Combine the results into a data frame
na_table = data.frame(Row = rownames(AllDroughtVars)[nans[, 1]], Column = colnames(AllDroughtVars)[nans[, 2]])

# Store the table
write.csv(na_table, "../Rout/tables/skipped_measurements/skippedMeasurementsDrought.csv")

# View the result
print(na_table)

# Re-impute
nTree = 20
# I'll add back a genotype column for grouping first
AllControlVars$Genotype = as.factor(sub("_C_[0-9]+", "", rownames(AllControlVars)))
AllDroughtVars$Genotype = as.factor(sub("_D_[0-9]+", "", rownames(AllDroughtVars)))

# Do imputation for control data
impute = groupWiseMissForest(AllControlVars,groupVars = c("Genotype"), ignoreVars = c(), ntree = nTree)
AllControlVars = impute$ImputedData
print("Control imputation OOB error:")
print(impute$OOBError)

impute = groupWiseMissForest(AllDroughtVars,groupVars = c("Genotype"), ignoreVars = c(), ntree = nTree)
AllDroughtVars = impute$ImputedData
print("Drought imputation OOB error:")
print(impute$OOBError)

#remove genotype column
AllControlVars$Genotype = NULL
AllDroughtVars$Genotype = NULL

# convert to data.frame (was tibble before due to some functions from tidyr)
AllControlVars = as.data.frame(AllControlVars)
AllDroughtVars = as.data.frame(AllDroughtVars)
harvest = as.data.frame(harvest)

if(transformMethod == "zscore"){
  save(list = c("AllControlVars","AllDroughtVars","harvest","harvest_xray"), file = "../Rout/Rdata/standardized_normalized_widened.Rdata")
}
if(transformMethod == "minmax"){
  save(list = c("AllControlVars","AllDroughtVars","harvest","harvest_xray"), file = "../Rout/Rdata/minmaxTransformed_normalized_widened.Rdata")
}


###################################################################
#-------------- POOLED TREATMENT PROCESSING ----------------------#
###################################################################
# For TPC and Pooled treatment TPP, Box-Cox transformation and standardisation
# needs to be done with pooled treatments. We repeat all steps from here.

# Load data from pre-processing up to these steps
load("../Rout/Rdata/reimputed_outliers.Rdata")

# Remove redundant traits again
traits_to_remove = read_csv("../PublData/Tier3_Traitstoremove.csv", col_names = FALSE)$X1

# Initialize a counter for the total number of removed variables
total_removed_count = 0

# Loop over each data frame and remove the specified columns
for (df_name in DFs) {
  df = get(df_name)
  
  # Identify columns to be removed
  cols_to_remove = names(df) %in% traits_to_remove
  
  # Count the number of columns to be removed
  removed_count = sum(cols_to_remove)
  
  # Update the total removed count
  total_removed_count = total_removed_count + removed_count
  
  # Remove the columns
  df = df[, !cols_to_remove]
  
  # Assign the modified data frame back to its original name
  assign(df_name, df)
  
  # Print the number of removed columns for the current data frame
  print(paste("Removed", removed_count, "variables from", df_name))
}

# Print the total number of removed variables
print(paste("Total number of removed variables:", total_removed_count))


###################################################################
#------------------- POOLED NON-NORMALITY ------------------------#
###################################################################
non_normality_counts = data.frame(DataFrame = character(), Count = numeric(), stringsAsFactors = FALSE)

# Initialize the list to store data frames with non-normality information
non_normality_list = list()

p.adjustment = "bonferroni"

for (df_name in DFs) {
  df = get(df_name)
  if(df_name %in% c("harvest","harvest_xray")){
    groupVars = NULL
  }
  else{
    groupVars = c("DAT")
  }
  ignoreVars = c("DAT", "Treatment", "Plant.ID", "Plant.Name","Weight","Weight.After.Watering","Plant_Temp","Env_Temp")
  
  shapiroTest = groupWiseNormalityTest(df, groupVars = groupVars, p.adjust.method = p.adjustment, ignoreVars = ignoreVars,test = "shapiro")
  
  # Check if any numerical column in adjusted_table has a value below 0.05
  if (any(shapiroTest$adjusted_table %>% dplyr::select(-all_of(groupVars)) < 0.05, na.rm = TRUE)) {
    rows_below_0_05 = shapiroTest$adjusted_table %>%
      pivot_longer(cols = -all_of(groupVars), names_to = "variable", values_to = "p.value") %>%
      filter(p.value < 0.05) %>%
      dplyr::select(all_of(groupVars), variable, p.value)
    
    # Store the result in the list
    non_normality_list[[df_name]] = rows_below_0_05
  }
  else{
     # State if all groups were normal
    print(paste(df_name,"all gaussian!"))
  }
}

# Print the non-normality counts table
non_normality_counts = data.frame(
  DataFrame = names(non_normality_list),
  Count = sapply(non_normality_list, nrow)
)

print(non_normality_counts)

# Box-Cox transform
non_normal_dfNames = names(non_normality_list)

for (df_name in DFs) {
  df = get(df_name)
  
  if (df_name %in% non_normal_dfNames){
    if(df_name %in% c("harvest","harvest_xray")){
      redistributed_df = boxcoxRedistributeNongrouping(df, non_normality_list[[df_name]])
    }
    else{
      redistributed_df = boxcoxRedistribute(df, non_normality_list[[df_name]], GroupVars = c("DAT"))
    }
    assign(df_name, redistributed_df)
  }
}


# Count non-normality again
non_normality_counts = data.frame(DataFrame = character(), Count = numeric(), stringsAsFactors = FALSE)

# Initialize the list to store data frames with non-normality information
non_normality_list = list()

p.adjustment = "bonferroni"

for (df_name in DFs) {
  df = get(df_name)
  if(df_name %in% c("harvest","harvest_xray")){
    groupVars = NULL
  }
  else{
    groupVars = c("DAT")
  }
  ignoreVars = c("DAT", "Treatment", "Plant.ID", "Plant.Name","Weight","Weight.After.Watering","Plant_Temp","Env_Temp")
  
  shapiroTest = groupWiseNormalityTest(df, groupVars = groupVars, p.adjust.method = p.adjustment, ignoreVars = ignoreVars,test = "shapiro")
  
  # Check if any numerical column in adjusted_table has a value below 0.05
  if (any(shapiroTest$adjusted_table %>% dplyr::select(-all_of(groupVars)) < 0.05, na.rm = TRUE)) {
    rows_below_0_05 = shapiroTest$adjusted_table %>%
      pivot_longer(cols = -all_of(groupVars), names_to = "variable", values_to = "p.value") %>%
      filter(p.value < 0.05) %>%
      dplyr::select(all_of(groupVars), variable, p.value)
    
    # Store the result in the list
    non_normality_list[[df_name]] = rows_below_0_05
  }
  else{
     # State if all groups were normal
    print(paste(df_name,"all gaussian!"))
  }
}

# Print the non-normality counts table
non_normality_counts = data.frame(
  DataFrame = names(non_normality_list),
  Count = sapply(non_normality_list, nrow)
)

print(non_normality_counts)

# Here we count the number of groups so we can calulate the non-normality rate to report
for (df_name in DFs) {
  df = get(df_name)
  
  if (df_name %in% c("harvest", "harvest_xray")) {
    groupVars = c("Plant.Name")
    ignoreN = 3
  } else {
    groupVars = c("DAT", "Plant.Name")
    ignoreN = 4
  }
  
  # Count the number of distinct group combinations
  n_groups = df %>%
    dplyr::distinct(across(all_of(groupVars))) %>%
    nrow()
  
  # Print the result
  print(paste(df_name, ": Number of groups =", n_groups*(ncol(df)-ignoreN)))
}



###################################################################
#----------------- POOLED STANDARDISATION ------------------------#
###################################################################
# standardize by time point
for (df_name in DFs){
  df = get(df_name)
  if(df_name %in% c("harvest","harvest_xray")){
    groupVars = NULL
  }
  else{
    groupVars = c("DAT")
  }
  if(transformMethod == "zscore"){
    assign(df_name,zScore(df,groupVars = groupVars,ignoreVars = c("Plant.Name","Plant.ID","DAT","Treatment")))
  }
  if(transformMethod == "minmax"){
      assign(df_name,minMaxNormalize(df,groupVars = groupVars,ignoreVars = c("Plant.Name","Plant.ID","DAT","Treatment")))
  }
}

###################################################################
#---------------- POOLED WIDENING AND MERGING --------------------#
###################################################################
# Get a list of all the predictor data frames
PredictorDFs = mget(setdiff(DFs,c("harvest","weight","harvest_xray")))

AllPredictorVars = widen_and_merge(PredictorDFs)

# We'll remove the plant.ID column and use these as row names instead, just leaving the predictor variables 
rownames(AllPredictorVars) = AllPredictorVars$Plant.ID
AllPredictorVars$Plant.ID = NULL


nTree = 20
# I'll add back a genotype and treatment column for grouping first
AllPredictorVars$Genotype = as.factor(sub("_(C|D)_[0-9]+", "", rownames(AllPredictorVars)))
# Create the Treatment column based on the row names
AllPredictorVars$Treatment = as.factor(ifelse(grepl("L[0-9]+_C_[0-9]+", rownames(AllPredictorVars)),
                                               "Control",
                                               ifelse(grepl("L[0-9]+_D_[0-9]+", rownames(AllPredictorVars)),
                                                      "Drought",
                                                      NA)))  # Use NA for any unmatched rows (optional)


# Impute missing values
impute = groupWiseMissForest(AllPredictorVars,groupVars = c("Genotype","Treatment"), ignoreVars = c(), ntree = nTree)
AllPredictorVars = impute$ImputedData

print(impute$OOBError)

#remove genotype column
AllPredictorVars$Genotype = NULL


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(transformMethod == "zscore"){
  save(list = c("AllPredictorVars","harvest","harvest_xray"), file = "../Rout/Rdata/StandardizedByDATWidened.Rdata")
}
if(transformMethod == "minmax"){
  save(list = c("AllPredictorVars","harvest","harvest_xray"), file = "../Rout/Rdata/minMaxTransformedByDATWidened.Rdata")
}

