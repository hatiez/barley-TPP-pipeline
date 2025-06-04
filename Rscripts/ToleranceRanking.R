source("ToleranceRankingFunctions.R")
load("../Rout/Rdata/StandardizedByDATWidened.Rdata")

harvest = as.data.frame(harvest)

if (!dir.exists("../Rout/ToleranceRanking/")) {
  dir.create("../Rout/ToleranceRanking/", recursive = TRUE)
}

################################################################
#--------------------- REMOVING SAMPLES -----------------------#
################################################################

# L7, L8 and L9 are genetically heterogeneous and need to be removed, as this
# method assumes genetic homogeneity
harvest = harvest[!(harvest$Plant.Name %in% c("L7","L8","L9")),]
harvest_xray = harvest_xray[!(harvest_xray$Plant.Name %in% c("L7","L8","L9")),]
AllPredictorVars = AllPredictorVars[!grepl("^L7|^L8|^L9", rownames(AllPredictorVars)), ]


################################################################
#----------------------- HARVEST TRAITS -----------------------#
################################################################

# Univariate drought effect size on harvest traits
effectSizes_harvest = harvest_testByTreatment(harvest, p_adjust_method = "bonferroni",test_method = "t.test")

# Open PDF device
pdf("../Rout/ToleranceRanking/harvest_univariate_barplots.pdf")

# Loop through traits and plot
for (trait in setdiff(colnames(harvest), c("Plant.ID", "Treatment", "Plant.Name"))) {
  plot_effect_size_barplot(effectSizes_harvest$effect_sizes, effectSizes_harvest$test_statistics, trait = trait, include_significance = TRUE)
}

# Close PDF device
dev.off()


# Multivariate drought effect size on harvest traits
# Open PDF device
pdf("../Rout/ToleranceRanking/harvest_permanova_barplots.pdf", height = 5, width = 8)

# using PERMANOVA
harvest_PERMANOVA = PERMANOVA_per_Genotype(harvest, genotypes = unique(harvest$Plant.Name))

plot_permanova_effect_size_barplot(harvest_PERMANOVA$effect_sizes,harvest_PERMANOVA$test_statistics, main = "Harvest Trait Treatment PERMANOVA by Line")

dev.off()


################################################################
#---------------------- TEMPORAL TRAITS -----------------------#
################################################################
# We need to add a genotype column
AllPredictorVars$Plant.Name = as.factor(sub("_(C|D)_[0-9]+", "", rownames(AllPredictorVars)))

# Applying to all plants
predictorEffectsize = PERMANOVA_per_Genotype(AllPredictorVars, permutations = 3000, genotypes = unique(AllPredictorVars$Plant.Name))

# Open PDF device
pdf("../Rout/ToleranceRanking/temporal_permanova_barplots.pdf", height = 5, width = 8)

plot_permanova_effect_size_barplot(predictorEffectsize$effect_sizes, predictorEffectsize$test_statistics, main = "Temporal Traits Treatment PERMANOVA by Genotype")

dev.off()



