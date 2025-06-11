# --- Install and Load Required Packages ---

# Define required packages
required_packages <- c(
  "dplyr", "factoextra", "mixOmics", "ggfortify", "PCAtools",
  "ggplot2", "ggalt", "missForest", "cluster", "ggsci", "rlang", "dplyr"
)

# Step 2: Identify missing packages
missing_packages <- required_packages[!(required_packages %in% rownames(installed.packages()))]

# Step 3: Prompt user to install missing packages
if (length(missing_packages) > 0) {
  message("The following packages are missing:\n", paste(missing_packages, collapse = ", "))
  install <- readline(prompt = "Do you want to install them now? [y/n]: ")
  if (tolower(install) == "y") {
    install.packages(missing_packages, dependencies = TRUE)
    rm(install)
  } else {
    stop("Required packages not installed. Script cannot continue.")
  }
}

# I had trouble with PCAtools and mixOmics. You can try:
# BiocManager::install('PCAtools')
# BiocManager::install('mixOmics')

# Step 4: Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Step 5: Clean up temporary variables
rm(required_packages, missing_packages)





#-------------------------------------------------------------------------------
# Load the dataset
load("../Rout/Rdata/StandardizedByDATWidened.Rdata")

if (!dir.exists("../Rout/")) {
  dir.create("../Rout/", recursive = TRUE)
}
# Output directory
if (!dir.exists("../Rout/Clustering_PCA")) {
  dir.create("../Rout/Clustering_PCA", recursive = TRUE)
}
folder_name = "../Rout/Clustering_PCA"

#Add plant ID and extract plant name from it
AllPredictorVars$Plant.ID<-rownames(AllPredictorVars)
AllPredictorVars$Plant.Name<-NA
AllPredictorVars <- AllPredictorVars %>%
  mutate(Plant.Name = sub("_.*", "", Plant.ID))


# Assign genotype categories (wild vs cultivated)
# This was not used in publication.
# The categories are based on the alleles Hs / Hv in locus HsDry2.2
AllPredictorVars <- AllPredictorVars %>% mutate(
  Plant.Name = as.character(Plant.Name), 
  Genotype = ifelse(Plant.Name %in% c("L2", "L4", "L6"), "wild",
                    ifelse(Plant.Name %in% c("L3", "L5", "L1"), "cultivated", NA))
)

# Remove rows with missing genotype info (genetically heterozygous lines)
AllPredictorVars <- AllPredictorVars %>%
  filter(!is.na(Genotype))


# 5) Build `df_save`: the data frame we'll actually feed into PCA
common_column   <- c("Treatment","Plant.Name","Genotype","Plant.ID")
final_predictor <- setdiff(colnames(AllPredictorVars), common_column)

df_save <- AllPredictorVars %>%
  dplyr::select(all_of(common_column), all_of(final_predictor))


# 6) Perform a standard PCA with `prcomp`
rownames(df_save) <- df_save$Plant.ID

pca.res <- prcomp(
  df_save %>%
    dplyr::select(-Treatment, -Plant.ID, -Genotype, -Plant.Name),
  center = TRUE, scale. = TRUE
)

# Also run PCAtools version if you need loadings etc.
pca_res <- PCAtools::pca(
  t(df_save %>% dplyr::select(-Treatment, -Plant.ID, -Genotype, -Plant.Name)),
  metadata = df_save %>% dplyr::select(Plant.ID, Treatment, Genotype, Plant.Name),
  center   = TRUE, scale = TRUE
)


# 7) Decide how many PCs to keep (e.g. 95% cumulative variance)
explained_variance  <- pca.res$sdev^2 / sum(pca.res$sdev^2)
cumulative_variance <- cumsum(explained_variance)
num_components      <- which(cumulative_variance >= 0.95)[1]


# 8) Extract the PCA scores and re-attach metadata via a join
df_pca <- as.data.frame(pca.res$x)[, 1:num_components]
df_pca$Plant.ID <- rownames(df_pca)

# Build a small metadata tibble
meta <- df_save %>%
  dplyr::select(Plant.ID, Treatment, Genotype, Plant.Name)

# Join so that every row of df_pca has its metadata
df_pca <- df_pca %>%
  left_join(meta, by = "Plant.ID") %>%
  mutate(
    Treatment = factor(Treatment),
    Genotype  = factor(Genotype),
    Plant.Name= factor(Plant.Name)
  )


# 9) PCA visualization using ggfortify
#    We pass the original data to `autoplot()` via the `data=` argument,
#    and then specify aesthetics by column names.
g1 <- autoplot(
  pca.res,
  data   = df_pca,
  colour = "Genotype",
  shape  = "Treatment",
  label  = FALSE,
  loadings = FALSE
) +
  ggtitle("PCA by Genotype and Treatment") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

# Save the plot
ggsave("PCA_genotype_treatment.png", plot = g1,
       path = folder_name, width = 8, height = 6, dpi = 300)


# 10) PCAtools biplot or loadings plot if you want
if (nrow(pca_res$loadings) > 1) {
  # Scree plot of explained variance
  scree <- PCAtools::screeplot(pca_res)
  ggsave("PCA_scree.png", plot = scree,
         path = folder_name, width = 20, height = 15, dpi = 300)
  
  # Loadings plot
  load <- plotloadings(pca_res, labSize = 3, rangeRetain = 0.05)
  ggsave("PCA_loadings.png", plot = load,
         path = folder_name, width = 20, height = 20, dpi = 300)
}


# 11) k-means clustering on the standardized data
#     First, build a trait-only matrix, standardized
trait_matrix <- df_save %>%
  dplyr::select(-all_of(common_column)) %>%
  scale(center = TRUE, scale = TRUE)

# Determine optimal k via silhouette
sil <- fviz_nbclust(trait_matrix, kmeans, method = "silhouette", k.max = 10)
ggsave("nbclust_silhouette.png", plot = sil,
       path = folder_name, width = 6, height = 4, dpi = 300)

best_k <- which.max(sil$data$y)

# Run k-means
set.seed(42)
km <- kmeans(trait_matrix, centers = best_k)

# Visualize clusters
cluster_plot <- fviz_cluster(km, data = trait_matrix) +
  geom_point(aes(color = df_save$Genotype), size = 2) +
  theme_minimal(base_size = 16)

ggsave("kmeans_cluster.png", plot = cluster_plot,
       path = folder_name, width = 8, height = 6, dpi = 300)
