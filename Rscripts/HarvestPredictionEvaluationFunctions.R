# Define required packages
.required_pkgs = c(
  "dplyr", "tidyr", "broom", "purrr", "ggplot2", "randomForest",
  "multcompView", "caret", "openxlsx", "ggpubr", "effsize",
  "rstatix", "DescTools", "rlang", "tools", "xtable"
)

# Identify missing packages
.missing_pkgs = .required_pkgs[!(.required_pkgs %in% installed.packages()[, "Package"])]

# Prompt user to install missing packages
if (length(.missing_pkgs) > 0) {
  message("The following packages are missing: ", paste(.missing_pkgs, collapse = ", "))
  .user_input = readline(prompt = "Do you want to install them now? [y/n]: ")
  if (tolower(.user_input) == "y") {
    install.packages(.missing_pkgs, dependencies = TRUE)
    remove(.user_input)
  } else {
    stop("Required packages not installed. Script cannot continue.")
  }
}

# Load all packages
invisible(lapply(.required_pkgs, library, character.only = TRUE))

# Clean up helper variables
rm(.required_pkgs, .missing_pkgs)



# Function for making shorter harvest trait names
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
  
  # Apply replacements only to matching values while preserving other labels
  x <- dplyr::recode(x, !!!replacements)
  
  return(x)
}

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

#' Compare Model Accuracies by Factor Groupings
#'
#' This function creates boxplots to compare model accuracy metrics (e.g., R², RMSE, MAE) 
#' across different groups defined by categorical variables in cross-validation results.
#' It supports statistical comparisons and faceting for multi-panel plots.
#'
#' @param allAccs A data frame containing cross-validation accuracy results and grouping variables.
#' @param accuracy_measure The accuracy metric to plot. Options: "Rsquared", "RMSE", or "MAE".
#' @param data_group_col The column name used to group data on the x-axis (e.g., model type or data setting).
#' @param color_by Optional column name to color boxplots; defaults to `data_group_col`.
#' @param pageBy Optional column name to separate plots by facet/page (e.g., response variable).
#' @param color_palette A vector of colors used to distinguish groups in the plot.
#' @param test_type The statistical test for comparison: "wilcox" (default) or "t".
#' @param p_adjust_method Method for p-value adjustment in pairwise tests (e.g., "holm", "fdr").
#' @param compare_scope Scope for pairwise comparisons: "all" (across all combinations) or "within" (within `data_group_col`).
#' @param show_comparisons Logical; whether to compute and include pairwise comparisons.
#' @param annotate_comparisons Logical; whether to annotate pairwise comparisons with significance labels.
#' @param annotate_global_test Logical; whether to include a global test p-value (e.g., ANOVA or Kruskal-Wallis).
#' @param main_title Optional base title for the plots.
#'
#' @return A named list where each element corresponds to a `pageBy` group. Each element contains:
#'         - `global_p`: Global test p-value (if applicable)
#'         - `pairwise`: Pairwise p-values matrix or list
#' 
#' @details
#' The function generates one plot per level of the `pageBy` factor. For two-group comparisons, 
#' it uses a direct test (t-test or Wilcoxon). For >2 groups, it performs a global test and optionally pairwise comparisons.
#' 
#' @examples
#' # Compare accuracy grouped by model and dataset:
#' compare_model_accuracies_byFactor(allAccs = my_results,
#'                                   accuracy_measure = "Rsquared",
#'                                   data_group_col = "Model",
#'                                   color_by = "Dataset",
#'                                   pageBy = "Target",
#'                                   compare_scope = "within")
compare_model_accuracies_byFactor = function(allAccs, accuracy_measure = "Rsquared",
                                              data_group_col,      # Factor column for x-axis grouping
                                              color_by = NULL,     # Factor column for box colors; if NULL, use data_group_col.
                                              pageBy = NULL,
                                              color_palette = c("orange", "lightblue", "lightgreen", "pink"),
                                              test_type = "wilcox",        # Either "wilcox" or "t"
                                              p_adjust_method = "holm",     # e.g., "fdr", "bonferroni", "holm", etc.
                                              compare_scope = "all",       # "all" or "within"
                                              show_comparisons = TRUE,     # Whether to compute and return pairwise p-values
                                              annotate_comparisons = TRUE,  # Whether to add annotations using stat_compare_means()
                                              annotate_global_test = FALSE,  #whether to add global test annotation
                                              main_title = NULL
) {
  ### --- Basic Checks ---
  if (!accuracy_measure %in% c("RMSE", "Rsquared", "MAE"))
    stop("Invalid accuracy measure. Choose from 'RMSE', 'Rsquared', or 'MAE'.")
  if (!pageBy %in% names(allAccs))
    stop("The data frame must have an 'pageBy' column.")
  if (is.null(color_by))
    color_by = data_group_col
  if (!data_group_col %in% names(allAccs))
    stop("The data grouping column is not found in the data frame.")
  if (!color_by %in% names(allAccs))
    stop("The coloring column is not found in the data frame.")
  
  plotPages = unique(allAccs[[pageBy]])
  results_list = list()
  plot_list = list()
  
  for (thisPage in plotPages) {
    dat = allAccs[allAccs[[pageBy]] == thisPage, ]
    
    # Decide the plot title depending on main_title
    # Otherwise, use main_title as-is.
    if (is.null(main_title)) {
      plot_title = paste("Model Accuracy Comparison -", thisPage)
    } else {
      plot_title = paste(main_title, thisPage)
    }
    
    ## ========== ONE-FACTOR CASE (data_group_col == color_by) ==========
    if (data_group_col == color_by) {
      dat[[data_group_col]] = factor(dat[[data_group_col]])
      group_var = data_group_col
      groups = levels(dat[[group_var]])
      
      # -- If exactly 2 groups, skip global test entirely --
      if (length(groups) == 2) {
        # Normality check if test_type == "t"
        if (test_type == "t") {
          # By default assume "t"
          pairwise_method = "t"
          for (g in groups) {
            vals = dat[dat[[group_var]] == g, accuracy_measure]
            if (length(vals) >= 3 && shapiro.test(vals)$p.value < 0.05) {
              pairwise_method = "wilcox"
              warning(paste("For pageBy:", thisPage,
                            "- Group:", g,
                            "non-normal data; switching to Wilcoxon pairwise test."))
              break
            }
          }
        } else {
          pairwise_method = "wilcox"
        }
        
        # Direct pairwise test (2-group scenario)
        if (pairwise_method == "t") {
          t_res = t.test(as.formula(paste(accuracy_measure, "~", group_var)), data = dat)
          pairwise_pvals = t_res$p.value
        } else {
          w_res = wilcox.test(as.formula(paste(accuracy_measure, "~", group_var)), data = dat)
          pairwise_pvals = w_res$p.value
        }
        
        p_global = NA  # no global test
        my_comparisons = list(groups)  # single comparison: 2 groups
        
        # Plot
        color_levels = levels(dat[[color_by]])
        num_levels = length(color_levels)
        color_palette_subset = rep(color_palette, length.out = num_levels)
        
        p = ggboxplot(dat, x = group_var, y = accuracy_measure,
                       fill = color_by, palette = color_palette_subset, add = "none") +
          scale_fill_manual(values = color_palette_subset) +
          labs(title = plot_title,
               y = if (accuracy_measure == "Rsquared") expression(R^2) else paste(accuracy_measure))
        
        # If we want to annotate the pairwise result:
        if (show_comparisons) {
          p = p + stat_compare_means(method = pairwise_method,
                                      comparisons = my_comparisons,
                                      label = "p.signif")
        }
        
        plot_list[[as.character(thisPage)]] = p
        results_list[[as.character(thisPage)]] = list(global_p = p_global, pairwise = pairwise_pvals)
        
      } else {
        # ----------  if length(groups) > 2 ----------
        
        ## Check normality only if the chosen test is "t"
        if (test_type == "t") {
          non_normal_found = FALSE
          for (g in groups) {
            vals = dat[dat[[group_var]] == g, accuracy_measure]
            if (length(vals) >= 3 && shapiro.test(vals)$p.value < 0.05) {
              non_normal_found = TRUE
              warning(paste("For pageBy:", thisPage, 
                            "- Group:", g, 
                            "data is not normal; switching global test to 'kruskal.test' and pairwise test to 'wilcox'."))
            }
          }
          if (non_normal_found) {
            global_method = "kruskal.test"
            pairwise_method = "wilcox"
          } else {
            global_method = "anova"
            pairwise_method = "t"
          }
        } else {
          global_method = "kruskal.test"
          pairwise_method = "wilcox"
        }
        
        ## Global test
        if (global_method == "anova") {
          anova_res = aov(as.formula(paste(accuracy_measure, "~", group_var)), data = dat)
          p_global = summary(anova_res)[[1]][["Pr(>F)"]][1]
        } else {
          kw = kruskal.test(as.formula(paste(accuracy_measure, "~", group_var)), data = dat)
          p_global = kw$p.value
        }
        
        ## Define pairwise comparisons
        my_comparisons = combn(groups, 2, simplify = FALSE)
        if (p_global < 0.05) {
          if (pairwise_method == "t") {
            p_test = pairwise.t.test(dat[[accuracy_measure]], dat[[group_var]], p.adjust.method = p_adjust_method)
          } else {
            p_test = pairwise.wilcox.test(dat[[accuracy_measure]], dat[[group_var]], p.adjust.method = p_adjust_method)
          }
          pairwise_pvals = p_test$p.value
        } else {
          pairwise_pvals = NA
        }
        
        ## Define color palette
        color_levels = levels(dat[[color_by]])
        num_levels = length(color_levels)
        color_palette_subset = rep(color_palette, length.out = num_levels)
        
        ## Create the one-factor plot
        p = ggboxplot(dat, x = group_var, y = accuracy_measure,
                       fill = color_by, palette = color_palette_subset, add = "none") +
          scale_fill_manual(values = color_palette_subset) +
          labs(title = plot_title,
               y = if (accuracy_measure == "Rsquared") expression(R^2) else paste(accuracy_measure))
        
        ## Global test annotation (annotate with the global test method)
        if (annotate_global_test) {
          p = p + stat_compare_means(method = global_method,
                                      label.y = max(dat[[accuracy_measure]], na.rm = TRUE) * 1.1)
        }
        
        ## Pairwise annotations, if requested
        if (show_comparisons) {
          if (p_global < 0.05) {
            p = p + stat_compare_means(comparisons = my_comparisons,
                                        method = pairwise_method, label = "p.signif")
          } else {
            p = p + stat_compare_means(label = "ns",
                                        label.y = max(dat[[accuracy_measure]], na.rm = TRUE) * 1.1)
          }
        }
        
        plot_list[[as.character(thisPage)]] = p
        results_list[[as.character(thisPage)]] = list(global_p = p_global, pairwise = pairwise_pvals)
      }
      
    } else {
      ## ========== TWO-FACTOR CASE ==========
      dat[[data_group_col]] = factor(dat[[data_group_col]])
      dat[[color_by]] = factor(dat[[color_by]])
      
      if (compare_scope == "all") {
        ## Create a combined factor for all-vs-all comparisons
        dat$Group = factor(paste(dat[[data_group_col]], dat[[color_by]], sep = " - "))
        groups = levels(dat$Group)
        
        # If exactly 2 combined groups, skip global test:
        if (length(groups) == 2) {
          # Normality check if test_type == "t"
          if (test_type == "t") {
            pairwise_method = "t"
            for (g in groups) {
              vals = dat[dat$Group == g, accuracy_measure]
              if (length(vals) >= 3 && shapiro.test(vals)$p.value < 0.05) {
                pairwise_method = "wilcox"
                warning(paste("For pageBy:", thisPage,
                              "- Group:", g,
                              "non-normal data; switching to Wilcoxon pairwise test."))
                break
              }
            }
          } else {
            pairwise_method = "wilcox"
          }
          # Direct pairwise test
          # We can do "t.test" or "wilcox.test" for 2 groups
          levels_g = levels(dat$Group)
          if (pairwise_method == "t") {
            t_res = t.test(as.formula(paste(accuracy_measure, "~ Group")), data = dat)
            pairwise_pvals = t_res$p.value
          } else {
            w_res = wilcox.test(as.formula(paste(accuracy_measure, "~ Group")), data = dat)
            pairwise_pvals = w_res$p.value
          }
          
          p_global = NA
          my_comparisons = list(levels_g)
          
          # Plot
          color_levels = levels(dat[[color_by]])
          num_levels = length(color_levels)
          color_palette_subset = rep(color_palette, length.out = num_levels)
          
          p = ggboxplot(dat, x = "Group", y = accuracy_measure,
                         fill = color_by, palette = color_palette_subset, add = "none") +
            scale_fill_manual(values = color_palette_subset) +
            scale_x_discrete(labels = function(x) sapply(x, function(group) {
              strsplit(group, " - ")[[1]][1]  # show only data_group_col
            })) +
            labs(title = plot_title,
                 y = if (accuracy_measure == "Rsquared") expression(R^2) else paste(accuracy_measure),
                 x = NULL) +
            theme(axis.text.x = element_text(size = 12, color = "black"),
                  axis.ticks.x = element_line(),
                  axis.title.x = element_blank())
          
          if (show_comparisons) {
            p = p + stat_compare_means(method = pairwise_method,
                                        comparisons = my_comparisons,
                                        label = "p.signif")
          }
          
          plot_list[[as.character(thisPage)]] = p
          results_list[[as.character(thisPage)]] = list(global_p = p_global, pairwise = pairwise_pvals)
          
        } else {
          # > 2 groups => do global + pairwise
          
          ## Check normality on the combined groups (if t-test chosen)
          if (test_type == "t") {
            non_normal_found = FALSE
            for (g in groups) {
              vals = dat[dat$Group == g, accuracy_measure]
              if (length(vals) >= 3 && shapiro.test(vals)$p.value < 0.05) {
                non_normal_found = TRUE
                warning(paste("For pageBy:", thisPage, 
                              "- Group:", g, 
                              "data is not normal; switching global test to 'kruskal.test' and pairwise test to 'wilcox'."))
                break
              }
            }
            if (non_normal_found) {
              global_method = "kruskal.test"
              pairwise_method = "wilcox"
            } else {
              global_method = "anova"
              pairwise_method = "t"
            }
          } else {
            global_method = "kruskal.test"
            pairwise_method = "wilcox"
          }
          
          if (global_method == "anova") {
            anova_res = aov(as.formula(paste(accuracy_measure, "~ Group")), data = dat)
            p_global = summary(anova_res)[[1]][["Pr(>F)"]][1]
          } else {
            kw = kruskal.test(as.formula(paste(accuracy_measure, "~ Group")), data = dat)
            p_global = kw$p.value
          }
          
          my_comparisons = combn(groups, 2, simplify = FALSE)
          if (p_global < 0.05) {
            if (pairwise_method == "t") {
              p_test = pairwise.t.test(dat[[accuracy_measure]], dat$Group, p.adjust.method = p_adjust_method)
            } else {
              p_test = pairwise.wilcox.test(dat[[accuracy_measure]], dat$Group, p.adjust.method = p_adjust_method)
            }
            pairwise_pvals = p_test$p.value
          } else {
            pairwise_pvals = NA
          }
          
          # For plotting
          color_levels = levels(dat[[color_by]])
          num_levels = length(color_levels)
          color_palette_subset = rep(color_palette, length.out = num_levels)
          
          p = ggboxplot(dat, x = "Group", y = accuracy_measure,
                         fill = color_by, palette = color_palette_subset, add = "none") +
            scale_fill_manual(values = color_palette_subset) +
            scale_x_discrete(labels = function(x) sapply(x, function(group) {
              strsplit(group, " - ")[[1]][1]
            })) +
            labs(title = plot_title,
                 y = if (accuracy_measure == "Rsquared") expression(R^2) else paste(accuracy_measure),
                 x = NULL) +
            theme(axis.text.x = element_text(size = 12, color = "black"),
                  axis.ticks.x = element_line(),
                  axis.title.x = element_blank())
          
          if (annotate_global_test) {
            p = p + stat_compare_means(method = global_method,
                                        label.y = max(dat[[accuracy_measure]], na.rm = TRUE) * 1.1)
          }
          
          if (show_comparisons) {
            if (p_global < 0.05) {
              p = p + stat_compare_means(comparisons = my_comparisons,
                                          method = pairwise_method, label = "p.signif")
            } else {
              p = p + stat_compare_means(label = "ns",
                                          label.y = max(dat[[accuracy_measure]], na.rm = TRUE) * 1.05)
            }
          }
          
          plot_list[[as.character(thisPage)]] = p
          results_list[[as.character(thisPage)]] = list(global_p = p_global, pairwise = pairwise_pvals)
        }
        
      } else if (compare_scope == "within") {
        ## Two-Factor Case: "within" comparisons (faceted by data_group_col)
        levels_dg = levels(dat[[data_group_col]])
        color_levels = levels(dat[[color_by]])
        num_levels = length(color_levels)
        color_palette_subset = rep(color_palette, length.out = num_levels)
        
        p = ggboxplot(dat, x = color_by, y = accuracy_measure,
                       fill = color_by, palette = color_palette_subset, add = "none") +
          scale_fill_manual(values = color_palette_subset) +
          facet_wrap(as.formula(paste("~", data_group_col)), nrow = 1, strip.position = "bottom")+
          labs(title = plot_title,
               y = if (accuracy_measure == "Rsquared") expression(R^2) else paste(accuracy_measure),
               x = NULL) +
          theme(strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(size = 12),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank())
        
        global_p_values = list()
        pairwise_list = list()
        
        for (lev in levels_dg) {
          subdat = dat[dat[[data_group_col]] == lev, ]
          groups_sub = levels(subdat[[color_by]])
          
          # If exactly 2 subgroups => skip global test
          if (length(groups_sub) == 2) {
            if (test_type == "t") {
              pairwise_method = "t"
              for (g in groups_sub) {
                vals = subdat[subdat[[color_by]] == g, accuracy_measure]
                if (length(vals) >= 3 && shapiro.test(vals)$p.value < 0.05) {
                  pairwise_method = "wilcox"
                  warning(paste("For pageBy:", thisPage, 
                                "- Group:", g, "in facet", lev,
                                "non-normal data; switching to Wilcoxon pairwise test."))
                  break
                }
              }
            } else {
              pairwise_method = "wilcox"
            }
            # Do the pairwise test directly
            if (pairwise_method == "t") {
              t_res = t.test(as.formula(paste(accuracy_measure, "~", color_by)), data = subdat)
              pairwise_list[[lev]] = t_res$p.value
            } else {
              w_res = wilcox.test(as.formula(paste(accuracy_measure, "~", color_by)), data = subdat)
              pairwise_list[[lev]] = w_res$p.value
            }
            global_p_values[[lev]] = NA  # no global test in 2-group scenario
            
          } else {
            # > 2 subgroups => do global test
            if (test_type == "t") {
              non_normal_found = FALSE
              for (g in groups_sub) {
                vals = subdat[subdat[[color_by]] == g, accuracy_measure]
                if (length(vals) >= 3 && shapiro.test(vals)$p.value < 0.05) {
                  non_normal_found = TRUE
                  warning(paste("For pageBy:", thisPage, 
                                "- Group:", g, "in facet", lev,
                                "data is not normal; switching global test to 'kruskal.test' and pairwise test to 'wilcox'."))
                  break
                }
              }
              if (non_normal_found) {
                facet_global_method = "kruskal.test"
                facet_pairwise_method = "wilcox"
              } else {
                facet_global_method = "anova"
                facet_pairwise_method = "t"
              }
            } else {
              facet_global_method = "kruskal.test"
              facet_pairwise_method = "wilcox"
            }
            
            if (facet_global_method == "anova") {
              anova_res = aov(as.formula(paste(accuracy_measure, "~", color_by)), data = subdat)
              p_val = summary(anova_res)[[1]][["Pr(>F)"]][1]
            } else {
              kw = kruskal.test(as.formula(paste(accuracy_measure, "~", color_by)), data = subdat)
              p_val = kw$p.value
            }
            global_p_values[[lev]] = p_val
            
            if (p_val < 0.05) {
              if (facet_pairwise_method == "t") {
                p_test = pairwise.t.test(subdat[[accuracy_measure]], subdat[[color_by]], p.adjust.method = p_adjust_method)
              } else {
                p_test = pairwise.wilcox.test(subdat[[accuracy_measure]], subdat[[color_by]], p.adjust.method = p_adjust_method)
              }
              pairwise_list[[lev]] = p_test$p.value
            } else {
              pairwise_list[[lev]] = NA
            }
          }
        }
        
        # Add global annotation - uses the last facet's method as a placeholder
        # (ggpubr won't annotate per facet automatically, but let's do a single line for completeness)
        if (length(color_levels) > 2) {
          # If more than 2 subgroups in any facet, pick the last calculated facet_global_method
          # (this is a known limitation: a single global annotation won't reflect different facets)
          if (annotate_global_test) {
            p = p + stat_compare_means(method = global_method,
                                        label.y = max(dat[[accuracy_measure]], na.rm = TRUE) * 1.1)
          }
        }
        
        if (show_comparisons) {
          comps = combn(levels(dat[[color_by]]), 2, simplify = FALSE)
          p = p + stat_compare_means(comparisons = comps,
                                      method = if (test_type == "t") "t" else "wilcox",
                                      label = "p.signif")
        }
        
        plot_list[[as.character(thisPage)]] = p
        results_list[[as.character(thisPage)]] = list(global_p = global_p_values,
                                                       pairwise = pairwise_list)
      }
    }
  }  # End thisPage loop
  
  # Print each plot
  for (n in names(plot_list)) {
    print(plot_list[[n]])
  }
  
  return(results_list)
}

#' Plot Random Forest Variable Importance
#'
#' Generates a horizontal barplot of the top `k` most important predictors
#' from a caret-trained random forest model. Trait names are cleaned using
#' `clean_label_TemporalTraits()` for readability.
#'
#' @param rf_model A caret-trained random forest model (`method = "rf"`).
#' @param top_k Integer. Number of top predictors to display. Defaults to 20.
#' @param main_title Plot title.
#' @param show_mean_r2 Logical. Whether to show mean R² and standard deviation from cross-validation.
#'
#' @return A base R barplot is displayed. No return value.
plot_rf_importance = function(rf_model, top_k = 20, main_title = "Top K Important Variables", show_mean_r2 = TRUE) {
  # Extract variable importance from the caret model
  var_importance = varImp(rf_model)$importance
  importance_df = data.frame(
    Variable = rownames(var_importance),
    Importance = var_importance$Overall,
    stringsAsFactors = FALSE
  )
  
  # Sort and select top k
  top_variables = head(importance_df[order(-importance_df$Importance), ], top_k)
  
  # Clean predictor names
  top_variables$Cleaned = sapply(top_variables$Variable, clean_label_TemporalTraits)
  
  # Preserve order in factor levels
  top_variables$Cleaned = factor(top_variables$Cleaned, levels = top_variables$Cleaned[order(-top_variables$Importance)])
  
  # Set up the plot
  par(mar = c(5, 12, 4, 4), oma = c(1, 3, 1, 3))
  
  bp = barplot(top_variables$Importance,
                names.arg = top_variables$Cleaned,
                las = 2, horiz = TRUE,
                main = main_title,
                xlab = "Importance",
                col = "lightblue",
                border = "black",
                xlim = c(0, max(top_variables$Importance) * 1.19))
  
  if (show_mean_r2) {
    mean_r2 = mean(rf_model$resample$Rsquared, na.rm = TRUE)
    sd_r2 = sd(rf_model$resample$Rsquared, na.rm = TRUE)
    r2_label = sprintf("Mean R^2: %.2f ± %.2f", mean_r2, sd_r2)
    mtext(side = 3, line = -6, r2_label, adj = 1, cex = 1)
  }
}

#' Plot LASSO Coefficients
#'
#' Creates a horizontal barplot of non-zero LASSO coefficients from a caret-trained
#' `glmnet` model. Predictor names are cleaned with `clean_label_TemporalTraits()`.
#'
#' @param lasso_model A caret-trained LASSO model (`method = "glmnet"`).
#' @param minAbsCoef Numeric. Minimum absolute coefficient value to retain. Default is 0 (no filtering).
#' @param main_title Plot title.
#' @param show_mean_accuracy Logical. Whether to show mean R² and standard deviation from cross-validation.
#'
#' @return Invisibly returns a list with one element: `plotted_vars`, the uncleaned names of plotted predictors.
plot_lasso_importance = function(lasso_model,
                                  minAbsCoef = 0,
                                  main_title = "LASSO Coefficients",
                                  show_mean_accuracy = TRUE) {
  if (!inherits(lasso_model, "train") || lasso_model$method != "glmnet") {
    stop("This function only supports LASSO models trained using caret with method = 'glmnet'.")
  }
  
  lasso_coefs = as.matrix(coef(lasso_model$finalModel, s = lasso_model$bestTune$lambda))
  lasso_coefs = lasso_coefs[-1, , drop = FALSE]
  
  coef_df = data.frame(
    Variable = rownames(lasso_coefs),
    Coefficient = lasso_coefs[, 1],
    stringsAsFactors = FALSE
  )
  
  coef_df = subset(coef_df, Coefficient != 0)
  
  if (minAbsCoef > 0) {
    num_filtered = sum(abs(coef_df$Coefficient) < minAbsCoef)
    if (num_filtered > 0) {
      message(sprintf("Filtering out %d coefficients with absolute values below %.4f",
                      num_filtered, minAbsCoef))
    }
    coef_df = subset(coef_df, abs(Coefficient) >= minAbsCoef)
  }
  
  coef_df = coef_df[order(-coef_df$Coefficient), ]
  coef_df$Cleaned = sapply(coef_df$Variable, clean_label_TemporalTraits)
  
  old_par = par(mar = c(5, 12, 1, 1), oma = c(1, 3, 1, 3))
  on.exit(par(old_par), add = TRUE)
  
  y_positions = seq_len(nrow(coef_df))
  x_min = min(coef_df$Coefficient)
  x_max = max(coef_df$Coefficient) * 1.07
  
  plot(
    0, 0, type = "n",
    xlim = c(1.1 * x_min, 1.1 * x_max),
    ylim = c(0.5, nrow(coef_df) + 1.5),
    xlab = "Coefficient",
    ylab = "",
    yaxt = "n",
    main = main_title
  )
  
  abline(h = c(y_positions - 0.5, max(y_positions) + 0.5), col = "gray80", lty = 1)
  abline(v = 0, col = "gray80", lty = 1)
  axis(2, at = y_positions, labels = coef_df$Cleaned, las = 2)
  
  half_height = 0.4
  for (i in seq_len(nrow(coef_df))) {
    this_coef = coef_df$Coefficient[i]
    y_i = y_positions[i]
    rect(
      xleft = min(0, this_coef),
      ybottom = y_i - half_height,
      xright = max(0, this_coef),
      ytop = y_i + half_height,
      col = "lightblue",
      border = "black"
    )
  }
  
  if (show_mean_accuracy && "Rsquared" %in% colnames(lasso_model$resample)) {
    mean_accuracy = mean(lasso_model$resample$Rsquared, na.rm = TRUE)
    sd_accuracy = sd(lasso_model$resample$Rsquared, na.rm = TRUE)
    accuracy_label = bquote("Mean R"^2 ~ ": " ~ .(sprintf("%.2f ± %.2f", mean_accuracy, sd_accuracy)))
    mtext(side = 3, line = -0, accuracy_label, adj = 0, cex = 1)
  }
  
  invisible(list(plotted_vars = coef_df$Variable))
}



#' Re-standardize Predictor Data from Pooled to Treatment-specific Scale
#'
#' Converts z-scored predictor data (standardized on pooled stats) to a new z-score scale
#' based on treatment-specific mean and SD.
#'
#' @param z_data Data frame of pooled z-score standardized predictors. Column names must match `traitMeans$WideName`.
#' @param traitMeans Data frame with columns: `WideName`, `Treatment`, `Mean`, `SD`.
#' @param target_treatment Character string. Name of the treatment group to re-standardize to (e.g., `"Control"` or `"Drought"`).
#'
#' @return A data frame of z-scored predictors on the treatment-specific scale.
restandardize = function(z_data, traitMeans, target_treatment) {
  # z_data: data frame with predictors on the pooled z-score scale.
  # traitMeans: data frame of summary stats with columns: WideName, Treatment, Mean, SD.
  # target_treatment: The target group for re-transformation ("Control" or "Drought").
  
  # Assume the column names of z_data match the 'WideName' values.
  predictors = colnames(z_data)
  
  # Subset the summary stats for pooled and for the target treatment
  pooled_stats = traitMeans[ traitMeans$Treatment == "Pooled", ]
  target_stats = traitMeans[ traitMeans$Treatment == target_treatment, ]
  
  # Match the order of the predictor columns to the rows in the summary stats
  pooled_stats = pooled_stats[ match(predictors, pooled_stats$WideName), ]
  target_stats = target_stats[ match(predictors, target_stats$WideName), ]
  
  # Reverse the pooled z-transformation to retrieve raw values
  # raw = z_data * pooled_SD + pooled_Mean
  raw_data = sweep(z_data, 2, pooled_stats$SD, FUN = "*")
  raw_data = sweep(raw_data, 2, pooled_stats$Mean, FUN = "+")
  
  # Re-apply the target z-transformation: new_z = (raw - target_Mean) / target_SD
  new_z = sweep(raw_data, 2, target_stats$Mean, FUN = "-")
  new_z = sweep(new_z, 2, target_stats$SD, FUN = "/")
  
  return(new_z)
}




