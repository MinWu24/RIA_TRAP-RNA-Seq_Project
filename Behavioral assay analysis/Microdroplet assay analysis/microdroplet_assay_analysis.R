# --------------
# Droplet assay
# --------------

assay = read.table("results.csv",sep=",",header = TRUE)
str(assay)
assay$Strain = factor(assay$Strain, levels = c("N2","mutant"))
levels(assay$Strain)

# Load required libraries
library(lme4)
library(ggplot2)
library(dplyr)
library(car)

# Set Strain as factor with specific levels
assay$Strain <- factor(assay$Strain, levels = c("N2", "mutant"))

# linear-mixed model
model_final <- lmer(Learning_index ~ Strain + (1 | Date), 
                    data = assay,
                    control = lmerControl(optimizer = "bobyqa"))

# Check if Date has any variation
var_date <- VarCorr(model_final)[[1]][1,1]
cat("Date variance:", var_date, "\n")
#Date variance: 1.693934e-18 
cat("Residual variance:", attr(VarCorr(model_final), "sc")^2, "\n")
#Residual variance: 0.007949926 

# 1. Get model predictions
assay$predicted <- fitted(model_final)
assay$residuals <- residuals(model_final)

# 2. Calculate marginal means for each strain
strain_means <- data.frame(
  Strain = c("N2", "mutant"),
  Mean = c(fixef(model_final)[1], 
           fixef(model_final)[1] + fixef(model_final)[2]),
  SE = c(summary(model_final)$coefficients[1, 2],
         summary(model_final)$coefficients[2, 2])
)
print(strain_means)

# 3. Calculate effect size (Cohen's d)
residual_sd <- attr(VarCorr(model_final), "sc")
strain_effect <- fixef(model_final)[2]
cohens_d <- abs(strain_effect) / residual_sd

# 4. Variance partitioning (ICC)
var_comp <- as.data.frame(VarCorr(model_final))
icc <- var_comp[1, "vcov"] / (var_comp[1, "vcov"] + var_comp[2, "vcov"])

# Fix random effects column names
if ("grp" %in% names(var_comp)) {
  random_effects <- var_comp[, c("grp", "var1", "var2", "vcov", "sdcor")]
  colnames(random_effects) <- c("Group", "Variable1", "Variable2", "Variance", "Std_Dev")
} else if ("groups" %in% names(var_comp)) {
  random_effects <- var_comp[, c("groups", "var1", "var2", "vcov", "sdcor")]
  colnames(random_effects) <- c("Group", "Variable1", "Variable2", "Variance", "Std_Dev")
} else {
  random_effects <- var_comp
  names(random_effects) <- c("Group", "Variable1", "Variable2", "Variance", "Std_Dev")
}

# 5. Confidence intervals for fixed effects
conf_int <- confint(model_final, method = "Wald")

# ============================================
# SAVE ALL STATISTICS TO CSV FILES
# ============================================

# 1. Save model summary statistics
model_summary <- summary(model_final)
fixed_effects <- as.data.frame(model_summary$coefficients)
fixed_effects$Parameter <- rownames(fixed_effects)
fixed_effects <- fixed_effects[, c("Parameter", "Estimate", "Std. Error", "t value")]
colnames(fixed_effects) <- c("Parameter", "Estimate", "Std_Error", "t_value")
write.csv(fixed_effects, "Fixed_Effects.csv", row.names = FALSE)

# 2. Save random effects
write.csv(random_effects, "Random_Effects.csv", row.names = FALSE)

# 3. Save strain means
write.csv(strain_means, "Strain_Means.csv", row.names = FALSE)

# 4. Save confidence intervals
conf_int_df <- as.data.frame(conf_int)
conf_int_df$Parameter <- rownames(conf_int_df)
conf_int_df <- conf_int_df[, c("Parameter", "2.5 %", "97.5 %")]
colnames(conf_int_df) <- c("Parameter", "CI_2.5", "CI_97.5")
write.csv(conf_int_df, "Confidence_Intervals.csv", row.names = FALSE)

# 5. Save ANOVA results
anova_results <- Anova(model_final, type = 3, test.statistic = "F")
anova_df <- as.data.frame(anova_results)
anova_df$Parameter <- rownames(anova_df)

# Check which columns exist and select appropriately
if ("F" %in% names(anova_df) & "Df" %in% names(anova_df) & "Df.res" %in% names(anova_df) & "Pr(>F)" %in% names(anova_df)) {
  anova_df <- anova_df[, c("Parameter", "F", "Df", "Df.res", "Pr(>F)")]
  colnames(anova_df) <- c("Parameter", "F_value", "Df", "Df_residual", "p_value")
} else if ("F" %in% names(anova_df) & "Df" %in% names(anova_df) & "Pr(>F)" %in% names(anova_df)) {
  anova_df <- anova_df[, c("Parameter", "F", "Df", "Pr(>F)")]
  colnames(anova_df) <- c("Parameter", "F_value", "Df", "p_value")
} else {
  # Keep all columns if structure is different
  anova_df <- anova_df[, c("Parameter", names(anova_df)[1:min(4, ncol(anova_df)-1)])]
}
write.csv(anova_df, "ANOVA_Results.csv", row.names = FALSE)

# 6. Save variance components and ICC
variance_stats <- data.frame(
  Metric = c("Date_Variance", "Residual_Variance", "Total_Variance", "ICC", "Cohen_s_d"),
  Value = c(var_comp[1, "vcov"], 
            var_comp[2, "vcov"], 
            var_comp[1, "vcov"] + var_comp[2, "vcov"],
            icc,
            cohens_d)
)
write.csv(variance_stats, "Variance_Statistics.csv", row.names = FALSE)

# 7. Save model predictions with actual data
predictions_df <- data.frame(
  Observation = 1:nrow(assay),
  Strain = assay$Strain,
  date = assay$date,
  Learning_Index = assay$Learning_Index,
  Predicted = assay$predicted,
  Residuals = assay$residuals
)
write.csv(predictions_df, "Model_Predictions.csv", row.names = FALSE)

# 8. Save model fit statistics
model_fit <- data.frame(
  Metric = c("AIC", "BIC", "LogLik", "REML_criterion"),
  Value = c(AIC(model_final), 
            BIC(model_final), 
            as.numeric(logLik(model_final)), 
            REMLcrit(model_final))
)
write.csv(model_fit, "Model_Fit_Statistics.csv", row.names = FALSE)

# 9. Create a comprehensive results summary
summary_stats <- data.frame(
  Parameter = c(
    "Number_of_Observations",
    "Number_of_Dates",
    "Number_of_Strains",
    "Intercept_N2_Mean",
    "Intercept_N2_SE",
    "Strain_mutant_Effect",
    "Strain_mutant_SE",
    "Strain_mutant_t_value",
    "Strain_mutant_p_value",
    "Cohen_s_d",
    "ICC_Date",
    "Date_Variance",
    "Residual_Variance",
    "Total_Variance",
    "AIC",
    "BIC",
    "LogLik"
  ),
  Value = c(
    nrow(assay),
    length(unique(assay$date)),
    length(unique(assay$Strain)),
    fixef(model_final)[1],
    summary(model_final)$coefficients[1, 2],
    fixef(model_final)[2],
    summary(model_final)$coefficients[2, 2],
    summary(model_final)$coefficients[2, 3],
    anova_results["Strain", "Pr(>F)"],
    cohens_d,
    icc,
    var_comp[1, "vcov"],
    var_comp[2, "vcov"],
    var_comp[1, "vcov"] + var_comp[2, "vcov"],
    AIC(model_final),
    BIC(model_final),
    as.numeric(logLik(model_final))
  )
)
write.csv(summary_stats, "Summary_Statistics.csv", row.names = FALSE)

# 10. Save strain means with confidence intervals
strain_means_ci <- strain_means
strain_means_ci$CI_lower <- strain_means$Mean - 1.96 * strain_means$SE
strain_means_ci$CI_upper <- strain_means$Mean + 1.96 * strain_means$SE
write.csv(strain_means_ci, "Strain_Means_with_CI.csv", row.names = FALSE)

# ============================================
# CREATE A MASTER RESULTS FILE
# ============================================

# Create a single CSV with all metrics (transposed format)
master_table <- data.frame(
  Metric = c(
    "Number_of_Observations",
    "Number_of_Dates", 
    "Number_of_Strains",
    "Intercept_N2_Mean",
    "Intercept_N2_SE",
    "Strain_syb9641_Effect",
    "Strain_syb9641_SE",
    "Strain_syb9641_t_value",
    "Strain_syb9641_p_value",
    "Cohen_s_d",
    "ICC_Date",
    "Date_Variance",
    "Residual_Variance",
    "Total_Variance",
    "AIC",
    "BIC",
    "LogLik",
    "REML_criterion"
  ),
  Value = c(
    nrow(assay),
    length(unique(assay$date)),
    length(unique(assay$Strain)),
    fixef(model_final)[1],
    summary(model_final)$coefficients[1, 2],
    fixef(model_final)[2],
    summary(model_final)$coefficients[2, 2],
    summary(model_final)$coefficients[2, 3],
    tryCatch(anova_results["Strain", "Pr(>F)"], error = function(e) NA),
    cohens_d,
    icc,
    var_comp[1, "vcov"],
    var_comp[2, "vcov"],
    var_comp[1, "vcov"] + var_comp[2, "vcov"],
    AIC(model_final),
    BIC(model_final),
    as.numeric(logLik(model_final)),
    REMLcrit(model_final)
  )
)
write.csv(master_table, "Master_Results_Table.csv", row.names = FALSE)
