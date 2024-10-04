# The datasets of main subgroups in this study include:
# 1. ukb_validation_male5059
# 2. ukb_validation_male6069
# 3. ukb_validation_female5059
# 4. ukb_validation_female6069
# 5. es_all_new_male5059
# 6. es_all_new_male6069
# 7. es_all_new_female5059
# 8. es_all_new_female6069

# Discrimination and reclassification analyses of all-cause mortality, using the subgroup dataset 1. ukb_validation_male5059 as an example.

# Load necessary RData file
load("~/Allcause_Mortality/allcausemortality.RData") 

# Load necessary libraries
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(dplyr)
library(survival)
library(Hmisc)
library(ggplot2)

# Define Cox regression and C-index calculation function
cox_regression_and_cindex <- function(data, covariate_indices, metabolite_indices, outcome_column, time_column) {
  # Extract covariates and metabolites
  covariates <- names(data)[covariate_indices]
  metabolites <- names(data)[metabolite_indices]
  
  # Build formulas
  formula_covariates <- as.formula(paste("Surv(", time_column, ",", outcome_column, ") ~", paste(covariates, collapse = " + ")))
  formula_full <- as.formula(paste("Surv(", time_column, ",", outcome_column, ") ~", paste(c(covariates, metabolites), collapse = " + ")))
  
  # Cox proportional hazards models
  model_covariates <- coxph(formula_covariates, data = data)
  model_full <- coxph(formula_full, data = data)
  
  # Calculate linear predictors (risk scores)
  lp_covariates <- predict(model_covariates, newdata = data, type = "lp")
  lp_full <- predict(model_full, newdata = data, type = "lp")
  
  # Calculate C-index for covariates-only model and full model
  cindex_covariates <- rcorr.cens(lp_covariates, Surv(data[[time_column]], data[[outcome_column]]))["C Index"]
  cindex_full <- rcorr.cens(lp_full, Surv(data[[time_column]], data[[outcome_column]]))["C Index"]
  
  # Plot ROC-like curves
  roc_data <- data.frame(
    TPR = c(cindex_covariates, cindex_full),
    FPR = c(1 - cindex_covariates, 1 - cindex_full),
    Model = factor(c("Covariates-only Model, C-index", "Full Model, C-index"))
  )
  
  roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(size = 0.5) +
    geom_abline(linetype = "dashed", color = "black") +
    labs(title = "", x = "False Positive Ratio", y = "True Positive Ratio") + 
    theme_minimal() +
    scale_color_manual(values = c("red", "black")) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black"),
      legend.position = c(0.57, 0.10),
      axis.text = element_text(size = 10),          
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    )
  
  print(roc_plot)
  
  return(list(
    cindex_covariates = cindex_covariates,
    cindex_full = cindex_full
  ))
}

# Set indices for covariates and metabolic risk score (MRS)
covariate_indices <- c(2:14)  
metabolite_indices <- c(260)  # Metabolic risk score (MRS)

# Set outcome column and time column names
outcome_column_name <- "mortality10y"  
time_column_name <- "followup_year"

# Example: Analyze for the `ukb_validation_male5059` dataset
result <- cox_regression_and_cindex(ukb_validation_male5059, covariate_indices, metabolite_indices, outcome_column_name, time_column_name)

# Print C-index values
print(paste("C-index (Covariates):", result$cindex_covariates))
print(paste("C-index (Full Model):", result$cindex_full))


# Load necessary packages for IDI and NRI calculation
if (!requireNamespace("PredictABEL", quietly = TRUE)) {
  install.packages("PredictABEL")
}
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}

library(survival)
library(PredictABEL)

# Function to calculate IDI and NRI based on Cox models
calculate_idi_nri <- function(data_frame, covariate_indices, metabolite_indices, outcome_column_name, time_column_name) {
  # Subset covariates and metabolite (MRS)
  covariates <- data_frame[, covariate_indices]
  metabolites <- data_frame[, metabolite_indices]
  all_covariates <- data_frame[, c(covariate_indices, metabolite_indices)]

  # Ensure the outcome and followup columns are numeric
  data_frame[[time_column_name]] <- as.numeric(as.character(data_frame[[time_column_name]]))
  data_frame[[outcome_column_name]] <- as.numeric(as.character(data_frame[[outcome_column_name]]))

  # Baseline Cox model using only covariates
  baseline_model <- coxph(Surv(data_frame[[time_column_name]], data_frame[[outcome_column_name]]) ~ ., data = covariates)

  # Extended Cox model including both covariates and metabolites (MRS)
  extended_model <- coxph(Surv(data_frame[[time_column_name]], data_frame[[outcome_column_name]]) ~ ., data = all_covariates)

  # Get linear predictors (risk scores) for both models
  baseline_lp <- predict(baseline_model, newdata = data_frame, type = "lp")
  extended_lp <- predict(extended_model, newdata = data_frame, type = "lp")

  # Calculate survival probabilities from the linear predictors
  baseline_prob <- 1 - exp(-baseline_lp)  # Approximate survival probability for baseline model
  extended_prob <- 1 - exp(-extended_lp)  # Approximate survival probability for extended model

  # Prepare data for reclassification
  event <- data_frame[[outcome_column_name]]
  data_for_reclassification <- data.frame(event, baseline_prob, extended_prob)

  # Set cut-off values for risk categories
  cutoffs <- c(0, 0.125, 0.25, 0.5, 1)  # Adjust cutoffs based on context

  # Calculate IDI and NRI using PredictABEL's reclassification function
  idi_nri_result <- reclassification(
    data = data_for_reclassification, 
    cOutcome = 1,  # Column index of the outcome variable
    predrisk1 = baseline_prob, 
    predrisk2 = extended_prob, 
    cutoff = cutoffs
  )

  return(idi_nri_result)
}

# Define covariates and metabolic risk score (MRS)
covariate_indices <- c(2:14)  # Example covariates
metabolite_indices <- c(260)  # Metabolic risk score (MRS)

# Define outcome and time column names
outcome_column_name <- "mortality10y"
time_column_name <- "followup_year"

# Example: Analyze IDI and NRI for the `ukb_validation_male5059` dataset
idi_nri_result <- calculate_idi_nri(ukb_validation_male5059, covariate_indices, metabolite_indices, outcome_column_name, time_column_name)

# Print IDI and NRI results
print(idi_nri_result)
