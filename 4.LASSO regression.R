# Load necessary RData file
load("~/Allcause_Mortality/allcausemortality.RData") 

# Set seed for reproducibility
set.seed(234)

# Generate random indices for splitting data into training and validation sets
train_index <- sample(seq_len(nrow(ukb_all_new_imputed_logall_z)), size = 0.7 * nrow(ukb_all_new_imputed_logall_z))

# Split data into training set and validation set
ukb_train_set <- ukb_all_new_imputed_logall_z[train_index, ]
ukb_validation_set <- ukb_all_new_imputed_logall_z[-train_index, ]

# Print first few rows of training and validation sets
print(head(ukb_train_set))
print(head(ukb_validation_set))

# Extract 4 subgroups from training set
ukb_train_male5059 <- subset(ukb_train_set, sex == 1 & (age >= 50 & age <= 59))
ukb_train_male6069 <- subset(ukb_train_set, sex == 1 & (age >= 60 & age <= 69))
ukb_train_female5059 <- subset(ukb_train_set, sex == 0 & (age >= 50 & age <= 59))
ukb_train_female6069 <- subset(ukb_train_set, sex == 0 & (age >= 60 & age <= 69))

# Extract 4 subgroups from internal validation set
ukb_validation_male5059 <- subset(ukb_validation_set, sex == 1 & (age >= 50 & age <= 59))
ukb_validation_male6069 <- subset(ukb_validation_set, sex == 1 & (age >= 60 & age <= 69))
ukb_validation_female5059 <- subset(ukb_validation_set, sex == 0 & (age >= 50 & age <= 59))
ukb_validation_female6069 <- subset(ukb_validation_set, sex == 0 & (age >= 60 & age <= 69))

# Extract 4 subgroups from external validation set
es_all_new_male5059 <- subset(es_all_new_imputed_logall_z, sex == 2 & (age >= 50 & age <= 59))
es_all_new_male6069 <- subset(es_all_new_imputed_logall_z, sex == 2 & (age >= 60 & age <= 69))
es_all_new_female5059 <- subset(es_all_new_imputed_logall_z, sex == 1 & (age >= 50 & age <= 59))
es_all_new_female6069 <- subset(es_all_new_imputed_logall_z, sex == 1 & (age >= 60 & age <= 69))

# Check the number of rows for each subgroup
print(nrow(ukb_train_male5059))
print(nrow(ukb_train_male6069))
print(nrow(ukb_train_female5059))
print(nrow(ukb_train_female6069))

print(nrow(ukb_validation_male5059))
print(nrow(ukb_validation_male6069))
print(nrow(ukb_validation_female5059))
print(nrow(ukb_validation_female6069))

print(nrow(es_all_new_male5059))
print(nrow(es_all_new_male6069))
print(nrow(es_all_new_female5059))
print(nrow(es_all_new_female6069))


# Define a function to perform Bootstrap LASSO analysis for a given dataset
perform_bootstrap_lasso <- function(data, group_name) {
  # Set the dependent variable (y)
  data$followup_year <- round(data$followup_year)
  y <- Surv(data$followup_year, data$mortality10y)
  
  # Set the independent variables (x)
  x <- data[, 15:259]  # Columns 15 to 259 are the metabolic biomarkers
  
  # Combine independent variables and dependent variables into a data frame
  analysis_data <- as.data.frame(cbind(x, followup_year = data$followup_year, mortality10y = data$mortality10y))
  
  # Define the Bootstrap LASSO function
  bootstrap_lasso <- function(data, indices) {
    d <- data[indices, ]
    x_boot <- as.matrix(d[, -c(ncol(d)-1, ncol(d))])  # Exclude followup_year and mortality10y
    y_boot <- Surv(d$followup_year, d$mortality10y)
    
    cv_lasso <- cv.glmnet(x_boot, y_boot, alpha = 1, family = "cox", nfolds = 10)
    best_lambda <- cv_lasso$lambda.1se
    lasso_model <- glmnet(x_boot, y_boot, alpha = 1, family = "cox", lambda = best_lambda)
    
    return(as.vector(coef(lasso_model)) != 0)
  }
  
  # Perform 1000 Bootstrap iterations
  library(boot)
  set.seed(234)
  results <- boot(analysis_data, statistic = bootstrap_lasso, R = 1000)
  
  # Count how many times each metabolite was selected
  selected_counts <- apply(results$t, 2, sum)
  
  # Select metabolites that were chosen in at least 95% of the Bootstrap samples
  selected_metabolites <- colnames(x)[selected_counts >= 0.95 * 1000]
  
  # Print the selected metabolites
  print(paste("Selected metabolites for", group_name, "(chosen in at least 95% of the bootstrap samples):"))
  print(selected_metabolites)
  
  # Save the selected metabolites to a CSV file
  output_filename <- paste0("bootstrap_selected_metabolites_", group_name, ".csv")
  write.csv(selected_metabolites, output_filename, row.names = FALSE)
}

# Apply the function to all four subgroups

# For male aged 50-59
perform_bootstrap_lasso(ukb_train_male5059, "male5059")

# For male aged 60-69
perform_bootstrap_lasso(ukb_train_male6069, "male6069")

# For female aged 50-59
perform_bootstrap_lasso(ukb_train_female5059, "female5059")

# For female aged 60-69
perform_bootstrap_lasso(ukb_train_female6069, "female6069")


# Cox regression analysis and forest plots
# Load necessary R packages
library(survival)
library(ggplot2)
library(dplyr)
library(forestplot)

# Define a function to perform COX regression analysis and return the results
cox_analysis <- function(data, biomarkers, covariates, outcome, followup) {
  results <- data.frame()
  data <- data %>% mutate(Surv_obj = Surv(get(followup), get(outcome)))
  for (biomarker in biomarkers) {
    formula <- as.formula(paste("Surv_obj ~", biomarker, "+", paste(covariates, collapse = "+")))
    cox_model <- coxph(formula, data = data)
    summary_cox <- summary(cox_model)
    hr <- summary_cox$coefficients[biomarker, "exp(coef)"]
    ci_lower <- summary_cox$conf.int[biomarker, "lower .95"]
    ci_upper <- summary_cox$conf.int[biomarker, "upper .95"]
    p_value <- summary_cox$coefficients[biomarker, "Pr(>|z|)"]
    schoenfeld_test <- cox.zph(cox_model)
    schoenfeld_p <- schoenfeld_test$table[biomarker, "p"]
    
    results <- rbind(results, data.frame(Biomarker = biomarker, HR = hr, CI_Lower = ci_lower, CI_Upper = ci_upper, P_Value = p_value, Schoenfeld_P = schoenfeld_p))
  }
  return(results)
}

# Define a function to check if necessary columns exist in the dataset
check_columns <- function(data, outcomes, followup) {
  all_columns <- c(outcomes, followup)
  if (!all(all_columns %in% colnames(data))) {
    stop("Required outcome or follow-up columns are missing from the dataset")
  }
}

# Combine the results and export them to a CSV file
combine_and_export <- function(results_train, results_validation, results_external, filename) {
  results_train$Dataset <- "Training"
  results_validation$Dataset <- "Validation"
  results_external$Dataset <- "External"
  all_results <- rbind(results_train, results_validation, results_external)
  write.csv(all_results, file = filename, row.names = FALSE)
}

# Load data and plot forest plot
plot_forest_from_csv <- function(filename, biomarkers, x_labels, outcome_name, tiff_filename = NULL) {
  results <- read.csv(filename)
  
  results <- results %>%
    mutate(Biomarker = factor(Biomarker, levels = biomarkers),
           Dataset = factor(Dataset, levels = c("Training", "Validation", "External"))) %>%
    arrange(Biomarker, Dataset)
  
  background_data <- results %>%
    mutate(index = row_number()) %>%
    mutate(group = (index - 1) %/% 3) %>%
    group_by(group) %>%
    mutate(bg_color = ifelse(group %% 2 == 0, "white", "grey95")) %>%
    ungroup()
  
  p <- ggplot() +
    geom_rect(data = background_data,
              aes(xmin = as.numeric(Biomarker) - 0.5, xmax = as.numeric(Biomarker) + 0.5, 
                  ymin = -Inf, ymax = Inf, fill = bg_color),
              inherit.aes = FALSE, alpha = 0.5) +
    geom_pointrange(data = results, 
                    aes(x = Biomarker, y = HR, ymin = CI_Lower, ymax = CI_Upper, color = Dataset),
                    position = position_dodge(width = 0.5), size = 0.15) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_continuous(limits = c(0.2, 1.8), breaks = seq(0.2, 1.8, by = 0.2)) +
    scale_x_discrete(labels = x_labels) + 
    labs(title = outcome_name, x = "", y = "HR per 1-SD increment") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey95"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 12),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    scale_color_manual(values = c("Training" = "#000000", "Validation" = "#CD3333", "External" = "#1874CD"),
                       labels = c("UK Biobank (70%)", "UK Biobank (30%)", "ESTHER Study")) +
    scale_fill_identity()
  
  if (!is.null(tiff_filename)) {
    ggsave(filename = tiff_filename, plot = p, device = "tiff", width = 7, height = 5, units = "in", dpi = 600)
  } else {
    print(p)
  }
}

# Define biomarker and covariate information for different subgroups
group_info <- list(
  male5059 = list(
    biomarkers = c("HDL_size", "XL_HDL_FC", "Omega_6_by_Omega_3", "GlycA", "bOHbutyrate", "Citrate", 
                   "Acetone", "Acetate", "Glucose", "Lactate", "LA_pct", "Tyr", "L_LDL_CE_pct", 
                   "His", "IDL_CE_pct", "L_VLDL_TG", "Val", "VLDL_size", "Albumin", "PUFA", "S_LDL_CE"),
    x_labels = c("HDL-size", "XL-HDL-FC", "Omega-6/Omega-3", "GlycA", "bOHbutyrate", "Citrate", 
                 "Acetone", "Acetate", "Glucose", "Lactate", "LA-pct", "Tyr", "L-LDL-CE-pct", 
                 "His", "IDL-CE-pct", "L-VLDL-TG", "Val", "VLDL-size", "Albumin", "PUFA", "S-LDL-CE"),
    train_data = ukb_train_male5059,
    validation_data = ukb_validation_male5059,
    external_data = es_all_new_male5059
  ),
  male6069 = list(
    biomarkers = c("XL_HDL_FC", "XL_HDL_PL", "Omega_6_by_Omega_3", "GlycA", "bOHbutyrate", "Acetone", 
                   "Citrate", "Phe", "Acetate", "Glucose", "Tyr", "LA_pct", "L_LDL_CE_pct", "LA", 
                   "His", "Leu", "IDL_CE_pct", "Albumin", "Val", "PUFA_by_MUFA", "Omega_3", 
                   "VLDL_size", "S_HDL_CE", "S_LDL_CE"),
    x_labels = c("XL-HDL-FC", "XL-HDL-PL", "Omega-6/Omega-3", "GlycA", "bOHbutyrate", "Acetone", 
                 "Citrate", "Phe", "Acetate", "Glucose", "Tyr", "LA-pct", "L-LDL-CE-pct", "LA", 
                 "His", "Leu", "IDL-CE-pct", "Albumin", "Val", "PUFA/MUFA", "Omega-3", "VLDL-size", 
                 "S-HDL-CE", "S-LDL-CE"),
    train_data = ukb_train_male6069,
    validation_data = ukb_validation_male6069,
    external_data = es_all_new_male6069
  ),
  # Additional subgroups like female5059 and female6069...
)

# Define covariates, outcome variables, and follow-up time
covariates13 <- c("age", "sex", "bmi", "sbp", "diab", "cvd", "ca", "smoke", "alc", "Total_C", "HDL_C", "HDL_TG", "Creatinine")
outcomes3 <- c("mortality10y", "cvd_mort10y", "ca_mort10y")
followup_year <- "followup_year"

# Process and analyze each subgroup
for (group in names(group_info)) {
  info <- group_info[[group]]
  
  # Check for necessary columns in the datasets
  check_columns(info$train_data, outcomes3, followup_year)
  check_columns(info$validation_data, outcomes3, followup_year)
  check_columns(info$external_data, outcomes3, followup_year)
  
  # Perform the analysis and save the results
  for (i in 1:length(outcomes3)) {
    outcome <- outcomes3[i]
    result_train <- cox_analysis(info$train_data, info$biomarkers, covariates13, outcome, followup_year)
    result_validation <- cox_analysis(info$validation_data, info$biomarkers, covariates13, outcome, followup_year)
    result_external <- cox_analysis(info$external_data, info$biomarkers, covariates13, outcome, followup_year)
    
    filename <- paste0("cox_results_", group, "_", i, ".csv")
    combine_and_export(result_train, result_validation, result_external, filename)
    
    # Plot forest plot
    plot_forest_from_csv(filename, info$biomarkers, info$x_labels, outcome, paste0("forest_plot_", group, "_", outcome, ".tiff"))
  }
}



# Function to calculate metabolic risk score (MRS)
calculate_mrs <- function(data, cox_model, biomarkers) {
  # Extract coefficients from the Cox model
  coefficients <- coef(cox_model)[biomarkers]
  intercept <- as.numeric(cox_model$coefficients[1])   
  # Subset the biomarkers from the dataset
  biomarker_data <- data[, biomarkers]
  
  # Ensure the biomarkers in the dataset match the coefficients
  if (!all(colnames(biomarker_data) %in% names(coefficients))) {
    stop("Some biomarkers in the dataset do not match the model coefficients.")
  }
  
  # Calculate MRS: sum of biomarker values * respective coefficients + intercept
  mrs <- rowSums(biomarker_data * coefficients) + intercept
  
  return(mrs)
}

# Function to calculate and add MRS for Training, Validation, and External datasets
calculate_and_add_mrs <- function(info, model) {
  # Generate MRS for Training, Validation, and External datasets
  info$train_data$MRS <- calculate_mrs(info$train_data, model, info$biomarkers)
  info$validation_data$MRS <- calculate_mrs(info$validation_data, model, info$biomarkers)
  info$external_data$MRS <- calculate_mrs(info$external_data, model, info$biomarkers)
  
  return(info)
}

# Define covariates, biomarkers, and follow-up time for each group
covariates13 <- c("age", "sex", "bmi", "sbp", "diab", "cvd", "ca", "smoke", "alc", "Total_C", "HDL_C", "HDL_TG", "Creatinine")
outcome <- "mortality10y" 
followup_year <- "followup_year"

# Define biomarker and covariate information for different subgroups
group_info <- list(
  male5059 = list(
    biomarkers = c("HDL_size", "XL_HDL_FC", "Omega_6_by_Omega_3", "GlycA", "bOHbutyrate", "Citrate", 
                   "Acetone", "Acetate", "Glucose", "Lactate", "LA_pct", "Tyr", "L_LDL_CE_pct", 
                   "His", "IDL_CE_pct", "L_VLDL_TG", "Val", "VLDL_size", "Albumin", "PUFA", "S_LDL_CE"),
    train_data = ukb_train_male5059,
    validation_data = ukb_validation_male5059,
    external_data = es_all_new_male5059
  ),
  male6069 = list(
    biomarkers = c("XL_HDL_FC", "XL_HDL_PL", "Omega_6_by_Omega_3", "GlycA", "bOHbutyrate", "Acetone", 
                   "Citrate", "Phe", "Acetate", "Glucose", "Tyr", "LA_pct", "L_LDL_CE_pct", "LA", 
                   "His", "Leu", "IDL_CE_pct", "Albumin", "Val", "PUFA_by_MUFA", "Omega_3", 
                   "VLDL_size", "S_HDL_CE", "S_LDL_CE"),
    train_data = ukb_train_male6069,
    validation_data = ukb_validation_male6069,
    external_data = es_all_new_male6069
  ),
  female5059 = list(
    biomarkers = c("GlycA", "Omega_6_by_Omega_3", "S_LDL_PL_pct", "Acetate", "L_VLDL_C_pct", "Acetoacetate", 
                   "Citrate", "Glucose", "L_LDL_CE_pct", "M_LDL_C_pct", "LA_pct", "Leu", "S_LDL_C_pct", 
                   "VLDL_size", "Val", "Albumin", "IDL_CE_pct"),
    train_data = ukb_train_female5059,
    validation_data = ukb_validation_female5059,
    external_data = es_all_new_female5059
  ),
  female6069 = list(
    biomarkers = c("GlycA", "Omega_6_by_Omega_3", "XS_VLDL_PL_pct", "Acetoacetate", "Glucose", "LA_pct", 
                   "L_LDL_CE_pct", "His", "Albumin", "Val", "S_HDL_CE", "IDL_CE_pct", 
                   "Omega_3", "VLDL_size", "PUFA", "M_LDL_CE"),
    train_data = ukb_train_female6069,
    validation_data = ukb_validation_female6069,
    external_data = es_all_new_female6069
  )
)

# Process each subgroup for "mortality10y" outcome
for (group in names(group_info)) {
  info <- group_info[[group]]
  
  # Check necessary columns in the datasets
  check_columns(info$train_data, outcome, followup_year)
  check_columns(info$validation_data, outcome, followup_year)
  check_columns(info$external_data, outcome, followup_year)
  
  # Perform Cox analysis on the Training dataset for "mortality10y"
  cox_model_train <- coxph(Surv(info$train_data[[followup_year]], info$train_data[[outcome]]) ~ .,
                           data = info$train_data[, c(info$biomarkers, covariates13)])
  
  # Calculate and add MRS for Training, Validation, and External datasets
  info <- calculate_and_add_mrs(info, cox_model_train)
    group_info[[group]] <- info
}
