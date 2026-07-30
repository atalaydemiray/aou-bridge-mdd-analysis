# Statistical analysis template
#
# This file contains planned main and supplementary analyses. Review the
# covariate vectors before treating any output as final.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")
analytic_data <- load_component("07_analytic_data.rds")

outcome <- "mdd_case"
exposures <- c(
  "prs_trans_ancestry", "pss10_score", "eds9_score", "ace_score"
)
core_covariates <- c(
  "age_years", "sex_at_birth", paste0("PC", 1:10)
)

# Income and education are not included in the current model template.

clinical_covariates <- "chronic_disease_count"

# Standardized continuous exposures make coefficients comparable.
for (variable in exposures) {
  analytic_data[[paste0(variable, "_z")]] <- as.numeric(scale(
    analytic_data[[variable]]
  ))
}

# Modified Poisson regression estimates adjusted prevalence ratios.
fit_prevalence_ratio <- function(data, formula) {
  variables <- all.vars(formula)
  model_data <- data[stats::complete.cases(data[variables]), variables]
  model <- stats::glm(
    formula,
    data = model_data,
    family = stats::poisson(link = "log")
  )
  robust <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC0")
  )
  list(model = model, robust = robust, n = nrow(model_data))
}

make_formula <- function(terms) {
  stats::reformulate(terms, response = outcome)
}

# Main result 1: each exposure in its own core-adjusted model.
independent_models <- lapply(paste0(exposures, "_z"), function(exposure) {
  fit_prevalence_ratio(
    analytic_data,
    make_formula(c(exposure, core_covariates))
  )
})
names(independent_models) <- exposures

# Main result 2: all exposures in one joint main-effects model.
joint_model <- fit_prevalence_ratio(
  analytic_data,
  make_formula(c(paste0(exposures, "_z"), core_covariates))
)

# Main result 3: one PRS-by-psychosocial interaction at a time.
interaction_models <- lapply(
  paste0(c("pss10_score", "eds9_score", "ace_score"), "_z"),
  function(exposure) {
    interaction <- paste0("prs_trans_ancestry_z * ", exposure)
    fit_prevalence_ratio(
      analytic_data,
      make_formula(c(interaction, core_covariates))
    )
  }
)
names(interaction_models) <- c("PSS", "EDS", "ACE")

# Supplementary extended adjustment adds the 0-5 chronic disease count as one
# covariate. The five component conditions are not adjusted for separately.
extended_joint_model <- fit_prevalence_ratio(
  analytic_data,
  make_formula(c(
    paste0(exposures, "_z"),
    core_covariates,
    clinical_covariates
  ))
)

model_results <- list(
  independent = independent_models,
  joint = joint_model,
  interactions = interaction_models,
  extended_joint = extended_joint_model
)
save_component(model_results, "08_model_results_template.rds")

message(
  "Template models completed. Add final manuscript tables, probability-scale ",
  "interaction contrasts, ancestry-stratified models, and approved ",
  "sensitivity analyses after the model specification is finalized."
)
