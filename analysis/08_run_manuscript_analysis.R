# Run the prespecified modified-Poisson manuscript models.
#
# The pooled analysis uses the diverse-population GWAS MDD PRS. Genomic-
# ancestry-stratified analyses use the European GWAS PRS for the EUR stratum
# and the diverse-population GWAS PRS for the AFR and AMR strata. Each PRS is
# standardized after the relevant final analytic sample or ancestry stratum is
# selected. PSS-10, EDS-9, and ACE-11 are standardized in the same population.
#
# Primary adjustment includes age, sex at birth, and PC1-PC10. An additional
# joint model adds the 0-5 chronic disease count. Sex-at-birth responses other
# than Female or Male remain missing and are handled by each model's complete-
# case selection. Income and education are not included in the current
# specification.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")
analytic_data <- load_component("07_analytic_data.rds")

outcome <- "mdd_case"
psychosocial_raw <- c(
  pss10_z = "pss10_score",
  eds9_z = "eds9_score",
  ace_z = "ace_score"
)
primary_covariates <- c(
  "age_years",
  "sex_at_birth",
  paste0("PC", 1:10)
)
additional_clinical_covariate <- "chronic_disease_count"

required_columns <- c(
  outcome,
  "prs_mdd_div_adjusted",
  "prs_mdd_eur_adjusted",
  "genomic_ancestry",
  unname(psychosocial_raw),
  primary_covariates
)
missing_columns <- setdiff(required_columns, names(analytic_data))
if (length(missing_columns) > 0L) {
  stop(
    "The analytic data are missing model fields: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

analytic_data[[outcome]] <- as.integer(analytic_data[[outcome]])
if (anyNA(analytic_data[[outcome]]) ||
    any(!analytic_data[[outcome]] %in% c(0L, 1L))) {
  stop("mdd_case must be coded as nonmissing 0/1.", call. = FALSE)
}

analytic_data$sex_at_birth[
  !analytic_data$sex_at_birth %in% c("Female", "Male")
] <- NA_character_
analytic_data$sex_at_birth <- factor(
  analytic_data$sex_at_birth,
  levels = c("Female", "Male")
)

scale_with_missing <- function(x, variable, population) {
  x <- suppressWarnings(as.numeric(x))
  available <- is.finite(x)
  if (sum(available) < 2L) {
    stop(
      variable, " has fewer than two finite values in ", population, ".",
      call. = FALSE
    )
  }
  center <- mean(x[available])
  spread <- stats::sd(x[available])
  if (!is.finite(spread) || spread <= 0) {
    stop(variable, " has no usable variation in ", population, ".", call. = FALSE)
  }
  standardized <- rep(NA_real_, length(x))
  standardized[available] <- (x[available] - center) / spread
  list(
    value = standardized,
    parameter = data.frame(
      population = population,
      variable = variable,
      n_available = sum(available),
      center = center,
      scale = spread
    )
  )
}

prepare_population <- function(data, population, prs_column) {
  if (!prs_column %in% names(data)) {
    stop("Missing PRS column: ", prs_column, call. = FALSE)
  }
  mapping <- c(prs_z = prs_column, psychosocial_raw)
  parameters <- vector("list", length(mapping))
  for (i in seq_along(mapping)) {
    result <- scale_with_missing(data[[mapping[[i]]]], mapping[[i]], population)
    data[[names(mapping)[i]]] <- result$value
    parameters[[i]] <- result$parameter
  }
  list(data = data, parameters = do.call(rbind, parameters))
}

make_formula <- function(terms) {
  stats::reformulate(terms, response = outcome)
}

extract_robust_estimates <- function(
  robust, model_id, population, target_terms
) {
  estimates <- data.frame(
    population = population,
    model_id = model_id,
    term = rownames(robust),
    log_prevalence_ratio = robust[, 1L],
    robust_standard_error = robust[, 2L],
    statistic = robust[, 3L],
    p_value = robust[, 4L],
    row.names = NULL,
    check.names = FALSE
  )
  estimates$prevalence_ratio <- exp(estimates$log_prevalence_ratio)
  estimates$confidence_interval_lower <- exp(
    estimates$log_prevalence_ratio - 1.96 * estimates$robust_standard_error
  )
  estimates$confidence_interval_upper <- exp(
    estimates$log_prevalence_ratio + 1.96 * estimates$robust_standard_error
  )
  estimates$report_target <- estimates$term %in% target_terms
  estimates
}

fit_prevalence_ratio <- function(
  data, formula, model_id, population, target_terms
) {
  variables <- all.vars(formula)
  model_data <- data[stats::complete.cases(data[variables]), variables, drop = FALSE]
  model_data$sex_at_birth <- droplevels(model_data$sex_at_birth)
  case_n <- sum(model_data[[outcome]] == 1L)
  control_n <- sum(model_data[[outcome]] == 0L)
  if (nrow(model_data) == 0L || case_n < 20L || control_n < 20L) {
    stop(
      model_id, " in ", population,
      " does not meet the minimum case/control result threshold.",
      call. = FALSE
    )
  }
  if (nlevels(model_data$sex_at_birth) != 2L) {
    stop(
      model_id, " in ", population,
      " does not contain both modeled sex-at-birth categories.",
      call. = FALSE
    )
  }

  model <- stats::glm(
    formula,
    data = model_data,
    family = stats::poisson(link = "log")
  )
  if (!isTRUE(model$converged) || anyNA(stats::coef(model))) {
    stop(model_id, " in ", population, " did not converge.", call. = FALSE)
  }
  robust <- lmtest::coeftest(
    model,
    vcov. = sandwich::vcovHC(model, type = "HC0")
  )
  if (any(!is.finite(robust))) {
    stop(
      model_id, " in ", population,
      " produced a nonfinite robust estimate.",
      call. = FALSE
    )
  }

  samples <- data.frame(
    population = population,
    model_id = model_id,
    total_n = nrow(model_data),
    case_n = case_n,
    control_n = control_n,
    maximum_fitted_prevalence = max(stats::fitted(model)),
    pearson_dispersion = sum(stats::residuals(model, type = "pearson")^2) /
      stats::df.residual(model)
  )
  estimates <- extract_robust_estimates(
    robust, model_id, population, target_terms
  )
  list(model = model, robust = robust, estimates = estimates, samples = samples)
}

run_model_suite <- function(data, population, prs_column) {
  prepared <- prepare_population(data, population, prs_column)
  model_data <- prepared$data
  models <- list()

  independent_targets <- c(
    prs = "prs_z", pss = "pss10_z", eds = "eds9_z", ace = "ace_z"
  )
  for (name in names(independent_targets)) {
    target <- independent_targets[[name]]
    model_id <- paste0("independent_", name)
    models[[model_id]] <- fit_prevalence_ratio(
      model_data,
      make_formula(c(target, primary_covariates)),
      model_id,
      population,
      target
    )
  }

  joint_exposures <- c("prs_z", "pss10_z", "eds9_z", "ace_z")
  models$joint_main_effects <- fit_prevalence_ratio(
    model_data,
    make_formula(c(joint_exposures, primary_covariates)),
    "joint_main_effects",
    population,
    joint_exposures
  )

  # Additional clinical adjustment requested by the group. This is not part
  # of the primary covariate set.
  models$joint_plus_chronic_count <- fit_prevalence_ratio(
    model_data,
    make_formula(c(
      joint_exposures,
      primary_covariates,
      additional_clinical_covariate
    )),
    "joint_plus_chronic_count",
    population,
    joint_exposures
  )

  interaction_targets <- c(
    pss = "pss10_z", eds = "eds9_z", ace = "ace_z"
  )
  for (name in names(interaction_targets)) {
    psychosocial <- interaction_targets[[name]]
    other_psychosocial <- setdiff(
      c("pss10_z", "eds9_z", "ace_z"), psychosocial
    )
    interaction_term <- paste0("prs_z:", psychosocial)
    model_id <- paste0("interaction_", name)
    models[[model_id]] <- fit_prevalence_ratio(
      model_data,
      make_formula(c(
        paste0("prs_z * ", psychosocial),
        other_psychosocial,
        primary_covariates
      )),
      model_id,
      population,
      interaction_term
    )
  }

  list(
    models = models,
    estimates = do.call(rbind, lapply(models, `[[`, "estimates")),
    samples = do.call(rbind, lapply(models, `[[`, "samples")),
    standardization = prepared$parameters
  )
}

# The pooled primary analysis includes all available genomic ancestries and
# uses the diverse-population GWAS score.
analysis_results <- list(
  pooled = run_model_suite(
    analytic_data,
    "Pooled",
    "prs_mdd_div_adjusted"
  )
)

# Prespecified strata and PRS source. Other ancestry groups remain in the
# pooled analysis but are not modeled separately in this initial set.
ancestry_specification <- data.frame(
  genomic_ancestry = c("EUR", "AFR", "AMR"),
  prs_column = c(
    "prs_mdd_eur_adjusted",
    "prs_mdd_div_adjusted",
    "prs_mdd_div_adjusted"
  )
)

for (i in seq_len(nrow(ancestry_specification))) {
  ancestry <- ancestry_specification$genomic_ancestry[i]
  stratum <- analytic_data[
    analytic_data$genomic_ancestry == ancestry,
    ,
    drop = FALSE
  ]
  analysis_results[[ancestry]] <- run_model_suite(
    stratum,
    ancestry,
    ancestry_specification$prs_column[i]
  )
}

model_estimates <- do.call(
  rbind, lapply(analysis_results, `[[`, "estimates")
)
model_samples <- do.call(
  rbind, lapply(analysis_results, `[[`, "samples")
)
standardization_parameters <- do.call(
  rbind, lapply(analysis_results, `[[`, "standardization")
)

if (any(model_samples$maximum_fitted_prevalence > 1)) {
  warning(
    "At least one modified-Poisson model has fitted values above 1. ",
    "Review 08_model_samples.csv before interpreting that specification."
  )
}

model_objects <- lapply(analysis_results, `[[`, "models")
save_component(model_objects, "08_model_results.rds")
save_component(model_estimates, "08_model_estimates.rds")
save_component(model_samples, "08_model_samples.rds")
save_component(
  standardization_parameters,
  "08_standardization_parameters.rds"
)

# These aggregate tables remain in the authorized Workbench output directory.
# They contain no participant identifiers. Review disclosure requirements
# before downloading or sharing any result.
data.table::fwrite(
  model_estimates[model_estimates$report_target, ],
  file.path(WORK_DIR, "08_primary_estimates.csv")
)
data.table::fwrite(
  model_samples,
  file.path(WORK_DIR, "08_model_samples.csv")
)
data.table::fwrite(
  standardization_parameters,
  file.path(WORK_DIR, "08_standardization_parameters.csv")
)

message(
  "Manuscript models completed for the pooled cohort and EUR, AFR, and AMR ",
  "strata. Review sample sizes, fitted prevalence, and robust estimates in ",
  WORK_DIR, "."
)
