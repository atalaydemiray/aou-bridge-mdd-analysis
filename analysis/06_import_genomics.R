# Import genomic quality-control fields, MDD polygenic scores, ancestry, and
# principal components.
#
# Keep the approved participant-level file inside the authorized Workbench and
# outside this Git repository. The expected one-row-per-participant fields are:
#
#   person_id, prs_mdd_div, prs_mdd_eur, ancestry_pred,
#   PC1 through PC10, flag, related
#
# `flag == 1` identifies a participant who failed genomic quality control.
# `related == 1` identifies a related participant selected for removal. Both
# groups are excluded before the genomic component is joined to the MDD cohort.
#
# The two PRS columns are adjusted across ancestry using PC1-PC5. The supplied
# approach first models each raw PRS mean using the five PCs, then models the
# squared residuals using the same PCs, and divides the mean-model residual by
# the square root of the absolute predicted residual variance. Raw scores are
# retained. All PRS values are standardized later, after the final analytic
# sample or ancestry stratum has been selected.
#
# In the active R Console, provide only the local Workbench path:
#
# Sys.setenv(
#   BRIDGE_MDD_GENOMICS_PATH =
#     "/path/inside/authorized/workbench/prs_mdd_all.txt"
# )
#
# The defaults below match the approved component. Optional environment
# variables allow a renamed column to be mapped without editing this script.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

GENOMICS_PATH <- Sys.getenv("BRIDGE_MDD_GENOMICS_PATH", unset = "")
PERSON_ID_COLUMN <- env_or_default(
  "BRIDGE_MDD_PERSON_ID_COLUMN", "person_id"
)
DIVERSE_PRS_COLUMN <- env_or_default(
  "BRIDGE_MDD_DIVERSE_PRS_COLUMN", "prs_mdd_div"
)
EUROPEAN_PRS_COLUMN <- env_or_default(
  "BRIDGE_MDD_EUROPEAN_PRS_COLUMN", "prs_mdd_eur"
)
ANCESTRY_COLUMN <- env_or_default(
  "BRIDGE_MDD_ANCESTRY_COLUMN", "ancestry_pred"
)
GENOMIC_QC_COLUMN <- env_or_default(
  "BRIDGE_MDD_GENOMIC_QC_COLUMN", "flag"
)
RELATED_COLUMN <- env_or_default(
  "BRIDGE_MDD_RELATED_COLUMN", "related"
)

if (!nzchar(GENOMICS_PATH)) {
  stop(
    "Genomic input is not configured. Set BRIDGE_MDD_GENOMICS_PATH in ",
    "this R Console session, then rerun analysis/06_import_genomics.R. ",
    "See README Step 6.",
    call. = FALSE
  )
}

genomics_raw <- read_component_file(path.expand(GENOMICS_PATH))

required_columns <- c(
  PERSON_ID_COLUMN,
  DIVERSE_PRS_COLUMN,
  EUROPEAN_PRS_COLUMN,
  ANCESTRY_COLUMN,
  paste0("PC", 1:10),
  GENOMIC_QC_COLUMN,
  RELATED_COLUMN
)
missing_columns <- setdiff(required_columns, names(genomics_raw))
if (length(missing_columns) > 0L) {
  stop(
    "The genomic component is missing: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

# Validate the participant key before applying any exclusion.
person_id <- as.character(genomics_raw[[PERSON_ID_COLUMN]])
if (anyNA(person_id) || any(!nzchar(trimws(person_id))) ||
    anyDuplicated(person_id)) {
  stop(
    "The genomic component must contain one nonmissing row per person_id.",
    call. = FALSE
  )
}

as_binary_flag <- function(x, name) {
  value <- suppressWarnings(as.integer(as.character(x)))
  if (anyNA(value) || any(!value %in% c(0L, 1L))) {
    stop(name, " must contain only nonmissing 0/1 values.", call. = FALSE)
  }
  value
}

qc_flag <- as_binary_flag(
  genomics_raw[[GENOMIC_QC_COLUMN]], GENOMIC_QC_COLUMN
)
related_flag <- as_binary_flag(
  genomics_raw[[RELATED_COLUMN]], RELATED_COLUMN
)

genomic_qc_flow <- data.frame(
  step = c(
    "Genomic component received",
    "Passed genomic quality control",
    "Passed genomic quality control and unrelated"
  ),
  n = c(
    nrow(genomics_raw),
    sum(qc_flag == 0L),
    sum(qc_flag == 0L & related_flag == 0L)
  )
)

eligible <- qc_flag == 0L & related_flag == 0L
if (!any(eligible)) {
  stop("No participant passed genomic QC and relatedness filters.", call. = FALSE)
}

genomics_qc <- data.frame(
  person_id = person_id[eligible],
  prs_mdd_div_raw = suppressWarnings(as.numeric(
    genomics_raw[[DIVERSE_PRS_COLUMN]][eligible]
  )),
  prs_mdd_eur_raw = suppressWarnings(as.numeric(
    genomics_raw[[EUROPEAN_PRS_COLUMN]][eligible]
  )),
  genomic_ancestry = toupper(trimws(as.character(
    genomics_raw[[ANCESTRY_COLUMN]][eligible]
  ))),
  genomics_raw[eligible, paste0("PC", 1:10), drop = FALSE],
  check.names = FALSE
)

numeric_columns <- c(
  "prs_mdd_div_raw", "prs_mdd_eur_raw", paste0("PC", 1:10)
)
for (column in numeric_columns) {
  genomics_qc[[column]] <- suppressWarnings(as.numeric(genomics_qc[[column]]))
  if (any(!is.finite(genomics_qc[[column]]))) {
    stop(
      "Genomic QC-passed rows contain missing or nonfinite values in ",
      column, ".",
      call. = FALSE
    )
  }
}
if (anyNA(genomics_qc$genomic_ancestry) ||
    any(!nzchar(genomics_qc$genomic_ancestry))) {
  stop(
    "Genomic QC-passed rows contain missing ancestry assignments.",
    call. = FALSE
  )
}

# Reproduce the supplied cross-ancestry PRS distribution adjustment using
# PC1-PC5. The absolute value in the residual-variance prediction is part of
# the supplied method. The script stops rather than silently retaining an
# undefined or infinite adjusted value.
adjust_pgs_across_ancestry <- function(
  data, pgs_column, pc_columns = paste0("PC", 1:5)
) {
  model_data <- data[c(pgs_column, pc_columns)]
  names(model_data)[1L] <- "PGS"

  mean_model <- stats::lm(PGS ~ ., data = model_data)
  variance_model <- stats::lm(
    stats::residuals(mean_model)^2 ~ .,
    data = model_data[pc_columns]
  )

  predicted_mean <- stats::predict(mean_model, model_data)
  predicted_variance <- stats::predict(variance_model, model_data)
  adjusted <- (model_data$PGS - predicted_mean) /
    sqrt(abs(predicted_variance))

  nan_index <- is.nan(adjusted)
  if (any(nan_index)) {
    adjusted[nan_index] <-
      (model_data$PGS[nan_index] - stats::coef(mean_model)[1L]) /
      sqrt(abs(stats::coef(variance_model)[1L]))
  }
  if (any(!is.finite(adjusted)) || stats::sd(adjusted) <= 0) {
    stop(
      "The cross-ancestry adjustment produced an invalid ", pgs_column,
      " result.",
      call. = FALSE
    )
  }
  as.numeric(adjusted)
}

genomics_qc$prs_mdd_div_adjusted <- adjust_pgs_across_ancestry(
  genomics_qc, "prs_mdd_div_raw"
)
genomics_qc$prs_mdd_eur_adjusted <- adjust_pgs_across_ancestry(
  genomics_qc, "prs_mdd_eur_raw"
)

genomics <- genomics_qc[c(
  "person_id",
  "prs_mdd_div_raw",
  "prs_mdd_eur_raw",
  "prs_mdd_div_adjusted",
  "prs_mdd_eur_adjusted",
  "genomic_ancestry",
  paste0("PC", 1:10)
)]
assert_unique_person(genomics, "genomics")

save_component(genomics, "06_genomics.rds")
save_component(genomic_qc_flow, "06_genomics_qc_flow.rds")
message(
  "Genomic component imported: ", nrow(genomics_raw), " participants; ",
  nrow(genomics), " passed genomic QC and relatedness filters."
)
