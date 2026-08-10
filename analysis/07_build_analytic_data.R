# Combine all study components into the genomic-QC-eligible analytic cohort.
#
# The MDD cohort is the initial row anchor. Participants without an approved
# genomic record, participants who failed genomic quality control, and related
# participants selected for removal are then excluded. Demographic, survey,
# and clinical components are left-joined by person_id, so survey nonresponse
# remains missing rather than becoming zero.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")

mdd <- load_component("01_mdd_cohort.rds")
demographics <- load_component("02_demographics.rds")
pss_eds <- load_component("03_pss_eds.rds")
ace <- load_component("04_ace.rds")
chronic <- load_component("05_chronic_disease_count.rds")
genomics <- load_component("06_genomics.rds")

# Character keys avoid accidental integer/integer64 mismatches across files.
components <- list(mdd, demographics, pss_eds, ace, chronic, genomics)
components <- lapply(components, function(x) {
  x$person_id <- as.character(x$person_id)
  x
})
mdd <- components[[1L]]
demographics <- components[[2L]]
pss_eds <- components[[3L]]
ace <- components[[4L]]
chronic <- components[[5L]]
genomics <- components[[6L]]

required_chronic_columns <- c("person_id", "chronic_disease_count")
if (!all(required_chronic_columns %in% names(chronic))) {
  stop(
    "The chronic component must contain person_id and chronic_disease_count.",
    call. = FALSE
  )
}
chronic <- chronic[required_chronic_columns]

# Script 06 contains only participants who passed genomic QC and relatedness
# filters and who have both PRS values, ancestry, and PC1-PC10. Matching this
# component is therefore an explicit analytic-cohort eligibility criterion.
mdd_with_genomics <- join_component(mdd, genomics, "genomics")
genomic_eligible <- !is.na(mdd_with_genomics$prs_mdd_div_adjusted) &
  !is.na(mdd_with_genomics$prs_mdd_eur_adjusted)

if (!any(genomic_eligible)) {
  stop(
    "No MDD cohort participant passed the genomic eligibility merge. ",
    "Check the person_id key before proceeding.",
    call. = FALSE
  )
}

analytic_data <- mdd_with_genomics[genomic_eligible, , drop = FALSE]
analytic_data <- join_component(analytic_data, demographics, "demographics")
analytic_data <- join_component(analytic_data, pss_eds, "PSS/EDS")
analytic_data <- join_component(analytic_data, ace, "ACE")
analytic_data <- join_component(
  analytic_data, chronic, "chronic disease count"
)
assert_unique_person(analytic_data, "analytic data")

if (any(!is.finite(analytic_data$prs_mdd_div_adjusted)) ||
    any(!is.finite(analytic_data$prs_mdd_eur_adjusted)) ||
    anyNA(analytic_data[paste0("PC", 1:10)])) {
  stop(
    "The analytic cohort contains incomplete genomic model variables.",
    call. = FALSE
  )
}

count_cases <- function(x) sum(x$mdd_case == 1L)
count_controls <- function(x) sum(x$mdd_case == 0L)
analytic_cohort_flow <- data.frame(
  step = c(
    "EHR/WGS/survey-eligible MDD case-control cohort",
    "After genomic QC, relatedness, and PRS availability"
  ),
  total_n = c(nrow(mdd), nrow(analytic_data)),
  case_n = c(count_cases(mdd), count_cases(analytic_data)),
  control_n = c(count_controls(mdd), count_controls(analytic_data))
)

save_component(analytic_data, "07_analytic_data.rds")
save_component(analytic_cohort_flow, "07_analytic_cohort_flow.rds")

# Record only reproducibility metadata; participant data remain in Workbench.
git_commit <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
  error = function(e) NA_character_
)
run_metadata <- list(
  run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  cdr = CDR,
  git_commit = git_commit,
  r_version = R.version.string,
  session_info = utils::sessionInfo()
)
save_component(run_metadata, "07_run_metadata.rds")

message(
  "Analytic cohort built: ", nrow(analytic_data), " participants; ",
  count_cases(analytic_data), " MDD cases and ",
  count_controls(analytic_data), " controls after genomic exclusions."
)
