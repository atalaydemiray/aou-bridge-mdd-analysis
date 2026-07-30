# Combine all study components into one analytic data frame.
#
# The MDD cohort is the row anchor. All components are left-joined by
# person_id, so survey nonresponse remains missing rather than becoming zero.
# WGS eligibility is already enforced in 01_build_mdd_cohort.R. The approved
# genomic component adds PRS, ancestry, and PCs; it does not redefine the
# cohort. Missing values remain available for model-specific complete-case
# selection and reporting.

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
mdd <- components[[1]]
demographics <- components[[2]]
pss_eds <- components[[3]]
ace <- components[[4]]
chronic <- components[[5]]
genomics <- components[[6]]

# The clinical adjustment uses one 0-5 chronic disease count.
required_chronic_columns <- c("person_id", "chronic_disease_count")
if (!all(required_chronic_columns %in% names(chronic))) {
  stop(
    "The chronic component must contain person_id and chronic_disease_count.",
    call. = FALSE
  )
}
chronic <- chronic[required_chronic_columns]

analytic_data <- mdd
analytic_data <- join_component(analytic_data, demographics, "demographics")
analytic_data <- join_component(analytic_data, pss_eds, "PSS/EDS")
analytic_data <- join_component(analytic_data, ace, "ACE")
analytic_data <- join_component(
  analytic_data, chronic, "chronic disease count"
)
analytic_data <- join_component(analytic_data, genomics, "genomics")

assert_unique_person(analytic_data, "analytic data")

# Confirm how much of the fixed MDD cohort received the imported genomic
# component. Script 07 does not silently restrict the cohort to matched rows.
genomic_match_n <- sum(analytic_data$person_id %in% genomics$person_id)
prs_available_n <- sum(!is.na(analytic_data$prs_trans_ancestry))
if (genomic_match_n == 0L) {
  stop(
    "No MDD cohort participant matched the imported genomic component. ",
    "Check the ID bridge before proceeding.",
    call. = FALSE
  )
}

save_component(analytic_data, "07_analytic_data.rds")

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
  "Analytic data built with ", nrow(analytic_data), " participants; ",
  genomic_match_n, " matched the genomic component and ",
  prs_available_n, " had a nonmissing trans-ancestry PRS."
)
