# Build every study data component in order.
# Statistical models are intentionally not run by this script. When genomic
# inputs are not configured in the current R session, the core build ends
# successfully after script 05.

core_scripts <- c(
  "analysis/00_preflight.R",
  "analysis/01_build_mdd_cohort.R",
  "analysis/02_build_demographics.R",
  "analysis/03_build_pss_eds.R",
  "analysis/04_build_ace.R",
  "analysis/05_build_chronic_conditions.R"
)

for (script in core_scripts) {
  message("\n--- Running ", script, " ---")
  source(script, local = .GlobalEnv)
}

genomics_configured <- nzchar(Sys.getenv(
  "BRIDGE_MDD_GENOMICS_PATH", unset = ""
))

if (genomics_configured) {
  for (script in c(
    "analysis/06_import_genomics.R",
    "analysis/07_build_analytic_data.R"
  )) {
    message("\n--- Running ", script, " ---")
    source(script, local = .GlobalEnv)
  }
} else {
  message(
    "\nCore build completed through script 05. ",
    "Set BRIDGE_MDD_GENOMICS_PATH in the ",
    "current R session before running scripts 06 and 07. See README Step 6."
  )
}
