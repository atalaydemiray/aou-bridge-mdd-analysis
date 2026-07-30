# Calculate age and sex-at-birth covariates.
#
# Age is calculated from year of birth using the study reference year.
# Unrecorded sex-at-birth values remain missing.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")

sql <- paste0(
"SELECT
  p.person_id,
  EXTRACT(YEAR FROM DATE '2026-07-29') - p.year_of_birth AS age_years,
  sex.concept_name AS sex_at_birth
FROM ", bq_table("person"), " p
LEFT JOIN ", bq_table("concept"), " sex
  ON p.sex_at_birth_concept_id = sex.concept_id
WHERE p.year_of_birth IS NOT NULL"
)

demographics <- run_query(sql)
assert_unique_person(demographics, "demographics")

# Keep the recorded label. The model script will treat unrecorded values as
# missing rather than assuming Female or Male.
demographics$sex_at_birth[
  is.na(demographics$sex_at_birth) |
    demographics$sex_at_birth %in% c("", "No matching concept")
] <- NA_character_

save_component(demographics, "02_demographics.rds")
message("Demographic covariates built.")
