# Calculate the five-condition chronic disease count from ICD source codes.

# A condition requires two diagnosis dates separated by at least 30 days.
# The five indicators are summed into one 0-5 count for model adjustment; the
# individual indicators are not entered together as separate model covariates.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")

# Codes are grouped here so the publication definition is visible in one file.
code_groups <- list(
  diabetes = list(
    ICD10CM = c("E08", "E09", "E10", "E11", "E13"),
    ICD9CM = "250"
  ),
  heart_disease = list(
    ICD10CM = c(
      "I20", "I21", "I22", "I23", "I24", "I25",
      "I42", "I43", "I48", "I50"
    ),
    ICD9CM = c(
      "410", "411", "412", "413", "414", "425", "427.3", "428"
    )
  ),
  hypertension = list(
    ICD10CM = c("I10", "I11", "I12", "I13", "I15", "I1A"),
    ICD9CM = c("401", "402", "403", "404", "405")
  ),
  ckd = list(
    ICD10CM = c(
      "E08.2", "E09.2", "E10.2", "E11.2", "E13.2", "I12", "I13",
      "N01", "N02", "N03", "N04", "N05", "N06", "N07", "N08",
      "N14", "N15", "N16", "N18", "N25", "N26", "N99.0", "Q61"
    ),
    ICD9CM = c(
      "250.4", "403", "404", "580", "581", "582", "583",
      "585", "586", "587"
    )
  ),
  chronic_pain = list(
    ICD10CM = c(
      "G89.21", "G89.22", "G89.28", "G89.29", "G89.3", "G89.4",
      "M54.1", "M60", "M79.1", "M79.2", "M79.7", "R53.82"
    ),
    ICD9CM = c("338.2", "338.3", "338.4", "729.1", "729.2", "780.71")
  )
)

rules <- do.call(rbind, lapply(names(code_groups), function(condition) {
  do.call(rbind, lapply(names(code_groups[[condition]]), function(vocabulary) {
    data.frame(
      condition = condition,
      vocabulary_id = vocabulary,
      concept_code_prefix = code_groups[[condition]][[vocabulary]]
    )
  }))
}))

rule_rows <- vapply(seq_len(nrow(rules)), function(i) {
  paste0(
    "STRUCT(", sql_quote(rules$condition[i]), " AS condition, ",
    sql_quote(rules$vocabulary_id[i]), " AS vocabulary_id, ",
    sql_quote(rules$concept_code_prefix[i]), " AS code_prefix)"
  )
}, character(1))

sql <- paste0(
"WITH rules AS (
  SELECT * FROM UNNEST([
    ", paste(rule_rows, collapse = ",\n    "), "
  ])
),
adult_ehr AS (
  SELECT p.person_id
  FROM ", bq_table("person"), " p
  JOIN ", bq_table("cb_search_person"), " cb USING (person_id)
  WHERE cb.has_ehr_data = 1
),
condition_events AS (
  SELECT DISTINCT
    co.person_id,
    r.condition,
    COALESCE(DATE(co.condition_start_datetime), co.condition_start_date)
      AS diagnosis_date
  FROM ", bq_table("condition_occurrence"), " co
  JOIN adult_ehr USING (person_id)
  JOIN ", bq_table("concept"), " source_concept
    ON co.condition_source_concept_id = source_concept.concept_id
  JOIN rules r
    ON UPPER(source_concept.vocabulary_id) = UPPER(r.vocabulary_id)
   AND STARTS_WITH(
     UPPER(source_concept.concept_code),
     UPPER(r.code_prefix)
   )
  WHERE COALESCE(
    DATE(co.condition_start_datetime), co.condition_start_date
  ) IS NOT NULL
),
condition_summary AS (
  SELECT
    person_id,
    condition,
    COUNT(DISTINCT diagnosis_date) >= 2
      AND DATE_DIFF(MAX(diagnosis_date), MIN(diagnosis_date), DAY) >= 30
      AS condition_present
  FROM condition_events
  GROUP BY person_id, condition
),
person_flags AS (
  SELECT
    person_id,
    LOGICAL_OR(condition = 'diabetes' AND condition_present) AS diabetes,
    LOGICAL_OR(condition = 'heart_disease' AND condition_present)
      AS heart_disease,
    LOGICAL_OR(condition = 'hypertension' AND condition_present)
      AS hypertension,
    LOGICAL_OR(condition = 'ckd' AND condition_present) AS ckd,
    LOGICAL_OR(condition = 'chronic_pain' AND condition_present)
      AS chronic_pain
  FROM condition_summary
  GROUP BY person_id
)
SELECT
  a.person_id,
  COALESCE(p.diabetes, FALSE) AS diabetes,
  COALESCE(p.heart_disease, FALSE) AS heart_disease,
  COALESCE(p.hypertension, FALSE) AS hypertension,
  COALESCE(p.ckd, FALSE) AS ckd,
  COALESCE(p.chronic_pain, FALSE) AS chronic_pain
FROM adult_ehr a
LEFT JOIN person_flags p USING (person_id)"
)

chronic_conditions <- run_query(sql)
assert_unique_person(chronic_conditions, "chronic conditions")

# Count how many of the five chronic diseases each participant has. Only the
# resulting count is saved for analysis; the individual flags are intermediate
# values used to make the construction explicit and reproducible.
condition_columns <- c(
  "diabetes", "heart_disease", "hypertension", "ckd", "chronic_pain"
)
chronic_conditions$chronic_disease_count <- as.integer(rowSums(
  chronic_conditions[condition_columns]
))
if (anyNA(chronic_conditions$chronic_disease_count) ||
    any(!chronic_conditions$chronic_disease_count %in% 0:5)) {
  stop(
    "The chronic disease count must be an integer from 0 through 5.",
    call. = FALSE
  )
}

chronic_disease_count <- chronic_conditions[
  c("person_id", "chronic_disease_count")
]
assert_unique_person(chronic_disease_count, "chronic disease count")
save_component(chronic_disease_count, "05_chronic_disease_count.rds")
message("Five-condition chronic disease count built.")
