# Build the WGS-eligible EHR phenotype cohort and study cohort flow.
#
# Cases use the complete official 37-code U.S. ICD-9-CM/ICD-10-CM mapping for
# PhecodeX MB_286.2 and require at least two diagnosis dates separated by 30
# or more days. Controls have no evidence under that mapping and no positive
# current-treatment survey response. The psychiatric-exclusion concept set
# and the one-year EHR observability rule are applied equally to cases and
# controls.
# Short-read whole-genome sequencing availability is a cohort-entry
# requirement; PRS, genomic ancestry, and principal components are attached
# later.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")

# Official 37-code U.S. ICD-9-CM/ICD-10-CM PhecodeX MB_286.2 mapping.
mb2862 <- data.frame(
  vocabulary_id = c(
    rep("ICD10CM", 21),
    rep("ICD9CM", 16)
  ),
  concept_code = c(
    "F33", "F32.9", "F33.0", "F33.1", "F33.8", "F33.3", "F33.4",
    "F33.42", "F33.41", "F33.40", "F33.9", "F32.89", "F33.2", "F32",
    "F32.0", "F32.1", "F32.2", "F32.3", "F32.4", "F32.5", "F32.8",
    "296.34", "296.24", "296.31", "296.32", "296.2", "296.20",
    "296.33", "296.35", "296.30", "296.3", "296.26", "296.25",
    "296.36", "296.22", "296.21", "296.23"
  )
)

# Use the official PhecodeX phenotype without selectively removing source
# codes. The repeated-date rule supplies additional EHR phenotype specificity.
mb2862$qualifying_mdd <- TRUE

stopifnot(
  nrow(mb2862) == 37L,
  sum(mb2862$qualifying_mdd) == 37L
)

code_rows <- vapply(seq_len(nrow(mb2862)), function(i) {
  paste0(
    "STRUCT(",
    sql_quote(mb2862$vocabulary_id[i]), " AS vocabulary_id, ",
    sql_quote(mb2862$concept_code[i]), " AS concept_code, ",
    if (mb2862$qualifying_mdd[i]) "TRUE" else "FALSE",
    " AS qualifying_mdd)"
  )
}, character(1))

# Six standard OMOP ancestors define the psychiatric-exclusion concept set.
psychiatric_exclusion_ancestors <- c(
  435783,   # Schizophrenia
  436665,   # Bipolar disorder
  4182210,  # Dementia
  4295956,  # Mood disorder due to a general medical condition
  37016268, # Opioid-induced mood disorder due to opioid abuse
  37018689  # Opioid-induced mood disorder due to opioid dependence
)

# Cohort entry requires evidence of participation in all three survey modules.
# This is module availability, not a requirement for a valid answer to every
# PSS, EDS, or ACE item.
survey_gate_ctes <- paste0(
"emotional_health_participants AS (
  SELECT DISTINCT person_id
  FROM ", bq_index_table("T_ENT_surveyOccurrence"), "
  WHERE survey_item_id IN (
    SELECT descendant
    FROM ", bq_index_table("T_HAD_surveyEmotionalHealth_default"), "
    WHERE ancestor = 60000
    UNION ALL SELECT 60000
  )
),
sdoh_participants AS (
  SELECT DISTINCT person_id
  FROM ", bq_index_table("T_ENT_surveyOccurrence"), "
  WHERE survey_item_id IN (
    SELECT descendant
    FROM ",
    bq_index_table("T_HAD_surveySocialDeterminantsOfHealth_default"), "
    WHERE ancestor = 32000
    UNION ALL SELECT 32000
  )
),
basics_participants AS (
  SELECT DISTINCT person_id
  FROM ", bq_index_table("T_ENT_surveyOccurrence"), "
  WHERE survey_item_id IN (
    SELECT descendant
    FROM ", bq_index_table("T_HAD_surveyBasics_default"), "
    WHERE ancestor = 1000
    UNION ALL SELECT 1000
  )
)"
)

sql <- paste0(
"WITH mb2862 AS (
  SELECT * FROM UNNEST([
    ", paste(code_rows, collapse = ",\n    "), "
  ])
),
psychiatric_exclusion_set AS (
  SELECT DISTINCT descendant_concept_id AS concept_id
  FROM ", bq_table("concept_ancestor"), "
  WHERE ancestor_concept_id IN (
    ", paste(psychiatric_exclusion_ancestors, collapse = ", "), "
  )
  UNION DISTINCT
  SELECT concept_id
  FROM UNNEST([
    ", paste(psychiatric_exclusion_ancestors, collapse = ", "), "
  ]) AS concept_id
),
", survey_gate_ctes, ",
eligible_participants AS (
  SELECT p.person_id
  FROM ", bq_table("person"), " p
  JOIN ", bq_index_table("T_ENT_person"), " idx
    ON idx.id = p.person_id
  JOIN emotional_health_participants eh
    ON eh.person_id = p.person_id
  JOIN sdoh_participants sdoh
    ON sdoh.person_id = p.person_id
  JOIN basics_participants basics
    ON basics.person_id = p.person_id
  WHERE idx.has_ehr_data IS TRUE
    AND idx.has_whole_genome_variant IS TRUE
    AND CAST(FLOOR(
      TIMESTAMP_DIFF(
        TIMESTAMP('2026-07-23 18:00:00.000'), idx.age, DAY
      ) / 365.25
    ) AS INT64) BETWEEN 18 AND 120
),
depression_events AS (
  SELECT
    co.person_id,
    COALESCE(DATE(co.condition_start_datetime), co.condition_start_date)
      AS diagnosis_date,
    codes.qualifying_mdd
  FROM ", bq_table("condition_occurrence"), " co
  JOIN eligible_participants USING (person_id)
  JOIN ", bq_table("concept"), " source_concept
    ON co.condition_source_concept_id = source_concept.concept_id
  JOIN mb2862 codes
    ON UPPER(source_concept.vocabulary_id) = UPPER(codes.vocabulary_id)
   AND UPPER(source_concept.concept_code) = UPPER(codes.concept_code)
),
depression_summary AS (
  SELECT
    person_id,
    COUNT(*) > 0 AS any_mb2862_evidence,
    COUNT(DISTINCT IF(qualifying_mdd, diagnosis_date, NULL))
      AS qualifying_date_count,
    DATE_DIFF(
      MAX(IF(qualifying_mdd, diagnosis_date, NULL)),
      MIN(IF(qualifying_mdd, diagnosis_date, NULL)),
      DAY
    ) AS qualifying_span_days
  FROM depression_events
  GROUP BY person_id
),
psychiatric_exclusions AS (
  SELECT DISTINCT co.person_id
  FROM ", bq_table("condition_occurrence"), " co
  JOIN eligible_participants USING (person_id)
  JOIN psychiatric_exclusion_set e
    ON co.condition_concept_id = e.concept_id
),
survey_positive AS (
  SELECT DISTINCT o.person_id
  FROM ", bq_table("observation"), " o
  JOIN eligible_participants USING (person_id)
  WHERE (o.observation_source_concept_id = 43530358
         OR o.observation_concept_id = 43530358)
    AND (o.value_source_concept_id = 43530029
         OR o.value_as_concept_id = 43530029)
),
encounters AS (
  SELECT
    v.person_id,
    COUNT(DISTINCT COALESCE(
      DATE(v.visit_start_datetime), v.visit_start_date
    )) AS encounter_date_count
  FROM ", bq_table("visit_occurrence"), " v
  JOIN eligible_participants USING (person_id)
  GROUP BY v.person_id
),
ehr_activity AS (
  SELECT person_id, condition_start_date AS activity_date
  FROM ", bq_table("condition_occurrence"), "
  WHERE condition_start_date IS NOT NULL
  UNION ALL
  SELECT person_id, visit_start_date
  FROM ", bq_table("visit_occurrence"), "
  WHERE visit_start_date IS NOT NULL
  UNION ALL
  SELECT person_id, procedure_date
  FROM ", bq_table("procedure_occurrence"), "
  WHERE procedure_date IS NOT NULL
  UNION ALL
  SELECT person_id, drug_exposure_start_date
  FROM ", bq_table("drug_exposure"), "
  WHERE drug_exposure_start_date IS NOT NULL
),
ehr_span AS (
  SELECT
    activity.person_id,
    DATE_DIFF(MAX(activity_date), MIN(activity_date), DAY)
      AS ehr_span_days
  FROM ehr_activity activity
  JOIN eligible_participants USING (person_id)
  GROUP BY activity.person_id
)
SELECT
  a.person_id,
  CASE
    WHEN px.person_id IS NOT NULL THEN 'excluded_psychiatric'
    WHEN COALESCE(d.qualifying_date_count, 0) >= 2
      AND d.qualifying_span_days >= 30 THEN 'case'
    WHEN NOT COALESCE(d.any_mb2862_evidence, FALSE)
      AND sp.person_id IS NULL THEN 'control'
    ELSE 'indeterminate'
  END AS phenotype_status,
  COALESCE(e.encounter_date_count, 0) >= 2
    AND COALESCE(s.ehr_span_days, 0) >= 365 AS ehr_observability_eligible
FROM eligible_participants a
LEFT JOIN depression_summary d USING (person_id)
LEFT JOIN psychiatric_exclusions px USING (person_id)
LEFT JOIN survey_positive sp USING (person_id)
LEFT JOIN encounters e USING (person_id)
LEFT JOIN ehr_span s USING (person_id)"
)

classified <- run_query(sql)
assert_unique_person(classified, "classified MDD cohort")

# Retain the case/control contrast after the one-year EHR gate.
mdd_cohort <- classified[
  classified$phenotype_status %in% c("case", "control") &
    classified$ehr_observability_eligible,
  c("person_id", "phenotype_status"),
  drop = FALSE
]
mdd_cohort$mdd_case <- as.integer(mdd_cohort$phenotype_status == "case")
assert_unique_person(mdd_cohort, "MDD cohort")
save_component(mdd_cohort, "01_mdd_cohort.rds")

# Count the broad eligibility stages without downloading participant records.
eligibility_sql <- paste0(
"WITH ", survey_gate_ctes, "
SELECT
  COUNT(DISTINCT p.person_id) AS all_participants,
  COUNT(DISTINCT IF(
    CAST(FLOOR(
      TIMESTAMP_DIFF(
        TIMESTAMP('2026-07-23 18:00:00.000'), idx.age, DAY
      ) / 365.25
    ) AS INT64) BETWEEN 18 AND 120,
    p.person_id, NULL
  )) AS adults,
  COUNT(DISTINCT IF(
    CAST(FLOOR(
      TIMESTAMP_DIFF(
        TIMESTAMP('2026-07-23 18:00:00.000'), idx.age, DAY
      ) / 365.25
    ) AS INT64) BETWEEN 18 AND 120
      AND idx.has_ehr_data IS TRUE,
    p.person_id, NULL
  )) AS adults_with_ehr,
  COUNT(DISTINCT IF(
    CAST(FLOOR(
      TIMESTAMP_DIFF(
        TIMESTAMP('2026-07-23 18:00:00.000'), idx.age, DAY
      ) / 365.25
    ) AS INT64) BETWEEN 18 AND 120
      AND idx.has_ehr_data IS TRUE
      AND idx.has_whole_genome_variant IS TRUE,
    p.person_id, NULL
  )) AS adults_with_ehr_and_wgs
  ,
  COUNT(DISTINCT IF(
    CAST(FLOOR(
      TIMESTAMP_DIFF(
        TIMESTAMP('2026-07-23 18:00:00.000'), idx.age, DAY
      ) / 365.25
    ) AS INT64) BETWEEN 18 AND 120
      AND idx.has_ehr_data IS TRUE
      AND idx.has_whole_genome_variant IS TRUE
      AND eh.person_id IS NOT NULL,
    p.person_id, NULL
  )) AS after_emotional_health,
  COUNT(DISTINCT IF(
    CAST(FLOOR(
      TIMESTAMP_DIFF(
        TIMESTAMP('2026-07-23 18:00:00.000'), idx.age, DAY
      ) / 365.25
    ) AS INT64) BETWEEN 18 AND 120
      AND idx.has_ehr_data IS TRUE
      AND idx.has_whole_genome_variant IS TRUE
      AND eh.person_id IS NOT NULL
      AND sdoh.person_id IS NOT NULL,
    p.person_id, NULL
  )) AS after_sdoh,
  COUNT(DISTINCT IF(
    CAST(FLOOR(
      TIMESTAMP_DIFF(
        TIMESTAMP('2026-07-23 18:00:00.000'), idx.age, DAY
      ) / 365.25
    ) AS INT64) BETWEEN 18 AND 120
      AND idx.has_ehr_data IS TRUE
      AND idx.has_whole_genome_variant IS TRUE
      AND eh.person_id IS NOT NULL
      AND sdoh.person_id IS NOT NULL
      AND basics.person_id IS NOT NULL,
    p.person_id, NULL
  )) AS after_basics
FROM ", bq_table("person"), " p
LEFT JOIN ", bq_index_table("T_ENT_person"), " idx
  ON idx.id = p.person_id
LEFT JOIN emotional_health_participants eh
  ON eh.person_id = p.person_id
LEFT JOIN sdoh_participants sdoh
  ON sdoh.person_id = p.person_id
LEFT JOIN basics_participants basics
  ON basics.person_id = p.person_id"
)
eligibility_counts <- run_query(eligibility_sql)

# One aggregate row per sequential step supports the manuscript flow diagram.
# Case/control counts begin only after phenotype classification.
after_psychiatric_exclusions <- sum(
  classified$phenotype_status != "excluded_psychiatric"
)
classified_case_control <- classified$phenotype_status %in% c(
  "case", "control"
)
before_observability <- classified[classified_case_control, , drop = FALSE]

cohort_flow <- data.frame(
  step = c(
    "All All of Us participants",
    "Age 18 years or older",
    "Age 18 years or older with EHR data",
    "Age 18 years or older with EHR and WGS data",
    "Recorded Emotional Health survey participation",
    "Recorded Social Determinants of Health survey participation",
    "Recorded Basics survey participation",
    "After psychiatric exclusions",
    "Classified as an MDD case or clean control",
    "After one-year EHR observability requirement"
  ),
  total_n = c(
    as.numeric(eligibility_counts$all_participants),
    as.numeric(eligibility_counts$adults),
    as.numeric(eligibility_counts$adults_with_ehr),
    as.numeric(eligibility_counts$adults_with_ehr_and_wgs),
    as.numeric(eligibility_counts$after_emotional_health),
    as.numeric(eligibility_counts$after_sdoh),
    as.numeric(eligibility_counts$after_basics),
    after_psychiatric_exclusions,
    nrow(before_observability),
    nrow(mdd_cohort)
  ),
  case_n = c(
    NA, NA, NA, NA, NA, NA, NA, NA,
    sum(before_observability$phenotype_status == "case"),
    sum(mdd_cohort$mdd_case == 1L)
  ),
  control_n = c(
    NA, NA, NA, NA, NA, NA, NA, NA,
    sum(before_observability$phenotype_status == "control"),
    sum(mdd_cohort$mdd_case == 0L)
  ),
  stringsAsFactors = FALSE
)
cohort_flow$excluded_from_prior <- c(
  NA,
  head(cohort_flow$total_n, -1L) - tail(cohort_flow$total_n, -1L)
)
save_component(cohort_flow, "01_mdd_cohort_flow.rds")
print(cohort_flow, row.names = FALSE)

# Draw a reproducible manuscript flow diagram using only base R graphics.
flow_pdf <- file.path(WORK_DIR, "Figure_cohort_flow.pdf")
grDevices::pdf(flow_pdf, width = 8.5, height = 11)
graphics::par(mar = rep(0.5, 4))
graphics::plot.new()
graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))

flow_y <- seq(0.94, 0.28, length.out = nrow(cohort_flow))
for (i in seq_len(nrow(cohort_flow))) {
  graphics::rect(
    0.12, flow_y[i] - 0.035, 0.88, flow_y[i] + 0.035,
    border = "#2F4F4F", col = "#F4F7F7", lwd = 1.2
  )
  graphics::text(
    0.50, flow_y[i],
    labels = paste0(
      cohort_flow$step[i], "\nN = ",
      format(cohort_flow$total_n[i], big.mark = ",", scientific = FALSE)
    ),
    cex = 0.78
  )
  if (i < nrow(cohort_flow)) {
    graphics::arrows(
      0.50, flow_y[i] - 0.038, 0.50, flow_y[i + 1L] + 0.038,
      length = 0.07, lwd = 1
    )
    graphics::text(
      0.76, mean(flow_y[c(i, i + 1L)]),
      labels = paste0(
        "Excluded: ",
        format(
          cohort_flow$excluded_from_prior[i + 1L],
          big.mark = ",", scientific = FALSE
        )
      ),
      cex = 0.64, adj = 0
    )
  }
}

final_y <- 0.12
graphics::arrows(
  0.50, flow_y[nrow(cohort_flow)] - 0.038,
  0.28, final_y + 0.045, length = 0.07, lwd = 1
)
graphics::arrows(
  0.50, flow_y[nrow(cohort_flow)] - 0.038,
  0.72, final_y + 0.045, length = 0.07, lwd = 1
)
graphics::rect(
  0.08, final_y - 0.045, 0.46, final_y + 0.045,
  border = "#2F4F4F", col = "#EAF3F8", lwd = 1.2
)
graphics::rect(
  0.54, final_y - 0.045, 0.92, final_y + 0.045,
  border = "#2F4F4F", col = "#EEF6EA", lwd = 1.2
)
graphics::text(
  0.27, final_y,
  labels = paste0(
    "MDD cases\nN = ",
    format(sum(mdd_cohort$mdd_case == 1L), big.mark = ",")
  ),
  cex = 0.82
)
graphics::text(
  0.73, final_y,
  labels = paste0(
    "Clean controls\nN = ",
    format(sum(mdd_cohort$mdd_case == 0L), big.mark = ",")
  ),
  cex = 0.82
)
grDevices::dev.off()

message(
  "WGS-eligible MDD cohort built: ",
  sum(mdd_cohort$mdd_case == 1L), " cases and ",
  sum(mdd_cohort$mdd_case == 0L), " controls."
)
