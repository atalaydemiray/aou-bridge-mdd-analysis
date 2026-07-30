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
    "All of Us participants",
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
# Retained participants follow the central path. Exclusions appear in a
# separate right-hand column so labels do not overlap the flow boxes.
format_flow_n <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

flow_pdf <- file.path(WORK_DIR, "Figure_cohort_flow.pdf")
grDevices::pdf(
  flow_pdf,
  width = 8.5,
  height = 11,
  title = "Participant selection and primary MDD cohort"
)
graphics::par(mar = rep(0.25, 4), family = "sans", xpd = NA)
graphics::plot.new()
graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))

graphics::text(
  0.50, 0.975,
  labels = "Participant selection and primary MDD cohort",
  font = 2, cex = 1.05, col = "#18343B"
)

main_left <- 0.06
main_right <- 0.69
main_mid <- mean(c(main_left, main_right))
exclusion_left <- 0.76
exclusion_right <- 0.97
box_half_height <- 0.027
flow_y <- seq(0.91, 0.25, length.out = nrow(cohort_flow))

for (i in seq_len(nrow(cohort_flow))) {
  final_cohort_step <- i == nrow(cohort_flow)
  box_fill <- if (final_cohort_step) "#DCEBF3" else "#F3F7F8"
  box_border <- if (final_cohort_step) "#1F5A70" else "#355C64"
  box_line_width <- if (final_cohort_step) 1.7 else 1.1

  graphics::rect(
    main_left, flow_y[i] - box_half_height,
    main_right, flow_y[i] + box_half_height,
    border = box_border, col = box_fill, lwd = box_line_width
  )

  count_label <- paste0("N = ", format_flow_n(cohort_flow$total_n[i]))
  if (!is.na(cohort_flow$case_n[i])) {
    count_label <- paste0(
      count_label,
      "  |  cases = ", format_flow_n(cohort_flow$case_n[i]),
      "  |  controls = ", format_flow_n(cohort_flow$control_n[i])
    )
  }

  graphics::text(
    main_mid, flow_y[i] + 0.009,
    labels = if (final_cohort_step) {
      "Primary MDD cohort after one-year EHR observability"
    } else {
      cohort_flow$step[i]
    },
    cex = 0.67,
    font = if (final_cohort_step) 2 else 1,
    col = "#17282C"
  )
  graphics::text(
    main_mid, flow_y[i] - 0.011,
    labels = count_label,
    cex = 0.61,
    col = "#243A3F"
  )

  if (i < nrow(cohort_flow)) {
    connector_y <- mean(flow_y[c(i, i + 1L)])

    graphics::arrows(
      main_mid, flow_y[i] - box_half_height,
      main_mid, flow_y[i + 1L] + box_half_height,
      length = 0.055, lwd = 1.0, col = "#33464A"
    )
    graphics::segments(
      main_mid, connector_y, exclusion_left, connector_y,
      lwd = 0.9, col = "#6C4B4B"
    )
    graphics::arrows(
      exclusion_left - 0.02, connector_y,
      exclusion_left, connector_y,
      length = 0.045, lwd = 0.9, col = "#6C4B4B"
    )
    graphics::rect(
      exclusion_left, connector_y - 0.019,
      exclusion_right, connector_y + 0.019,
      border = "#966B6B", col = "#FBF2F1", lwd = 0.9
    )
    graphics::text(
      mean(c(exclusion_left, exclusion_right)), connector_y,
      labels = paste0(
        "Excluded, n = ",
        format_flow_n(cohort_flow$excluded_from_prior[i + 1L])
      ),
      cex = 0.58, col = "#583B3B"
    )
  }
}

final_y <- 0.085
case_left <- 0.07
case_right <- 0.43
control_left <- 0.57
control_right <- 0.93
final_half_height <- 0.040

graphics::arrows(
  main_mid, flow_y[nrow(cohort_flow)] - box_half_height,
  mean(c(case_left, case_right)), final_y + final_half_height,
  length = 0.060, lwd = 1.1, col = "#33464A"
)
graphics::arrows(
  main_mid, flow_y[nrow(cohort_flow)] - box_half_height,
  mean(c(control_left, control_right)), final_y + final_half_height,
  length = 0.060, lwd = 1.1, col = "#33464A"
)

graphics::rect(
  case_left, final_y - final_half_height,
  case_right, final_y + final_half_height,
  border = "#2D6176", col = "#E4F0F6", lwd = 1.5
)
graphics::rect(
  control_left, final_y - final_half_height,
  control_right, final_y + final_half_height,
  border = "#55704B", col = "#EDF4E8", lwd = 1.5
)
graphics::text(
  mean(c(case_left, case_right)), final_y + 0.013,
  labels = "MDD cases", font = 2, cex = 0.78, col = "#173B49"
)
graphics::text(
  mean(c(case_left, case_right)), final_y - 0.016,
  labels = paste0(
    "n = ", format_flow_n(cohort_flow$case_n[nrow(cohort_flow)])
  ),
  cex = 0.72, col = "#173B49"
)
graphics::text(
  mean(c(control_left, control_right)), final_y + 0.013,
  labels = "Clean controls", font = 2, cex = 0.78, col = "#34492E"
)
graphics::text(
  mean(c(control_left, control_right)), final_y - 0.016,
  labels = paste0(
    "n = ", format_flow_n(cohort_flow$control_n[nrow(cohort_flow)])
  ),
  cex = 0.72, col = "#34492E"
)
graphics::text(
  0.50, 0.018,
  labels = paste(
    "Model-specific sample sizes may be smaller because of missing exposure,",
    "genomic, or covariate data."
  ),
  cex = 0.47, col = "#4B5659"
)
grDevices::dev.off()

message(
  "WGS-eligible MDD cohort built: ",
  sum(mdd_cohort$mdd_case == 1L), " cases and ",
  sum(mdd_cohort$mdd_case == 0L), " controls."
)
