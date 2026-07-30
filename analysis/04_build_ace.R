# Calculate the modified 11-item Adverse Childhood Experiences (ACE-11) score.
#
# A valid item is coded 1 for Yes/any reported occurrence and 0 for No/Never.
# Other responses are missing. Participants with at least 8 valid items receive
# a prorated score: (number exposed / number answered) * 11.
#
# Scoring note: for parental separation/divorce item 1333208, answer
# concept 1704106 ("Parents not married") is currently treated as missing. It
# is neither coded 0 nor coded 1 and does not enter the answered-item
# denominator. Any future change to this rule should be documented explicitly.
#
# The questionnaire should provide one answer per participant and ACE item.
# The CDR can contain exact duplicate OMOP rows, so identical records are
# collapsed. The script still stops if duplicate rows contain different
# answer concepts.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")

ace_items <- data.frame(
  item_concept_id = c(
    1332946, 1332947, 1332948, 1332950, 1333208,
    1332857, 1333210, 1333212, 1332951, 1333202, 1333203
  ),
  response_type = c(rep("yes_no", 5), rep("frequency", 6))
)

item_rows <- vapply(seq_len(nrow(ace_items)), function(i) {
  paste0(
    "STRUCT(", ace_items$item_concept_id[i], " AS item_concept_id, ",
    sql_quote(ace_items$response_type[i]), " AS response_type)"
  )
}, character(1))

sql <- paste0(
"WITH items AS (
  SELECT * FROM UNNEST([
    ", paste(item_rows, collapse = ",\n    "), "
  ])
),
responses AS (
  SELECT
    o.person_id,
    i.item_concept_id,
    i.response_type,
    COALESCE(
      NULLIF(o.value_source_concept_id, 0),
      NULLIF(o.value_as_concept_id, 0)
    ) AS answer_concept_id,
    LOWER(TRIM(answer.concept_name)) AS answer_name
  FROM ", bq_table("observation"), " o
  JOIN items i
    ON o.observation_source_concept_id = i.item_concept_id
    OR o.observation_concept_id = i.item_concept_id
  LEFT JOIN ", bq_table("concept"), " answer
    ON answer.concept_id = COALESCE(
      NULLIF(o.value_source_concept_id, 0),
      NULLIF(o.value_as_concept_id, 0)
    )
),
scored_rows AS (
  SELECT
    person_id,
    item_concept_id,
    answer_concept_id,
    CASE
      WHEN response_type = 'yes_no' AND answer_name = 'yes' THEN 1
      WHEN response_type = 'yes_no' AND answer_name = 'no' THEN 0
      WHEN response_type = 'frequency'
        AND answer_name IN ('once', 'more than once') THEN 1
      WHEN response_type = 'frequency' AND answer_name = 'never' THEN 0
      -- Concept 1704106, 'Parents not married', remains unscored.
      -- Any future change to this rule should be documented explicitly.
      ELSE NULL
    END AS item_score
  FROM responses
),
person_item AS (
  SELECT
    person_id,
    item_concept_id,
    COUNT(DISTINCT COALESCE(
      CAST(answer_concept_id AS STRING), 'NULL'
    )) AS response_signature_count,
    MAX(item_score) AS item_score
  FROM scored_rows
  GROUP BY person_id, item_concept_id
),
person_score AS (
  SELECT
    person_id,
    COUNTIF(response_signature_count > 1) AS conflicting_ace_items,
    COUNT(item_score) AS ace_items_answered,
    SUM(item_score) AS ace_exposures_reported
  FROM person_item
  GROUP BY person_id
)
SELECT
  person_id,
  conflicting_ace_items,
  ace_items_answered,
  ace_exposures_reported,
  IF(
    ace_items_answered >= 8,
    SAFE_DIVIDE(ace_exposures_reported, ace_items_answered) * 11,
    NULL
  ) AS ace_score
FROM person_score"
)

ace <- run_query(sql)
assert_unique_person(ace, "ACE")
if (any(ace$conflicting_ace_items > 0, na.rm = TRUE)) {
  stop(
    "Conflicting ACE answers were found for at least one participant-item ",
    "pair. Resolve the extraction before scoring.",
    call. = FALSE
  )
}
ace$conflicting_ace_items <- NULL
save_component(ace, "04_ace.rds")
message("ACE scores built for participants with at least 8 valid items.")
