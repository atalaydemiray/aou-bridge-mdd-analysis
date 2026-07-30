# Calculate scores for the 10-item Perceived Stress Scale (PSS-10)
# and the 9-item Everyday Discrimination Scale (EDS-9).
#
# PSS-10 is available with 8-10 scored responses and EDS-9 with 7-9 scored
# responses. Each score is the completed-item mean multiplied by the total
# number of scale items. The exact item and answer mappings are stated below.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")

items <- data.frame(
  scale = c(rep("PSS10", 10), rep("EDS9", 9)),
  item_concept_id = c(
    40192452, 40192381, 40192491, 40192419, 40192525,
    40192506, 40192449, 40192445, 40192396, 40192462,
    40192466, 40192489, 40192416, 40192490, 40192380,
    40192395, 40192496, 40192519, 40192451
  ),
  reverse_scored = c(
    FALSE, FALSE, FALSE, TRUE, TRUE,
    FALSE, TRUE, TRUE, FALSE, FALSE,
    rep(FALSE, 9)
  )
)

answers <- data.frame(
  scale = c(rep("PSS10", 5), rep("EDS9", 6)),
  answer_concept_id = c(
    40192465, 40192430, 40192429, 40192477, 40192424,
    40192465, 40192464, 40192453, 40192461, 40192391, 40192421
  ),
  item_score = c(0:4, 0:4, NA)
)

item_rows <- vapply(seq_len(nrow(items)), function(i) {
  paste0(
    "STRUCT(", items$item_concept_id[i], " AS item_concept_id, ",
    sql_quote(items$scale[i]), " AS scale, ",
    if (items$reverse_scored[i]) "TRUE" else "FALSE",
    " AS reverse_scored)"
  )
}, character(1))

answer_rows <- vapply(seq_len(nrow(answers)), function(i) {
  score <- if (is.na(answers$item_score[i])) {
    "CAST(NULL AS INT64)"
  } else {
    as.character(answers$item_score[i])
  }
  paste0(
    "STRUCT(", answers$answer_concept_id[i], " AS answer_concept_id, ",
    sql_quote(answers$scale[i]), " AS scale, ", score, " AS item_score)"
  )
}, character(1))

sql <- paste0(
"WITH items AS (
  SELECT * FROM UNNEST([
    ", paste(item_rows, collapse = ",\n    "), "
  ])
),
answers AS (
  SELECT * FROM UNNEST([
    ", paste(answer_rows, collapse = ",\n    "), "
  ])
),
responses AS (
  SELECT
    o.person_id,
    i.scale,
    i.reverse_scored,
    a.item_score
  FROM ", bq_table("observation"), " o
  JOIN items i
    ON o.observation_source_concept_id = i.item_concept_id
    OR o.observation_concept_id = i.item_concept_id
  LEFT JOIN answers a
    ON a.scale = i.scale
   AND a.answer_concept_id = COALESCE(
     NULLIF(o.value_source_concept_id, 0),
     NULLIF(o.value_as_concept_id, 0)
   )
),
scored AS (
  SELECT
    person_id,
    scale,
    CASE
      WHEN item_score IS NULL THEN NULL
      WHEN reverse_scored THEN 4 - item_score
      ELSE item_score
    END AS item_score
  FROM responses
),
person_scores AS (
  SELECT
    person_id,
    COUNTIF(scale = 'PSS10' AND item_score IS NOT NULL)
      AS pss10_items_answered,
    SUM(IF(scale = 'PSS10', item_score, NULL)) AS pss10_sum,
    COUNTIF(scale = 'EDS9' AND item_score IS NOT NULL)
      AS eds9_items_answered,
    SUM(IF(scale = 'EDS9', item_score, NULL)) AS eds9_sum
  FROM scored
  GROUP BY person_id
)
SELECT
  person_id,
  pss10_items_answered,
  IF(
    pss10_items_answered BETWEEN 8 AND 10,
    SAFE_DIVIDE(pss10_sum, pss10_items_answered) * 10,
    NULL
  ) AS pss10_score,
  eds9_items_answered,
  IF(
    eds9_items_answered BETWEEN 7 AND 9,
    SAFE_DIVIDE(eds9_sum, eds9_items_answered) * 9,
    NULL
  ) AS eds9_score
FROM person_scores"
)

pss_eds <- run_query(sql)
assert_unique_person(pss_eds, "PSS/EDS")
save_component(pss_eds, "03_pss_eds.rds")
message("PSS-10 and EDS-9 scores built.")
