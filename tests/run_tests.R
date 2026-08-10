#!/usr/bin/env Rscript

# Public, data-free checks for the BRIDGE MDD R repository. These tests load
# selected production assignments without sourcing Workbench setup or queries.

if (!file.exists("README.md") || !dir.exists("analysis")) {
  stop("Run tests/run_tests.R from the repository root.", call. = FALSE)
}

passed <- 0L

check <- function(label, code) {
  tryCatch(
    code(),
    error = function(error) {
      stop("FAIL: ", label, ": ", conditionMessage(error), call. = FALSE)
    }
  )
  passed <<- passed + 1L
  cat("PASS: ", label, "\n", sep = "")
}

expect_error_matching <- function(code, pattern) {
  message <- tryCatch(
    {
      code()
      NA_character_
    },
    error = function(error) conditionMessage(error)
  )
  stopifnot(!is.na(message), grepl(pattern, message, ignore.case = TRUE))
}

assignment_name <- function(expression) {
  if (!is.call(expression) || length(expression) < 3L) {
    return(NA_character_)
  }
  operator <- expression[[1L]]
  target <- expression[[2L]]
  is_assignment <- identical(operator, as.name("<-")) ||
    identical(operator, as.name("="))
  if (!is_assignment || !is.symbol(target)) {
    return(NA_character_)
  }
  as.character(target)
}

load_assignments <- function(path, names, environment) {
  expressions <- parse(path)
  loaded <- character()
  for (expression in expressions) {
    target <- assignment_name(expression)
    if (!is.na(target) && target %in% names) {
      eval(expression, envir = environment)
      loaded <- c(loaded, target)
    }
  }
  missing <- setdiff(names, loaded)
  if (length(missing) > 0L) {
    stop(
      "Could not load assignments from ", path, ": ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(environment)
}

compact_source <- function(path) {
  gsub("[[:space:]]+", " ", paste(readLines(path, warn = FALSE), collapse = " "))
}

contains_all <- function(text, required) {
  all(vapply(required, function(value) grepl(value, text, fixed = TRUE), logical(1)))
}

check("all R files parse", function() {
  files <- sort(unique(unlist(lapply(
    c("analysis", "scripts", "tests"),
    function(path) {
      if (!dir.exists(path)) return(character())
      list.files(path, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
    }
  ))))
  stopifnot(length(files) > 0L)
  invisible(lapply(files, parse))
})

helpers <- new.env(parent = baseenv())
helpers$CDR <- "example-project.example_cdr"
helpers$CDR_INDEX <- "example-project.example_cdr_index"
load_assignments(
  "analysis/00_preflight.R",
  c(
    "bq_table", "bq_index_table", "sql_quote", "assert_unique_person",
    "join_component", "read_component_file"
  ),
  helpers
)

check("SQL identifiers and values are quoted safely", function() {
  stopifnot(
    identical(helpers$bq_table("person"), "`example-project.example_cdr.person`"),
    identical(
      helpers$bq_index_table("T_ENT_person"),
      "`example-project.example_cdr_index.T_ENT_person`"
    ),
    identical(helpers$sql_quote(c("plain", "O'Brien")), c("'plain'", "'O''Brien'"))
  )
  expect_error_matching(function() helpers$bq_table("person; DROP TABLE"), "unsafe")
  expect_error_matching(function() helpers$bq_index_table("index-name"), "unsafe")
})

check("participant keys must be present, unique, and nonmissing", function() {
  stopifnot(isTRUE(helpers$assert_unique_person(
    data.frame(person_id = c("p1", "p2")),
    "synthetic"
  )))
  expect_error_matching(
    function() helpers$assert_unique_person(data.frame(value = 1), "synthetic"),
    "no person_id"
  )
  expect_error_matching(
    function() helpers$assert_unique_person(
      data.frame(person_id = c("p1", "p1")),
      "synthetic"
    ),
    "one nonmissing row"
  )
  expect_error_matching(
    function() helpers$assert_unique_person(
      data.frame(person_id = c("p1", NA_character_)),
      "synthetic"
    ),
    "one nonmissing row"
  )
})

check("left joins preserve the cohort and missingness", function() {
  base <- data.frame(person_id = c("p2", "p1"), mdd_case = c(1L, 0L))
  component <- data.frame(person_id = c("p1", "p3"), score = c(7, 9))
  joined <- helpers$join_component(base, component, "synthetic component")
  stopifnot(
    identical(joined$person_id, base$person_id),
    identical(joined$mdd_case, base$mdd_case),
    is.na(joined$score[1L]),
    identical(joined$score[2L], 7)
  )
  expect_error_matching(
    function() helpers$join_component(
      base,
      data.frame(person_id = c("p1", "p2"), mdd_case = c(0L, 1L)),
      "overlap"
    ),
    "replace existing columns"
  )
})

check("delimited genomic component files are supported", function() {
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(c("person_id\tprs", "p1\t0.1", "p2\t0.2"), path)
  component <- helpers$read_component_file(path)
  stopifnot(
    identical(component$person_id, c("p1", "p2")),
    identical(component$prs, c(0.1, 0.2))
  )
})

phenotype <- new.env(parent = baseenv())
load_assignments(
  "analysis/01_build_mdd_cohort.R",
  c("mb2862", "psychiatric_exclusion_ancestors"),
  phenotype
)

check("MDD code and exclusion registries remain locked", function() {
  keys <- paste(phenotype$mb2862$vocabulary_id, phenotype$mb2862$concept_code)
  stopifnot(
    nrow(phenotype$mb2862) == 37L,
    sum(phenotype$mb2862$vocabulary_id == "ICD10CM") == 21L,
    sum(phenotype$mb2862$vocabulary_id == "ICD9CM") == 16L,
    !anyNA(keys),
    !anyDuplicated(keys),
    identical(
      phenotype$psychiatric_exclusion_ancestors,
      c(435783, 436665, 4182210, 4295956, 37016268, 37018689)
    )
  )
})

check("primary phenotype and observability rules remain explicit", function() {
  code <- compact_source("analysis/01_build_mdd_cohort.R")
  stopifnot(contains_all(code, c(
    "idx.has_ehr_data IS TRUE",
    "idx.has_whole_genome_variant IS TRUE",
    "COUNT(DISTINCT IF(qualifying_mdd, diagnosis_date, NULL)) AS qualifying_date_count",
    "COALESCE(d.qualifying_date_count, 0) >= 2 AND d.qualifying_span_days >= 30 THEN 'case'",
    "COALESCE(e.encounter_date_count, 0) >= 2 AND COALESCE(s.ehr_span_days, 0) >= 365 AS ehr_observability_eligible"
  )))
})

scores <- new.env(parent = baseenv())
load_assignments("analysis/03_build_pss_eds.R", c("items", "answers"), scores)

check("PSS-10 and EDS-9 mappings remain locked", function() {
  stopifnot(
    nrow(scores$items) == 19L,
    sum(scores$items$scale == "PSS10") == 10L,
    sum(scores$items$scale == "EDS9") == 9L,
    sum(scores$items$reverse_scored & scores$items$scale == "PSS10") == 4L,
    !anyDuplicated(scores$items$item_concept_id),
    identical(scores$answers$item_score[scores$answers$scale == "PSS10"], 0:4),
    identical(scores$answers$item_score[scores$answers$scale == "EDS9"], 0:5)
  )
  code <- compact_source("analysis/03_build_pss_eds.R")
  stopifnot(contains_all(code, c(
    "pss10_items_answered BETWEEN 8 AND 10",
    "SAFE_DIVIDE(pss10_sum, pss10_items_answered) * 10",
    "eds9_items_answered BETWEEN 7 AND 9",
    "SAFE_DIVIDE(eds9_sum, eds9_items_answered) * 9"
  )))
})

ace <- new.env(parent = baseenv())
load_assignments("analysis/04_build_ace.R", "ace_items", ace)

check("ACE-11 mapping and missing-answer rule remain locked", function() {
  stopifnot(
    nrow(ace$ace_items) == 11L,
    sum(ace$ace_items$response_type == "yes_no") == 5L,
    sum(ace$ace_items$response_type == "frequency") == 6L,
    !anyDuplicated(ace$ace_items$item_concept_id)
  )
  code <- compact_source("analysis/04_build_ace.R")
  stopifnot(contains_all(code, c(
    "Concept 1704106, 'Parents not married', remains unscored.",
    "ace_items_answered >= 8",
    "SAFE_DIVIDE(ace_exposures_reported, ace_items_answered) * 11"
  )))
})

genomics <- new.env(parent = baseenv())
load_assignments(
  "analysis/06_import_genomics.R",
  "adjust_pgs_across_ancestry",
  genomics
)

check("cross-ancestry PRS adjustment returns finite varying values", function() {
  set.seed(20260810)
  synthetic <- data.frame(
    prs = stats::rnorm(300),
    PC1 = stats::rnorm(300),
    PC2 = stats::rnorm(300),
    PC3 = stats::rnorm(300),
    PC4 = stats::rnorm(300),
    PC5 = stats::rnorm(300)
  )
  adjusted <- genomics$adjust_pgs_across_ancestry(synthetic, "prs")
  stopifnot(
    length(adjusted) == nrow(synthetic),
    all(is.finite(adjusted)),
    stats::sd(adjusted) > 0
  )
})

check("genomic QC and PRS source rules remain explicit", function() {
  import_code <- compact_source("analysis/06_import_genomics.R")
  analysis_code <- compact_source("analysis/08_run_manuscript_analysis.R")
  stopifnot(
    contains_all(import_code, c(
      "qc_flag == 0L & related_flag == 0L",
      "prs_mdd_div_raw",
      "prs_mdd_eur_raw",
      "prs_mdd_div_adjusted",
      "prs_mdd_eur_adjusted",
      "paste0(\"PC\", 1:5)"
    )),
    contains_all(analysis_code, c(
      "\"Pooled\", \"prs_mdd_div_adjusted\"",
      "\"prs_mdd_eur_adjusted\", \"prs_mdd_div_adjusted\", \"prs_mdd_div_adjusted\"",
      "\"EUR\", \"AFR\", \"AMR\""
    ))
  )
})

check("manuscript models retain modified Poisson with HC0 covariance", function() {
  code <- compact_source("analysis/08_run_manuscript_analysis.R")
  stopifnot(contains_all(code, c(
    "family = stats::poisson(link = \"log\")",
    "sandwich::vcovHC(model, type = \"HC0\")",
    "additional_clinical_covariate <- \"chronic_disease_count\"",
    "models$joint_plus_chronic_count",
    "primary_covariates"
  )))
})

cat("PASS: ", passed, " R code checks completed\n", sep = "")
