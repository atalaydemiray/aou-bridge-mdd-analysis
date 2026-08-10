# BRIDGE AoU MDD study: complete preflight
#
# Run this first from the repository root. It checks the Workbench
# environment, installs/loads the required R packages, and defines the small
# set of functions used by the remaining analysis scripts.

options(stringsAsFactors = FALSE)

if (!file.exists("README.md") || !dir.exists("analysis")) {
  stop("Run this script from the repository root.", call. = FALSE)
}

# Parse every analysis script before any participant-level query is submitted.
analysis_files <- sort(list.files(
  "analysis", pattern = "[.]R$", full.names = TRUE
))
invisible(lapply(analysis_files, parse))

required_packages <- c("bigrquery", "data.table", "lmtest", "sandwich")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    # Install only required dependencies. Suggested packages add a large,
    # unnecessary compilation burden in a clean Workbench app.
    dependencies = NA
  )
}
invisible(lapply(required_packages, library, character.only = TRUE))

# The RStudio Console may start without the environment variables that are
# already available in the app Terminal. If either required value is missing,
# use the Workbench-provided loader and import only these two values. Do not
# import the complete shell environment because it can contain credentials.
load_workbench_console_env <- function() {
  required <- c("WORKSPACE_CDR", "GOOGLE_PROJECT")
  current <- Sys.getenv(required, unset = "")
  if (all(nzchar(current))) {
    return(invisible(FALSE))
  }

  loader <- "/home/rstudio/load-env.sh"
  if (!file.exists(loader)) {
    return(invisible(FALSE))
  }

  command <- paste(
    "bash -lc",
    shQuote(paste(
      "source /home/rstudio/load-env.sh >/dev/null 2>&1;",
      "printf 'WORKSPACE_CDR=%s\\nGOOGLE_PROJECT=%s\\n'",
      "\"$WORKSPACE_CDR\" \"$GOOGLE_PROJECT\""
    ))
  )
  env_lines <- tryCatch(
    system(command, intern = TRUE, ignore.stderr = TRUE),
    error = function(e) character()
  )

  loaded <- setNames(rep("", length(required)), required)
  for (key in required) {
    prefix <- paste0(key, "=")
    match <- env_lines[startsWith(env_lines, prefix)]
    if (length(match) == 1L) {
      loaded[[key]] <- substring(match, nchar(prefix) + 1L)
    }
  }

  missing <- !nzchar(current)
  values_to_set <- loaded[missing & nzchar(loaded)]
  if (length(values_to_set) > 0L) {
    do.call(Sys.setenv, as.list(values_to_set))
  }
  invisible(length(values_to_set) > 0L)
}

load_workbench_console_env()

# Workbench supplies these variables when a CDR is attached.
CDR <- Sys.getenv("WORKSPACE_CDR", unset = "")
BILLING_PROJECT <- Sys.getenv("GOOGLE_PROJECT", unset = "")
if (!nzchar(CDR) || !nzchar(BILLING_PROJECT)) {
  stop(
    "WORKSPACE_CDR or GOOGLE_PROJECT is missing after attempting the ",
    "Workbench environment setup. Confirm that this is an R Analysis ",
    "Environment for AoU app with CDR C2025Q4R6 attached, then restart ",
    "the app and source analysis/00_preflight.R again.",
    call. = FALSE
  )
}
if (!grepl("C2025Q4R6", CDR, fixed = TRUE)) {
  stop(
    "This release targets CDR C2025Q4R6. Review and document any CDR update ",
    "before running the manuscript analysis.",
    call. = FALSE
  )
}

# Data Explorer availability flags are populated in the release-matched index
# dataset rather than the ordinary CDR cb_search_person table. The environment
# variable permits an explicit override if a future release changes the index
# suffix.
CDR_INDEX <- Sys.getenv(
  "AOU_CDR_INDEX",
  unset = paste0(CDR, "_index_061026")
)

# Participant-level products remain in Workbench and are ignored by Git.
WORK_DIR <- Sys.getenv(
  "AOU_BRIDGE_WORK_DIR",
  unset = path.expand("~/aou_bridge_mdd_work")
)
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

# Return a fully qualified CDR table reference.
bq_table <- function(name) {
  if (!grepl("^[A-Za-z0-9_]+$", name)) {
    stop("Unsafe table name: ", name, call. = FALSE)
  }
  paste0("`", CDR, ".", name, "`")
}

# Return a fully qualified Data Explorer index-table reference.
bq_index_table <- function(name) {
  if (!grepl("^[A-Za-z0-9_]+$", name)) {
    stop("Unsafe index table name: ", name, call. = FALSE)
  }
  paste0("`", CDR_INDEX, ".", name, "`")
}

# Submit Standard SQL and download the result into R.
run_query <- function(sql) {
  job <- bigrquery::bq_project_query(
    BILLING_PROJECT,
    sql,
    use_legacy_sql = FALSE,
    quiet = FALSE
  )
  as.data.frame(
    bigrquery::bq_table_download(job, bigint = "integer64", quiet = FALSE)
  )
}

# Quote a character vector for an inline SQL list.
sql_quote <- function(x) {
  paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
}

# Every manuscript component must contain one unique row per participant.
assert_unique_person <- function(x, name) {
  if (!"person_id" %in% names(x)) {
    stop(name, " has no person_id column.", call. = FALSE)
  }
  if (anyNA(x$person_id) || anyDuplicated(x$person_id)) {
    stop(name, " must contain one nonmissing row per person_id.", call. = FALSE)
  }
  invisible(TRUE)
}

# Left joins preserve the MDD cohort exactly and never replace existing fields.
join_component <- function(base, component, name) {
  assert_unique_person(base, "analytic base")
  assert_unique_person(component, name)
  overlap <- intersect(
    setdiff(names(base), "person_id"),
    setdiff(names(component), "person_id")
  )
  if (length(overlap) > 0L) {
    stop(
      name, " would replace existing columns: ",
      paste(overlap, collapse = ", "),
      call. = FALSE
    )
  }
  index <- match(base$person_id, component$person_id)
  joined <- cbind(
    base,
    component[index, setdiff(names(component), "person_id"), drop = FALSE]
  )
  if (nrow(joined) != nrow(base) ||
      !identical(joined$person_id, base$person_id)) {
    stop("The ", name, " join changed the analytic base.", call. = FALSE)
  }
  joined
}

save_component <- function(x, filename) {
  saveRDS(x, file.path(WORK_DIR, filename))
  invisible(x)
}

load_component <- function(filename) {
  path <- file.path(WORK_DIR, filename)
  if (!file.exists(path)) {
    stop("Missing ", filename, ". Run its build script first.", call. = FALSE)
  }
  readRDS(path)
}

read_component_file <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    stop("Component file does not exist: ", path, call. = FALSE)
  }
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("csv", "txt", "tsv", "tab")) {
    return(data.table::fread(path, data.table = FALSE))
  }
  if (extension == "rds") {
    return(readRDS(path))
  }
  stop(
    "Use a CSV, TXT/TSV, or RDS component file: ", path,
    call. = FALSE
  )
}

.preflight_passed <- TRUE
message(
  "Preflight passed for ", CDR, ". ",
  "Continue by sourcing the numbered R scripts in this Console."
)
