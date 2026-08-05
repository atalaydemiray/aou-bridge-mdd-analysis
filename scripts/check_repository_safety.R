#!/usr/bin/env Rscript

# Guard the public repository boundary. This check cannot certify privacy; it
# rejects common participant-data formats, protected paths, cloud identifiers,
# and credential signatures before they are merged into the public repository.

if (!file.exists("README.md") || !dir.exists("analysis")) {
  stop(
    "Run scripts/check_repository_safety.R from the repository root.",
    call. = FALSE
  )
}

git_files <- system2(
  "git",
  c("ls-files", "--cached", "--others", "--exclude-standard"),
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(git_files, "status")
if (!is.null(status) && status != 0L) {
  stop("Could not enumerate repository files with git ls-files.", call. = FALSE)
}

files <- sort(unique(git_files[nzchar(git_files)]))
files <- files[file.exists(files) & !dir.exists(files)]
if (length(files) == 0L) {
  stop("No repository files were found.", call. = FALSE)
}

fail_if_any <- function(values, label) {
  if (length(values) > 0L) {
    stop(label, ": ", paste(values, collapse = ", "), call. = FALSE)
  }
}

lower_files <- tolower(files)
blocked_extensions <- c(
  "csv", "tsv", "csv.gz", "tsv.gz", "rds", "rda", "rdata",
  "parquet", "feather", "fst", "sas7bdat", "dta", "sav", "xls",
  "xlsx", "db", "sqlite", "duckdb", "ipynb", "zip", "7z"
)
blocked_extension_files <- files[vapply(
  lower_files,
  function(path) any(endsWith(path, paste0(".", blocked_extensions))),
  logical(1)
)]
fail_if_any(blocked_extension_files, "Protected or generated file type found")

blocked_paths <- files[grepl(
  "^(data|work|cache|results|outputs)/",
  lower_files,
  perl = TRUE
)]
fail_if_any(blocked_paths, "Protected output directory found")

blocked_names <- files[grepl(
  "(^|/)([.]renviron|[.]env([.]|$)|credentials?([._-]|$)|service[-_]?account([._-]|$))",
  lower_files,
  perl = TRUE
)]
fail_if_any(blocked_names, "Credential-like filename found")

read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

text <- vapply(files, read_text, character(1), USE.NAMES = TRUE)
content_patterns <- c(
  private_key = paste0("-----BEGIN ", "(?:RSA |EC |OPENSSH )?", "PRIVATE KEY-----"),
  google_api_key = "AIza[0-9A-Za-z_-]{30,}",
  github_token = "(github_pat_[0-9A-Za-z_]{20,}|gh[pousr]_[0-9A-Za-z]{20,})",
  aws_access_key = "AKIA[0-9A-Z]{16}",
  service_account_key = paste0('"private_', 'key"[[:space:]]*:'),
  cloud_bucket = "gs://[a-z0-9][a-z0-9._-]{2,}",
  workbench_project = "wb-[a-z][a-z0-9-]*-[0-9]{3,}"
)

for (label in names(content_patterns)) {
  matched <- names(text)[vapply(
    text,
    function(value) grepl(content_patterns[[label]], value, perl = TRUE),
    logical(1)
  )]
  fail_if_any(matched, paste("Potential", label, "content found"))
}

ignore_lines <- trimws(readLines(".gitignore", warn = FALSE))
required_ignores <- c(
  ".Renviron", ".env", "data/", "work/", "cache/", "results/",
  "outputs/", "*.csv", "*.tsv", "*.rds", "*.parquet", "*.ipynb"
)
missing_ignores <- setdiff(required_ignores, ignore_lines)
fail_if_any(missing_ignores, "Required .gitignore rule missing")

readme <- gsub(
  "[[:space:]]+",
  " ",
  read_text("README.md")
)
required_notices <- c(
  "No participant-level data",
  "remain inside an authorized Researcher Workbench",
  "do not rerun the live CDR or genomic analysis"
)
missing_notices <- required_notices[!vapply(
  required_notices,
  function(value) grepl(value, readme, fixed = TRUE),
  logical(1)
)]
fail_if_any(missing_notices, "Required README data-safety notice missing")

cat(
  "PASS: repository data-safety checks completed for ",
  length(files),
  " files\n",
  sep = ""
)
