# Import the PRS, genomic ancestry, and principal-component results.
#
# IMPORTANT: Script 01 confirms that participants have short-read WGS data,
# but WGS availability is not a PRS result. Before running the full analysis,
# obtain these two approved files inside the authorized Workbench:
#
# 1. A one-row-per-participant genomic-results CSV or RDS file containing:
#      research_id
#      prs_trans_ancestry
#      genomic_ancestry
#      PC1 through PC10
#
# 2. A one-to-one CSV or RDS ID-bridge file containing:
#      research_id
#      person_id
#
# The genomic-results file must contain the finalized trans-ancestry MDD PRS,
# assigned genomic ancestry, and the ancestry principal components produced by
# the genomic analysis. This script imports those results; it does not
# calculate a PRS or infer ancestry.
#
# Keep both files inside the authorized Workbench. Never place participant
# data, identifiers, file paths, or bucket names in this Git repository.
# This script reads local CSV or RDS files. If an approved input is supplied
# as a gs:// URI, first use gsutil cp in the Workbench Terminal to copy it to a
# directory outside this Git repository, such as ~/aou_bridge_mdd_inputs/.
#
# In the SAME R Console session used to run this script, point to the files:
#
# Sys.setenv(
#   BRIDGE_MDD_GENOMICS_PATH =
#     "/path/inside/authorized/workbench/genomics_results.csv",
#   BRIDGE_MDD_ID_MAP_PATH =
#     "/path/inside/authorized/workbench/genomic_id_bridge.csv"
# )
#
# The default column names are those shown above. If the supplied files use
# different names, also set one or more of:
#
# Sys.setenv(
#   BRIDGE_MDD_GENOMICS_RESEARCH_ID = "genomic_file_id_column",
#   BRIDGE_MDD_MAP_RESEARCH_ID = "id_bridge_research_id_column",
#   BRIDGE_MDD_MAP_PERSON_ID = "id_bridge_person_id_column"
# )
#
# Sys.setenv() lasts only for the current R session. Repeat it after restarting
# the R app or Console. See README Step 6 before running scripts 06 through 08.

if (!exists(".preflight_passed")) source("analysis/00_preflight.R")

GENOMICS_PATH <- Sys.getenv("BRIDGE_MDD_GENOMICS_PATH", unset = "")
GENOMICS_ID_MAP_PATH <- Sys.getenv("BRIDGE_MDD_ID_MAP_PATH", unset = "")

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

GENOMICS_RESEARCH_ID <- env_or_default(
  "BRIDGE_MDD_GENOMICS_RESEARCH_ID", "research_id"
)
MAP_RESEARCH_ID <- env_or_default(
  "BRIDGE_MDD_MAP_RESEARCH_ID", "research_id"
)
MAP_PERSON_ID <- env_or_default(
  "BRIDGE_MDD_MAP_PERSON_ID", "person_id"
)

if (!nzchar(GENOMICS_PATH) || !nzchar(GENOMICS_ID_MAP_PATH)) {
  stop(
    "Genomic inputs are not configured. In this R Console session, set ",
    "BRIDGE_MDD_GENOMICS_PATH and BRIDGE_MDD_ID_MAP_PATH with Sys.setenv(), ",
    "then rerun analysis/06_import_genomics.R. See README Step 6.",
    call. = FALSE
  )
}

genomics_raw <- read_component_file(path.expand(GENOMICS_PATH))
id_map <- read_component_file(path.expand(GENOMICS_ID_MAP_PATH))

required_genomic <- c(
  GENOMICS_RESEARCH_ID,
  "prs_trans_ancestry",
  "genomic_ancestry",
  paste0("PC", 1:10)
)
missing_genomic <- setdiff(required_genomic, names(genomics_raw))
if (length(missing_genomic) > 0L) {
  stop(
    "The genomics file is missing: ",
    paste(missing_genomic, collapse = ", "),
    call. = FALSE
  )
}

required_map <- c(MAP_RESEARCH_ID, MAP_PERSON_ID)
missing_map <- setdiff(required_map, names(id_map))
if (length(missing_map) > 0L) {
  stop(
    "The genomic ID map is missing: ",
    paste(missing_map, collapse = ", "),
    call. = FALSE
  )
}

# Convert join keys to character before validation and matching. This avoids
# silent failure when one CSV stores an identifier as text and the other as a
# numeric or integer64 field.
genomics_key <- as.character(genomics_raw[[GENOMICS_RESEARCH_ID]])
map_research_key <- as.character(id_map[[MAP_RESEARCH_ID]])
map_person_id <- as.character(id_map[[MAP_PERSON_ID]])

invalid_key <- function(x) anyNA(x) || any(!nzchar(trimws(x)))
if (invalid_key(genomics_key) ||
    invalid_key(map_research_key) ||
    invalid_key(map_person_id)) {
  stop("Genomic identifiers must be nonmissing and nonblank.", call. = FALSE)
}
if (anyDuplicated(genomics_key) ||
    anyDuplicated(map_research_key) ||
    anyDuplicated(map_person_id)) {
  stop("Genomic identifiers must map one-to-one.", call. = FALSE)
}

map_index <- match(genomics_key, map_research_key)
person_id <- map_person_id[map_index]
if (anyNA(person_id)) {
  stop(
    "Some genomic research identifiers do not map to person_id.",
    call. = FALSE
  )
}

genomics <- data.frame(
  person_id = person_id,
  genomics_raw[
    c("prs_trans_ancestry", "genomic_ancestry", paste0("PC", 1:10))
  ],
  check.names = FALSE
)
assert_unique_person(genomics, "genomics")
save_component(genomics, "06_genomics.rds")
message(
  "Genomic component imported and mapped to person_id for ",
  nrow(genomics), " participants."
)
