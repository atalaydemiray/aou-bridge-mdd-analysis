# All of Us MDD (BRIDGE) manuscript analysis

[![R code checks](https://github.com/atalaydemiray/aou-bridge-mdd-analysis/actions/workflows/r-code.yml/badge.svg?branch=main&event=push)](https://github.com/atalaydemiray/aou-bridge-mdd-analysis/actions/workflows/r-code.yml)
[![Data-safety checks](https://github.com/atalaydemiray/aou-bridge-mdd-analysis/actions/workflows/data-safety.yml/badge.svg?branch=main&event=push)](https://github.com/atalaydemiray/aou-bridge-mdd-analysis/actions/workflows/data-safety.yml)

This public R repository reproduces the data construction and planned
statistical analysis for the BRIDGE group study of genetic liability,
psychosocial exposures, and lifetime recorded major depressive disorder
(MDD).

No participant-level data, identifiers, credentials, bucket names, or
workspace-specific paths may be committed.

## Automated checks

GitHub Actions runs two public-repository workflows after each push and pull
request. The R code workflow parses every R script, exercises the production
join and validation helpers with synthetic data, and verifies locked phenotype
and scoring definitions. The data-safety workflow separately rejects tracked
participant-data formats, credentials, cloud resource identifiers, and
protected output directories.

These checks use no All of Us participant data and do not rerun the live CDR or
genomic analysis. Each passing badge confirms only its named public-repository
check; final scientific validation still occurs inside the authorized
Researcher Workbench.

## Primary analysis specification

- Adult All of Us participants with EHR and short-read whole-genome
  sequencing data, identified by `T_ENT_person.has_whole_genome_variant`
  in the release-matched Data Explorer index.
- Participation in the Emotional Health, Social Determinants of Health, and
  Basics survey modules.
- MDD cases: the complete official 37-code U.S. ICD-9-CM/ICD-10-CM
  PhecodeX `MB_286.2` MDD mapping, with at least two diagnosis dates
  separated by at least 30 days and no upper time limit.
- Controls: no evidence under the full `MB_286.2` mapping and no positive
  current-depression-treatment survey response.
- Six psychiatric-exclusion concept families apply to both groups.
- Cases and controls require at least two encounter dates and at least 365
  days between first and last EHR activity.
- Primary adjustment uses age, sex at birth, and genomic principal components.
- PSS-10 scores each of 10 items from 0 to 4, reverse-scores four positively
  worded items, and requires at least 8 of 10 scored items. The completed-item
  mean is multiplied by 10, producing a prorated score from 0 to 40.
- EDS-9 scores each of 9 items from 0 to 5 and requires at least 7 of 9 scored
  items. The completed-item mean is multiplied by 9, producing a prorated
  score from 0 to 45.
- ACE scores each of 11 items as 1 for an adverse experience and 0 for no
  reported adverse experience, and requires at least 8 of 11 scored items.
  The completed-item mean is multiplied by 11, producing a prorated score
  from 0 to 11. Other responses, including “Parents not married,” remain
  unscored.
- The clinical adjustment is a 0-5 count of diabetes, heart disease,
  hypertension, chronic kidney disease, and chronic pain.
- The primary effect measure is the adjusted prevalence ratio from modified
  Poisson regression with robust standard errors.

See [METHODS.md](METHODS.md) for the operational definitions.

## Repository structure

```text
analysis/   all ordered, self-contained R scripts
```

There is no separate helper-code folder. Data-extraction SQL is embedded in
the R script that creates each research component. Genomic file locations are
provided only in the active R session and are never stored in the repository.

| Script | Purpose |
|---|---|
| `00_preflight.R` | Prepare and verify the Workbench environment |
| `01_build_mdd_cohort.R` | Build the WGS-eligible EHR case/control cohort and cohort flow |
| `02_build_demographics.R` | Build age and sex-at-birth covariates |
| `03_build_pss_eds.R` | Build prorated PSS-10 (8+ items) and EDS-9 (7+ items) scores |
| `04_build_ace.R` | Build the prorated ACE-11 score for participants with 8+ scored items |
| `05_build_chronic_conditions.R` | Build the five-condition chronic disease count |
| `06_import_genomics.R` | Import PRS, ancestry, and PCs |
| `07_build_analytic_data.R` | Combine all components by `person_id` |
| `08_manuscript_analysis_template.R` | Main and supplementary model template |
| `99_build_data.R` | Run scripts 00 through 07 in order |

## Before starting

Each researcher needs:

1. Authorized access to an All of Us Controlled Tier workspace.
2. CDR `C2025Q4R6` attached to that workspace.
3. Permission to create an R Analysis Environment for AoU app.

This repository is public. No GitHub login, repository invitation, or SSH key
is required.

## Step 1: Attach the repository to the workspace

Before creating the R app, open the intended Workbench workspace and select
**Resources** → **Git repositories**.

1. Select **Add repository**.
2. Use the resource name `bridge_mdd_analysis`.
3. Enter this public repository URL:

```text
https://github.com/atalaydemiray/aou-bridge-mdd-analysis.git
```

4. Select **Add**.

After the R app is created and opened, the repository should appear at:

```text
~/repos/bridge_mdd_analysis
```

Confirm this in the RStudio Terminal:

```sh
ls ~/repos/bridge_mdd_analysis
```

If the app was already running when the repository was added, mount the new
resource and check again:

```sh
wb resource mount
ls ~/repos/bridge_mdd_analysis
```

If the repository still does not appear, clone the public repository directly
from the RStudio Terminal:

```sh
mkdir -p ~/repos
cd ~/repos
git clone https://github.com/atalaydemiray/aou-bridge-mdd-analysis.git bridge_mdd_analysis
```

## Step 2: Create and open the R app

In the workspace, select **Apps**, then **New app instance**.

1. Choose **R Analysis Environment for AoU**.
2. Give the app a recognizable name.
3. A configuration with 4 CPUs, approximately 26 GB memory, and a 50 GB disk
   was sufficient for the validation run.
4. Enable an automatic stop time, such as four hours.
5. Create the app, wait until its status is **Running**, and launch it.

The app opens RStudio in the browser. RStudio has two different command areas:

- **Terminal** accepts commands such as `cd`, `git pull`, `ls`, and `Rscript`.
- **Console** is an interactive R session.

Use the **Terminal** for Git commands and the **Console** for the R analysis.
The preflight script automatically loads the two required Workbench variables
when the Console does not inherit them. It imports only the CDR dataset
identifier and billing project, not the complete shell environment.

## Step 3: Update the repository in the Terminal

Open the **Terminal** tab in RStudio and run:

```sh
cd ~/repos/bridge_mdd_analysis
git pull --ff-only
pwd
git rev-parse --short HEAD
```

`pwd` should end in `/repos/bridge_mdd_analysis`. If `cd` fails, return to
Step 1 and confirm that the repository is attached or cloned.

## Step 4: Run scripts 00 through 05 in the R Console

Select the **Console** tab. Set the repository as the R working directory:

```r
setwd("~/repos/bridge_mdd_analysis")
```

Then run each command separately and wait for its completion message and a new
`>` prompt before continuing:

```r
source("analysis/00_preflight.R")
source("analysis/01_build_mdd_cohort.R")
source("analysis/02_build_demographics.R")
source("analysis/03_build_pss_eds.R")
source("analysis/04_build_ace.R")
source("analysis/05_build_chronic_conditions.R")
```

Opening a numbered script in the Source pane and selecting **Source** is
equivalent, provided the working directory is the repository root. Scripts 01
through 08 automatically source preflight when it has not already passed in
the current Console session.

The first run may install several R packages. Allow installation to finish.
Messages about packages masking functions are expected and are not failures.

Script 01 prints the complete cohort flow and the final case/control counts.
Record the CDR version and Git commit with those results. Counts generated
with an earlier phenotype implementation should not be used as a benchmark
for the current complete PhecodeX case definition.

Scripts 02 through 05 intentionally build source-wide, one-row-per-person
components. Their raw row counts can exceed the final analytical cohort
because script 07 later anchors every join to the MDD cohort.

## Step 5: Check the outputs

All participant-level products remain inside Workbench. By default, scripts
write to:

```text
~/aou_bridge_mdd_work
```

Check the output files from the Console:

```r
list.files("~/aou_bridge_mdd_work")
readRDS("~/aou_bridge_mdd_work/01_mdd_cohort_flow.rds")
```

After scripts 01 through 05, the folder should contain:

```text
01_mdd_cohort.rds
01_mdd_cohort_flow.rds
02_demographics.rds
03_pss_eds.rds
04_ace.rds
05_chronic_disease_count.rds
Figure_cohort_flow.pdf
```

Do not download, commit, or move participant-level RDS files outside the
authorized Workbench.

## Step 6: Add genomics when the approved files are available

Script 01 uses the release-matched WGS availability flag for cohort
eligibility. That flag confirms that WGS data are available; it is not a PRS,
genomic ancestry assignment, or principal-component result. Script 06 does
not calculate those genomic measures. Before running it, obtain two approved
files inside the authorized Workbench:

1. A one-row-per-participant genomic-results CSV or RDS file containing
   `research_id`, `prs_trans_ancestry`, `genomic_ancestry`, and `PC1` through
   `PC10`.
2. A one-to-one CSV or RDS ID-bridge file containing `research_id` and the
   corresponding AoU `person_id`.

The genomic-results file must contain the finalized trans-ancestry MDD PRS,
genomic ancestry assignment, and ancestry principal components produced by
the genomic analysis. The ID bridge is needed because the genomic research
identifier is not assumed to be the AoU `person_id`. The two identifiers must
map one-to-one.

Keep both files inside the authorized Workbench. Do not copy them into this
repository or place their paths, identifiers, or bucket names in Git.
Script 06 reads local CSV or RDS files. If an approved file is supplied as a
`gs://` URI, first use `gsutil cp` in the Workbench Terminal to copy it into a
directory outside the repository, such as `~/aou_bridge_mdd_inputs/`.

In the R Console, set the two file paths for the current session:

```r
Sys.setenv(
  BRIDGE_MDD_GENOMICS_PATH =
    "/path/inside/authorized/workbench/genomics_results.csv",
  BRIDGE_MDD_ID_MAP_PATH =
    "/path/inside/authorized/workbench/genomic_id_bridge.csv"
)
```

These are placeholders. Replace them with the approved local Workbench paths.
Both inputs may be CSV or RDS files. `Sys.setenv()` changes only the active R
session; it does not copy a file or write a path into Git. Repeat this command
after restarting the R app or Console.

The default column names are:

| File | Required column |
|---|---|
| Genomic results | `research_id` |
| ID bridge | `research_id` |
| ID bridge | `person_id` |

If the supplied identifier columns have different names, set the relevant
optional variables in the same R Console:

```r
Sys.setenv(
  BRIDGE_MDD_GENOMICS_RESEARCH_ID = "genomic_file_id_column",
  BRIDGE_MDD_MAP_RESEARCH_ID = "id_bridge_research_id_column",
  BRIDGE_MDD_MAP_PERSON_ID = "id_bridge_person_id_column"
)
```

Then run:

```r
source("analysis/06_import_genomics.R")
source("analysis/07_build_analytic_data.R")
```

Script 06 verifies the required genomic columns, enforces one-to-one
identifier mapping, converts the approved results to `person_id`, and writes
`06_genomics.rds` inside the Workbench output folder. Script 07 left-joins
every component to the MDD cohort by `person_id` and preserves missing values.
Its completion message reports the fixed cohort size, the number matched to
the genomic component, and the number with a nonmissing trans-ancestry PRS.
Review those three numbers before running models.

After setting both required path variables, the complete data build can
instead be rerun with:

```r
source("analysis/99_build_data.R")
```

If the two required path variables are not set, script 99 ends successfully
after script 05 and explains that genomics is still pending. It does not run
script 07 or create the final analytic data file.

## Step 7: Run statistical models only when finalized

The model template is deliberately separate:

```r
source("analysis/08_manuscript_analysis_template.R")
```

Run it only after the final analytical specification and genomic component
are approved and scripts 06 and 07 have completed successfully.

## Step 8: Stop the app

When finished, return to the workspace **Apps** page and stop the app to avoid
unnecessary compute charges. Files on the app disk remain available when the
same app is restarted.

## Beginner beta-test checklist

The beta tester should confirm:

- the repository appears under `~/repos/bridge_mdd_analysis`;
- `00_preflight.R` reports the expected CDR and passes;
- each script from 01 through 05 reaches its completion message;
- script 01 prints a complete cohort flow and nonzero final case/control
  counts;
- all seven expected files appear in `~/aou_bridge_mdd_work`;
- no participant-level file is committed to Git; and
- the Workbench app is stopped after testing.

After the approved genomic files become available, also confirm:

- script 06 reports that the genomic component was imported and mapped;
- script 07 reports a nonzero genomic match and nonmissing PRS count;
- `06_genomics.rds`, `07_analytic_data.rds`, and `07_run_metadata.rds` appear
  in `~/aou_bridge_mdd_work`; and
- script 08 is not run until the imported genomic component and model
  specification are approved.

Record the Git commit shown by `git rev-parse --short HEAD` and report the
exact script and full error message if any step fails.

## Collaboration and disclosure

This repository contains manuscript code only. Participant-level All of Us
data, genomic identifiers, credentials, workspace paths, and restricted output
must remain inside an authorized Researcher Workbench.

To propose a change:

1. Create a branch:

```sh
git switch -c initials-short-description
```

2. Edit the relevant numbered script or documentation file.
3. Run `analysis/00_preflight.R` and the affected build script in Workbench.
4. Review `git status` and `git diff`.
5. Commit the specific reviewed files and open a pull request.

Do not use `git add .` without first reviewing every listed file.

- Do not commit participant-level output or counts that violate the current
  All of Us dissemination policy.
- Preserve missing values during joins; absence of a survey response is not
  zero exposure.
- Record the Git commit and CDR release for every manuscript run.

GitHub distributes shared code only. It does not transfer All of Us data or
collaborator genomic files between workspaces.
