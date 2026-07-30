# Operational methods

## Study design

This is a cross-sectional prevalence analysis among adult All of Us
participants with EHR and genomic data. The outcome is lifetime
recorded MDD in available EHR data. Associations are not interpreted as
incident risk or causal effects.

Genomic-data eligibility is defined by short-read whole-genome sequencing
availability in the release-matched Data Explorer index field
`T_ENT_person.has_whole_genome_variant`.

The cohort also requires recorded participation in each of the Emotional
Health, Social Determinants of Health, and Basics survey modules, reproducing
the required survey-module availability gates. Module participation does not
imply that every survey item has a valid response; score-specific completeness
rules are applied subsequently.

## MDD cases and controls

The primary phenotype is derived from the official 37-code U.S.
ICD-9-CM/ICD-10-CM mapping for PhecodeX `MB_286.2`, major depressive
disorder. All 37 mapped source codes qualify as case evidence; the official
mapping is used without selectively removing individual source codes.

PhecodeX `MB_286.2` was selected instead of traditional Phecode v1.2
`296.22` because the newer mapping includes the ICD-9-CM mild-MDD codes
`296.21` and `296.31`, while assigning ICD-10-CM `F32.81` (premenstrual
dysphoric disorder) to the separate PhecodeX phenotype `MB_286.3`.
PhecodeX was developed through clinical curation and aligns its phenotype
structure more directly with ICD-10-CM.

A case has at least two distinct qualifying MDD diagnosis dates separated by
at least 30 days. There is no upper separation limit. A control has no
condition record under the full `MB_286.2` map and no available positive
response to the current-depression-treatment survey item.

The following standard OMOP ancestors and descendants are excluded from both
groups:

- 435783, Schizophrenia;
- 436665, Bipolar disorder;
- 4182210, Dementia;
- 4295956, Mood disorder due to a general medical condition;
- 37016268, Opioid-induced mood disorder due to opioid abuse; and
- 37018689, Opioid-induced mood disorder due to opioid dependence.

Both cases and controls require at least two distinct encounter dates and at
least 365 days between first and last recorded EHR activity.

The cohort-flow output reports sequential counts for all participants,
adulthood, EHR availability, WGS availability, the three required survey
modules, psychiatric exclusions, case/control classification, and the
one-year EHR observability requirement. Valid score and model-variable
completeness are reported separately for each model.

## Psychosocial exposures

PSS-10 uses the ten reviewed question concept IDs and reverses items
40192419, 40192525, 40192449, and 40192445. Responses are scored 0 through 4.
With 8-10 scored responses, the completed-item mean is multiplied by 10.

EDS-9 uses the nine reviewed question concept IDs. To reproduce the supplied
component, response concept 40192421 remains unscored. With 7-9 scored
responses, the completed-item mean is multiplied by 9.

The modified ACE score uses the 11 reviewed childhood-adversity question
concept IDs. Yes or any reported occurrence is scored 1, while No or Never is
scored 0. Other responses remain unscored; for the parental-separation item,
“Parents not married” also remains unscored because it does not establish
separation or divorce. With 8-11 scored responses, the completed-item mean is
multiplied by 11.

## Demographic, clinical, and genomic covariates

Primary adjustment uses age, sex at birth, and genomic principal components.

The clinical component contains diabetes, heart disease, hypertension,
chronic kidney disease, and chronic pain. Each condition requires two diagnosis
dates separated by at least 30 days. The five indicators are summed into
a chronic disease count ranging from 0 to 5. The count, rather than the five
individual conditions, is used as the clinical adjustment variable.

The genomic component supplies the trans-ancestry MDD PRS, genomic
ancestry, and PC1-PC10. PRS construction is not duplicated in this repository.

## Statistical analysis

The primary model is modified Poisson regression with a log link and robust
standard errors, producing adjusted prevalence ratios. Planned analyses are:

1. one core-adjusted model for each exposure;
2. a joint main-effects model with PRS, PSS, EDS, and ACE;
3. one PRS-by-psychosocial interaction at a time;
4. an extended model adding the 0-5 chronic disease count; and
5. ancestry-stratified and approved sensitivity analyses when the model
   specification is finalized.

Continuous exposures are standardized before modeling. Each model uses its
own complete-case sample, and the sample size must accompany its result.

## Data protection

Participant-level inputs and outputs remain in Researcher Workbench. Public
tables and figures must follow the current All of Us dissemination policy,
including suppression of restricted small cells and derived values.
