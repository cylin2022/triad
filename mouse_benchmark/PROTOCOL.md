# Mouse immunogenicity retrospective benchmark protocol

This is repository documentation for a computational benchmark, not a claim of
prospective or independent experimental validation. Fixed before inspecting model
performance. Scope: H-2-Kb / H-2-Db; other alleles require separate evaluation.

## Data and endpoint

Download all IEDB T-cell records with mouse host ancestry (NCBITaxon:10090) and
explicit H2-Kb or H2-Db annotation. Retain original assay IDs, reference IDs, host,
restriction evidence, antigen, immunization description, dates, and endpoint.
Check pagination completeness, unique IDs, and SHA-256. Do not balance classes.

Primary cohort: C57BL/6, C57BL/6J or C57BL/6N hosts explicitly annotated in IEDB;
unmodified canonical linear peptides, lengths 8–11; IFNg release assays; known
reference ID; Positive (including graded positives) or Negative outcomes.
Tetramer binding, MHC binding alone and tumor control are not endpoint labels.
Explicit transgenic/knockout hosts are excluded by the host allowlist; experiments
using transferred T cells or other interventions can still occur in otherwise
normal hosts. These require subsequent manual curation before model training.

Collapse eligible assays to unique peptide:allele pairs. Exclude a pair entirely
if eligible records contain both positive and negative labels. A negative means
negative under reported experimental conditions, not universally nonimmunogenic.
Retain all supporting references; do not infer an allele from mouse strain.
IEDB allele annotations based on predictions are inventoried, not treated as
experimentally proven restriction. Primary analysis includes these annotations;
a sensitivity cohort requires binding-assay, single-allele or T-cell restriction
evidence. This stronger-evidence cohort is not proof of TCR-specific restriction.

## Frozen baselines

No fitting, score-direction selection, threshold optimization, or new model
training. Use the installed netMHCpan 4.2 without pathogen/neo context, with BA.
Run all eligible unique peptides against both Kb and Db; join using peptide AND
allele, requiring complete, unique and finite predictions.

1. Current production Calis implementation, unchanged (higher is better).
2. netMHCpan negative EL percentile rank (higher is better).
3. netMHCpan negative BA percentile rank (higher is better).
4. netMHCpan EL score (higher is better).
5. Mouse presentation/binding base used by the platform, without Calis:
   `s_pres = clip(.5*Score_EL + .5*clip(1-Rank_EL/2,0,1),0,1)`;
   `s_aff = clip(1-log10(clip(Affinity_nM,1,50000))/log10(50000),0,1)`;
   `base = .55*s_pres + .45*s_aff`.
6. Current mouse composite: `round(clip(base*clip(1+.35*Calis,.5,1.5),0,1),4)`.

The current Calis implementation uses its existing length and anchor behavior;
9-mer-only sensitivity analysis tests a narrower domain. Scores are not calibrated
probabilities. Source and binary hashes are recorded to identify evaluated code.

## Evaluation

Primary metrics: average precision (AP; stepwise PR area, not trapezoidal PR-AUC),
ROC-AUC, tie-aware expected precision/recall at K=20 and K=50. Random-ranking AP
reference is prevalence. Boundary ties receive proportional credit, not favorable
input ordering. No fake random negatives. No random train/test split.

Report per allele; 9-mers; stronger restriction evidence; EL rank <=2 binders;
post-2019 references only with no eligible pre-2020 observation of that pair;
and removal of the largest contributing study. Also evaluate removal of studies
whose titles explicitly mention cross-reactivity/alloreactivity, removal of the
largest reference/sequence component, and stronger evidence within each allele.
These additional audits were chosen after inspecting data composition but before
inspecting model scores. Title screening is not comprehensive manual curation.
The date analysis is a sensitivity
analysis, not guaranteed model-independent validation. Audit pooled study and
antigen-source imbalance. Report within-study metrics for studies containing at
least five positive and five negative pairs, and unweighted macro summaries.

For the primary cohort, construct connected components joining records sharing
any reference or peptide, including equal-length peptides differing by <=1 residue
(at least 87.5% identity for 8–11-mers). Use 1,000 paired component-bootstrap
replicates, seed 20261009, for AP intervals and composite-minus-base AP difference.
Skip single-class replicates and report valid counts. Report group concentration:
intervals with few effective groups / a dominant component are unstable, and do
not establish generalization. This graph also inventories whether leakage-resistant
future training splits are feasible. This does not capture all protein-family
homology or relationships among different publications from the same experiment.
No train/test split is claimed in this run.

## Decision

Assess added value of Calis versus the unmodulated base, across studies, alleles
and sensitivity cohorts; do not declare a mouse model necessary solely from a low
pooled AUC. If scores are weak or inconsistent, recommend a staged mouse-specific
training experiment only after context/annotation cleanup and locked external
validation design. If evidence is insufficient, say so. Public records may overlap
Calis/netMHCpan training data; full training overlap is unknown. This benchmark
cannot certify independent external validation or prospective mouse accuracy.

The published `vincentlaboratories/neoag` model requires mutant/reference peptide
pairs and mutation positions, absent from the general IEDB extract. Do not invent
these inputs or apply that SNV-specific model to arbitrary microbial peptides.
