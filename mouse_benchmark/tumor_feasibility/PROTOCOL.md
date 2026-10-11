# Tumor neoantigen feasibility: frozen analysis plan

Date: 2026-10-10 Asia/Taipei. Scope selected explicitly by the user: tumor
neoantigens. No production integration. No entropy or free-energy model.

## Biological endpoint and sources

Estimate IFN-gamma reactivity to a mutant peptide after immunization in standard
C57BL/6 mice. This is vaccine-induced reactivity, NOT spontaneous recognition of
untreated tumors, natural presentation, tumor rejection, or human immunogenicity.
Vaccine formulations differ; report each external study separately.

Use the existing IEDB snapshot with freshly retrieved `tcell_export` details.
Require a canonical 8–11mer, explicit Kb/Db, IFNg assay, no modification, standard
B6 host, no adoptive transfer, and actual assay antigen equal to mutant epitope.
The last rule is crucial: an epitope record can contain a NEGATIVE assay to its
wild-type counterpart. That negative is not a negative mutant response.
Require an identified, canonical same-length WT counterpart differing at exactly
one residue (SNV-focused v1); preserve mutant/WT, locus and assay IDs.
Exclude untreated, tumor-only, in-vitro-only and post-vaccination tumor-challenge
records. Conflicting outcomes within the retained same-study setting are excluded,
not majority-voted. No random negatives and no balancing/downsampling.

Manually reviewed primary studies (full texts and IEDB assay comments):

| PMID | Role | Tumor source | Evidence/setting |
|---|---|---|---|
|31515258 (2019)|train|B16F10, BBN963, MB49, UPPL1541|Smith supplementary table supplies mutant/WT and cell line; peptide vaccination; exclude P815 mHA, T11 and class II|
|33879241 (2021)|validation|E0771|polyepitope DNA/SLP immunization; retain pre-challenge assays only|
|36950115 (2023)|locked test|LLC1|neoantigen-pulsed BMDC vaccination with in-vitro restimulation|
|40586324 (2025)|locked test|CMT93|pooled peptide immunization; mutant responses only, remove WT challenge controls|

Other studies remain a review queue, not implicitly validated training data.
The protocol selects these studies for verified mutation identity, vaccination
endpoint and distinct tumor models, not according to model performance.
Allele labels are as reported and may be prediction-derived. No claim of proven
allele-specific T-cell restriction is made. The available Kb/Db snapshot omits
other alleles and lengths; this is not a complete neoantigen literature census.

## Isolation and locking, before fitting

Hold out all retained 2023/2025 study records. Remove development records sharing
any source-protein identifier/name or mutant/WT sequence with a holdout, or having
normalized Levenshtein distance <=0.20 to either holdout mutant OR WT sequence.
Apply the same purge to training versus retained validation. Do not remove test
examples to improve class balance. Check zero cross-split study and tumor-model
overlap as well. Protein identifiers are incomplete; this is not whole-protein
family separation or proof of independence from upstream netMHCpan training.

Write train/validation/test tables and their SHA-256, source checksums, source-code
checksums and UTC timestamp to `split_lock.json` before training. Refuse to
overwrite a lock. This is an audited research holdout, not a security boundary:
labels were inspected for curation, and some records were used in the preceding
baseline analysis. It is independent of fitting the NEW models, not a previously
unseen prospective external cohort.

## Frozen models and features

No hyperparameter search or test-driven revisions. Seed 20261010.
Use local netMHCpan 4.2 BA/EL predictions for BOTH mutant and WT, without labels.
Primary model: StandardScaler + L2 LogisticRegression(C=0.1, max_iter=5000).
Binding-only logistic comparator uses log10(1+EL rank), log10(1+BA rank),
log10(clipped affinity), length and Kb indicator.
Primary adds WT-minus-mutant log affinity, WT-minus-mutant log EL rank,
mutant hydrophobic fraction (AVILMFWY), aromatic fraction (FWY), signed charged
fraction ((KR)-(DE))/length, proline fraction, normalized mutation position,
and mutation hydrophobicity-indicator change. These are descriptive features,
not energies. Scaling is fitted only on training rows.
Exploratory GradientBoostingClassifier(n_estimators=100, learning_rate=0.05,
max_depth=2, min_samples_leaf=10, subsample=1) uses the same extended features.
No class weights; no calibrated-probability claims. Predict_proba outputs are
research ranking scores. Do not refit on validation or select a different primary
model after inspecting results. Validation is a diagnostic only.

Frozen baselines: Calis, negative mutant EL/BA rank, current mouse composite and
mouse base. Training features exclude reference/study/cell-line identifiers,
experimental binding measurements, assay outcomes, original GBM predictions and
publication-derived ELISpot scores.

## Evaluation and decision

Evaluate validation first, then test once, after serialized model artifacts and
their hashes are written. Report AP, ROC-AUC, expected tie-aware Precision@5/10/20,
per allele, per held-out study and pooled test. Always show prevalence and N.
Single-class AUC/AP are unavailable, never coerced to a favorable value.
Bootstrap by mutant/WT similarity/source groups WITHIN each held-out study (1000
paired replicates) for primary-minus-base AP differences. This gives conditional
candidate-sampling uncertainty, NOT uncertainty over new laboratories/studies;
two external studies cannot characterize the latter reliably. Treat all confidence
intervals as exploratory, especially if groups or positives/negatives are few.

Run fixed leave-one-tumor-model-out diagnostics within training with cross-fold
sequence/source purging; no tuning and no resplitting the locked external test.
If a training fold is single-class, report unavailable. A viable deployment claim
would require reproducible improvement over strong baselines within both studies
and both supported alleles, adequate independent data, calibration and further
external validation. This run is a feasibility experiment only. Failure to improve
is a valid result; do not iterate on the same test set to manufacture success.
