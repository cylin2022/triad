# Research model card: mouse tumor neoantigen expansion v2

**Status: not promoted; no production integration.** Primary model:
`run/LR_tumor_primary.joblib`. Comparators: `LR_binding.joblib`, exploratory gradient boosting `GB_tumor_exploratory.joblib`, and expanded cohort model `LR_tumor_expanded.joblib`. The existing production NetMHCpan/Calis composite pipeline remains unchanged.

Intended research endpoint: IFN-gamma response to an SNV-derived mutant peptide after vaccination in standard C57BL/6 mice. Supported in this benchmark: H-2-Kb, H-2-Db, canonical unmodified 8–11mers with a same-length wild-type counterpart differing by exactly one amino acid residue. No claim is made for spontaneous antitumor recognition, whole-tumor rejection efficacy, human neoantigens, indels/frameshifts, or other mouse MHC haplotypes.

## 1. Data Provenance & Split Architecture

- **Training (Core + Expansion):** 113 pairs (79 positive, 34 negative; 69.9% prevalence). Sources: Smith et al. 2019 (PMID: 31515258, 100 pairs across B16F10, BBN963, MB49, UPPL1541) and Fernando et al. 2019 (PMID: 31292145, 13 pairs from MC-38 peptide vaccination).
- **Validation:** 13 pairs (3 positive, 10 negative) from E0771 (PMID: 33879241).
- **Holdout Test 1 (Prior Benchmark):** 69 pairs (49 positive, 20 negative; 71.0% prevalence) from LLC1 (PMID: 36950115, N=49) and CMT93 (PMID: 40586324, N=20). Evaluated under frozen lock.
- **Holdout Test 2 (New Prospective Screening Benchmark):** 216 pairs (22 positive, 194 negative; 10.2% prevalence) from Capietto et al. 2020 (JEM, PMC7144530 / PMID: 31940002) in MC-38 and TRAMP-C1. Zero sequence overlap, zero Levenshtein distance <= 0.20 to training (6 identical peptides shared with Fernando 2019 were strictly purged).
- **Exploratory Cohort (Cross-Strain):** 185 pairs (18 positive, 167 negative) from BALB/c (CT-26, EMT-6, H-2-Kd, H-2-Ld, H-2-Dd).

## 2. Feature Architecture

All features are computed with zero cohort/label leakage:
- 5 Binding features: `log_EL_rank`, `log_BA_rank`, `log_affinity`, `length`, `is_Kb`.
- 8 Extended sequence/mutation features: `WT_minus_mutant_log_affinity`, `WT_minus_mutant_log_EL_rank`, `hydrophobic_fraction`, `aromatic_fraction`, `charge_fraction`, `proline_fraction`, `mutation_position_fraction`, `hydrophobic_change`.
- Standard scaling fitted exclusively on training rows. No class weighting, hyperparameter tuning, or calibration fitting.

## 3. Empirical Holdout Performance

| Cohort | Metric | NetMHCpan BA Rank | NetMHCpan EL Rank | Mouse Composite | LR Binding | LR Tumor Primary (v2) | LR Tumor Expanded |
|---|---|---:|---:|---:|---:|---:|---:|
| **Prior Holdout** (N=69, Prev=0.710) | **AP ↑** | **0.8815** | 0.8273 | 0.8384 | 0.8787 | 0.8094 | 0.8117 |
| | **ROC-AUC ↑** | 0.7628 | 0.6459 | 0.6694 | **0.7653** | 0.7204 | 0.6418 |
| | **Precision@10 ↑** | **0.90** | **0.90** | **0.90** | **0.90** | 0.60 | **0.90** |
| **Capietto Holdout** (N=216, Prev=0.102) | **AP ↑** | **0.2999** | 0.2287 | 0.2025 | 0.2447 | 0.2429 | 0.2694 |
| | **ROC-AUC ↑** | **0.6922** | 0.6059 | 0.6326 | 0.6860 | 0.6251 | 0.6371 |
| | **Precision@10 ↑** | **0.50** | 0.40 | 0.20 | 0.40 | 0.40 | 0.30 |
| **Combined Holdout** (N=285) | **AP ↑** | 0.3319 | 0.2664 | 0.2584 | 0.3986 | 0.3713 | **0.4607** |
| | **ROC-AUC ↑** | 0.5793 | 0.4777 | 0.4880 | **0.6999** | 0.6261 | 0.6869 |
| | **Precision@10 ↑** | 0.50 | 0.40 | 0.20 | 0.40 | 0.50 | **0.70** |

## 4. Scientific Insights & Decision

1. **Within each independent study, NetMHCpan BA Rank remains the single strongest point estimate predictor.** On both the prior holdout (AP 0.8815) and the new prospective low-prevalence Capietto holdout (AP 0.2999), BA rank outperforms models trained on empirical sequence features.
2. **Sequence/mutation features do not generalize across distinct vaccine platforms.** Features that appear descriptive in one study (e.g. SLP + poly(I:C) vs BMDC vs peptide pools) exhibit covariate shift and fail to provide robust gains in external laboratories.
3. **Deployment Decision:** **Do not promote `LR_tumor_primary` or `LR_tumor_expanded` to production.** The production system remains unchanged. The expanded dataset and benchmark artifacts are retained as an immutable, reproducible evaluation suite.
