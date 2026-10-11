# Mouse Tumor Neoantigen Expansion Benchmark (v2): Protocol

**Date:** 2026-10-10 Asia/Taipei  
**Scope:** Expanded mouse tumor neoantigen immunogenicity benchmark and model re-evaluation.  
**Constraint:** Strict data isolation; zero test leakage; no production integration.

## 1. Scientific Background & Motivation

In feasibility study v1 (`mouse_benchmark/tumor_feasibility/`), machine learning models (L2 Logistic Regression and Gradient Boosting) trained on 100 neoantigen pairs from Smith et al. (2019, PMID: 31515258) failed to outperform the baseline NetMHCpan BA rank on the independent holdout test set (LLC1 and CMT93, N=69, AP 0.7567 vs 0.8815).

A critical limitation identified in v1 was:
1. Training data came from a single publication (Smith et al. 2019) with a high positive prevalence (73.0% positive).
2. Sample size was small (100 training pairs; 27 negatives).
3. The evaluation cohort had varying prevalence across laboratories (LLC1 was 85.7% positive, CMT93 was 35.0% positive).

In this v2 expansion, we systematically incorporate newly curated, reliable mouse tumor neoantigen datasets to address two key research questions:
- **Question 1 (Holdout Generalization):** Does sequence-derived feature modeling improve ranking over NetMHCpan binding when evaluated on a large, prospective, systematically screened tumor neoantigen cohort with low prevalence?
- **Question 2 (Expanded Training):** Does training with expanded multi-study data and authentic negative neoepitopes improve out-of-study ranking performance on external cohorts?

## 2. Cohort Curation & Data Sources

All cohorts represent vaccine-induced CD8+ T-cell responses (IFN-gamma ELISpot / multimer) to SNV-derived mutant peptides with verified single amino acid substitution wild-type counterparts in mice.

| Cohort | PMID / PMC | Tumor Model | Assay / Platform | N (Kb/Db) | Pos / Neg | Role |
|---|---|---|---|---|---|---|
| **Smith 2019** | 31515258 | B16F10, BBN963, MB49, UPPL1541 | Peptide vaccine + IFNg ELISpot | 100 | 73 / 27 | Training (Core) |
| **Fernando 2019** | 31292145 | MC-38 | Neoepitope pool vaccine + individual 9mer IFNg ELISpot | 13 | 6 / 7 | Training (Expansion) |
| **Validation E0771** | 33879241 | E0771 | Polyepitope DNA/SLP vaccine + IFNg ELISpot | 10 | 7 / 3 | Validation |
| **Frozen Test 1** | 36950115, 40586324 | LLC1, CMT93 | BMDC / pooled peptide vaccine | 69 | 49 / 20 | External Holdout 1 (v1 benchmark) |
| **Capietto 2020 (C57BL/6)** | 31940002 / PMC7144530 | MC-38, TRAMP-C1 | SLP vaccine + poly(I:C) + anti-CD40, IFNg ELISpot / multimer | 224 (218 purged) | 22 / 202 (21 / 197 purged) | External Holdout 2 (Large prospective test) |
| **Capietto 2020 (BALB/c)** | 31940002 / PMC7144530 | CT-26, EMT-6 | SLP vaccine + poly(I:C) + anti-CD40, Kd/Ld/Dd | 185 | 18 / 167 | Cross-strain exploratory cohort |

### Inclusion & Exclusion Rules
- **Inclusion:** Canonical 8–11mers, standard 20 amino acids, exact single-residue substitution with WT counterpart, explicit H-2 MHC restriction, in vivo vaccination followed by ex vivo IFN-gamma restimulation.
- **Exclusion:** In vitro only priming, untreated spontaneous recognition, post-challenge tumor-only assays, viral escape variants, class II (I-Ab), non-canonical lengths, ambiguous alleles.
- **Overlap & Leakage Purge:** Any test candidate sharing mutant or WT sequence, or having normalized Levenshtein distance <= 0.20 to any training sequence, is strictly purged from the holdout split. 6 identical peptides shared between Fernando 2019 and Capietto MC-38 are purged from the Capietto holdout test split (224 -> 218 pairs).

## 3. Predeclared Models & Features

All features are computed identically to v1 without using labels or cohort IDs:
- **Binding Features (5):** `log_EL_rank`, `log_BA_rank`, `log_affinity`, `length`, `is_Kb`.
- **Extended Sequence & Mutation Features (8):**
  - `WT_minus_mutant_log_affinity`
  - `WT_minus_mutant_log_EL_rank`
  - `hydrophobic_fraction` (AVILMFWY)
  - `aromatic_fraction` (FWY)
  - `charge_fraction` ((KR - DE) / length)
  - `proline_fraction`
  - `mutation_position_fraction`
  - `hydrophobic_change`
- **Models Fitted:**
  - `LR_binding`: Pipeline(StandardScaler, LogisticRegression(C=0.1, max_iter=5000, random_state=20261010))
  - `LR_tumor_primary`: Pipeline(StandardScaler, LogisticRegression(C=0.1, max_iter=5000, random_state=20261010))
  - `GB_tumor_exploratory`: GradientBoostingClassifier(n_estimators=100, learning_rate=0.05, max_depth=2, min_samples_leaf=10, random_state=20261010)
  - `LR_tumor_expanded`: Fitted on combined training data (Smith + Fernando + Capietto MC-38) to evaluate whether large negative training cohorts improve generalization.
- **Baselines:** Calis immunogenicity score, NetMHCpan EL Rank Score (-Rank_EL), NetMHCpan BA Rank Score (-Rank_BA), Mouse Base, Mouse Composite.

## 4. Evaluation Metrics & Locking

- Metric reporting: Average Precision (AP), ROC-AUC, Precision@5, Precision@10, Precision@20, Recall@5/10/20.
- Stratified reporting: Per study (LLC1, CMT93, Capietto), per allele (H-2-Kb, H-2-Db), and pooled holdouts.
- Uncertainty estimation: Paired bootstrap (1,000 replicates) for AP differences between primary models and baselines.
- Lock mechanism: Source hashes, data split SHA-256 hashes, and code hashes locked in `split_lock.json` before training.
