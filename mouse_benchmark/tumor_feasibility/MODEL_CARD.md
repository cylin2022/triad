# Research model card: mouse tumor neoantigen feasibility v1

**Status: not promoted; no production integration.** Primary model:
`run/LR_tumor_primary.joblib`. Comparators: binding-only logistic regression and
exploratory gradient boosting. The existing netMHCpan/Calis pipeline is unchanged.

Intended research endpoint: IFN-gamma response to an SNV-derived mutant peptide
after vaccination in C57BL/6 mice. Supported in this experiment: H-2-Kb, H-2-Db,
canonical unmodified 8–11mers with a same-length WT counterpart differing by one
residue. No support claim for other alleles, CD4, indels/frameshifts, human tumors,
natural tumor recognition, treatment benefit, or prospective screening accuracy.

Training: 100 pairs, 73 positive, 27 negative, from Smith et al. 2019; four tumor
cell lines but only one publication. All source and outcome records are in the
locked tables. Validation: 10 E0771 pairs. Test: 69 LLC1/CMT93 pairs from two other
papers. Study, tumor-model, source identifier and mutant/WT similarity separation
are audited. Upstream netMHCpan training overlap and whole-protein-family homology
are not fully known. Prior broad benchmark exposure is explicitly acknowledged.

Architecture: L2 logistic regression, C=0.1, StandardScaler fitted only on training
data; 13 predeclared features from mutant/WT sequences and netMHCpan BA/EL output.
No study ID, tumor label, experimental binding, ELISpot score or original published
GBM score enters training. No hyperparameter search, class balancing or calibration.
`predict_proba` values must be treated as uncalibrated ranking scores, not biological
activation probabilities. No entropy/energy calculation is included.

Fixed test performance: primary AP 0.7567, ROC-AUC 0.6224, Precision@10 0.60.
Test prevalence is 0.7101. netMHCpan BA rank AP 0.8815; existing composite AP 0.8384.
The primary is worse by point estimate within both external studies. This does not
establish that sequence/mutation information is universally unhelpful; small data,
single-study training and vaccination-context shift limit the conclusion.

Use the report and per-study metrics, not the pooled AP alone. The LLC1 test is
85.7% positive, CMT93 35.0%; a high pooled AP can largely reflect study composition.
Conditional sequence-group bootstrap intervals are not confidence intervals for
performance across new laboratories or tumor models.

Decision: retain as a reproducible feasibility experiment. Add independently
curated tumor vaccination cohorts and plan a new external holdout before future
model selection. The current test has been evaluated and must not be reused as an
untouched test for subsequent tuning.
