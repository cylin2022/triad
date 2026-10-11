"""Regression and data integrity tests for v2 tumor neoantigen expansion."""
import json
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parent
from prepare import sha, overlap, close
from train import verify_lock, feature_frame, FEATURES, BINDING


class ExpansionV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lock = verify_lock()
        cls.df_train = pd.read_csv(ROOT / 'locked/train.csv')
        cls.df_val = pd.read_csv(ROOT / 'locked/validation.csv')
        cls.df_prior = pd.read_csv(ROOT / 'locked/test_frozen_prior.csv')
        cls.df_cap = pd.read_csv(ROOT / 'locked/test_capietto_holdout.csv')
        cls.scored = pd.read_csv(ROOT / 'run/scored_cohort.csv')
        cls.metrics = pd.read_csv(ROOT / 'run/external_metrics.csv')

    def test_lock_integrity(self):
        """All files in locked/ must match split_lock.json exactly."""
        for name, expected in self.lock['files'].items():
            self.assertEqual(sha(ROOT / 'locked' / name), expected, f'Hash mismatch on {name}')

    def test_zero_leakage_between_train_and_holdouts(self):
        """Zero sequence or Levenshtein distance <= 0.20 between train and holdouts."""
        train_rows = self.df_train.to_dict('records')
        for holdout_name, holdout_df in [('test_frozen_prior', self.df_prior), ('test_capietto_holdout', self.df_cap)]:
            hits = [
                (x['Peptide'], y['Peptide'], overlap(x, y))
                for x in train_rows
                for y in holdout_df.to_dict('records')
                if overlap(x, y)
            ]
            self.assertEqual(len(hits), 0, f'Leakage detected between train and {holdout_name}: {hits}')

    def test_purged_overlap_tracked(self):
        """Purged CSV must record exactly the purged peptides from Capietto."""
        purged = pd.read_csv(ROOT / 'locked/purged_capietto_from_holdout.csv')
        self.assertGreater(len(purged), 0)
        # All purged peptides must not appear in test_capietto_holdout
        self.assertEqual(len(set(purged.Peptide) & set(self.df_cap.Peptide)), 0)

    def test_features_independent_of_labels(self):
        """Feature extraction must not access labels or cohort information."""
        sample = self.df_train.head(5).copy()
        sample['Rank_EL'] = 0.5
        sample['Rank_BA'] = 0.8
        sample['Affinity_nM'] = 120.0
        sample['WT_Rank_EL'] = 1.5
        sample['WT_Rank_BA'] = 2.0
        sample['WT_Affinity_nM'] = 800.0

        sample['Label'] = 1
        f1 = feature_frame(sample)
        sample['Label'] = 0
        f2 = feature_frame(sample)
        np.testing.assert_array_equal(f1.to_numpy(), f2.to_numpy())
        self.assertEqual(list(f1.columns), FEATURES)

    def test_scaler_fitted_only_on_training(self):
        """StandardScaler must reflect only training data distribution."""
        import joblib
        model = joblib.load(ROOT / 'run/LR_tumor_primary.joblib')
        scaler = model.named_steps['standardscaler']
        train_features = feature_frame(self.scored[self.scored.Split == 'train'])
        np.testing.assert_allclose(scaler.mean_, train_features.mean().to_numpy(), rtol=1e-5)

    def test_metrics_reproducibility(self):
        """External metrics in CSV must match freshly computed AP on scored cohort."""
        prior_part = self.scored[self.scored.Split == 'test_frozen_prior']
        y_prior = prior_part.Label.to_numpy()
        s_ba = prior_part.BA_Rank_Score.to_numpy()
        computed_ap = average_precision_score(y_prior, s_ba)
        
        row = self.metrics[
            (self.metrics.Split == 'test_frozen_prior') &
            (self.metrics.Cohort == 'pooled') &
            (self.metrics.Model == 'BA_Rank_Score')
        ]
        self.assertEqual(len(row), 1)
        self.assertAlmostEqual(computed_ap, row.iloc[0]['AP'], places=4)


if __name__ == '__main__':
    unittest.main()
