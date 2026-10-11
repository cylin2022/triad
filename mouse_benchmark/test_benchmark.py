"""Tests for endpoint filtering, leakage prevention, metrics and production parity."""
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from curate import group_components, rejection, strong_evidence
from evaluate import metrics, top_k
from predict import add_scores


class BenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.record = {"host_organism_name": "Mus musculus C57BL/6", "linear_sequence": "SIINFEKL",
                       "structure_type": "Linear peptide", "e_modification": None,
                       "assay_names": "IFNg release|ELISPOT", "qualitative_measure": "Positive-Low", "reference_id": 1}

    def test_binding_is_not_functional_label(self):
        self.assertEqual(rejection(self.record), "eligible")
        self.record["assay_names"] = "qualitative binding|multimer/tetramer"
        self.assertEqual(rejection(self.record), "not_IFNg_functional_assay")

    def test_transgenic_and_modified_peptides_are_excluded(self):
        self.record["host_organism_name"] = "Mus musculus OT-1 TCR Tg"
        self.assertNotEqual(rejection(self.record), "eligible")
        self.record["host_organism_name"] = "Mus musculus C57BL/6"
        self.record["e_modification"] = "Phosphorylation|PHOS"
        self.assertEqual(rejection(self.record), "modified_peptide")

    def test_ambiguous_amino_acids_are_excluded(self):
        self.record["linear_sequence"] = "SIINFEKX"
        self.assertNotEqual(rejection(self.record), "eligible")

    def test_predicted_restriction_is_not_strong_evidence(self):
        self.assertFalse(strong_evidence("MHC binding prediction"))
        self.assertFalse(strong_evidence("Cited reference"))
        self.assertTrue(strong_evidence("T cell assay -Mismatched MHC molecules"))

    def test_transitive_study_and_sequence_leakage_groups(self):
        df = pd.DataFrame({"Peptide": ["SIINFEKL", "SIINFEKM", "AAAAAAAA", "WWWWWWWW", "SIINFEKL"],
                           "Reference_IDs": ["1", "2", "2;3", "4", "5"]})
        groups = group_components(df)
        self.assertEqual(groups[0], groups[1])  # one substitution
        self.assertEqual(groups[1], groups[2])  # shared study
        self.assertEqual(groups[0], groups[4])  # same peptide, other allele/study
        self.assertNotEqual(groups[0], groups[3])

    def test_boundary_ties_are_order_invariant(self):
        # Top one clear positive; remaining two tie, one positive, one negative.
        p, r = top_k([1, 0, 1], [3, 2, 2], 2)
        self.assertEqual(p, .75)
        self.assertEqual(r, .75)
        self.assertEqual((p, r), top_k([1, 1, 0], [2, 3, 2], 2))
        self.assertEqual(top_k([1], [1], 20), (None, None))

    def test_constant_scores_and_single_class(self):
        df = pd.DataFrame({"Label": [1, 0, 0, 0], "score": [1, 1, 1, 1]})
        self.assertEqual(metrics(df, "score")["AP"], .25)
        self.assertEqual(metrics(df, "score")["ROC_AUC"], .5)
        self.assertIsNone(metrics(df.iloc[:1], "score")["AP"])

    def test_composite_matches_actual_production_pipeline(self):
        from batch_runner import BatchJob
        from unittest.mock import Mock
        # Cover both alleles, clipping bounds, varied ranks and peptide lengths.
        data = pd.DataFrame({"Peptide": ["SIINFEKL", "ASNENMETM", "AAAAAAAAAA", "WWWWWWWWWWW"],
                             "MHC": ["H-2-Kb", "H-2-Db", "H-2-Kb", "H-2-Db"],
                             "Score_EL": [.9, .5, .0, .1], "Rank_EL": [.01, 1., 98., 3.],
                             "Rank_BA": [.1, 2, 50, 3], "Affinity_nM": [.5, 100, 60000, 1500]})
        expected = add_scores(data).set_index(["Peptide", "MHC"])
        job = BatchJob.__new__(BatchJob)
        job.alleles, job.meta, job.logger, job.job_id = ["H-2-Kb", "H-2-Db"], {}, Mock(), "parity"
        job._finalize(data.copy())
        actual = job.results.set_index(["Peptide", "MHC"]).loc[expected.index]
        np.testing.assert_array_equal(actual.Composite_Score, expected.Mouse_Composite)
        np.testing.assert_array_equal(actual.Immunogenicity_Score, expected.Calis)


if __name__ == "__main__":
    unittest.main()
