"""Curation, isolation and no-label-feature regression checks."""
import json
import unittest
import numpy as np
import pandas as pd
from prepare import ROOT, close, overlap, select_record
from train import feature_frame, verify_lock, FEATURES


class TumorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw={r['tcell_id']:r for r in [json.loads(l) for l in (ROOT.parent/'data/iedb_mouse_kb_db.jsonl').read_text().splitlines()]}
        cls.details={r['assay_id']:r for r in json.loads((ROOT/'sources/iedb_tumor_details.json').read_text())}
        cls.smith=pd.read_excel(ROOT/'sources/smith_supp.xlsx')

    def test_wild_type_negative_is_not_mutant_negative(self):
        good,item=select_record(self.raw[23044448],self.details[23044448],self.smith)
        bad,_=select_record(self.raw[23044468],self.details[23044468],self.smith)
        self.assertEqual(good,'eligible');self.assertEqual(item['Label'],1)
        self.assertEqual(bad,'assay_tests_WT_or_other_antigen')

    def test_untreated_is_not_vaccinated_negative(self):
        reason,_=select_record(self.raw[16007693],self.details[16007693],self.smith)
        self.assertEqual(reason,'no_vaccination')

    def test_tumor_challenge_kept_out_of_vaccine_endpoint(self):
        reason,_=select_record(self.raw[16007666],self.details[16007666],self.smith)
        self.assertEqual(reason,'second_in_vivo_process')

    def test_similarity_includes_indels_and_wild_type(self):
        self.assertTrue(close('SIINFEKL','SIINFEKLA'))
        self.assertTrue(close('SIINFEKL','SIINFEKM'))
        self.assertFalse(close('AAAAAAAA','WWWWWWWW'))
        a={'Peptide':'SIINFEKL','WT':'SIINFEKM','Source_Tokens':'x'}
        b={'Peptide':'SIINFEKV','WT':'SIINFEKM','Source_Tokens':'y'}
        self.assertIsNotNone(overlap(a,b))
        b={'Peptide':'AAAAAAAA','WT':'AAAAAAAV','Source_Tokens':'x'}
        self.assertEqual(overlap(a,b),'same_source_protein')

    def test_actual_locked_isolation(self):
        verify_lock()
        splits={s:pd.read_csv(ROOT/'locked'/(s+'.csv')).to_dict('records') for s in ['train','validation','test']}
        for a,b in [('train','test'),('train','validation'),('validation','test')]:
            self.assertFalse(any(overlap(x,y) for x in splits[a] for y in splits[b]))
            self.assertFalse({r['PMID'] for r in splits[a]} & {r['PMID'] for r in splits[b]})
            self.assertFalse({r['Tumor_Model'] for r in splits[a]} & {r['Tumor_Model'] for r in splits[b]})

    def test_detail_identity_matches_snapshot(self):
        for assay,d in self.details.items():
            r=self.raw[assay]
            self.assertEqual(d['epitope__name'],r['linear_sequence'])
            self.assertEqual(d['mhc_restriction__name'],r['mhc_allele_name'])
            self.assertEqual(d['host__name'],r['host_organism_name'])

    def test_feature_values_do_not_depend_on_labels_or_study(self):
        df=pd.DataFrame({'Peptide':['SIINFEKL'],'WT':['SIINFEKM'],'MHC':['H-2-Kb'],
                         'Rank_EL':[.1],'Rank_BA':[.2],'Affinity_nM':[50.],
                         'WT_Rank_EL':[.4],'WT_Affinity_nM':[500.], 'Label':[1],'PMID':['train']})
        before=feature_frame(df)
        df['Label']=0;df['PMID']='test';after=feature_frame(df)
        np.testing.assert_array_equal(before,after)
        self.assertEqual(list(before),FEATURES)
        self.assertAlmostEqual(before.WT_minus_mutant_log_affinity.iloc[0],1.)

    def test_scaler_fit_only_training(self):
        import joblib
        if not (ROOT/'run/model_manifest.json').exists():self.skipTest('Training not yet completed')
        scored=pd.read_csv(ROOT/'run/scored_cohort.csv')
        model=joblib.load(ROOT/'run/LR_tumor_primary.joblib')
        np.testing.assert_allclose(model.named_steps['standardscaler'].mean_,feature_frame(scored[scored.Split=='train']).mean().to_numpy())


if __name__=='__main__':unittest.main()
