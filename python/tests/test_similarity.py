"""Port of tests/test_similarity.cpp."""

from __future__ import annotations

import unittest

from _support import LLKATestCase, make_structure

import llka

EXPECTED_TORSIONS = {
    llka.DinucleotideTorsion.DELTA_1: 2.11367168175336,
    llka.DinucleotideTorsion.EPSILON_1: 3.03091618163782,
    llka.DinucleotideTorsion.ZETA_1: -1.54423502253339,
    llka.DinucleotideTorsion.ALPHA_2: -0.98817404111552,
    llka.DinucleotideTorsion.BETA_2: -3.12713306928336,
    llka.DinucleotideTorsion.GAMMA_2: 0.91145018867645,
    llka.DinucleotideTorsion.DELTA_2: 1.72629995189555,
    llka.DinucleotideTorsion.CHI_1: -2.13194310890551,
    llka.DinucleotideTorsion.CHI_2: -2.22249417305504,
}


class TestSimilarity(LLKATestCase):
    def setUp(self) -> None:
        self.stru = make_structure('REAL_1BNA_A_6_7_ATOMS')

    def assertMetrics(self, metrics: llka.StepMetrics, expected: dict[str, float]) -> None:
        for name, value in expected.items():
            self.assertFlt(getattr(metrics, name), value, name)

    def test_should_fail(self) -> None:
        dinucleotide = llka.extract_extended_backbone(self.stru)
        with self.assertRaises(llka.InvalidArgumentError):
            llka.dinucleotide_torsion(llka.DinucleotideTorsion.CHI_1, dinucleotide)

    def test_all_metrics(self) -> None:
        all_metrics = llka.extract_metrics_structure(self.stru)

        for torsion, expected in EXPECTED_TORSIONS.items():
            with self.subTest(torsion=llka.dinucleotide_torsion_name(torsion)):
                atoms = llka.dinucleotide_torsion(torsion, all_metrics)
                self.assertEqual(len(atoms), 4)
                self.assertFlt(llka.measure_dihedral(*atoms), expected, llka.dinucleotide_torsion_name(torsion))

        cc = llka.cross_residue_metric(llka.CrossResidueMetric.DIST_CC, all_metrics)
        self.assertFlt(llka.measure_distance(cc[0], cc[1]), 4.78625636588764, 'CC')

        nn = llka.cross_residue_metric(llka.CrossResidueMetric.DIST_NN, all_metrics)
        self.assertFlt(llka.measure_distance(nn[0], nn[1]), 4.22325786567669, 'NN')

        mu = llka.cross_residue_metric(llka.CrossResidueMetric.TOR_MU, all_metrics)
        self.assertFlt(llka.measure_dihedral(*mu), 0.45071211146595, 'mu')

    def test_all_metrics_2(self) -> None:
        metrics = llka.calculate_step_metrics(make_structure('REAL_1BNA_B_16_17_ATOMS'))
        self.assertMetrics(metrics, {
            'delta_1': 2.37223786208374,
            'epsilon_1': 3.03929251506922,
            'zeta_1': -1.71720034870345,
            'alpha_2': -0.98818512153652,
            'beta_2': -2.95903498082587,
            'gamma_2': 0.93968862228876,
            'delta_2': 2.55870427626711,
            'chi_1': -2.00305820661614,
            'chi_2': -1.85775524738115,
            'CC': 5.04421341737243,
            'NN': 4.45444407305783,
            'mu': 0.58047587191140,
        })

    def test_all_metrics_difference_against_reference(self) -> None:
        metrics = llka.calculate_step_metrics_difference_against_reference(
            make_structure('REAL_1BNA_B_16_17_ATOMS'), llka.NtC.BB00
        )
        self.assertMetrics(metrics, {
            'delta_1': -0.033174913014841,
            'epsilon_1': -0.156579878257593,
            'zeta_1': 0.059893895677181,
            'alpha_2': -0.005913818514107,
            'beta_2': 0.189713522697103,
            'gamma_2': 0.167729494131665,
            'delta_2': 0.147706447562145,
            'chi_1': -0.127701925348230,
            'chi_2': -0.081882733476915,
            'CC': 0.097213417372426,
            'NN': 0.081444073057826,
            'mu': 0.134718780952045,
        })

    def test_similarity(self) -> None:
        stru = make_structure('REAL_1BNA_A_3_4_ATOMS')
        similar = llka.measure_step_similarity_ntc(stru, llka.NtC.AB01)
        self.assertFlt(similar.rmsd, 0.176588293714412, 'wrong RMSD')
        self.assertFlt(similar.euclidean_distance, 25.0592315110376, 'wrong euclidean distance')

    def test_similarity_multiple(self) -> None:
        stru = make_structure('REAL_1BNA_A_3_4_ATOMS')
        results = llka.measure_step_similarity_ntc_multiple(stru, [llka.NtC.AB01, llka.NtC.AB03, llka.NtC.BB00])

        self.assertEqual(len(results), 3)
        self.assertFlt(results[0].rmsd, 0.176588293714412, 'wrong RMSD')
        self.assertFlt(results[0].euclidean_distance, 25.0592315110376, 'wrong euclidean distance')
        self.assertFlt(results[1].rmsd, 0.24827518856867, 'wrong RMSD')
        self.assertFlt(results[1].euclidean_distance, 46.3910103838844, 'wrong euclidean distance')
        self.assertFlt(results[2].rmsd, 0.348261210522226, 'wrong RMSD')
        self.assertFlt(results[2].euclidean_distance, 56.5297526013971, 'wrong euclidean distance')

    def test_similarity_structure(self) -> None:
        stru = make_structure('REAL_1BNA_A_3_4_ATOMS')
        similar = llka.measure_step_similarity_structure(stru, stru)
        self.assertFlt(similar.euclidean_distance, 0.0, 'step is not similar to itself')


if __name__ == '__main__':
    unittest.main()
