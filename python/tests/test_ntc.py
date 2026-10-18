"""Port of tests/test_ntc.cpp."""

from __future__ import annotations

import unittest

from _support import LLKATestCase, make_structure

import llka
from llka import CrossResidueMetric as XR
from llka import DinucleotideTorsion as Tor

# Expected atom names for a DC-DG step and an A-PSU step. Cross-residue metric
# names are the values that libLLKA currently returns. Lists of pairs are used
# because DinucleotideTorsion and CrossResidueMetric values compare equal as ints.
DC_DG = [
    (Tor.DELTA_1, ("C5'", "C4'", "C3'", "O3'")),
    (Tor.EPSILON_1, ("C4'", "C3'", "O3'", 'P')),
    (Tor.ZETA_1, ("C3'", "O3'", 'P', "O5'")),
    (Tor.ALPHA_2, ("O3'", 'P', "O5'", "C5'")),
    (Tor.BETA_2, ('P', "O5'", "C5'", "C4'")),
    (Tor.GAMMA_2, ("O5'", "C5'", "C4'", "C3'")),
    (Tor.DELTA_2, ("C5'", "C4'", "C3'", "O3'")),
    (Tor.CHI_1, ("O4'", "C1'", 'N1', 'C2')),
    (Tor.CHI_2, ("O4'", "C1'", 'N9', 'C4')),
    (XR.DIST_CC, ('C2', 'C4', None, None)),
    (XR.DIST_NN, ('N1', 'N9', None, None)),
    (XR.TOR_MU, ('C2', 'N1', 'N9', 'C4')),
]

A_PSU = [
    (Tor.DELTA_1, ("C5'", "C4'", "C3'", "O3'")),
    (Tor.EPSILON_1, ("C4'", "C3'", "O3'", 'P')),
    (Tor.ZETA_1, ("C3'", "O3'", 'P', "O5'")),
    (Tor.ALPHA_2, ("O3'", 'P', "O5'", "C5'")),
    (Tor.BETA_2, ('P', "O5'", "C5'", "C4'")),
    (Tor.GAMMA_2, ("O5'", "C5'", "C4'", "C3'")),
    (Tor.DELTA_2, ("C5'", "C4'", "C3'", "O3'")),
    (Tor.CHI_1, ("O4'", "C1'", 'N9', 'C4')),
    (Tor.CHI_2, ("O4'", "C1'", 'C5', 'C4')),
    (XR.DIST_CC, ('C4', 'C4', None, None)),
    (XR.DIST_NN, ('N9', 'C5', None, None)),
    (XR.TOR_MU, ('C4', 'N9', 'C5', 'C4')),
]


class TestNtC(LLKATestCase):
    def test_cana(self) -> None:
        self.assertEqual(llka.name_to_cana('AAw'), llka.CANA.AAw, 'wrong conversion from name to CANA item')
        with self.assertRaises(ValueError):
            llka.name_to_cana('this is wrong')

        self.assertEqual(llka.cana_to_name(llka.CANA.miB), 'miB', 'wrong conversion from CANA item to name')
        self.assertEqual(llka.cana_to_name(llka.CANA.INVALID), 'NAN', 'wrong conversion from invalid CANA item to name')

    def test_ntc(self) -> None:
        self.assertEqual(llka.name_to_ntc('BA17'), llka.NtC.BA17, 'wrong conversion from name to NtC item')
        with self.assertRaises(ValueError):
            llka.name_to_ntc('this is wrong')

        self.assertEqual(llka.ntc_to_name(llka.NtC.BB1S), 'BB1S', 'wrong conversion from NtC item to name')
        self.assertEqual(llka.ntc_to_name(llka.NtC.INVALID), 'NANT', 'wrong conversion from invalid NtC item to name')

    def test_ntc_names_match_enum(self) -> None:
        for ntc in llka.NtC:
            if ntc != llka.NtC.INVALID:
                self.assertEqual(llka.ntc_to_name(ntc), ntc.name)
                self.assertEqual(llka.name_to_ntc(ntc.name), ntc)

    def test_difference_against_reference(self) -> None:
        stru = make_structure('REAL_1BNA_A_1_2_ATOMS')
        metrics = llka.calculate_step_metrics_difference_against_reference(stru, llka.NtC.BB04)

        self.assertFlt(metrics.delta_1, 0.291304749269346, 'delta_1 metric is wrong')
        self.assertFlt(metrics.epsilon_1, 0.303056944684855, 'epsilon_1 metric is wrong')
        self.assertFlt(metrics.zeta_1, 0.032871663591650, 'zeta_1 metric is wrong')
        self.assertFlt(metrics.alpha_2, -0.355335537090472, 'alpha_2 metric is wrong')
        self.assertFlt(metrics.beta_2, 0.299443860928775, 'beta_2 metric is wrong')
        self.assertFlt(metrics.gamma_2, -0.104152469812837, 'gamma_2 metric is wrong')
        self.assertFlt(metrics.delta_2, -0.206873055232490, 'delta_2 metric is wrong')
        self.assertFlt(metrics.chi_1, -0.131009933521641, 'chi_1 metric is wrong')
        self.assertFlt(metrics.chi_2, -0.053816152662312, 'chi_2 metric is wrong')
        self.assertFlt(metrics.CC, 0.228744666804201, 'CC metric is wrong')
        self.assertFlt(metrics.NN, 0.003706138611939, 'NN metric is wrong')
        self.assertFlt(metrics.mu, 0.064645999883097, 'mu metric is wrong')

        with self.assertRaises(llka.InvalidArgumentError):
            llka.calculate_step_metrics_difference_against_reference(stru, llka.NtC.INVALID)

    def _check_names(self, get_quad, expected) -> None:
        for metric, names in expected:
            with self.subTest(metric=metric.name):
                self.assertEqual(tuple(get_quad(metric)), names)

    def _quad_from_bases(self, first: str, second: str):
        def get(metric):
            if isinstance(metric, Tor):
                return llka.dinucleotide_torsion_atoms_from_bases(first, second, metric)
            return llka.cross_residue_metric_atoms_from_bases(first, second, metric)
        return get

    def _quad_from_structure(self, stru: llka.Structure):
        def get(metric):
            if isinstance(metric, Tor):
                return llka.dinucleotide_torsion_atoms_from_structure(stru, metric)
            return llka.cross_residue_metric_atoms_from_structure(stru, metric)
        return get

    def test_torsion_names_bases(self) -> None:
        self._check_names(self._quad_from_bases('DC', 'DG'), DC_DG)
        self._check_names(self._quad_from_bases('A', 'PSU'), A_PSU)

        with self.assertRaises(llka.InvalidArgumentError):
            llka.dinucleotide_torsion_atoms_from_bases('Organic', 'A', Tor.DELTA_1)
        with self.assertRaises(llka.InvalidArgumentError):
            llka.dinucleotide_torsion_atoms_from_bases('A', 'Unicorn', Tor.DELTA_1)

    def test_torsion_names_structure(self) -> None:
        self._check_names(self._quad_from_structure(make_structure('REAL_1BNA_A_1_2_ATOMS')), DC_DG)
        self._check_names(self._quad_from_structure(make_structure('REAL_1BZT_A_12_13_ATOMS')), A_PSU)

    def test_structure_is_step(self) -> None:
        self.assertTrue(llka.is_step(make_structure('REAL_1BNA_A_1_2_ATOMS')))
        self.assertFalse(llka.is_step(make_structure('REAL_1DK1_B_27_28_ATOMS')))
        with self.assertRaises(llka.MissingAtomsError):
            llka.structure_is_step(llka.Structure())

    def test_ntc_structure(self) -> None:
        stru = llka.ntc_structure(llka.NtC.AB01)
        self.assertTrue(llka.is_step(stru))
        with self.assertRaises(llka.InvalidArgumentError):
            llka.ntc_structure(llka.NtC.INVALID)

    def test_backbone_atom_index(self) -> None:
        backbone = llka.extract_backbone(make_structure('REF_AB01_ATOMS'))
        self.assertEqual(llka.backbone_atom_index(llka.BackboneAtom(llka.ResidueInStep.FIRST_RESIDUE, "C5'"), backbone), 0)
        self.assertEqual(llka.backbone_atom_index(llka.BackboneAtom(llka.ResidueInStep.SECOND_RESIDUE, 'P'), backbone), 5)
        with self.assertRaises(llka.InvalidArgumentError):
            llka.backbone_atom_index(llka.BackboneAtom(llka.ResidueInStep.FIRST_RESIDUE, 'XX'), backbone)


if __name__ == '__main__':
    unittest.main()
