"""Port of tests/test_nucleotide.cpp."""

from __future__ import annotations

import unittest

from _support import LLKATestCase, make_structure

import llka


class TestNucleotide(LLKATestCase):
    def setUp(self) -> None:
        self.stru = make_structure('REAL_1BNA_A_3_4_ATOMS')

    def test_extract_nucleotide(self) -> None:
        first = llka.extract_nucleotide(self.stru, 1, 'A', 3)
        self.assertEqual(len(first), 19, 'wrong number of atoms in nucleotide')
        self.assertEqual(first[0].label_seq_id, first[-1].label_seq_id, 'sequence ids of first and last atom do not match')

        second = llka.extract_nucleotide(self.stru, 1, 'A', 4)
        self.assertEqual(len(second), 22, 'wrong number of atoms in nucleotide')
        self.assertEqual(second[0].label_seq_id, second[-1].label_seq_id, 'sequence ids of first and last atom do not match')

        self.assertEqual(first[0].label_seq_id + 1, second[0].label_seq_id, 'nucleotides have the same sequence id')

        self.assertEqual(len(llka.extract_nucleotide(self.stru, 1, 'A', 99)), 0, 'nonexistent nucleotide')

    def test_extract_ribose(self) -> None:
        nucl = llka.extract_nucleotide(self.stru, 1, 'A', 3)

        ribose = llka.extract_ribose(nucl)
        self.assertEqual([a.label_atom_id for a in ribose], ["C4'", "O4'", "C1'", "C2'", "C3'"], 'wrong ribose atoms')

        ribose_view = llka.extract_ribose_view(nucl)
        self.assertEqual(list(ribose_view), list(ribose), 'ribose view does not match')

        with self.assertRaises(llka.MissingAtomsError):
            llka.extract_ribose(llka.Structure())

    def test_sugar_pucker(self) -> None:
        nucl = llka.extract_nucleotide(self.stru, 1, 'A', 3)
        self.assertEqual(llka.sugar_pucker(nucl), llka.SugarPucker.O4_ENDO, 'Wrong sugar pucker')

        nucl_view = llka.extract_nucleotide_view(self.stru, 1, 'A', 4)
        self.assertEqual(llka.sugar_pucker(nucl_view), llka.SugarPucker.C2_ENDO, 'Wrong sugar pucker')

    def test_ribose_metrics(self) -> None:
        nucl = llka.extract_nucleotide(self.stru, 1, 'A', 3)
        metrics = llka.ribose_metrics(nucl)
        self.assertEqual(metrics.pucker, llka.SugarPucker.O4_ENDO, 'Wrong sugar pucker')

        with self.assertRaises(llka.InvalidArgumentError):
            llka.ribose_metrics(llka.Structure())

    def test_ribose_metrics_view(self) -> None:
        nucl_view = llka.extract_nucleotide_view(self.stru, 1, 'A', 4)
        metrics = llka.ribose_metrics(nucl_view)
        self.assertEqual(metrics.pucker, llka.SugarPucker.C2_ENDO, 'Wrong sugar pucker')
        self.assertEqual(metrics, llka.ribose_metrics(nucl_view.to_structure()), 'View and structure metrics differ')

        empty_view = llka.extract_nucleotide_view(self.stru, 1, 'A', 99)
        with self.assertRaises(llka.InvalidArgumentError):
            llka.ribose_metrics(empty_view)

    def test_is_nucleotide_compound(self) -> None:
        self.assertTrue(llka.is_nucleotide_compound('DA'))
        self.assertTrue(llka.is_nucleotide_compound('PSU'))
        self.assertFalse(llka.is_nucleotide_compound('ALA'))


if __name__ == '__main__':
    unittest.main()
