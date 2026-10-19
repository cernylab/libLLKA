"""Port of tests/test_measurements.cpp."""

from __future__ import annotations

import unittest

from _support import LLKATestCase, make_structure

import llka


class TestMeasurements(LLKATestCase):
    def setUp(self) -> None:
        self.stru = make_structure('REAL_1BNA_A_1_2_ATOMS')

    def atom(self, name: str, seq_id: int = 1) -> llka.Atom:
        atom = self.stru.find_atom(name, label_seq_id=seq_id)
        self.assertIsNotNone(atom, f'atom {name} {seq_id} not found')
        return atom

    def test_distances(self) -> None:
        self.assertFlt(llka.measure_distance(self.atom("O5'"), self.atom("C5'")), 1.4378821231241452,
                       "unexpected distance between O5' -> C5' atoms")
        self.assertFlt(llka.measure_distance(self.atom("C1'"), self.atom('C2')), 2.5047702489450008,
                       "unexpected distance between C1' -> C2 atoms")
        self.assertFlt(llka.measure_distance(self.atom("O5'"), self.atom("O5'", 2)), 7.2086283022500197,
                       "unexpected distance between O5'_1 -> O5'_2 atoms")

    def test_angles(self) -> None:
        self.assertFlt(llka.measure_angle(self.atom("O5'"), self.atom("C5'"), self.atom("C4'")), 1.885722272657443,
                       "unexpected angle between O5' - C5' - C4' atoms")
        self.assertFlt(llka.measure_angle(self.atom("C5'"), self.atom("C4'"), self.atom("C3'")), 2.0342231990866746,
                       "unexpected angle between C5' - C4' - C3' atoms")
        self.assertFlt(llka.measure_angle(self.atom('O2'), self.atom('C2'), self.atom('C5')), 3.1036907363579922,
                       'unexpected angle between O2 - C2 - C5 atoms')

    def test_dihedrals(self) -> None:
        self.assertFlt(
            llka.measure_dihedral(self.atom("O5'"), self.atom("C5'"), self.atom("C4'"), self.atom("O4'")),
            0.93724348718692196,
            "unexpected dihedral angle between O5' - C5' - C4' - O4' atoms",
        )
        self.assertFlt(
            llka.measure_dihedral(self.atom("O5'"), self.atom("C5'"), self.atom("C4'"), self.atom("C3'")),
            3.0409263486359723,
            "unexpected dihedral angle between O5' - C5' - C4' - C3' atoms",
        )
        self.assertFlt(
            llka.measure_dihedral(self.atom('C4'), self.atom('N3'), self.atom('C2'), self.atom('N1')),
            0.0075379644566811376,
            'unexpected dihedral angle between C4 - N3 - C2 - N1 atoms',
        )

    def test_single_precision(self) -> None:
        a, b, c, d = self.atom("O5'"), self.atom("C5'"), self.atom("C4'"), self.atom("O4'")
        tolerance = 1.0e-5
        self.assertAlmostEqual(llka.measure_distance(a, b, single_precision=True), 1.4378821231241452, delta=tolerance)
        self.assertAlmostEqual(llka.measure_angle(a, b, c, single_precision=True), 1.885722272657443, delta=tolerance)
        self.assertAlmostEqual(llka.measure_dihedral(a, b, c, d, single_precision=True), 0.93724348718692196, delta=tolerance)


if __name__ == '__main__':
    unittest.main()
