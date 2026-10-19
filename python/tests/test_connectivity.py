"""Port of tests/test_connectivity.cpp."""

from __future__ import annotations

import unittest

from _support import LLKATestCase, make_structure

import llka


class TestConnectivity(LLKATestCase):
    def assertConnectivity(self, conn: llka.Connectivity, c5: float, o3: float) -> None:
        self.assertFlt(conn.C5_prime_distance, c5, "wrong C5' distance")
        self.assertFlt(conn.O3_prime_distance, o3, "wrong O3' distance")

    def test_connectivity(self) -> None:
        # Test based on results for 1BNA molecule
        real_first = make_structure('REAL_1BNA_A_2_3_ATOMS')
        real_second = make_structure('REAL_1BNA_A_3_4_ATOMS')
        ref_first = make_structure('REF_BA01_ATOMS')
        ref_second = make_structure('REF_AB01_ATOMS')

        conn = llka.measure_step_connectivity_structures(real_first, ref_first, real_second, ref_second)
        self.assertConnectivity(conn, 0.0778497140143741, 0.110550489844076)

    def test_connectivity_multiple(self) -> None:
        real_first = make_structure('REAL_1BNA_A_3_4_ATOMS')
        real_second = make_structure('REAL_1BNA_A_4_5_ATOMS')
        ref_first = make_structure('REF_AB01_ATOMS')
        refs_second = [make_structure('REF_BA01_ATOMS'), make_structure('REF_BB00_ATOMS'), make_structure('REF_BB04_ATOMS')]

        results = llka.measure_step_connectivity_structures_multiple(real_first, ref_first, real_second, refs_second)

        self.assertEqual(len(results), 3)
        self.assertConnectivity(results[0], 0.222853768522876, 0.332800460520272)
        self.assertConnectivity(results[1], 0.398108806100796, 0.29940718941947)
        self.assertConnectivity(results[2], 0.379975441318778, 0.372926953223463)

    def test_connectivity_ntcs(self) -> None:
        real_first = make_structure('REAL_1BNA_A_2_3_ATOMS')
        real_second = make_structure('REAL_1BNA_A_3_4_ATOMS')

        conn = llka.measure_step_connectivity_ntcs(real_first, llka.NtC.BA01, real_second, llka.NtC.AB01)
        self.assertConnectivity(conn, 0.0778497140143741, 0.110550489844076)

        with self.assertRaises(llka.InvalidArgumentError):
            llka.measure_step_connectivity_ntcs(real_first, llka.NtC.INVALID, real_second, llka.NtC.AB01)

    def test_connectivity_multiple_ntcs_first(self) -> None:
        real_first = make_structure('REAL_1BNA_A_3_4_ATOMS')
        real_second = make_structure('REAL_1BNA_A_4_5_ATOMS')

        results = llka.measure_step_connectivity_ntcs_multiple_first(
            real_first, [llka.NtC.AB02, llka.NtC.BBS1, llka.NtC.BB07], real_second, llka.NtC.BB04
        )

        self.assertEqual(len(results), 3)
        self.assertConnectivity(results[0], 0.172379441824832, 0.218599268902842)
        self.assertConnectivity(results[1], 0.22807852040331, 0.118886962245301)
        self.assertConnectivity(results[2], 0.289601415171216, 0.106908891049413)

    def test_connectivity_multiple_ntcs_second(self) -> None:
        real_first = make_structure('REAL_1BNA_A_3_4_ATOMS')
        real_second = make_structure('REAL_1BNA_A_4_5_ATOMS')

        results = llka.measure_step_connectivity_ntcs_multiple_second(
            real_first, llka.NtC.AB01, real_second, [llka.NtC.BA01, llka.NtC.BB00, llka.NtC.BB04]
        )

        self.assertEqual(len(results), 3)
        self.assertConnectivity(results[0], 0.222853768522876, 0.332800460520272)
        self.assertConnectivity(results[1], 0.398108806100796, 0.29940718941947)
        self.assertConnectivity(results[2], 0.379975441318778, 0.372926953223463)


if __name__ == '__main__':
    unittest.main()
