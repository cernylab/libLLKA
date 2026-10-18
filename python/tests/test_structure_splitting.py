"""Port of tests/test_structure_splitting.cpp."""

from __future__ import annotations

import unittest

from _support import LLKATestCase, make_structure


class TestStructureSplitting(LLKATestCase):
    def test_split_alt_ids(self) -> None:
        stru = make_structure('REAL_1DK1_B_26_28_ATOMS')

        splitted = stru.split_by_alt_ids()

        self.assertEqual(len(splitted), 2, 'wrong number of alt-id-splitted structures')
        alt_ids = [alt_id for alt_id, _ in splitted]
        strus = [s for _, s in splitted]

        self.assertEqual(len(strus[0]), 68, 'wrong number of atoms in splitted structure')
        self.assertEqual(len(strus[1]), 68, 'wrong number of atoms in splitted structure')

        self.assertEqual(alt_ids, ['A', 'B'], 'wrong alternate position ids')

        self.assertEqual(strus[0][0].label_alt_id, '', 'wrong alternate position id')
        self.assertEqual(strus[0][67].label_alt_id, 'A', 'wrong alternate position id')
        self.assertEqual(strus[1][0].label_alt_id, '', 'wrong alternate position id')
        self.assertEqual(strus[1][67].label_alt_id, 'B', 'wrong alternate position id')

    def test_split_without_alt_ids(self) -> None:
        stru = make_structure('REAL_1BNA_A_1_2_ATOMS')

        splitted = stru.split_by_alt_ids()

        self.assertEqual(len(splitted), 1)
        self.assertEqual(splitted[0][0], '')
        self.assertEqual(splitted[0][1], stru)

    def test_split_dinucleotides(self) -> None:
        stru = make_structure('REAL_1DK1_B_26_28_ATOMS')

        steps = stru.split_to_dinucleotide_steps()
        self.assertEqual(len(steps), 3, 'wrong number of steps')

        self.assertEqual(steps[0][0].label_seq_id, 26, 'wrong sequence id')
        self.assertEqual(steps[0][-1].label_seq_id, 27, 'wrong sequence id')

        self.assertEqual(steps[2][0].label_seq_id, 27, 'wrong sequence id')
        self.assertEqual(steps[2][0].label_alt_id, 'B', 'wrong alternate position id')
        self.assertEqual(steps[2][-1].label_seq_id, 28, 'wrong sequence id')
        self.assertEqual(steps[2][-1].label_alt_id, 'B', 'wrong alternate position id')

    def test_steps_outlive_each_other(self) -> None:
        stru = make_structure('REAL_1DK1_B_26_28_ATOMS')
        steps = stru.split_to_dinucleotide_steps()
        last = steps[-1]
        del steps, stru

        # The remaining step must still be valid after the other steps were released
        self.assertEqual(last[0].label_seq_id, 27)
        last.append(last[0].replace(id=100000))
        self.assertEqual(last[-1].id, 100000)


if __name__ == '__main__':
    unittest.main()
