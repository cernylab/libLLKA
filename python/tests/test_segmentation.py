"""Port of tests/test_segmentation.cpp."""

from __future__ import annotations

import unittest

from _support import LLKATestCase, make_structure

import llka


class TestSegmentation(LLKATestCase):
    def test_segmentation(self) -> None:
        stru = make_structure('REAL_1DK1_B_26_28_ATOMS')
        segs = llka.structure_segments(stru)

        self.assertEqual(len(segs), 1, 'Wrong number of models')
        model = segs[0]
        self.assertEqual(len(model.atoms), 113, 'Wrong number of atoms in model')
        self.assertEqual(model.atoms[0].id, 541, 'Wrong atom id of the first atom in model')

        self.assertEqual(len(model.chains), 1, 'Wrong number of chains in model')
        chain = model.chains[0]
        self.assertEqual(len(chain.residues), 3, 'Wrong number of residues in chain')
        self.assertEqual(len(chain.atoms), 113, 'Wrong number of atoms in chain')
        self.assertEqual(chain.atoms[0].id, 541, 'Wrong atom id of the first atom in chain')

        residue_a = chain.residues[0]
        self.assertEqual(len(residue_a.atoms), 23, 'Wrong number of atoms in residue')
        self.assertEqual(residue_a.atoms[0].id, 541, 'Wrong atom id of first atom in residue')
        self.assertEqual(residue_a.atoms[0].label_seq_id, 26, 'Wrong label_seq_id of first atom in residue')
        self.assertEqual(residue_a.atoms[22].id, 563, 'Wrong atom id of last atom in residue')
        self.assertEqual(residue_a.atoms[22].label_seq_id, 26, 'Wrong label_seq_id of last atom in residue')

        residue_b = chain.residues[2]
        self.assertEqual(len(residue_b.atoms), 46, 'Wrong number of atoms in residue')
        self.assertEqual(residue_b.atoms[0].id, 608, 'Wrong atom id of first atom in residue')
        self.assertEqual(residue_b.atoms[0].label_seq_id, 28, 'Wrong label_seq_id of first atom in residue')
        self.assertEqual(residue_b.atoms[0].label_alt_id, 'A', 'Wrong label_alt_id of first atom in residue')
        self.assertEqual(residue_b.atoms[45].id, 653, 'Wrong atom id of last atom in residue')
        self.assertEqual(residue_b.atoms[45].label_seq_id, 28, 'Wrong label_seq_id of last atom in residue')
        self.assertEqual(residue_b.atoms[45].label_alt_id, 'B', 'Wrong label_alt_id of last atom in residue')

    def test_segmentation_multiple_chains(self) -> None:
        stru = make_structure('REAL_1CJG_A_B_1_2_M1')
        segs = llka.structure_segments(stru)

        self.assertEqual(len(segs), 1, 'Wrong number of models')

        model = segs[0]
        self.assertEqual(len(model.atoms), 82, 'Wrong number of atoms in model')
        self.assertEqual(model.atoms[0].id, 1397, 'Wrong atom id of the first atom in model')

        self.assertEqual(len(model.chains), 2, 'Wrong number of chains in model')

        expected = {
            'A': ((19, 1397), (22, 1416)),
            'B': ((19, 2363), (22, 2382)),
        }
        for chain, (asym_id, residues) in zip(model.chains, expected.items()):
            with self.subTest(chain=asym_id):
                self.assertEqual(len(chain.residues), 2, 'Wrong number of residues in chain')
                self.assertEqual(len(chain.atoms), 41, 'Wrong number of atoms in chain')
                for residue, (n_atoms, first_id) in zip(chain.residues, residues):
                    self.assertEqual(len(residue.atoms), n_atoms, 'Wrong number of atoms in residue')
                    self.assertEqual(residue.atoms[0].id, first_id, 'Wrong atom id on residue')
                    self.assertEqual(residue.atoms[0].label_asym_id, asym_id, 'Wrong label_asym_id on residue')

    def test_segmentation_multiple_models(self) -> None:
        stru = make_structure('REAL_1CJG_A_1_6_M6_M1_M2_M3')
        segs = llka.structure_segments(stru)

        self.assertEqual(len(segs), 4, 'Wrong number of models')

        model_numbers = [1, 2, 3, 6]
        residue_seq_ids = [1, 6]

        for model, model_num in zip(segs, model_numbers):
            with self.subTest(model=model_num):
                self.assertEqual(len(model.atoms), 38, 'Wrong number of atoms in model')
                self.assertEqual(model.atoms[0].label_atom_id, 'N', 'Wrong atom id of the first atom in model')
                self.assertEqual(model.atoms[0].pdbx_PDB_model_num, model_num, 'Wrong pdbx_PDB_model_num')

                self.assertEqual(len(model.chains), 1, 'Wrong number of chains in model')
                chain = model.chains[0]
                self.assertEqual(len(chain.residues), 2, 'Wrong number of residues in chain')
                self.assertEqual(len(chain.atoms), 38, 'Wrong number of atoms in chain')
                self.assertEqual(chain.atoms[0].label_asym_id, 'A', 'Wrong label_asym_id on first atom in chain')

                first, second = chain.residues
                self.assertEqual(len(first.atoms), 19, 'Wrong number of atoms in first residue')
                self.assertEqual(first.atoms[0].label_atom_id, 'N', 'Wrong label_atom_id on first atom in residue')
                self.assertEqual(first.atoms[18].label_atom_id, 'HE3', 'Wrong label_atom_id on last atom in residue')
                self.assertEqual(first.atoms[0].label_seq_id, residue_seq_ids[0], 'Wrong label_seq_id')
                self.assertEqual(first.atoms[18].label_seq_id, residue_seq_ids[0], 'Wrong label_seq_id')

                self.assertEqual(len(second.atoms), 19, 'Wrong number of atoms in second residue')
                self.assertEqual(second.atoms[0].label_atom_id, 'N', 'Wrong label_atom_id on first atom in residue')
                self.assertEqual(second.atoms[18].label_atom_id, 'HD23', 'Wrong label_atom_id on last atom in residue')
                self.assertEqual(second.atoms[0].label_seq_id, residue_seq_ids[1], 'Wrong label_seq_id')
                self.assertEqual(second.atoms[18].label_seq_id, residue_seq_ids[1], 'Wrong label_seq_id')

    def test_segments_invalidated_by_structure_change(self) -> None:
        stru = make_structure('REAL_1DK1_B_26_28_ATOMS')
        segs = llka.structure_segments(stru)
        residue = segs[0].chains[0].residues[0]

        stru.append(stru[0].replace(id=100000))

        with self.assertRaises(RuntimeError):
            residue.atoms[0]


if __name__ == '__main__':
    unittest.main()
