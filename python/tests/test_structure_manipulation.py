"""Port of tests/test_structure_manipulation.cpp."""

from __future__ import annotations

import copy
import unittest

from _support import LLKATestCase, make_structure

import llka


def _atom(atom_id: int, name: str = "C5'", comp: str = 'DC', type_symbol: str = 'C',
          coords=(1.059, -3.973, 1.922), seq_id: int = 1, model: int = 1) -> llka.Atom:
    return llka.Atom(
        id=atom_id, type_symbol=type_symbol, label_atom_id=name, label_entity_id='1', label_comp_id=comp,
        label_asym_id='A', coords=coords, label_seq_id=seq_id, auth_seq_id=seq_id, pdbx_PDB_model_num=model,
    )


class TestStructureManipulation(LLKATestCase):
    def test_comparison(self) -> None:
        atom_a = _atom(1, comp='DA', coords=(10, 11, 12))
        atom_aa = _atom(2, comp='DA', coords=(10, 11, 12))
        atom_b = _atom(3, name="C4'", comp='DA', coords=(10, 11, 12), seq_id=2)
        atom_b2 = _atom(4, name="C4'", comp='DA', coords=(10, 11, 12), seq_id=2, model=2)

        self.assertFalse(llka.compare_atoms(atom_a, atom_aa), 'compare identical atoms')
        self.assertTrue(llka.compare_atoms(atom_a, atom_aa, ignore_id=True), 'compare identical atoms (ignore ID)')
        self.assertFalse(llka.compare_atoms(atom_a, atom_b, ignore_id=True), 'compare different atoms')
        self.assertFalse(llka.compare_atoms(atom_b, atom_b2), 'compare atoms from different models')

    def test_duplication(self) -> None:
        atom_a = _atom(1, comp='DA', coords=(10, 11, 12))
        atom_b = copy.copy(atom_a)
        self.assertTrue(llka.compare_atoms(atom_a, atom_b), 'duplicate atom')

        stru = llka.Structure([atom_a, _atom(2)])
        dup = stru.copy()
        self.assertEqual(dup, stru, 'duplicate structure')
        self.assertIsNot(dup, stru)

    def test_append_atom(self) -> None:
        stru = llka.Structure([_atom(1), _atom(2), _atom(3), _atom(4)])
        to_append = _atom(5, name="C4'", coords=(0.551, -3.732, 0.531))

        stru.append(to_append)

        self.assertEqual(len(stru), 5, 'append atom, count')
        self.assertTrue(llka.compare_atoms(stru[4], to_append), 'append atom, compare appended')

        stru.append(llka.Atom(
            id=6, type_symbol='C', label_atom_id="C4'", label_entity_id='1', label_comp_id='DA', label_asym_id='A',
            label_seq_id=6, auth_seq_id=6, coords=llka.Point(10.0, 11.0, 12.0),
        ))

        self.assertEqual(len(stru), 6, 'append atom by params, count')
        self.assertEqual(stru[5].type_symbol, 'C', 'append atom by params, type_symbol')
        self.assertEqual(stru[5].label_comp_id, 'DA', 'append atom by params, label_comp_id')
        self.assertEqual(stru[5].coords.x, 10.0, 'append atom by params, X coordinate')
        self.assertEqual(stru[5].auth_comp_id, 'DA', 'auth_comp_id defaults to label_comp_id')

    def test_append_structure(self) -> None:
        first = llka.Structure([_atom(idx) for idx in range(1, 5)])
        second = llka.Structure([_atom(idx, name='N9', comp='DA', type_symbol='N') for idx in range(5, 10)])

        first.extend(second)

        self.assertEqual(len(first), 9, 'wrong number of atoms')
        self.assertEqual(first[0].auth_atom_id, "C5'", 'wrong atom name')
        self.assertEqual(first[8].auth_atom_id, 'N9', 'wrong atom name')
        self.assertEqual(first[8].id, 9, 'wrong atom id')
        self.assertEqual(len(second), 5, 'appended structure must not be modified')

    def test_removal(self) -> None:
        atoms = [_atom(idx) for idx in range(1, 5)]

        for removed, expected in ((1, [2, 3, 4]), (4, [1, 2, 3]), (2, [1, 3, 4])):
            with self.subTest(removed=removed):
                stru = llka.Structure(atoms)
                stru.remove_atom_by_id(removed)
                self.assertEqual(len(stru), 3, 'remove atom, count')
                self.assertEqual([a.id for a in stru], expected, 'remove atom, remaining atoms')

        stru = llka.Structure(atoms)
        with self.assertRaises(ValueError):
            stru.remove_atom_by_id(42)

    def test_find_atoms(self) -> None:
        stru = make_structure('REAL_1DK1_B_27_28_ATOMS')

        self.assertEqual(stru.find_atom_index('OP2', 'A', 'B', 27, 'B'), 5, 'expected atom not found')
        self.assertEqual(stru.find_atom('OP2', 'A', 'B', 27, 'B'), stru[5], 'expected atom not found')

        self.assertEqual(stru.find_atom_index('OP2', 'A', 'B', 27, 'A'), 4, 'expected atom not found')
        self.assertEqual(stru.find_atom_index('OP2', 'A', 'B', 27), 4, 'expected atom not found')

        self.assertEqual(stru.find_atom_index("O5'"), 6, 'expected atom not found')
        self.assertIsNone(stru.find_atom("O5'", label_seq_id=29), 'did not expect to find an atom')

    def test_find_atoms_view(self) -> None:
        stru = make_structure('REAL_1DK1_B_27_28_ATOMS')
        view = llka.extract_nucleotide_view(stru, 1, 'B', 27)
        self.assertGreater(len(view), 0)

        idx = view.find_atom_index('OP2', 'A', 'B', 27, 'B')
        self.assertIsNotNone(idx, 'expected atom not found')
        self.assertEqual(view[idx], stru[5], 'wrong atom found')
        self.assertEqual(view.find_atom('OP2', 'A', 'B', 27, 'B'), stru[5], 'wrong atom found')

        self.assertEqual(view.find_atom('OP2', 'A', 'B', 27, 'A'), stru[4], 'wrong atom found')
        self.assertEqual(view.find_atom('OP2', 'A', 'B', 27), stru[4], 'wrong atom found')
        self.assertEqual(view.find_atom("O5'"), stru[6], 'wrong atom found')

        # The view contains only residue 27
        self.assertIsNotNone(stru.find_atom("O5'", label_seq_id=28))
        self.assertIsNone(view.find_atom("O5'", label_seq_id=28), 'did not expect to find an atom')

        self.assertIsNone(llka.extract_nucleotide_view(stru, 1, 'B', 99).find_atom("O5'"), 'empty view')

        stru.append(stru[0].replace(id=100000))
        with self.assertRaises(RuntimeError):
            view.find_atom("O5'")

    def test_sequence_protocol(self) -> None:
        stru = llka.Structure([_atom(idx) for idx in range(1, 6)])

        self.assertEqual([a.id for a in stru[1:4]], [2, 3, 4], 'slicing')
        self.assertEqual(stru[-1].id, 5, 'negative index')

        stru[0] = stru[0].replace(id=10)
        self.assertEqual(stru[0].id, 10, 'item assignment')

        del stru[1]
        self.assertEqual([a.id for a in stru], [10, 3, 4, 5], 'item deletion')

        stru.insert(1, _atom(2))
        self.assertEqual([a.id for a in stru], [10, 2, 3, 4, 5], 'insertion')

        with self.assertRaises(IndexError):
            stru[5]


if __name__ == '__main__':
    unittest.main()
