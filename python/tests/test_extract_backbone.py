"""Port of tests/test_extract_backbone.cpp."""

from __future__ import annotations

import unittest

from _support import LLKATestCase, make_structure

import llka

EXPECTED_ORDER = [
    ("C5'", 1), ("C4'", 1), ("O4'", 1), ("C3'", 1), ("O3'", 1),
    ('P', 2), ("O5'", 2), ("C5'", 2), ("C4'", 2), ("O4'", 2), ("C3'", 2), ("O3'", 2),
]


def _base_atoms(kind: llka.BaseKind) -> tuple[str, str]:
    return ('N9', 'C4') if kind == llka.BaseKind.PURINE else ('N1', 'C2')


def _extended_expected_order(info: llka.StepInfo, first_seq_id: int = 1, second_seq_id: int = 2):
    first = _base_atoms(info.first_base_kind)
    second = _base_atoms(info.second_base_kind)
    return [
        ("C5'", first_seq_id), ("C4'", first_seq_id), ("O4'", first_seq_id), ("C3'", first_seq_id),
        ("O3'", first_seq_id), ("C1'", first_seq_id), (first[0], first_seq_id), (first[1], first_seq_id),
        ('P', second_seq_id), ("O5'", second_seq_id), ("C5'", second_seq_id), ("C4'", second_seq_id),
        ("O4'", second_seq_id), ("C3'", second_seq_id), ("O3'", second_seq_id), ("C1'", second_seq_id),
        (second[0], second_seq_id), (second[1], second_seq_id),
    ]


class TestExtractBackbone(LLKATestCase):
    def assertOrder(self, atoms, expected) -> None:
        self.assertEqual(len(atoms), len(expected), 'number of backbone atoms is wrong')
        for atom, (name, seq_id) in zip(atoms, expected):
            self.assertEqual(atom.label_atom_id, name, 'label_atom_id is wrong')
            self.assertEqual(atom.label_seq_id, seq_id, 'label_seq_id is wrong')

    def test_plain_backbone(self) -> None:
        stru = make_structure('REF_AB01_ATOMS')
        backbone = llka.extract_backbone(stru)
        self.assertIsInstance(backbone, llka.Structure)
        self.assertOrder(backbone, EXPECTED_ORDER)

    def test_plain_backbone_view(self) -> None:
        stru = make_structure('REF_AB01_ATOMS')
        backbone = llka.extract_backbone_view(stru)
        self.assertIsInstance(backbone, llka.StructureView)
        self.assertIs(backbone.structure, stru)
        self.assertOrder(backbone, EXPECTED_ORDER)

    def test_plain_backbone_alt_pos(self) -> None:
        stru = make_structure('REAL_1DK1_B_27_28_ATOMS')
        with self.assertRaises(llka.MultipleAltIdsError) as cm:
            llka.extract_backbone(stru)
        self.assertEqual(cm.exception.retcode, llka.RetCode.E_MULTIPLE_ALT_IDS)

        stru = make_structure('REAL_1DK1_B_27_28_ALT_A_ONLY_ATOMS')
        llka.extract_backbone(stru)

    def test_extended_backbone(self) -> None:
        stru = make_structure('REF_AB01_ATOMS')
        info = llka.structure_is_step(stru)
        backbone = llka.extract_extended_backbone(stru)
        self.assertOrder(backbone, _extended_expected_order(info))

    def test_extended_backbone_view(self) -> None:
        stru = make_structure('REF_AB01_ATOMS')
        info = llka.structure_is_step(stru)
        backbone = llka.extract_extended_backbone_view(stru)
        self.assertOrder(backbone, _extended_expected_order(info))

    def test_extended_backbone_pseudouridine(self) -> None:
        stru = make_structure('REAL_1BZT_A_12_13_ATOMS')
        info = llka.structure_is_step(stru)
        self.assertEqual(info.second_base_kind, llka.BaseKind.NON_STANDARD_BASE, 'Wrong kind of second base')

        expected = _extended_expected_order(info, 12, 13)
        expected[16] = ('C5', 13)
        expected[17] = ('C4', 13)

        backbone = llka.extract_extended_backbone(stru)
        self.assertOrder(backbone, expected)

    def test_view_invalidated_by_structure_change(self) -> None:
        stru = make_structure('REF_AB01_ATOMS')
        backbone = llka.extract_backbone_view(stru)
        self.assertTrue(backbone.is_valid)

        stru.remove_atom_by_id(stru[0].id)

        self.assertFalse(backbone.is_valid)
        with self.assertRaises(RuntimeError):
            list(backbone)


if __name__ == '__main__':
    unittest.main()
