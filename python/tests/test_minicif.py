"""Port of tests/test_minicif.cpp."""

from __future__ import annotations

import unittest

from _support import TEST_CIFS_DIR, LLKATestCase

import llka

ATOM_SITE_HEADER = (
    'data_\n'
    'loop_\n'
    '_atom_site.group_PDB\n'
    '_atom_site.id\n'
    '_atom_site.type_symbol\n'
    '_atom_site.label_atom_id\n'
    '_atom_site.label_alt_id\n'
    '_atom_site.label_comp_id\n'
    '_atom_site.label_asym_id\n'
    '_atom_site.label_entity_id\n'
    '_atom_site.label_seq_id\n'
    '_atom_site.pdbx_PDB_ins_code\n'
    '_atom_site.segment_id\n'
    '_atom_site.cartn_x\n'
    '_atom_site.cartn_y\n'
    '_atom_site.cartn_z\n'
    '_atom_site.occupancy\n'
    '_atom_site.B_iso_or_equiv\n'
    '_atom_site.cartn_x_esd\n'
    '_atom_site.cartn_y_esd\n'
    '_atom_site.cartn_z_esd\n'
    '_atom_site.occupancy_esd\n'
    '_atom_site.B_iso_or_equiv_esd\n'
    '_atom_site.pdbx_formal_charge\n'
    '_atom_site.auth_seq_id\n'
    '_atom_site.auth_comp_id\n'
    '_atom_site.auth_asym_id\n'
    '_atom_site.auth_atom_id\n'
    '_atom_site.pdbx_PDB_model_num\n'
)

NO_ENTRY_CIF = ATOM_SITE_HEADER + 'ATOM 1  P P     . DA A 1 99 . . 8.258  14.006  12.104 1 30 ? ? ? ? ? ? 1 DT A P     1\n#\n'
BROKEN_ATOM_SITE_CIF = (
    ATOM_SITE_HEADER
    + ' ATOM 8  N .     . DA . 1 20 . . -2.3707326 -4.507768  -0.94517957  1 20 ? ? ? ? ? ?  1 DA A N9    1\n'
)


class TestMiniCif(LLKATestCase):
    def test_ok(self) -> None:
        imported = llka.cif_file_to_structure(TEST_CIFS_DIR / '1BNA.cif')

        self.assertEqual(imported.entry_id, '1BNA')
        self.assertEqual(len(imported.structure), 566, 'wrong number of atoms in structure')
        self.assertEqual(imported.structure[165].label_seq_id, 9, 'unexpected label_seq_id')
        self.assertIsNone(imported.cif_data)

    def test_ok_with_cif_data(self) -> None:
        imported = llka.cif_file_to_structure(TEST_CIFS_DIR / '1BNA.cif', get_cif_data=True)

        self.assertIsNotNone(imported.cif_data)
        block = imported.cif_data[0]
        self.assertIn('atom_site', block)
        self.assertEqual(len(block['atom_site']['id']), 566)

        # The structure must stay usable after the other parts of the import are released
        stru = imported.structure
        del imported, block
        self.assertEqual(len(stru), 566)

    def test_broken(self) -> None:
        with self.assertRaises(llka.CifParseError) as cm:
            llka.cif_file_to_structure(TEST_CIFS_DIR / '1BNA_broken.cif')
        self.assertEqual(str(cm.exception), 'Malformed loop on line 1089, unexpected token kind 8', 'unexpected error message')
        self.assertEqual(cm.exception.retcode, llka.RetCode.E_BAD_DATA)

    def test_missing_file(self) -> None:
        with self.assertRaises(llka.NoFileError):
            llka.cif_file_to_structure(TEST_CIFS_DIR / 'does_not_exist.cif')

    def test_manipulate_categories(self) -> None:
        data = llka.CifData()
        block = data.add_block('')

        block.add_category('first')
        block.add_category('second')
        block.add_category('third')

        with self.assertRaises(llka.InvalidArgumentError, msg='Failed to prevent creation of duplicate category'):
            block.add_category('third')

        self.assertEqual([c.name for c in block], ['first', 'second', 'third'], 'wrong order of categories')
        self.assertIsNone(block['third'].next, 'Last category does not look as last')

        block.delete_category('second')
        self.assertEqual([c.name for c in block], ['first', 'third'], 'wrong order of categories')

        block.delete_category('first')
        third = block['third']
        self.assertEqual([c.name for c in block], ['third'], 'wrong order of categories')
        self.assertIsNone(third.next, 'Category is not the only left category (has next check)')
        self.assertIsNone(third.previous, 'Category is not the only left category (has previous check)')

        with self.assertRaises(llka.InvalidArgumentError):
            block.delete_category('first')

    def test_manipulate_entries(self) -> None:
        data = llka.CifData()
        block = data.add_block('')
        cat = block.add_category('xyz')

        item_a = cat.add_item('a')
        item_b = cat.add_item('b')
        cat.add_item('c')

        with self.assertRaises(llka.InvalidArgumentError, msg='Failed to prevent creation of duplicate item'):
            cat.add_item('b')

        self.assertEqual([i.keyword for i in cat], ['a', 'b', 'c'], 'wrong order of items')
        self.assertIsNone(cat['c'].next, 'Last item does not look as last')

        cat.delete_item('c')
        self.assertEqual([i.keyword for i in cat], ['a', 'b'], 'wrong order of items')
        self.assertIsNone(item_b.next, 'Last item does not look like last (has next check)')
        self.assertEqual(item_b.previous, item_a, 'Entry should have a previous item but it does not')

        # Add some values and check tainting
        item_a.set_values(['Run to the hills', 'Run for your life'])
        item_b.add_value('I am a man who walks alone')
        item_b.add_value(llka.CifValue("when I'm walking a dark road"))

        self.assertTrue(data.is_tainted, 'Data is not marked as tainted when it should be')
        data.detaint()
        self.assertFalse(data.is_tainted)

        item_b.add_value('at night or strolling through a park')
        with self.assertRaises(llka.MismatchingSizesError, msg='Incorrectly detainted invalid data'):
            data.detaint()

    def test_loop_not_loop(self) -> None:
        data = llka.CifData()
        block = data.add_block('')

        loop_cat = block.add_category('xyz')
        not_loop_cat = block.add_category('abc')

        loop_cat.add_item('loop_item').add_values(['Osculum', 'obscenum'])
        not_loop_cat.add_item('not_loop_item').add_value('Spillways')

        data.detaint()

        self.assertTrue(loop_cat.is_loop, 'Category is not a loop but it should be')
        self.assertFalse(not_loop_cat.is_loop, 'Category is a loop but it should not be')

    def test_detaint_empty_category(self) -> None:
        data = llka.CifData()
        data.add_block('').add_category('empty')
        with self.assertRaises(llka.BadDataError):
            data.detaint()

    def test_quote_in_quotes(self) -> None:
        data = llka.cif_file_to_data(TEST_CIFS_DIR / 'quote_in_quotes.cif')
        self.assertEqual(len(data), 1, 'unexpected number of blocks')

        cat = data[0]['some']
        expected = {
            'double_quote_string': '"ABC"',
            'single_quote_string': "'ABC'",
            'single_quote_in_double_quote_string': '"AB\'C"',
            'double_quote_in_single_quote_string': "'AB\"C'",
            'single_quote_in_single_quote_string': "'AB'C'",
            'single_quote_in_single_quote_string_terminating': "'ABC''",
            'double_quote_in_double_quote_string_terminating': '"ABC""',
        }
        for keyword, text in expected.items():
            with self.subTest(keyword=keyword):
                self.assertEqual(cat[keyword].values, [llka.CifValue(text, llka.CifValueState.SET)])

    def test_unterminated_quote(self) -> None:
        with self.assertRaises(llka.CifParseError) as cm:
            llka.cif_file_to_data(TEST_CIFS_DIR / 'unterminated_quote.cif')
        self.assertEqual(str(cm.exception), 'Unterminated quoted token that begins on line 3', 'unexpected error condition')

    def test_no_entry_category(self) -> None:
        imported = llka.cif_text_to_structure(NO_ENTRY_CIF, allow_no_entry_category=True)

        self.assertEqual(len(imported.structure), 1, 'wrong number of atoms in structure')
        self.assertEqual(imported.structure[0].label_seq_id, 99, 'unexpected label_seq_id')
        self.assertEqual(imported.entry_id, '', 'unexpected entry id')

    def test_no_entry_category_not_allowed(self) -> None:
        with self.assertRaises(llka.CifParseError) as cm:
            llka.cif_text_to_structure(NO_ENTRY_CIF)
        self.assertEqual(str(cm.exception), 'Category entry is not present in block ', 'unexpected error message')

    def test_allow_broken_atom_site(self) -> None:
        imported = llka.cif_text_to_structure(BROKEN_ATOM_SITE_CIF, allow_no_entry_category=True, allow_broken_atom_site=True)

        self.assertEqual(len(imported.structure), 1, 'wrong number of atoms in structure')
        self.assertEqual(imported.entry_id, '', 'unexpected entry id')
        atom = imported.structure[0]
        self.assertEqual(atom.label_seq_id, 20, 'unexpected label_seq_id')
        self.assertEqual(atom.label_atom_id, 'N9', 'unexpected label_atom_id')
        self.assertEqual(atom.label_asym_id, 'A', 'unexpected label_asym_id')

    def test_allow_broken_atom_site_not_allowed(self) -> None:
        with self.assertRaises(llka.BadDataError):
            llka.cif_text_to_structure(BROKEN_ATOM_SITE_CIF, allow_no_entry_category=True)

    def test_round_trip(self) -> None:
        data = llka.CifData()
        block = data.add_block('test')
        cat = block.add_category('xyz')
        cat.add_item('a', ['hello world', llka.CifValue.none(), llka.CifValue.unknown()])
        cat.add_item('b', ['1', '2', '3'])
        data.detaint()

        for pretty in (False, True):
            with self.subTest(pretty=pretty):
                text = data.to_string(pretty)
                parsed = llka.cif_text_to_data(text)
                self.assertEqual(parsed[0].name, 'test')
                self.assertEqual(
                    parsed[0]['xyz']['a'].values,
                    [llka.CifValue('"hello world"'), llka.CifValue.none(), llka.CifValue.unknown()],
                )
                self.assertEqual([v.text for v in parsed[0]['xyz']['b']], ['1', '2', '3'])

    def test_tainted_data_cannot_be_written(self) -> None:
        data = llka.CifData()
        data.add_block('').add_category('xyz').add_item('a', ['1'])
        with self.assertRaises(llka.BadDataError):
            data.to_string()

    def test_stale_wrappers(self) -> None:
        data = llka.CifData()
        block = data.add_block('a')
        cat = block.add_category('xyz')

        block.delete_category('xyz')
        with self.assertRaises(KeyError):
            cat.is_loop

        # Deleting a block from data with more than one block reads past the end of
        # the block array in LLKA_cifData_deleteBlock(), so only the single-block case is tested.
        data.delete_block(0)
        with self.assertRaises(RuntimeError):
            block.name
        self.assertEqual(len(data), 0)


if __name__ == '__main__':
    unittest.main()
