"""Runs the Python example programs in python/examples.

If the ``LLKA_C_EXAMPLES_DIR`` environment variable points to the built C
examples, the output of each Python example is compared byte by byte with the
output of the corresponding C example. CTest sets the variable when the C
examples are built.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from _support import ASSETS_DIR, TEST_CIFS_DIR, REPO_ROOT, testing_atoms

import llka

EXAMPLES_DIR = REPO_ROOT / 'python' / 'examples'

# Python example script -> C example program
EXAMPLES = {
    'mirror_cif.py': 'mirror_cif',
    'simple_ntc_assignment.py': 'simple_NtC_assignment',
    'classify_and_write_cif.py': 'classify_and_write_cif',
}

_ATOM_SITE_COLUMNS = [
    'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id',
    'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z',
    'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num',
]


def _write_cif(path: Path, entry_id: str, atoms) -> None:
    """Writes atoms as a minimal mmCIF file."""
    none = llka.CifValue.none()
    data = llka.CifData()
    block = data.add_block(entry_id)
    block.add_category('entry').add_item('id', [entry_id])
    rows = [
        [
            'ATOM', str(a.id), a.type_symbol, a.label_atom_id, a.label_alt_id or none, a.label_comp_id,
            a.label_asym_id, a.label_entity_id or '1', str(a.label_seq_id), a.pdbx_PDB_ins_code or none,
            f'{a.coords.x:.3f}', f'{a.coords.y:.3f}', f'{a.coords.z:.3f}', str(a.auth_seq_id), a.auth_comp_id,
            a.auth_asym_id, a.auth_atom_id, str(a.pdbx_PDB_model_num),
        ]
        for a in atoms
    ]
    cat = block.add_category('atom_site')
    for column, keyword in enumerate(_ATOM_SITE_COLUMNS):
        cat.add_item(keyword, [row[column] for row in rows])
    data.detaint()
    path.write_text(data.to_string(pretty=True))


def _make_inputs(directory: Path) -> list[Path]:
    """Creates input files that exercise alternate positions, insertion codes, multiple models etc."""
    atoms = testing_atoms()
    base = atoms['REAL_1BNA_A_1_2_ATOMS']
    inputs = {
        'altids.cif': ('1DK1', atoms['REAL_1DK1_B_26_28_ATOMS']),
        'violations.cif': ('1DK1', atoms['REAL_1DK1_B_5_6_ATOMS']),
        'syn.cif': ('3VOK', atoms['REAL_3VOK_U_1_2_ATOMS']),
        'pseudouridine.cif': ('1BZT', atoms['REAL_1BZT_A_12_13_ATOMS']),
        'inscode.cif': ('1BNA', [a.replace(pdbx_PDB_ins_code='A') if a.label_seq_id == 2 else a for a in base]),
        'multimodel.cif': ('1BNA', list(base) + [
            a.replace(id=a.id + 10000, pdbx_PDB_model_num=2, coords=llka.Point(a.coords.x + 0.05, a.coords.y, a.coords.z))
            for a in base
        ]),
    }
    paths = [TEST_CIFS_DIR / '1BNA.cif', TEST_CIFS_DIR / '1BNA_broken.cif']
    for name, (entry_id, stru_atoms) in inputs.items():
        _write_cif(directory / name, entry_id, stru_atoms)
        paths.append(directory / name)
    return paths


def _run(command: list[str]) -> subprocess.CompletedProcess:
    # The classification parameter files are expected in the working directory
    return subprocess.run(command, cwd=ASSETS_DIR, capture_output=True, timeout=600)


class TestExamples(unittest.TestCase):
    tmp: tempfile.TemporaryDirectory
    inputs: list[Path]

    @classmethod
    def setUpClass(cls) -> None:
        cls.tmp = tempfile.TemporaryDirectory()
        cls.inputs = _make_inputs(Path(cls.tmp.name))

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmp.cleanup()

    def run_example(self, script: str, *args: str) -> subprocess.CompletedProcess:
        return _run([sys.executable, str(EXAMPLES_DIR / script), *args])

    def test_simple_ntc_assignment(self) -> None:
        result = self.run_example('simple_ntc_assignment.py', str(TEST_CIFS_DIR / '1BNA.cif'))
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        out = result.stdout.decode()
        self.assertIn('* Step 1\n  DNATCO name:                       1bna_A_DC_1_DG_2\n', out)
        self.assertIn('  Assigned NtC:  BB04, Closest NtC:  BB04\n', out)

    def test_classify_and_write_cif(self) -> None:
        result = self.run_example('classify_and_write_cif.py', str(TEST_CIFS_DIR / '1BNA.cif'))
        self.assertEqual(result.returncode, 0, result.stderr.decode())

        # The output ends with a comment line followed by an empty line. The libLLKA Cif parser
        # reads past the end of such input, so strip the trailing whitespace before parsing.
        data = llka.cif_text_to_data(result.stdout.decode().rstrip() + '\n')
        block = data[0]
        self.assertEqual(block['ndb_struct_ntc_overall']['entry_id'].values, [llka.CifValue('1BNA')])
        self.assertEqual(block['ndb_struct_ntc_step']['name'][0], llka.CifValue('1bna_A_DC_1_DG_2'))
        self.assertEqual(block['ndb_struct_ntc_step_summary']['assigned_NtC'][0], llka.CifValue('BB04'))

    def test_mirror_cif(self) -> None:
        result = self.run_example('mirror_cif.py', str(TEST_CIFS_DIR / '1BNA.cif'))
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        self.assertTrue(result.stdout.decode().startswith('Data frame 1BNA\n_entry\n\tid 1BNA\n'))

    def test_errors(self) -> None:
        for script in EXAMPLES:
            with self.subTest(script=script):
                self.assertNotEqual(self.run_example(script).returncode, 0, 'missing argument')
                self.assertNotEqual(self.run_example(script, str(TEST_CIFS_DIR / '1BNA_broken.cif')).returncode, 0)

    def test_output_matches_c_examples(self) -> None:
        c_examples_dir = os.environ.get('LLKA_C_EXAMPLES_DIR')
        if not c_examples_dir:
            self.skipTest('LLKA_C_EXAMPLES_DIR is not set')

        for script, program in EXAMPLES.items():
            c_program = Path(c_examples_dir) / program / program
            for args in [[str(p)] for p in self.inputs] + [[], ['/nonexistent.cif']]:
                with self.subTest(example=script, args=args):
                    c_result = _run([str(c_program), *args])
                    py_result = self.run_example(script, *args)
                    self.assertEqual(py_result.returncode, c_result.returncode, 'exit codes differ')
                    self.assertEqual(py_result.stdout, c_result.stdout, 'standard outputs differ')
                    self.assertEqual(py_result.stderr, c_result.stderr, 'standard error outputs differ')

        result = self.run_example('classify_and_write_cif.py', '--version')
        c_result = _run([str(Path(c_examples_dir) / 'classify_and_write_cif' / 'classify_and_write_cif'), '--version'])
        self.assertEqual(result.stdout, c_result.stdout, '--version outputs differ')


if __name__ == '__main__':
    unittest.main()
