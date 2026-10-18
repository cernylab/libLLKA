"""Shared helpers of the Python test suite."""

from __future__ import annotations

import functools
import re
import unittest
from pathlib import Path

import llka

REPO_ROOT = Path(__file__).resolve().parents[2]
ASSETS_DIR = REPO_ROOT / 'assets'
TEST_CIFS_DIR = ASSETS_DIR / 'test_cifs'
TESTING_STRUCTURES_HEADER = REPO_ROOT / 'tests' / 'testing_structures.h'

# Same tolerance as EffedUp::cmpFlt() in the C++ test suite
FLOAT_TOLERANCE = 1.0e-11

_TOKEN_RE = re.compile(
    r"""
    \s*(?:
        (?P<string>"(?:[^"\\]|\\.)*")
      | (?P<char>'(?:[^'\\]|\\.)')
      | (?P<number>[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?[uUlL]*)
      | (?P<ident>[A-Za-z_]\w*)
      | (?P<punct>[{},;=\[\]])
    )
    """,
    re.VERBOSE,
)

_ARRAY_RE = re.compile(r'inline\s+constinit\s+LLKA_Atom\s+(\w+)\s*\[\]\s*=\s*\{')


def _tokenize(text: str):
    pos = 0
    while pos < len(text):
        m = _TOKEN_RE.match(text, pos)
        if m is None or m.end() == pos:
            raise ValueError(f'Cannot tokenize test data near: {text[pos:pos + 40]!r}')
        pos = m.end()
        kind = m.lastgroup
        yield kind, m.group(kind)
        if kind == 'punct' and m.group(kind) == ';':
            return


def _value(kind: str, text: str):
    if kind == 'string':
        return bytes(text[1:-1], 'utf-8').decode('unicode_escape')
    if kind == 'char':
        return bytes(text[1:-1], 'utf-8').decode('unicode_escape')
    if kind == 'number':
        text = text.rstrip('uUlL')
        return float(text) if any(ch in text for ch in '.eE') else int(text)
    if kind == 'ident':
        return {'LLKA_NO_INSCODE': '', 'LLKA_NO_ALTID': ''}[text]
    raise ValueError(f'Unexpected token {text!r}')


def _parse_initializer(tokens) -> list:
    """Parses a brace-enclosed C initializer list into nested Python lists."""
    out: list = []
    for kind, text in tokens:
        if kind == 'punct':
            if text == '{':
                out.append(_parse_initializer(tokens))
            elif text == '}':
                return out
            elif text == ',':
                continue
            else:
                raise ValueError(f'Unexpected {text!r} in initializer')
        else:
            out.append(_value(kind, text))
    raise ValueError('Unterminated initializer')


def _to_atom(fields: list) -> llka.Atom:
    (type_symbol, label_atom_id, label_entity_id, label_comp_id, label_asym_id,
     auth_atom_id, auth_comp_id, auth_asym_id, coords, atom_id, label_seq_id, auth_seq_id,
     model_num, ins_code, alt_id) = fields
    return llka.Atom(
        id=atom_id,
        type_symbol=type_symbol,
        label_atom_id=label_atom_id,
        label_entity_id=label_entity_id,
        label_comp_id=label_comp_id,
        label_asym_id=label_asym_id,
        auth_atom_id=auth_atom_id,
        auth_comp_id=auth_comp_id,
        auth_asym_id=auth_asym_id,
        coords=llka.Point(*(float(c) for c in coords)),
        label_seq_id=label_seq_id,
        auth_seq_id=auth_seq_id,
        pdbx_PDB_model_num=model_num,
        pdbx_PDB_ins_code=ins_code,
        label_alt_id=alt_id,
    )


@functools.cache
def testing_atoms() -> dict[str, tuple[llka.Atom, ...]]:
    """Atoms of all testing structures defined in ``tests/testing_structures.h``, keyed by array name."""
    text = TESTING_STRUCTURES_HEADER.read_text()
    structures: dict[str, tuple[llka.Atom, ...]] = {}
    for m in _ARRAY_RE.finditer(text):
        tokens = _tokenize(text[m.end():])
        atoms = _parse_initializer(tokens)
        structures[m.group(1)] = tuple(_to_atom(a) for a in atoms)

    # Cross-check against the _LEN constants
    for name, length in re.findall(r'inline\s+constinit\s+size_t\s+(\w+)_LEN\s*=\s*(\d+)\s*;', text):
        if len(structures[name]) != int(length):
            raise ValueError(f'Parsed {len(structures[name])} atoms of {name} but the header says {length}')

    return structures


def make_structure(name: str) -> llka.Structure:
    """Creates a structure from one of the testing structures, e.g. ``"REAL_1BNA_A_1_2_ATOMS"``."""
    return llka.Structure(testing_atoms()[name])


class LLKATestCase(unittest.TestCase):
    def assertFlt(self, actual: float, expected: float, msg: str | None = None) -> None:
        self.assertAlmostEqual(actual, expected, delta=FLOAT_TOLERANCE, msg=msg)

    def assertPoint(self, actual: llka.Point, expected: tuple[float, float, float], msg: str | None = None) -> None:
        for a, e in zip(actual, expected):
            self.assertFlt(a, e, msg)
