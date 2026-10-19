"""Verifies that the ctypes mirrors of the C structures and enums match the C headers.

The test generates a C program that prints the sizes and field offsets of all
structures and the values of all enumerators, compiles it and compares the
output with the ctypes definitions and the Python enums.

The test needs a C compiler. It uses the ``LLKA_TEST_C_COMPILER`` and
``LLKA_TEST_C_FLAGS`` environment variables which CTest sets up. Without them,
the test is skipped.
"""

from __future__ import annotations

import ctypes
import enum
import os
import shlex
import subprocess
import tempfile
import unittest
from pathlib import Path

from _support import REPO_ROOT

import llka
from llka import _capi, minicif

# Fields whose names differ between the ctypes mirror and the C definition
_C_FIELD_NAMES = {
    ('LLKA_Points', 'points'): 'data',  # Anonymous union in C++, named "data" in C
}

_HEADERS = [
    'llka_main.h',
    'llka_structure.h',
    'llka_ntc.h',
    'llka_nucleotide.h',
    'llka_superposition.h',
    'llka_segmentation.h',
    'llka_connectivity_similarity.h',
    'llka_classification.h',
    'llka_resource_loaders.h',
    'llka_minicif.h',
    'llka_tracing.h',
]


def _enum_constants() -> list[tuple[str, str, int]]:
    """Returns (description, C enumerator, Python value) triples."""
    def members(cls: type[enum.Enum], name_of) -> list[tuple[str, str, int]]:
        out = []
        for member in cls:
            c_name = name_of(member.name)
            if c_name is not None:
                out.append((f'{cls.__name__}.{member.name}', c_name, int(member.value)))
        return out

    def invalid(suffix: str):
        return lambda name: f'LLKA_INVALID_{suffix}' if name == 'INVALID' else f'LLKA_{name}'

    constants = []
    constants += members(llka.RetCode, lambda n: f'LLKA_{n}')
    constants += members(llka.BaseKind, lambda n: f'LLKA_{n}')
    constants += members(llka.CrossResidueMetric, lambda n: f'LLKA_XR_{n}')
    constants += members(llka.DinucleotideTorsion, lambda n: f'LLKA_TOR_{n}')
    constants += members(llka.ResidueInStep, lambda n: f'LLKA_{n}')
    constants += members(llka.CANA, invalid('CANA'))
    constants += members(llka.NtC, invalid('NTC'))
    constants += members(llka.SugarPucker, invalid('SUGAR_PUCKER'))
    constants += members(llka.SugarPuckerNameBrevity, lambda n: f'LLKA_SPN_{n}')
    constants += members(llka.ResourceType, lambda n: f'LLKA_RES_{n}')
    constants += members(
        llka.ClassificationViolation,
        lambda n: 'LLKA_CLASSIFICATION_OK' if n == 'OK' else f'LLKA_CLASSIFICATION_E_{n}',
    )
    constants += members(llka.ViolatingTorsion, lambda n: None if n == 'NONE' else f'LLKA_CLASSIFICATION_E_TORSION_{n}')
    constants += members(
        llka.CifValueState,
        lambda n: {'SET': 'LLKA_MINICIF_VALUE_SET', 'NONE': 'LLKA_MINICIF_VALUE_NONE', 'UNKNOWN': 'LLKA_MINICIF_VALUE_UNKW'}[n],
    )
    constants += [
        ('minicif._GET_CIFDATA', 'LLKA_MINICIF_GET_CIFDATA', minicif._GET_CIFDATA),
        ('minicif._ALLOW_NO_ENTRY_CATEGORY', 'LLKA_MINICIF_ALLOW_NO_ENTRY_CATEGORY', minicif._ALLOW_NO_ENTRY_CATEGORY),
        ('minicif._ALLOW_BROKEN_ATOMSITE', 'LLKA_MINICIF_ALLOW_BROKEN_ATOMSITE', minicif._ALLOW_BROKEN_ATOMSITE),
        ('_capi.LLKA_INVALID_BKBN_ATOM_INDEX', 'LLKA_INVALID_BKBN_ATOM_INDEX', _capi.LLKA_INVALID_BKBN_ATOM_INDEX),
    ]
    return constants


_ENUM_TYPES = [
    'LLKA_RetCode', 'LLKA_BaseKind', 'LLKA_CrossResidueMetric', 'LLKA_DinucleotideTorsion', 'LLKA_ResidueInStep',
    'LLKA_CANA', 'LLKA_NtC', 'LLKA_SugarPucker', 'LLKA_SugarPuckerNameBrevity', 'LLKA_ResourceType',
    'LLKA_CifDataValueState',
]


def _generate_program() -> tuple[str, dict[str, int]]:
    """Generates the C program and the expected output keyed by the printed label."""
    expected: dict[str, int] = {}
    lines = ['#include <stddef.h>', '#include <stdio.h>']
    lines += [f'#include <{h}>' for h in _HEADERS]
    lines += ['int main(void)', '{']

    def emit(label: str, c_expr: str, value: int) -> None:
        expected[label] = value
        lines.append(f'    printf("%s %llu\\n", "{label}", (unsigned long long)({c_expr}));')

    for name, cls in _capi.STRUCTURES.items():
        emit(f'sizeof:{name}', f'sizeof({name})', ctypes.sizeof(cls))
        for field_name, *_ in cls._fields_:
            c_field = _C_FIELD_NAMES.get((name, field_name), field_name)
            emit(f'offsetof:{name}.{field_name}', f'offsetof({name}, {c_field})', getattr(cls, field_name).offset)

    for enum_type in _ENUM_TYPES:
        emit(f'sizeof:{enum_type}', f'sizeof({enum_type})', ctypes.sizeof(_capi.c_enum))

    for label, c_name, value in _enum_constants():
        # Compare as unsigned 64-bit values so that both sides agree on wrapping
        emit(f'value:{label}', f'(unsigned long long){c_name}', value & 0xFFFFFFFFFFFFFFFF)

    lines += ['    return 0;', '}']
    return '\n'.join(lines) + '\n', expected


class TestAbi(unittest.TestCase):
    def test_layouts_and_constants_match(self) -> None:
        compiler = os.environ.get('LLKA_TEST_C_COMPILER')
        if not compiler:
            self.skipTest('LLKA_TEST_C_COMPILER is not set')
        flags = shlex.split(os.environ.get('LLKA_TEST_C_FLAGS', ''))

        source, expected = _generate_program()
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / 'abi_check.c'
            exe = Path(tmp) / 'abi_check'
            src.write_text(source)
            subprocess.run(
                [compiler, *flags, f'-I{REPO_ROOT / "include"}', str(src), '-o', str(exe)],
                check=True, capture_output=True, text=True,
            )
            output = subprocess.run([str(exe)], check=True, capture_output=True, text=True).stdout

        actual = {}
        for line in output.splitlines():
            label, value = line.rsplit(' ', 1)
            actual[label] = int(value)

        mismatches = {label: (value, actual.get(label)) for label, value in expected.items() if actual.get(label) != value}
        self.assertEqual(mismatches, {}, 'Mismatches between Python and C (python value, C value)')


if __name__ == '__main__':
    unittest.main()
