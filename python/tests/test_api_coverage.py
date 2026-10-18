"""Checks that the bindings cover the whole callable C API of libLLKA."""

from __future__ import annotations

import re
import unittest
from pathlib import Path

from _support import REPO_ROOT

from llka import _capi

INCLUDE_DIR = REPO_ROOT / 'include'
PACKAGE_DIR = Path(_capi.__file__).parent

_DECL_RE = re.compile(r'^LLKA_API\s+[^;{]*?\b(LLKA_\w+)\s*\(', re.S | re.M)


def declared_functions() -> set[str]:
    names: set[str] = set()
    for header in INCLUDE_DIR.glob('llka_*.h'):
        if header.name == 'llka_cpp.h':
            continue
        names.update(_DECL_RE.findall(header.read_text()))
    return names


class TestApiCoverage(unittest.TestCase):
    def test_all_functions_have_prototypes(self) -> None:
        declared = declared_functions()
        self.assertGreater(len(declared), 100, 'header parsing seems broken')

        missing = declared - set(_capi.PROTOTYPES) - _capi.BROKEN_FUNCTIONS
        self.assertEqual(missing, set(), 'C functions without a ctypes prototype')

        unknown = (set(_capi.PROTOTYPES) | _capi.BROKEN_FUNCTIONS) - declared
        self.assertEqual(unknown, set(), 'prototypes of functions that are not declared in the headers')

        self.assertEqual(set(_capi.PROTOTYPES) & _capi.BROKEN_FUNCTIONS, set(), 'broken functions must not be bound')

    def test_all_prototypes_are_used(self) -> None:
        sources = '\n'.join(p.read_text() for p in PACKAGE_DIR.glob('*.py') if p.name != '_capi.py')
        # Helpers in _capi.py use some functions too; exclude the tables that merely list the functions
        capi = (PACKAGE_DIR / '_capi.py').read_text()
        capi = re.sub(r'\n(?:PROTOTYPES|UNWRAPPED_FUNCTIONS|BROKEN_FUNCTIONS|FILESYSTEM_FUNCTIONS)\b.*?\n[})]', '', capi, flags=re.S)
        sources += capi
        unused = {name for name in _capi.PROTOTYPES if not re.search(rf'\b{name}\b', sources)}
        self.assertEqual(unused, set(_capi.UNWRAPPED_FUNCTIONS), 'C functions not used by the Python API')

    def test_broken_functions_are_not_exported(self) -> None:
        # If any of these becomes callable, it should be wrapped and removed from BROKEN_FUNCTIONS
        for name in _capi.BROKEN_FUNCTIONS:
            with self.subTest(function=name):
                self.assertFalse(_capi.has_function(name))


if __name__ == '__main__':
    unittest.main()
