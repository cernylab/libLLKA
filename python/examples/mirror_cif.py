#!/usr/bin/env python3
"""Loads and parses a mmCIF file and prints its contents back to standard output.

This is the Python counterpart of the ``mirror_cif`` C example. It demonstrates
how to load and read mmCIF files with libLLKA.

Usage: mirror_cif.py <input_file.cif>
"""

from __future__ import annotations

import sys

import llka


def format_value(value: llka.CifValue) -> str:
    if value.state == llka.CifValueState.SET:
        return value.text
    if value.state == llka.CifValueState.NONE:
        return '.'
    return '?'


def mirror(data: llka.CifData) -> None:
    out = sys.stdout

    for block in data:
        out.write(f'Data frame {block.name}\n')

        for category in block:
            out.write(f'_{category.name}\n')

            for item in category:
                out.write('\t' + item.keyword)
                for value in item.values:
                    out.write(' ' + format_value(value))
                out.write('\n')

            out.write('\n')

        out.write('\n')


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print('No input file', file=sys.stderr)
        return 1

    try:
        data = llka.cif_file_to_data(argv[1])
    except llka.LLKAError as ex:
        print(f'Failed to read Cif file: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
        if isinstance(ex, llka.CifParseError):
            print(ex, file=sys.stderr)
        return 1

    mirror(data)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
