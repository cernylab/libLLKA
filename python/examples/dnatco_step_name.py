"""Construction of DNATCO step names.

This is the Python counterpart of ``dnatco_step_name.c`` used by the C examples.
The logic is a bit lengthy so it is kept in a separate module to keep the
example code brief and easy to follow.
"""

from __future__ import annotations

from collections.abc import Sequence

import llka


def _lower_entry_id(entry_id: str) -> str:
    # Same as the C example: the first character of a PDB ID is a digit,
    # only the remaining three characters are converted to lowercase
    return entry_id[:1] + entry_id[1:4].lower() + entry_id[4:]


def step_alt_ids(step: llka.Structure) -> tuple[str, str]:
    """Returns alternate position IDs of the first and the second nucleotide of a step.

    An empty string means that the nucleotide has no alternate positions.
    """
    atoms = list(step)
    if not atoms:
        raise ValueError('Step has no atoms')

    seq_id_first = atoms[0].label_seq_id
    alt_id_first = atoms[0].label_alt_id

    # Check if the first nucleotide in the step has any altId
    idx = 1
    while idx < len(atoms):
        if alt_id_first or atoms[idx].label_seq_id != seq_id_first:
            break
        alt_id_first = atoms[idx].label_alt_id
        idx += 1

    # Scroll to the second nucleotide
    while idx < len(atoms) and atoms[idx].label_seq_id == seq_id_first:
        idx += 1

    if idx == len(atoms):
        raise ValueError('Step seems to have only one nucleotide')

    alt_id_second = atoms[idx].label_alt_id
    idx += 1

    # Check if the second nucleotide has any altId
    while idx < len(atoms) and not alt_id_second:
        alt_id_second = atoms[idx].label_alt_id
        idx += 1

    return alt_id_first, alt_id_second


def derive_dnatco_step_name(entry_id: str, has_multiple_models: bool, step: llka.Structure) -> str:
    """Returns the DNATCO name of a step, e.g. ``1bna_A_DC_1_DG_2``."""
    alt_id_first, alt_id_second = step_alt_ids(step)
    first_atom = step[0]
    last_atom = step[-1]

    name = _lower_entry_id(entry_id)
    if has_multiple_models:
        name += f'-m{first_atom.pdbx_PDB_model_num}'

    name += f'_{first_atom.auth_asym_id}_{first_atom.auth_comp_id}'
    if alt_id_first:
        name += f'.{alt_id_first}'
    name += f'_{first_atom.auth_seq_id}'
    if first_atom.pdbx_PDB_ins_code:
        name += f'.{first_atom.pdbx_PDB_ins_code}'

    name += f'_{last_atom.auth_comp_id}'
    if alt_id_second:
        name += f'.{alt_id_second}'
    name += f'_{last_atom.auth_seq_id}'
    if last_atom.pdbx_PDB_ins_code:
        name += f'.{last_atom.pdbx_PDB_ins_code}'

    return name


def structure_has_multiple_models(steps: Sequence[llka.Structure]) -> bool:
    """Returns whether the steps come from more than one model."""
    model_nums = {step[0].pdbx_PDB_model_num for step in steps if len(step) > 0}
    return len(model_nums) > 1
