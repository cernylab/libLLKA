#!/usr/bin/env python3
"""Assigns NtC classes to a structure and writes the results as mmCIF.

This is the Python counterpart of the ``classify_and_write_cif`` C example.
It reads a structure from a mmCIF file, performs the NtC assignment and prints
the original mmCIF data extended with categories that contain the results of
the assignment. The output mostly matches the output of DNATCO 4.1.

The classification parameter files are looked up in the working directory,
in ``${CCP4}/share/dnatco`` and in the directory given by the
``DNATCO_ASSETS_PATH`` environment variable, in this order.

Usage: classify_and_write_cif.py [--version] <input_file.cif>
"""

from __future__ import annotations

import math
import os
import sys
from collections.abc import Sequence

import llka
from dnatco_step_name import derive_dnatco_step_name, structure_has_multiple_models

VERSION = '0.0.2'
MAX_CLOSE_ENOUGH_RMSD = 0.5

OVERALL_KEYWORDS = [
    'entry_id', 'confal_score', 'confal_percentile', 'ntc_version', 'cana_version',
    'num_steps', 'num_classified', 'num_unclassified', 'num_unclassified_rmsd_close',
]
NTC_KEYWORDS = [
    'id', 'name', 'PDB_model_number',
    'label_entity_id_1', 'label_asym_id_1', 'label_seq_id_1', 'label_comp_id_1', 'label_alt_id_1',
    'label_entity_id_2', 'label_asym_id_2', 'label_seq_id_2', 'label_comp_id_2', 'label_alt_id_2',
    'auth_asym_id_1', 'auth_seq_id_1', 'auth_asym_id_2', 'auth_seq_id_2', 'PDB_ins_code_1', 'PDB_ins_code_2',
]
SUMMARY_KEYWORDS = [
    'step_id', 'assigned_CANA', 'assigned_NtC', 'confal_score', 'euclidean_distance_NtC_ideal',
    'cartesian_rmsd_closest_NtC_representative', 'closest_CANA', 'closest_NtC', 'closest_step_golden',
]
PARAMETERS_KEYWORDS = [
    'step_id', 'tor_delta_1', 'tor_epsilon_1', 'tor_zeta_1', 'tor_alpha_2', 'tor_beta_2', 'tor_gamma_2',
    'tor_delta_2', 'tor_chi_1', 'tor_chi_2', 'dist_NN', 'dist_CC', 'tor_NCCN',
    'diff_tor_delta_1', 'diff_tor_epsilon_1', 'diff_tor_zeta_1', 'diff_tor_alpha_2', 'diff_tor_beta_2',
    'diff_tor_gamma_2', 'diff_tor_delta_2', 'diff_tor_chi_1', 'diff_tor_chi_2', 'diff_dist_NN', 'diff_dist_CC',
    'diff_tor_NCCN', 'confal_tor_delta_1', 'confal_tor_epsilon_1', 'confal_tor_zeta_1', 'confal_tor_alpha_2',
    'confal_tor_beta_2', 'confal_tor_gamma_2', 'confal_tor_delta_2', 'confal_tor_chi_1', 'confal_tor_chi_2',
    'confal_dist_NN', 'confal_dist_CC', 'confal_tor_NCCN', 'details',
]
SUGAR_KEYWORDS = [
    'step_id', 'P_1', 'tau_1', 'Pn_1', 'P_2', 'tau_2', 'Pn_2',
    'nu_1_1', 'nu_1_2', 'nu_1_3', 'nu_1_4', 'nu_1_5', 'nu_2_1', 'nu_2_2', 'nu_2_3', 'nu_2_4', 'nu_2_5',
    'diff_nu_1_1', 'diff_nu_1_2', 'diff_nu_1_3', 'diff_nu_1_4', 'diff_nu_1_5',
    'diff_nu_2_1', 'diff_nu_2_2', 'diff_nu_2_3', 'diff_nu_2_4', 'diff_nu_2_5',
]

TORSIONS = ['delta_1', 'epsilon_1', 'zeta_1', 'alpha_2', 'beta_2', 'gamma_2', 'delta_2', 'chi_1', 'chi_2']
NU_ANGLES = ['nu_0', 'nu_1', 'nu_2', 'nu_3', 'nu_4']
DETAILS_TORSION_NAMES = ['d1', 'e1', 'z1', 'a2', 'b2', 'g2', 'd2', 'ch1', 'ch2']

NONE = llka.CifValue.none()


# Angle conversions with the same operation order as the inline functions in llka_util.h
# so that the output is identical to the output of the C example

def deg2rad(deg: float) -> float:
    return math.pi * deg / 180.0


def rad2deg(rad: float) -> float:
    return 180.0 * rad / math.pi


def full_angle(rad: float) -> float:
    """Converts an angle from the [-PI; PI] to the [0; 2PI] range."""
    return (2.0 * math.pi + rad) if rad < 0.0 else rad + 0.0  # + 0.0 turns -0.0 into 0.0 like the C code does


def deg1(rad: float) -> str:
    """Angle in degrees in the [0; 360] range with one decimal place."""
    return f'{rad2deg(full_angle(rad)):.1f}'


def diff1(rad: float) -> str:
    """Angle difference in degrees with one decimal place."""
    return f'{rad2deg(rad):.1f}'


def rounded(value: float) -> str:
    """Rounds a positive value to an integer the way the C example does."""
    return str(int(value + 0.5))


def find_asset(filename: str) -> str:
    """Returns the path of a classification parameters file. See the module documentation for the search order."""
    candidates = [filename]
    ccp4 = os.environ.get('CCP4')
    if ccp4:
        candidates.append(f'{ccp4}/share/dnatco/{filename}')
    assets_path = os.environ.get('DNATCO_ASSETS_PATH')
    if assets_path:
        candidates.append(f'{assets_path}/{filename}')

    for path in candidates:
        if os.path.isfile(path) and os.access(path, os.R_OK):
            return path

    # Not found, the attempt to load the file will fail with an appropriate error
    return filename


def initialize_classification_context() -> llka.ClassificationContext:
    limits = llka.ClassificationLimits(
        minimum_nearest_neighbors=7,
        number_of_used_nearest_neighbors=11,
        average_neighbors_torsion_cutoff=deg2rad(28.0),
        nearest_neighbor_torsions_cutoff=deg2rad(28.0),
        total_distance_cutoff=deg2rad(60.0),
        pseudorotation_cutoff=deg2rad(72.0),
        minimum_cluster_votes=0.001111,
    )

    resources = {}
    for key, filename, kind, description in (
        ('golden_steps', 'golden_steps.csv', llka.ResourceType.GOLDEN_STEPS, 'golden steps'),
        ('clusters', 'clusters.csv', llka.ResourceType.CLUSTERS, 'clusters'),
        ('confals', 'confals.csv', llka.ResourceType.CONFALS, 'confals'),
        ('cluster_nu_angles', 'nu_angles.csv', llka.ResourceType.AVERAGE_NU_ANGLES, 'average nu angles'),
        ('confal_percentiles', 'confal_percentiles.csv', llka.ResourceType.CONFAL_PERCENTILES, 'confal percentiles'),
    ):
        try:
            resources[key] = llka.load_resource_file(find_asset(filename), kind)
        except llka.LLKAError as ex:
            print(
                f'Failed to load {description}: {llka.error_to_string(ex.retcode)}\n'
                'The file is expected in cwd, ${CCP4}/share/dnatco, or you can set DNATCO_ASSETS_PATH environment variable.',
                file=sys.stderr,
            )
            sys.exit(1)

    try:
        return llka.ClassificationContext(**resources, limits=limits, max_close_enough_rmsd=MAX_CLOSE_ENOUGH_RMSD)
    except llka.LLKAError as ex:
        print(f'Failed to initialize classification context: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
        sys.exit(1)


def residue_alt_id(stru: llka.Structure, seq_id: int) -> str:
    """Returns the alternate position ID of a residue in a step.

    Some atoms of a residue may have alternate positions while others don't, e.g. backbone
    atoms may have alt_id 'A' or 'B' while base atoms have none. Use the first alt_id found.
    """
    for atom in stru:
        if atom.label_seq_id == seq_id and atom.label_alt_id:
            return atom.label_alt_id
    return ''


def ntc_name(ntc: llka.NtC | None) -> str:
    return llka.ntc_to_name(llka.NtC.INVALID if ntc is None else ntc)


def cana_name(cana: llka.CANA | None) -> str:
    return llka.cana_to_name(llka.CANA.INVALID if cana is None else cana)


def details(cs: llka.ClassifiedStep) -> str:
    """Short codes describing why a step could not be classified."""
    V = llka.ClassificationViolation
    codes = []

    if cs.violations & V.AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT:
        codes += [f'cAn{name}' for bit, name in enumerate(DETAILS_TORSION_NAMES) if cs.violating_torsions_average & (1 << bit)]
    if cs.violations & V.NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT:
        # The C example reports the torsions from violatingTorsionsAverage here as well; kept for identical output
        codes += [f'cNn{name}' for bit, name in enumerate(DETAILS_TORSION_NAMES) if cs.violating_torsions_average & (1 << bit)]
    if cs.violations & (V.NN_TOO_LOW | V.NN_TOO_HIGH):
        codes.append('cNN')
    if cs.violations & (V.CC_TOO_LOW | V.CC_TOO_HIGH):
        codes.append('cCC')
    if cs.violations & (V.MU_TOO_LOW | V.MU_TOO_HIGH):
        codes.append('cmu')
    if cs.violations & V.TOTAL_DISTANCE_TOO_HIGH:
        codes.append('cMB')
    if cs.violations & V.FIRST_PSEUDOROTATION_TOO_DIFFERENT:
        codes.append('cP')
    if cs.violations & V.SECOND_PSEUDOROTATION_TOO_DIFFERENT:
        codes.append('cP1')
    if cs.violations & V.DELTA_TORSION_ANGLE_REJECTED:
        codes.append('DELTA')

    return ';'.join(codes)


def add_category(block: llka.CifDataBlock, name: str, keywords: list[str], rows: list[list]) -> None:
    """Adds a category whose items are filled column by column from the rows."""
    category = block.add_category(name)
    for column, keyword in enumerate(keywords):
        category.add_item(keyword, [row[column] for row in rows])


def extend_cif_data(
    cif_data: llka.CifData,
    attempted_steps: Sequence[llka.AttemptedClassifiedStep],
    avg_confal: llka.AverageConfal,
    entry_id: str,
    steps: Sequence[llka.Structure],
    has_multiple_models: bool,
) -> None:
    # We expect the structure to be in the first block and therefore we extend the first block
    block = cif_data[0]

    summary_rows = []
    ntc_rows = []
    parameters_rows = []
    sugar_rows = []

    assigned_count = 0
    unassigned_but_close_enough_count = 0
    total_count = 0  # Steps that could not be classified at all are not counted

    for stru, attempted in zip(steps, attempted_steps):
        # Attempted step could not have been classified, ignore it
        if attempted.error is not None:
            continue
        cs = attempted.step

        first_atom = stru[0]
        last_atom = stru[-1]
        first_alt_id = residue_alt_id(stru, first_atom.label_seq_id)
        second_alt_id = residue_alt_id(stru, last_atom.label_seq_id)

        total_count += 1
        if cs.violations == llka.ClassificationViolation.OK:
            assigned_count += 1
        elif cs.violations & llka.ClassificationViolation.UNASSIGNED_BUT_CLOSE_ENOUGH:
            unassigned_but_close_enough_count += 1

        step_id = str(total_count)

        summary_rows.append([
            step_id,
            cana_name(cs.assigned_cana),
            ntc_name(cs.assigned_ntc),
            rounded(cs.confal_score.total),
            f'{cs.euclidean_distance_ntc_ideal:.1f}',
            f'{cs.rmsd_to_closest_ntc:.3f}',
            cana_name(cs.closest_cana),
            ntc_name(cs.closest_ntc),
            cs.closest_golden_step or NONE,
        ])

        ntc_rows.append([
            step_id,
            derive_dnatco_step_name(entry_id, has_multiple_models, stru),
            str(first_atom.pdbx_PDB_model_num),
            first_atom.label_entity_id,
            first_atom.label_asym_id,
            str(first_atom.label_seq_id),
            first_atom.auth_comp_id,  # The C example writes auth_comp_id to label_comp_id_1 as well
            first_alt_id or NONE,
            last_atom.label_entity_id,
            last_atom.label_asym_id,
            str(last_atom.label_seq_id),
            last_atom.auth_comp_id,
            second_alt_id or NONE,
            first_atom.auth_asym_id,
            str(first_atom.auth_seq_id),
            last_atom.auth_asym_id,
            str(last_atom.auth_seq_id),
            first_atom.pdbx_PDB_ins_code or NONE,
            last_atom.pdbx_PDB_ins_code or NONE,
        ])

        metrics = cs.metrics
        diffs = cs.differences_from_ntc_averages
        confal = cs.confal_score
        parameters_rows.append([
            step_id,
            # Actual step metrics
            *[deg1(getattr(metrics, t)) for t in TORSIONS],
            f'{metrics.NN:.2f}',
            f'{metrics.CC:.2f}',
            deg1(metrics.mu),
            # Metrics differences from averages of the closest NtC class
            *[diff1(getattr(diffs, t)) for t in TORSIONS],
            f'{diffs.NN:.2f}',
            f'{diffs.CC:.2f}',
            diff1(diffs.mu),
            # Confals
            *[f'{getattr(confal, t):.1f}' for t in TORSIONS],
            f'{confal.NN:.1f}',
            f'{confal.CC:.1f}',
            f'{confal.mu:.1f}',
            details(cs) or NONE,
        ])

        sugar_rows.append([
            step_id,
            deg1(cs.ribose_pseudorotation_1),
            deg1(cs.tau_1),
            llka.sugar_pucker_to_name(cs.sugar_pucker_1, llka.SugarPuckerNameBrevity.VERY_TERSE),
            deg1(cs.ribose_pseudorotation_2),
            deg1(cs.tau_2),
            llka.sugar_pucker_to_name(cs.sugar_pucker_2, llka.SugarPuckerNameBrevity.VERY_TERSE),
            *[deg1(getattr(cs.nu_angles_1, nu)) for nu in NU_ANGLES],
            *[deg1(getattr(cs.nu_angles_2, nu)) for nu in NU_ANGLES],
            *[diff1(getattr(cs.nu_angle_differences_1, nu)) for nu in NU_ANGLES],
            *[diff1(getattr(cs.nu_angle_differences_2, nu)) for nu in NU_ANGLES],
        ])

    overall_rows = [[
        entry_id,
        f'{avg_confal.score:.1f}',
        str(int(avg_confal.percentile)),
        llka.INTERNAL_NTC_VERSION,
        llka.INTERNAL_CANA_VERSION,
        str(total_count),
        str(assigned_count),
        str(total_count - assigned_count),
        str(unassigned_but_close_enough_count),
    ]]

    add_category(block, 'ndb_struct_ntc_overall', OVERALL_KEYWORDS, overall_rows)
    add_category(block, 'ndb_struct_ntc_step', NTC_KEYWORDS, ntc_rows)
    add_category(block, 'ndb_struct_ntc_step_summary', SUMMARY_KEYWORDS, summary_rows)
    add_category(block, 'ndb_struct_ntc_step_parameters', PARAMETERS_KEYWORDS, parameters_rows)
    add_category(block, 'ndb_struct_sugar_step_parameters', SUGAR_KEYWORDS, sugar_rows)

    # Adding items and values marks the Cif data as "tainted". Before the data can be
    # written out, it must pass a validation check that clears the flag.
    try:
        cif_data.detaint()
    except llka.LLKAError as ex:
        print(f'Extended Cif data is invalid: {ex}', file=sys.stderr)
        sys.exit(1)


def main(argv: list[str]) -> int:
    if len(argv) > 1 and argv[1] == '--version':
        print(f'classify_and_write_cif version {VERSION}')
        print(f'DNATCO v5.0, NtC v{llka.INTERNAL_NTC_VERSION}, CANA v{llka.INTERNAL_CANA_VERSION}')
        return 0

    if len(argv) < 2:
        print('No input CIF file', file=sys.stderr)
        return 1

    # Read the structure together with the complete Cif data that will be extended
    try:
        imported = llka.cif_file_to_structure(argv[1], get_cif_data=True)
    except llka.LLKAError as ex:
        print(f'Cannot read CIF file: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
        if isinstance(ex, llka.CifParseError):
            print(f'CIF parsing error: {ex}', file=sys.stderr)
        return 1

    try:
        steps = imported.structure.split_to_dinucleotide_steps()
    except llka.LLKAError as ex:
        print(f'Failed to split the structure into steps: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
        return 1

    ctx = initialize_classification_context()

    try:
        attempted_steps = ctx.classify_steps(steps)
    except llka.LLKAError as ex:
        print(f'Failed to classify steps: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
        return 1
    avg_confal = ctx.average_confal(attempted_steps)

    # NtC assignment is done. Extend Cif data with NtC assignment results
    extend_cif_data(
        imported.cif_data, attempted_steps, avg_confal, imported.entry_id, steps, structure_has_multiple_models(steps)
    )

    try:
        print(imported.cif_data.to_string(pretty=True))
    except llka.LLKAError as ex:
        print(f'Cannot write out CifData as Cif string: {llka.error_to_string(ex.retcode)}', file=sys.stderr)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
