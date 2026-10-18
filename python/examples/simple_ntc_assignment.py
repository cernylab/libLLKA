#!/usr/bin/env python3
"""Assigns NtC classes to a nucleic acid structure and prints the results.

This is the Python counterpart of the ``simple_NtC_assignment`` C example.
It accepts a mmCIF file with the structure as a command line parameter and:

1) Reads the mmCIF file and creates an internal representation of the structure.
2) Splits the structure into a list of dinucleotide steps.
3) Classifies each step, meaning that it attempts to assign NtC class to each
   step and calculates the accompanying information.
4) Prints out the result.

The classification parameter files (clusters.csv, confals.csv, golden_steps.csv,
nu_angles.csv and confal_percentiles.csv) must be present in the working directory.

Usage: simple_ntc_assignment.py <input_file.cif>
"""

from __future__ import annotations

import math
import sys

import llka
from dnatco_step_name import derive_dnatco_step_name, structure_has_multiple_models

MAX_CLOSE_ENOUGH_RMSD = 0.5


def deg2rad(deg: float) -> float:
    # Same operation order as LLKA_deg2rad() so that the limits are bit-identical to the C example
    return math.pi * deg / 180.0


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

    # Load classification parameters from files
    resources = {}
    for key, path, kind, description in (
        ('golden_steps', 'golden_steps.csv', llka.ResourceType.GOLDEN_STEPS, 'golden steps'),
        ('clusters', 'clusters.csv', llka.ResourceType.CLUSTERS, 'clusters'),
        ('confals', 'confals.csv', llka.ResourceType.CONFALS, 'confals'),
        ('cluster_nu_angles', 'nu_angles.csv', llka.ResourceType.AVERAGE_NU_ANGLES, 'average nu angles'),
        ('confal_percentiles', 'confal_percentiles.csv', llka.ResourceType.CONFAL_PERCENTILES, 'confal percentiles'),
    ):
        try:
            resources[key] = llka.load_resource_file(path, kind)
        except llka.LLKAError as ex:
            print(f'Failed to load {description}: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
            sys.exit(1)

    try:
        # The context keeps internal copies of the parameters, the loaded resources are released when they go out of scope
        return llka.ClassificationContext(**resources, limits=limits, max_close_enough_rmsd=MAX_CLOSE_ENOUGH_RMSD)
    except llka.LLKAError as ex:
        print(f'Failed to initialize classification context: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
        sys.exit(1)


def ntc_name(ntc: llka.NtC | None) -> str:
    return llka.ntc_to_name(llka.NtC.INVALID if ntc is None else ntc)


def cana_name(cana: llka.CANA | None) -> str:
    return llka.cana_to_name(llka.CANA.INVALID if cana is None else cana)


def ljust_bytes(s: str, width: int) -> str:
    """Pads a string to a width in bytes the way printf("%-*s") does with UTF-8 strings."""
    return s + ' ' * max(0, width - len(s.encode('utf-8')))


def annotate_structure(entry_id: str, steps: list[llka.Structure]) -> None:
    has_multiple_models = structure_has_multiple_models(steps)  # Needed to derive DNATCO step name

    # The classification process is highly parametrized and requires external data to operate
    ctx = initialize_classification_context()

    try:
        attempted_steps = ctx.classify_steps(steps)
    except llka.LLKAError as ex:
        print(f'Failed to classify steps: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
        sys.exit(1)

    # The order of classified steps is the same as the order of the input steps
    for idx, (stru, attempted) in enumerate(zip(steps, attempted_steps)):
        # The step was rejected by the classification routine.
        # The structure may not be a valid NtC step or it has unacceptable geometry
        if attempted.error is not None:
            print(f'WARN: Could not classify step {idx}: {llka.error_to_string(attempted.error.retcode)}')
            continue

        step = attempted.step
        dnatco_name = derive_dnatco_step_name(entry_id, has_multiple_models, stru)

        # Print basic step classification info
        print(
            f'* Step {idx + 1}\n'
            f'  DNATCO name:                       {dnatco_name}\n'
            f'  Assigned NtC:  {ntc_name(step.assigned_ntc):<4}, Closest NtC:  {ntc_name(step.closest_ntc):<4}\n'
            f'  Assigned CANA: {cana_name(step.assigned_cana):<4}, Closest CANA: {cana_name(step.closest_cana):<4}\n'
            f'  Golden step DNATCO name:           {step.closest_golden_step}\n'
            f'  RMSD to closest NtC:               {step.rmsd_to_closest_ntc:g}'
        )

        if not step.is_assigned:
            # If a step could not have been fully classified, show why.
            names = ''.join(
                llka.classification_violation_to_name(1 << bit) + ' '
                for bit in range(32)
                if step.violations & (1 << bit)
            )
            print('  Classification violations:         ' + names)

        # Show confal scores. The total score is zero if the step has not been fully classified.
        score = step.confal_score
        lines = ['  Confal scores:']
        for torsion in llka.DinucleotideTorsion:
            value = getattr(score, torsion.name.lower())
            lines.append(f'  {ljust_bytes(llka.dinucleotide_torsion_name(torsion, True), 6)} {int(value + 0.5)}')
        lines.append(f'  {ljust_bytes(llka.cross_residue_metric_name(llka.CrossResidueMetric.DIST_CC, True), 5)} {int(score.CC + 0.5)}')
        lines.append(f'  {ljust_bytes(llka.cross_residue_metric_name(llka.CrossResidueMetric.DIST_NN, True), 5)} {int(score.NN + 0.5)}')
        lines.append(f'  {ljust_bytes(llka.cross_residue_metric_name(llka.CrossResidueMetric.TOR_MU, True), 6)} {int(score.mu + 0.5)}')
        lines.append(f'  total {int(score.total + 0.5)}')
        print('\n'.join(lines))
        print()


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print('No input CIF file', file=sys.stderr)
        return 1

    # Read CIF file and create a structure from its contents
    try:
        imported = llka.cif_file_to_structure(argv[1])
    except llka.LLKAError as ex:
        print(f'Cannot read CIF file: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
        # If the CIF file looks malformed, display more information about why the parsing failed
        if isinstance(ex, llka.CifParseError):
            print(f'CIF parsing error: {ex}', file=sys.stderr)
        return 1

    # The imported structure contains the entire structure. In order to do the NtC assignment,
    # the structure has to be split up to dinucleotide steps.
    try:
        steps = imported.structure.split_to_dinucleotide_steps()
    except llka.LLKAError as ex:
        print(f'Failed to split the structure into steps: {llka.error_to_string(ex.retcode)}', file=sys.stderr)
        return 1

    annotate_structure(imported.entry_id, steps)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
