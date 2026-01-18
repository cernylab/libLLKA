#!/usr/bin/env python3
"""
Example: Basic structure classification (without trajectory).

This example demonstrates how to use pyllka to classify nucleic acid
conformations in a single structure from an mmCIF file.
"""

import sys
import math
import pyllka

def main():
    if len(sys.argv) < 2:
        print("Usage: python basic_classification.py <structure.cif>")
        sys.exit(1)

    structure_path = sys.argv[1]

    # Load structure
    print(f"Loading structure from {structure_path}")
    structure = pyllka.load_structure(structure_path)
    print(f"Loaded {structure.n_atoms} atoms")

    # Print some atom information
    print("\nFirst 5 atoms:")
    for i in range(min(5, structure.n_atoms)):
        atom = structure[i]
        print(f"  {i}: {atom.label_atom_id} {atom.label_comp_id} {atom.label_seq_id} "
              f"({atom.coords.x:.2f}, {atom.coords.y:.2f}, {atom.coords.z:.2f})")

    # Get coordinates as numpy array
    coords = structure.get_coordinates()
    print(f"\nCoordinates shape: {coords.shape}")
    print(f"Coordinate range: x=[{coords[:,0].min():.2f}, {coords[:,0].max():.2f}], "
          f"y=[{coords[:,1].min():.2f}, {coords[:,1].max():.2f}], "
          f"z=[{coords[:,2].min():.2f}, {coords[:,2].max():.2f}]")

    # Initialize classification context
    print("\nInitializing classification context...")
    try:
        ctx = pyllka.create_classification_context()
        print("Classification context created")
    except NotImplementedError as e:
        print(f"Note: {e}")
        print("\nClassification not yet fully implemented - need to extend bindings.cpp")
        return

    # Split into dinucleotides
    print("\nSplitting structure into dinucleotide steps...")
    try:
        steps = pyllka.split_to_dinucleotides(structure)
        print(f"Found {len(steps)} dinucleotide steps")
    except NotImplementedError as e:
        print(f"Note: {e}")
        return

    # Classify all steps
    print("\nClassifying dinucleotide steps...")
    results = pyllka.classify_structure(structure, ctx)

    # Print results
    print(f"\n=== Classification Results ===")
    for i, result in enumerate(results):
        print(f"Step {i+1}:")
        print(f"  Assigned NtC: {pyllka.ntc_to_name(result.assigned_ntc)}")
        print(f"  Assigned CANA: {pyllka.cana_to_name(result.assigned_cana)}")
        print(f"  Closest NtC: {pyllka.ntc_to_name(result.closest_ntc)}")
        print(f"  Closest CANA: {pyllka.cana_to_name(result.closest_cana)}")
        print(f"  RMSD to closest: {result.rmsd_to_closest_ntc:.3f} Å")
        print(f"  Euclidean distance: {result.euclidean_distance_ntc_ideal:.3f}")

        # Print pseudorotation values (convert from radians to degrees)
        print(f"  Ribose pseudorotation 1: {math.degrees(result.ribose_pseudorotation_1):.2f}°")
        print(f"  Ribose pseudorotation 2: {math.degrees(result.ribose_pseudorotation_2):.2f}°")
        print(f"  Sugar pucker 1: {pyllka.sugar_pucker_to_name(result.sugar_pucker_1)}")
        print(f"  Sugar pucker 2: {pyllka.sugar_pucker_to_name(result.sugar_pucker_2)}")
        print(f"  Tau 1: {math.degrees(result.tau_1):.2f}°")
        print(f"  Tau 2: {math.degrees(result.tau_2):.2f}°")

        # Check if fully assigned (no violations)
        if result.violations == pyllka.CLASSIFICATION_OK:
            print(f"  ✓ Fully assigned (no violations)")
            print(f"  Confal score: {result.confal_score.total:.3f}")
        else:
            print(f"  ⚠ Has violations (code: {result.violations}):")
            # Decode violation flags
            for flag_name in dir(pyllka):
                if flag_name.startswith('CLASSIFICATION_E_'):
                    flag_value = getattr(pyllka, flag_name)
                    if result.violations & flag_value:
                        print(f"    - {flag_name}")
        print()  # blank line between steps


if __name__ == '__main__':
    main()
