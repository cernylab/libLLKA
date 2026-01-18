#!/usr/bin/env python3
"""
Example: Calculate similarity and connectivity for nucleic acid structures.

This example demonstrates how to:
1. Load a structure and classify dinucleotide steps
2. Calculate similarity between steps and reference NtCs
3. Calculate connectivity between consecutive steps
4. Output results as JSON

This is the Python equivalent of the C++ similarity_connectivity example.
"""

import sys
import json
import argparse
import pyllka


# All NtCs ordered as in the C++ example
ALL_NTCS = [
    pyllka.NtC.AA00, pyllka.NtC.AA02, pyllka.NtC.AA03, pyllka.NtC.AA04, pyllka.NtC.AA08,
    pyllka.NtC.AA09, pyllka.NtC.AA01, pyllka.NtC.AA05, pyllka.NtC.AA06, pyllka.NtC.AA10,
    pyllka.NtC.AA11, pyllka.NtC.AA07, pyllka.NtC.AA12, pyllka.NtC.AA13,
    pyllka.NtC.AB01, pyllka.NtC.AB02, pyllka.NtC.AB03, pyllka.NtC.AB04, pyllka.NtC.AB05,
    pyllka.NtC.BA01, pyllka.NtC.BA05, pyllka.NtC.BA09, pyllka.NtC.BA08, pyllka.NtC.BA10,
    pyllka.NtC.BA13, pyllka.NtC.BA16, pyllka.NtC.BA17,
    pyllka.NtC.BB00, pyllka.NtC.BB01, pyllka.NtC.BB17, pyllka.NtC.BB02, pyllka.NtC.BB03,
    pyllka.NtC.BB11, pyllka.NtC.BB16, pyllka.NtC.BB04, pyllka.NtC.BB05, pyllka.NtC.BB07,
    pyllka.NtC.BB08, pyllka.NtC.BB10, pyllka.NtC.BB12, pyllka.NtC.BB13, pyllka.NtC.BB14,
    pyllka.NtC.BB15, pyllka.NtC.BB20,
    pyllka.NtC.IC01, pyllka.NtC.IC02, pyllka.NtC.IC03, pyllka.NtC.IC04, pyllka.NtC.IC05,
    pyllka.NtC.IC06, pyllka.NtC.IC07,
    pyllka.NtC.OP01, pyllka.NtC.OP02, pyllka.NtC.OP03, pyllka.NtC.OP04, pyllka.NtC.OP05,
    pyllka.NtC.OP06, pyllka.NtC.OP07, pyllka.NtC.OP08, pyllka.NtC.OP09, pyllka.NtC.OP10,
    pyllka.NtC.OP11, pyllka.NtC.OP12, pyllka.NtC.OP13, pyllka.NtC.OP14, pyllka.NtC.OP15,
    pyllka.NtC.OP16, pyllka.NtC.OP17, pyllka.NtC.OP18, pyllka.NtC.OP19, pyllka.NtC.OP20,
    pyllka.NtC.OP21, pyllka.NtC.OP22, pyllka.NtC.OP23, pyllka.NtC.OP24, pyllka.NtC.OP25,
    pyllka.NtC.OP26, pyllka.NtC.OP27, pyllka.NtC.OP28, pyllka.NtC.OP29, pyllka.NtC.OP30,
    pyllka.NtC.OP31, pyllka.NtC.OPS1, pyllka.NtC.OP1S,
    pyllka.NtC.AAS1, pyllka.NtC.AB1S, pyllka.NtC.AB2S,
    pyllka.NtC.BB1S, pyllka.NtC.BB2S, pyllka.NtC.BBS1,
    pyllka.NtC.ZZ01, pyllka.NtC.ZZ02, pyllka.NtC.ZZ1S, pyllka.NtC.ZZ2S, pyllka.NtC.ZZS1, pyllka.NtC.ZZS2,
]


def calculate_similarity(steps, classified_steps, step_id=None, ntc_filter=None, rmsd_cutoff=None):
    """
    Calculate similarity between steps and reference NtCs.

    Args:
        steps: List of Structure objects (dinucleotide steps)
        classified_steps: List of ClassifiedStep objects
        step_id: Optional specific step ID to analyze (1-indexed), None for all steps
        ntc_filter: Optional NtC to filter results
        rmsd_cutoff: Optional RMSD cutoff to filter results

    Returns:
        Dictionary with similarity data
    """
    results = {}

    # Determine which steps to process
    if step_id is not None:
        step_indices = [step_id - 1]  # Convert to 0-indexed
    else:
        step_indices = range(len(steps))

    for idx in step_indices:
        step = steps[idx]
        cs = classified_steps[idx]
        step_key = f"step_{idx + 1}"

        # If filtering by NtC, check if this step matches
        if ntc_filter is not None and cs.assigned_ntc != ntc_filter:
            continue

        step_data = {
            "step_id": idx + 1,
            "assigned_ntc": pyllka.ntc_to_name(cs.assigned_ntc),
            "similarities": []
        }

        # Calculate similarity against all NtCs
        # Note: We calculate one at a time to avoid issues with batch processing
        for ntc in ALL_NTCS:
            try:
                sim = pyllka.measure_step_similarity_ntc(step, ntc)

                # Apply RMSD cutoff if specified
                if rmsd_cutoff is not None and sim.rmsd > rmsd_cutoff:
                    continue

                step_data["similarities"].append({
                    "ntc": pyllka.ntc_to_name(ntc),
                    "rmsd": round(sim.rmsd, 3),
                    "euclidean_distance": round(sim.euclidean_distance, 1)
                })

            except pyllka.LLKAError as e:
                # Skip NtCs that can't be measured (e.g., wrong base composition)
                continue

        results[step_key] = step_data

    return results


def calculate_connectivity(steps, classified_steps, step_id=None, ntc_filter=None,
                          prev_ntc_filter=None, next_ntc_filter=None, distance_cutoff=None):
    """
    Calculate connectivity between consecutive steps.

    This function follows the C++ example's approach: for each step, iterate through
    all possible NtCs that the current step could be, and for each, calculate connectivity
    with all possible NtCs for the neighboring steps.

    Args:
        steps: List of Structure objects (dinucleotide steps)
        classified_steps: List of ClassifiedStep objects
        step_id: Optional specific step ID to analyze (1-indexed), None for all steps
        ntc_filter: Optional NtC to filter current step (test only this NtC for current step)
        prev_ntc_filter: Optional NtC to filter previous step connections
        next_ntc_filter: Optional NtC to filter next step connections
        distance_cutoff: Optional distance cutoff to filter results

    Returns:
        Dictionary with connectivity data
    """
    results = {}

    # Determine which steps to process
    if step_id is not None:
        step_indices = [step_id - 1]  # Convert to 0-indexed
    else:
        step_indices = range(len(steps))

    # Determine which NtCs to test for the current step
    if ntc_filter is not None:
        current_ntcs_to_check = [ntc_filter]
    else:
        current_ntcs_to_check = ALL_NTCS

    for idx in step_indices:
        cs = classified_steps[idx]

        # For each possible NtC that the current step could be
        for current_ntc in current_ntcs_to_check:
            step_key = f"step_{idx + 1}_as_{pyllka.ntc_to_name(current_ntc)}"

            step_data = {
                "step_id": idx + 1,
                "tested_as_ntc": pyllka.ntc_to_name(current_ntc),
                "assigned_ntc": pyllka.ntc_to_name(cs.assigned_ntc),
                "prev_step_connectivity": [],
                "next_step_connectivity": []
            }

            # Calculate connectivity with previous step
            if idx > 0:
                # Determine which NtCs to test for previous step
                if prev_ntc_filter is not None:
                    prev_ntcs_to_check = [prev_ntc_filter]
                else:
                    prev_ntcs_to_check = ALL_NTCS

                for prev_ntc in prev_ntcs_to_check:
                    try:
                        conn = pyllka.measure_step_connectivity_ntcs(
                            steps[idx - 1], prev_ntc,
                            steps[idx], current_ntc
                        )

                        max_dist = max(conn.c5_prime_distance, conn.o3_prime_distance)

                        # Apply distance cutoff if specified
                        if distance_cutoff is not None and max_dist > distance_cutoff:
                            continue

                        step_data["prev_step_connectivity"].append({
                            "ntc": pyllka.ntc_to_name(prev_ntc),
                            "c5_prime_distance": round(conn.c5_prime_distance, 3),
                            "o3_prime_distance": round(conn.o3_prime_distance, 3)
                        })

                    except pyllka.LLKAError as e:
                        # Skip NtCs that can't be measured
                        continue

            # Calculate connectivity with next step
            if idx < len(steps) - 1:
                # Determine which NtCs to test for next step
                if next_ntc_filter is not None:
                    next_ntcs_to_check = [next_ntc_filter]
                else:
                    next_ntcs_to_check = ALL_NTCS

                for next_ntc in next_ntcs_to_check:
                    try:
                        conn = pyllka.measure_step_connectivity_ntcs(
                            steps[idx], current_ntc,
                            steps[idx + 1], next_ntc
                        )

                        max_dist = max(conn.c5_prime_distance, conn.o3_prime_distance)

                        # Apply distance cutoff if specified
                        if distance_cutoff is not None and max_dist > distance_cutoff:
                            continue

                        step_data["next_step_connectivity"].append({
                            "ntc": pyllka.ntc_to_name(next_ntc),
                            "c5_prime_distance": round(conn.c5_prime_distance, 3),
                            "o3_prime_distance": round(conn.o3_prime_distance, 3)
                        })

                    except pyllka.LLKAError as e:
                        # Skip NtCs that can't be measured
                        continue

            # Only add to results if there's connectivity data
            if step_data["prev_step_connectivity"] or step_data["next_step_connectivity"]:
                results[step_key] = step_data

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Calculate similarity and connectivity for nucleic acid structures"
    )
    parser.add_argument("input", metavar="INPUT", help="Input mmCIF file")
    parser.add_argument("-s", "--similarity", metavar="FILE",
                       help="Calculate similarity and write to FILE (JSON format)")
    parser.add_argument("-c", "--connectivity", metavar="FILE",
                       help="Calculate connectivity and write to FILE (JSON format)")
    parser.add_argument("-x", "--step", type=int, metavar="N",
                       help="Analyze only step N (1-indexed)")
    parser.add_argument("-n", "--ntc", metavar="NTC",
                       help="Filter by NtC type (e.g., AA01)")
    parser.add_argument("-b", "--prev-ntc", metavar="NTC",
                       help="Filter previous step connectivity by NtC type")
    parser.add_argument("-a", "--next-ntc", metavar="NTC",
                       help="Filter next step connectivity by NtC type")
    parser.add_argument("-r", "--rmsd", type=float, metavar="RMSD",
                       help="RMSD cutoff for similarity (Angstroms)")
    parser.add_argument("-d", "--distance-cutoff", type=float, metavar="DIST",
                       help="Distance cutoff for connectivity (Angstroms)")
    parser.add_argument("-v", "--version", action="store_true",
                       help="Show version information")

    args = parser.parse_args()

    if args.version:
        print(f"Similarity & Connectivity (Python version)")
        print(f"DNATCO v5.0, NtC v{pyllka.NTC_VERSION}, CANA v{pyllka.CANA_VERSION}")
        return 0

    if not args.similarity and not args.connectivity:
        parser.error("At least one of --similarity or --connectivity must be specified")

    # Parse NtC filters if provided
    ntc_filter = None
    if args.ntc:
        ntc_filter = pyllka.name_to_ntc(args.ntc)
        if ntc_filter == pyllka.NtC.INVALID_NTC:
            print(f"Error: Invalid NtC '{args.ntc}'", file=sys.stderr)
            return 1

    prev_ntc_filter = None
    if args.prev_ntc:
        prev_ntc_filter = pyllka.name_to_ntc(args.prev_ntc)
        if prev_ntc_filter == pyllka.NtC.INVALID_NTC:
            print(f"Error: Invalid prev NtC '{args.prev_ntc}'", file=sys.stderr)
            return 1

    next_ntc_filter = None
    if args.next_ntc:
        next_ntc_filter = pyllka.name_to_ntc(args.next_ntc)
        if next_ntc_filter == pyllka.NtC.INVALID_NTC:
            print(f"Error: Invalid next NtC '{args.next_ntc}'", file=sys.stderr)
            return 1

    # Validate step ID
    if args.step is not None and args.step <= 0:
        print("Error: Step ID must be positive", file=sys.stderr)
        return 1

    # Load structure
    print(f"Loading structure from {args.input}...", file=sys.stderr)
    try:
        imported = pyllka.load_structure_with_cif(args.input)
        structure = imported.get_structure()
    except Exception as e:
        print(f"Error loading structure: {e}", file=sys.stderr)
        return 1

    print(f"Number of atoms: {len(structure)}", file=sys.stderr)

    # Split to dinucleotide steps
    print("Splitting structure to dinucleotide steps...", file=sys.stderr)
    steps = pyllka.split_to_dinucleotides(structure)
    print(f"Number of steps: {len(steps)}", file=sys.stderr)

    if args.step is not None and args.step > len(steps):
        print(f"Error: Step ID {args.step} is larger than number of steps ({len(steps)})", file=sys.stderr)
        return 1

    # Create classification context and classify
    print("Loading classification data...", file=sys.stderr)
    ctx = pyllka.create_classification_context()

    print("Classifying steps...", file=sys.stderr)
    classified_steps = pyllka.classify_structure(structure, ctx)

    assigned = sum(1 for cs in classified_steps if cs.violations == pyllka.CLASSIFICATION_OK)
    print(f"Assigned steps: {assigned}/{len(classified_steps)}", file=sys.stderr)

    # Calculate similarity if requested
    if args.similarity:
        print("Calculating similarity...", file=sys.stderr)
        similarity_results = calculate_similarity(
            steps, classified_steps,
            step_id=args.step,
            ntc_filter=ntc_filter,
            rmsd_cutoff=args.rmsd
        )

        with open(args.similarity, 'w') as f:
            json.dump(similarity_results, f, indent=2)
        print(f"Similarity results written to {args.similarity}", file=sys.stderr)

    # Calculate connectivity if requested
    if args.connectivity:
        print("Calculating connectivity...", file=sys.stderr)
        connectivity_results = calculate_connectivity(
            steps, classified_steps,
            step_id=args.step,
            ntc_filter=ntc_filter,
            prev_ntc_filter=prev_ntc_filter,
            next_ntc_filter=next_ntc_filter,
            distance_cutoff=args.distance_cutoff
        )

        with open(args.connectivity, 'w') as f:
            json.dump(connectivity_results, f, indent=2)
        print(f"Connectivity results written to {args.connectivity}", file=sys.stderr)

    print("Done!", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
