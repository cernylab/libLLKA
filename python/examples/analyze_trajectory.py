#!/usr/bin/env python3
"""
Example: Analyze DNA/RNA conformations in an MD trajectory.

This example demonstrates how to use pyllka with MDAnalysis to analyze
nucleic acid conformations (NtC/CANA classifications) across an MD trajectory.

Requirements:
- pyllka (with MDAnalysis support)
- A reference mmCIF file (topology)
- A trajectory file (DCD, XTC, TRR, etc.)

The mmCIF and trajectory must have matching atom counts and ordering.
"""

import sys
import argparse
from pathlib import Path
import numpy as np

try:
    import pyllka
    from pyllka import TrajectoryAnalyzer
    from pyllka import ntc_to_name, cana_to_name
except ImportError:
    print("Error: pyllka not installed or not built with MDAnalysis support")
    print("Install with: pip install pyllka[mdanalysis]")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze nucleic acid conformations in MD trajectory"
    )
    parser.add_argument(
        "topology",
        help="Path to mmCIF topology file (reference structure)"
    )
    parser.add_argument(
        "trajectory",
        help="Path to trajectory file (DCD, XTC, TRR, etc.)"
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="Starting frame index (default: 0)"
    )
    parser.add_argument(
        "--stop",
        type=int,
        default=None,
        help="Stopping frame index (default: all frames)"
    )
    parser.add_argument(
        "--step",
        type=int,
        default=1,
        help="Frame step size (default: 1)"
    )
    parser.add_argument(
        "--assets-dir",
        help="Directory containing DNATCO classification assets"
    )
    parser.add_argument(
        "--output",
        help="Output file for results (CSV format)"
    )
    parser.add_argument(
        "--timeseries",
        type=int,
        help="Extract timeseries for specific dinucleotide step index"
    )

    args = parser.parse_args()

    # Validate input files
    if not Path(args.topology).exists():
        print(f"Error: Topology file not found: {args.topology}")
        sys.exit(1)

    if not Path(args.trajectory).exists():
        print(f"Error: Trajectory file not found: {args.trajectory}")
        sys.exit(1)

    # Create analyzer
    print(f"Loading topology: {args.topology}")
    print(f"Loading trajectory: {args.trajectory}")

    try:
        analyzer = TrajectoryAnalyzer(
            args.topology,
            args.trajectory,
            assets_dir=args.assets_dir
        )
    except Exception as e:
        print(f"Error initializing analyzer: {e}")
        sys.exit(1)

    print(f"Universe contains {analyzer.n_atoms} atoms")
    print(f"Trajectory contains {len(analyzer.universe.trajectory)} frames")

    # Analyze trajectory
    if args.timeseries is not None:
        print(f"\nExtracting timeseries for step {args.timeseries}")
        analyze_timeseries(analyzer, args)
    else:
        print(f"\nAnalyzing frames {args.start} to {args.stop or 'end'} "
              f"(step={args.step})")
        analyze_all_frames(analyzer, args)


def analyze_all_frames(analyzer, args):
    """Analyze all frames and print/save results."""
    try:
        results = analyzer.analyze_trajectory(
            start=args.start,
            stop=args.stop,
            step=args.step,
            verbose=True
        )
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)

    print(f"\nAnalyzed {len(results)} frames")

    # Print summary statistics
    print("\n=== Summary ===")
    if results:
        n_steps = len(results[0])
        print(f"Number of dinucleotide steps per frame: {n_steps}")

        # Count NtC occurrences across all frames
        ntc_counts = {}
        cana_counts = {}

        for frame_results in results:
            for step_result in frame_results:
                ntc = ntc_to_name(step_result.assigned_ntc)
                cana = cana_to_name(step_result.assigned_cana)

                ntc_counts[ntc] = ntc_counts.get(ntc, 0) + 1
                cana_counts[cana] = cana_counts.get(cana, 0) + 1

        print("\nMost common NtC classes:")
        for ntc, count in sorted(ntc_counts.items(), key=lambda x: -x[1])[:10]:
            print(f"  {ntc}: {count}")

        print("\nCANA distribution:")
        for cana, count in sorted(cana_counts.items(), key=lambda x: -x[1]):
            print(f"  {cana}: {count}")

    # Save results if requested
    if args.output:
        save_results_csv(results, args.output)
        print(f"\nResults saved to {args.output}")


def analyze_timeseries(analyzer, args):
    """Extract and analyze timeseries for a specific step."""
    try:
        timeseries = analyzer.get_step_timeseries(
            step_index=args.timeseries,
            start=args.start,
            stop=args.stop,
            step=args.step
        )
    except Exception as e:
        print(f"Error extracting timeseries: {e}")
        sys.exit(1)

    print(f"Extracted timeseries for {len(timeseries['ntc'])} frames")

    # Print statistics
    print("\n=== Timeseries Statistics ===")
    print(f"Mean RMSD: {np.mean(timeseries['rmsd']):.3f} Å")
    print(f"Std RMSD: {np.std(timeseries['rmsd']):.3f} Å")

    # Count transitions
    ntc_names = [ntc_to_name(ntc) for ntc in timeseries['ntc']]
    unique_ntcs = set(ntc_names)
    print(f"\nNtC classes observed: {len(unique_ntcs)}")
    for ntc in sorted(unique_ntcs):
        count = ntc_names.count(ntc)
        print(f"  {ntc}: {count} ({100*count/len(ntc_names):.1f}%)")

    # Count conformational transitions
    transitions = 0
    for i in range(1, len(ntc_names)):
        if ntc_names[i] != ntc_names[i-1]:
            transitions += 1
    print(f"\nConformational transitions: {transitions}")

    # Save timeseries if requested
    if args.output:
        save_timeseries_csv(timeseries, args.output)
        print(f"\nTimeseries saved to {args.output}")


def save_results_csv(results, output_path):
    """Save all frame results to CSV."""
    import csv

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'frame', 'step_index', 'ntc', 'cana',
            'rmsd', 'violations'
        ])

        for frame_idx, frame_results in enumerate(results):
            for step_idx, step_result in enumerate(frame_results):
                writer.writerow([
                    frame_idx,
                    step_idx,
                    ntc_to_name(step_result.assigned_ntc),
                    cana_to_name(step_result.assigned_cana),
                    f"{step_result.rmsd_to_closest_ntc:.4f}",
                    step_result.violations
                ])


def save_timeseries_csv(timeseries, output_path):
    """Save timeseries to CSV."""
    import csv

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['frame', 'ntc', 'cana', 'rmsd', 'violations'])

        for i in range(len(timeseries['ntc'])):
            writer.writerow([
                i,
                ntc_to_name(timeseries['ntc'][i]),
                cana_to_name(timeseries['cana'][i]),
                f"{timeseries['rmsd'][i]:.4f}",
                timeseries['violations'][i]
            ])


if __name__ == '__main__':
    main()
