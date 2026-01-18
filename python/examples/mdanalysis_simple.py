#!/usr/bin/env python3
"""
Example: Simple MDAnalysis integration demo.

This example shows the basic workflow of updating libLLKA structures
with coordinates from MDAnalysis trajectory frames.
"""

import sys
import numpy as np

try:
    import MDAnalysis as mda
except ImportError:
    print("Error: MDAnalysis not installed")
    print("Install with: pip install MDAnalysis")
    sys.exit(1)

import pyllka


def main():
    if len(sys.argv) < 3:
        print("Usage: python mdanalysis_simple.py <topology.cif> <trajectory.dcd>")
        sys.exit(1)

    topology_path = sys.argv[1]
    trajectory_path = sys.argv[2]

    # Load with MDAnalysis
    print(f"Loading with MDAnalysis:")
    print(f"  Topology: {topology_path}")
    print(f"  Trajectory: {trajectory_path}")
    u = mda.Universe(topology_path, trajectory_path)
    print(f"  {u.atoms.n_atoms} atoms, {len(u.trajectory)} frames")

    # Load with pyllka
    print(f"\nLoading with pyllka:")
    structure = pyllka.load_structure(topology_path)
    print(f"  {structure.n_atoms} atoms")

    # Verify atom counts match
    if u.atoms.n_atoms != structure.n_atoms:
        print(f"\nERROR: Atom count mismatch!")
        print(f"  MDAnalysis: {u.atoms.n_atoms}")
        print(f"  pyllka: {structure.n_atoms}")
        sys.exit(1)

    print(f"\n✓ Atom counts match: {structure.n_atoms} atoms")

    # Demonstrate coordinate updates
    print(f"\n=== Demonstrating coordinate updates ===")

    # Get initial coordinates from pyllka structure
    initial_coords = structure.get_coordinates()
    print(f"Initial pyllka coords (first atom): {initial_coords[0]}")

    # Iterate through first few frames
    n_frames_to_show = min(5, len(u.trajectory))
    print(f"\nIterating through first {n_frames_to_show} frames:")

    for i, ts in enumerate(u.trajectory[:n_frames_to_show]):
        # Get coordinates from MDAnalysis
        mda_coords = u.atoms.positions

        # Update pyllka structure
        pyllka.update_coordinates(structure, mda_coords)

        # Verify update worked
        updated_coords = structure.get_coordinates()

        # Check that coordinates match
        if not np.allclose(updated_coords, mda_coords):
            print(f"  Frame {i}: ERROR - coordinates don't match!")
        else:
            print(f"  Frame {i}: ✓ Coordinates updated successfully")
            print(f"    First atom: {updated_coords[0]}")

    # Demonstrate that we can access atom metadata
    print(f"\n=== Atom metadata (from pyllka structure) ===")
    print("First 3 atoms:")
    for i in range(min(3, structure.n_atoms)):
        atom = structure[i]
        print(f"  [{i}] {atom.label_atom_id:4s} {atom.label_comp_id:3s} "
              f"{atom.label_seq_id:4d} chain {atom.label_asym_id}")

    print("\n✓ Demo completed successfully!")
    print("\nNext steps:")
    print("  1. Implement classification bindings in bindings.cpp")
    print("  2. Use TrajectoryAnalyzer for full trajectory analysis")
    print("  3. See analyze_trajectory.py for complete example")


if __name__ == '__main__':
    main()
