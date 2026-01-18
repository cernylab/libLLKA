"""Utility functions for structure manipulation."""

import numpy as np


def update_coordinates(structure, coords):
    """
    Update structure coordinates from a numpy array.

    This is the key function for trajectory analysis - it efficiently updates
    the coordinates of an existing structure without recreating it.

    Parameters
    ----------
    structure : Structure
        The structure to update
    coords : array_like
        Coordinates array of shape (n_atoms, 3)
        Can be from MDAnalysis: universe.atoms.positions

    Raises
    ------
    ValueError
        If coords array has wrong shape
    LLKAError
        If update fails

    Examples
    --------
    >>> import MDAnalysis as mda
    >>> u = mda.Universe("topology.cif", "trajectory.dcd")
    >>> structure = load_structure("topology.cif")
    >>> for ts in u.trajectory:
    ...     update_coordinates(structure, u.atoms.positions)
    ...     # Now structure has coordinates from current frame
    """
    coords = np.asarray(coords, dtype=np.float64)

    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(
            f"Coordinates must have shape (n_atoms, 3), got {coords.shape}"
        )

    if coords.shape[0] != structure.n_atoms:
        raise ValueError(
            f"Coordinate count mismatch: structure has {structure.n_atoms} atoms, "
            f"but coords has {coords.shape[0]} rows"
        )

    # Use the C++ implementation
    structure.update_coordinates(coords)


def split_to_dinucleotides(structure):
    """
    Split a structure into dinucleotide steps.

    Parameters
    ----------
    structure : Structure
        The structure to split

    Returns
    -------
    list of Structure
        List of dinucleotide structures

    Raises
    ------
    LLKAError
        If splitting fails

    Examples
    --------
    >>> structure = load_structure("dna.cif")
    >>> steps = split_to_dinucleotides(structure)
    >>> print(f"Found {len(steps)} dinucleotide steps")
    """
    from ._pyllka_core import split_to_dinucleotides_impl

    return split_to_dinucleotides_impl(structure)


def get_coordinates(structure):
    """
    Get structure coordinates as a numpy array.

    Parameters
    ----------
    structure : Structure
        The structure

    Returns
    -------
    np.ndarray
        Coordinates array of shape (n_atoms, 3)

    Examples
    --------
    >>> structure = load_structure("structure.cif")
    >>> coords = get_coordinates(structure)
    >>> print(coords.shape)  # (n_atoms, 3)
    """
    return structure.get_coordinates()
