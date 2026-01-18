"""I/O functions for loading and saving structures."""

from ._pyllka_core import load_structure_from_file


def load_structure(path):
    """
    Load a structure from an mmCIF file.

    Parameters
    ----------
    path : str
        Path to the mmCIF file

    Returns
    -------
    Structure
        The loaded structure

    Raises
    ------
    LLKAError
        If the file cannot be loaded or parsed

    Examples
    --------
    >>> structure = load_structure("structure.cif")
    >>> print(f"Loaded {structure.n_atoms} atoms")
    """
    return load_structure_from_file(path)


def save_structure(structure, path):
    """
    Save a structure to an mmCIF file.

    Parameters
    ----------
    structure : Structure
        The structure to save
    path : str
        Output file path

    Raises
    ------
    LLKAError
        If the file cannot be written

    Examples
    --------
    >>> save_structure(structure, "output.cif")
    """
    raise NotImplementedError("save_structure not yet implemented")
