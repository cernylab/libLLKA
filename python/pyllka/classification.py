"""Classification functions for nucleic acid structures."""

import os
from pathlib import Path
from ._pyllka_core import (
    ClassificationContext,
    load_resource_file,
)


def get_assets_dir():
    """Get the directory containing DNATCO classification assets."""
    # First, try package-bundled assets
    package_assets = Path(__file__).parent / "assets"
    if package_assets.exists():
        return str(package_assets)

    # Try environment variable
    env_path = os.environ.get("DNATCO_ASSETS_PATH")
    if env_path and Path(env_path).exists():
        return env_path

    # Try CCP4 installation
    ccp4 = os.environ.get("CCP4")
    if ccp4:
        ccp4_assets = Path(ccp4) / "share" / "dnatco"
        if ccp4_assets.exists():
            return str(ccp4_assets)

    raise FileNotFoundError(
        "Cannot find DNATCO classification assets. "
        "Set DNATCO_ASSETS_PATH environment variable or ensure assets are installed."
    )


def create_classification_context(assets_dir=None):
    """
    Create a classification context for NtC/CANA classification.

    This loads the required CSV parameter files (clusters, golden steps,
    confals, nu angles, and confal percentiles).

    Parameters
    ----------
    assets_dir : str, optional
        Directory containing the classification CSV files.
        If not provided, will search in:
        1. Package bundled assets
        2. DNATCO_ASSETS_PATH environment variable
        3. $CCP4/share/dnatco

    Returns
    -------
    ClassificationContext
        The initialized classification context

    Raises
    ------
    FileNotFoundError
        If asset files cannot be found
    LLKAError
        If assets cannot be loaded

    Examples
    --------
    >>> ctx = create_classification_context()
    >>> results = classify_structure(structure, ctx)
    """
    if assets_dir is None:
        assets_dir = get_assets_dir()

    assets_path = Path(assets_dir)

    # Load required CSV files
    clusters_path = assets_path / "clusters.csv"
    golden_steps_path = assets_path / "golden_steps.csv"
    confals_path = assets_path / "confals.csv"
    nu_angles_path = assets_path / "nu_angles.csv"
    confal_percentiles_path = assets_path / "confal_percentiles.csv"

    # Check all files exist
    for path in [clusters_path, golden_steps_path, confals_path,
                 nu_angles_path, confal_percentiles_path]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required asset file not found: {path}\n"
                f"Asset directory: {assets_dir}"
            )

    # Import the C++ implementation
    from ._pyllka_core import create_classification_context_impl

    # Create context using C++ binding
    return create_classification_context_impl(
        str(clusters_path),
        str(golden_steps_path),
        str(confals_path),
        str(nu_angles_path),
        str(confal_percentiles_path)
    )


def classify_structure(structure, context):
    """
    Classify all dinucleotide steps in a structure.

    Parameters
    ----------
    structure : Structure
        The structure to classify (will be split into dinucleotide steps)
    context : ClassificationContext
        The classification context

    Returns
    -------
    list of ClassifiedStep
        Classification results for each dinucleotide step

    Raises
    ------
    LLKAError
        If classification fails

    Examples
    --------
    >>> ctx = create_classification_context()
    >>> structure = load_structure("structure.cif")
    >>> results = classify_structure(structure, ctx)
    >>> for result in results:
    ...     print(f"NtC: {result.assigned_ntc}, CANA: {result.assigned_cana}")
    """
    from ._pyllka_core import split_to_dinucleotides_impl, classify_steps_multiple_impl

    # Split structure into dinucleotide steps
    steps = split_to_dinucleotides_impl(structure)

    # Classify all steps
    results = classify_steps_multiple_impl(steps, context)

    return results


def classify_step(step, context):
    """
    Classify a single dinucleotide step.

    Parameters
    ----------
    step : Structure
        A dinucleotide step structure (must be a valid dinucleotide step)
    context : ClassificationContext
        The classification context

    Returns
    -------
    ClassifiedStep
        Classification result

    Raises
    ------
    LLKAError
        If classification fails or step is not a valid dinucleotide

    Examples
    --------
    >>> ctx = create_classification_context()
    >>> structure = load_structure("structure.cif")
    >>> steps = split_to_dinucleotides(structure)
    >>> result = classify_step(steps[0], ctx)
    >>> print(f"NtC: {result.assigned_ntc}, RMSD: {result.rmsd_to_closest_ntc:.2f}")
    """
    from ._pyllka_core import classify_step_impl

    return classify_step_impl(step, context)
