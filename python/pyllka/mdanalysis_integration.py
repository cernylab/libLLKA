"""MDAnalysis integration for trajectory analysis."""

try:
    import MDAnalysis as mda
except ImportError:
    raise ImportError(
        "MDAnalysis is required for trajectory analysis. "
        "Install with: pip install MDAnalysis"
    )

import numpy as np
from .io import load_structure
from .classification import create_classification_context, classify_structure
from .utils import update_coordinates


class TrajectoryAnalyzer:
    """
    Analyze nucleic acid conformations in MD trajectories.

    This class provides a convenient interface for analyzing trajectories
    with MDAnalysis + libLLKA.

    Parameters
    ----------
    topology_path : str
        Path to mmCIF file serving as topology (must match trajectory atom order)
    trajectory_path : str, optional
        Path to trajectory file (DCD, XTC, TRR, etc.)
        If not provided, only the topology structure will be used
    assets_dir : str, optional
        Directory containing DNATCO classification assets

    Attributes
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis universe
    structure : Structure
        The libLLKA structure (topology)
    context : ClassificationContext
        The classification context
    n_atoms : int
        Number of atoms

    Examples
    --------
    >>> analyzer = TrajectoryAnalyzer("topology.cif", "trajectory.dcd")
    >>> results = analyzer.analyze_trajectory()
    >>> print(f"Analyzed {len(results)} frames")

    >>> # Analyze specific frames
    >>> for i, frame_results in enumerate(analyzer.iter_frames(start=0, stop=100)):
    ...     print(f"Frame {i}: {len(frame_results)} steps")
    """

    def __init__(self, topology_path, trajectory_path=None, assets_dir=None):
        """Initialize the trajectory analyzer."""
        # Load topology with MDAnalysis
        if trajectory_path:
            self.universe = mda.Universe(topology_path, trajectory_path)
        else:
            self.universe = mda.Universe(topology_path)

        # Load topology with libLLKA
        self.structure = load_structure(topology_path)

        # Validate atom counts match
        if self.universe.atoms.n_atoms != self.structure.n_atoms:
            raise ValueError(
                f"Atom count mismatch: MDAnalysis has {self.universe.atoms.n_atoms} atoms, "
                f"but libLLKA structure has {self.structure.n_atoms} atoms. "
                "Ensure the mmCIF topology matches the trajectory."
            )

        # Initialize classification context
        self.context = create_classification_context(assets_dir)

        self.n_atoms = self.structure.n_atoms

    def _update_structure_from_universe(self):
        """Update libLLKA structure with current MDAnalysis coordinates."""
        update_coordinates(self.structure, self.universe.atoms.positions)

    def classify_current_frame(self):
        """
        Classify the current trajectory frame.

        Returns
        -------
        list of ClassifiedStep
            Classification results for current frame

        Examples
        --------
        >>> analyzer.universe.trajectory[10]  # Go to frame 10
        >>> results = analyzer.classify_current_frame()
        """
        self._update_structure_from_universe()
        return classify_structure(self.structure, self.context)

    def iter_frames(self, start=None, stop=None, step=None):
        """
        Iterate over trajectory frames and yield classification results.

        Parameters
        ----------
        start : int, optional
            Starting frame index
        stop : int, optional
            Stopping frame index
        step : int, optional
            Step size

        Yields
        ------
        list of ClassifiedStep
            Classification results for each frame

        Examples
        --------
        >>> for frame_results in analyzer.iter_frames(start=0, stop=100, step=10):
        ...     print(f"Frame has {len(frame_results)} classified steps")
        """
        for ts in self.universe.trajectory[start:stop:step]:
            yield self.classify_current_frame()

    def analyze_trajectory(self, start=None, stop=None, step=None, verbose=True):
        """
        Analyze the entire trajectory (or a subset of frames).

        Parameters
        ----------
        start : int, optional
            Starting frame index
        stop : int, optional
            Stopping frame index
        step : int, optional
            Step size
        verbose : bool, default=True
            Print progress information

        Returns
        -------
        list of list of ClassifiedStep
            Results for each frame: results[frame_idx][step_idx]

        Examples
        --------
        >>> results = analyzer.analyze_trajectory(start=0, stop=100)
        >>> print(f"Analyzed {len(results)} frames")
        >>> # Access results for frame 0, step 0
        >>> print(results[0][0].assigned_ntc)
        """
        results = []
        frame_slice = self.universe.trajectory[start:stop:step]

        for i, ts in enumerate(frame_slice):
            if verbose and i % 100 == 0:
                print(f"Analyzing frame {i}/{len(frame_slice)}")

            frame_results = self.classify_current_frame()
            results.append(frame_results)

        return results

    def get_step_timeseries(self, step_index, start=None, stop=None, step=None):
        """
        Extract timeseries for a specific dinucleotide step across frames.

        Parameters
        ----------
        step_index : int
            Index of the dinucleotide step to track
        start : int, optional
            Starting frame index
        stop : int, optional
            Stopping frame index
        step : int, optional
            Step size

        Returns
        -------
        dict
            Dictionary with arrays for each property:
            - 'ntc': NtC classifications
            - 'cana': CANA classifications
            - 'rmsd': RMSD values
            - 'violations': Violation flags

        Examples
        --------
        >>> timeseries = analyzer.get_step_timeseries(step_index=0, start=0, stop=100)
        >>> print(timeseries['ntc'])  # NtC for each frame
        >>> print(timeseries['rmsd'])  # RMSD for each frame
        """
        ntcs = []
        canas = []
        rmsds = []
        violations = []

        for frame_results in self.iter_frames(start, stop, step):
            if step_index >= len(frame_results):
                raise IndexError(
                    f"Step index {step_index} out of range "
                    f"(frame has {len(frame_results)} steps)"
                )

            result = frame_results[step_index]
            ntcs.append(result.assigned_ntc)
            canas.append(result.assigned_cana)
            rmsds.append(result.rmsd_to_closest_ntc)
            violations.append(result.violations)

        return {
            'ntc': np.array(ntcs),
            'cana': np.array(canas),
            'rmsd': np.array(rmsds),
            'violations': np.array(violations),
        }

    def __repr__(self):
        return (
            f"TrajectoryAnalyzer(n_atoms={self.n_atoms}, "
            f"n_frames={len(self.universe.trajectory)})"
        )


def analyze_trajectory(topology_path, trajectory_path, start=None, stop=None,
                       step=None, assets_dir=None, verbose=True):
    """
    Convenience function to analyze a trajectory in one call.

    Parameters
    ----------
    topology_path : str
        Path to mmCIF topology file
    trajectory_path : str
        Path to trajectory file
    start : int, optional
        Starting frame index
    stop : int, optional
        Stopping frame index
    step : int, optional
        Step size
    assets_dir : str, optional
        Directory containing DNATCO assets
    verbose : bool, default=True
        Print progress

    Returns
    -------
    list of list of ClassifiedStep
        Results for each frame: results[frame_idx][step_idx]

    Examples
    --------
    >>> results = analyze_trajectory("topology.cif", "trajectory.dcd",
    ...                              start=0, stop=100)
    >>> print(f"Analyzed {len(results)} frames")
    """
    analyzer = TrajectoryAnalyzer(topology_path, trajectory_path, assets_dir)
    return analyzer.analyze_trajectory(start, stop, step, verbose)
