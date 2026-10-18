"""Segmentation of a structure into models, chains and residues (``llka_segmentation.h``)."""

from __future__ import annotations

import weakref
from collections.abc import Iterator, Sequence
from ctypes import byref
from typing import overload

from ._capi import LLKA_StructureSegments, lib
from .structure import Structure, StructureView, require_structure

__all__ = ['Chain', 'Model', 'Residue', 'StructureSegments', 'structure_segments']


class Residue:
    """A residue of a chain. ``atoms`` is a view of the atoms of the segmented structure."""

    __slots__ = ('atoms',)

    def __init__(self, atoms: StructureView) -> None:
        self.atoms = atoms

    def __repr__(self) -> str:
        return f'<Residue with {len(self.atoms)} atoms>'


class Chain:
    """A chain of a model. ``atoms`` is a view of all atoms of the chain."""

    __slots__ = ('atoms', 'residues')

    def __init__(self, residues: list[Residue], atoms: StructureView) -> None:
        self.residues = residues
        self.atoms = atoms

    def __repr__(self) -> str:
        return f'<Chain with {len(self.residues)} residues and {len(self.atoms)} atoms>'


class Model:
    """A model of a structure. ``atoms`` is a view of all atoms of the model."""

    __slots__ = ('atoms', 'chains')

    def __init__(self, chains: list[Chain], atoms: StructureView) -> None:
        self.chains = chains
        self.atoms = atoms

    def __repr__(self) -> str:
        return f'<Model with {len(self.chains)} chains and {len(self.atoms)} atoms>'


class _SegmentsOwner:
    __slots__ = ('_c', '_finalizer', '__weakref__')

    def __init__(self, c: LLKA_StructureSegments) -> None:
        self._c = c
        self._finalizer = weakref.finalize(self, lib.LLKA_destroyStructureSegments, byref(c))


class StructureSegments(Sequence[Model]):
    """Models of a structure, sorted by model number, each split into chains and residues.

    Chains are sorted by ``label_asym_id`` and residues by ``label_seq_id``.
    All atom lists are views of the segmented structure and become invalid
    when atoms are added to or removed from the structure.
    """

    __slots__ = ('_models',)

    def __init__(self, models: list[Model]) -> None:
        self._models = models

    @property
    def models(self) -> list[Model]:
        return list(self._models)

    def __len__(self) -> int:
        return len(self._models)

    @overload
    def __getitem__(self, idx: int) -> Model: ...
    @overload
    def __getitem__(self, idx: slice) -> list[Model]: ...

    def __getitem__(self, idx):
        return self._models[idx]

    def __iter__(self) -> Iterator[Model]:
        return iter(self._models)

    def __repr__(self) -> str:
        return f'<StructureSegments with {len(self)} models>'


def structure_segments(stru: Structure) -> StructureSegments:
    """Splits a structure into models, chains and residues."""
    require_structure(stru, 'structure_segments()')

    c = lib.LLKA_structureSegments(byref(stru._c))
    if c.nModels == 0:
        return StructureSegments([])

    owner = _SegmentsOwner(c)

    def view(atoms, n_atoms: int) -> StructureView:
        return StructureView._borrow(atoms, n_atoms, stru, owner)

    models = []
    for model_idx in range(c.nModels):
        c_model = c.models[model_idx]
        chains = []
        for chain_idx in range(c_model.nChains):
            c_chain = c_model.chains[chain_idx]
            residues = [
                Residue(view(c_chain.residues[res_idx].atoms, c_chain.residues[res_idx].nAtoms))
                for res_idx in range(c_chain.nResidues)
            ]
            chains.append(Chain(residues, view(c_chain.atoms, c_chain.nAtoms)))
        models.append(Model(chains, view(c_model.atoms, c_model.nAtoms)))

    return StructureSegments(models)
