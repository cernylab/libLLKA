/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKAINTERNAL_SEGMENTATION_H
#define _LLKAINTERNAL_SEGMENTATION_H

#include <map>
#include <string>

namespace LLKAInternal {

/*
 * Counts the number of atoms in a residue.
 * Mapping key is label_seq_id
 */
using ResidueCounts = std::map<int32_t, size_t>;

/*
 * Counts the number of residues in a chain.
 * Mapping key is label_asym_id
 */
struct _ChainCounts {
    ResidueCounts residues;
    size_t nAtoms;
};
using ChainCounts = std::map<std::string, _ChainCounts>;

/*
 * Counts the number of chains in a model.
 * Mapping key is pdbx_PDB_model_num
 */
struct _ModelCounts {
    ChainCounts chains;
    size_t nAtoms;
};
using ModelCounts = std::map<int32_t, _ModelCounts>;

template <typename AtomType>
static
auto countAtomsForSegmentation(const AtomType *atoms, const size_t nAtoms)
{
    ModelCounts models;

    for (size_t atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
        const auto &atom = atoms[atomIdx];

        if (!models.contains(atom.pdbx_PDB_model_num))
            models.emplace(atom.pdbx_PDB_model_num, _ModelCounts{});

        auto &model = models.at(atom.pdbx_PDB_model_num);
        model.nAtoms++;

        std::string labelAsymId{atom.label_asym_id};
        if (!model.chains.contains(labelAsymId))
            model.chains.emplace(labelAsymId, _ChainCounts{});

        auto &chain = model.chains.at(labelAsymId);
        chain.nAtoms++;

        if (!chain.residues.contains(atom.label_seq_id))
            chain.residues.emplace(atom.label_seq_id, 0);

        auto &residue = chain.residues.at(atom.label_seq_id);
        residue++;
    }

    return models;
}

} /* namespace LLKAInternal */

#endif /* _LLKAINTERNAL_SEGMENTATION_H */
