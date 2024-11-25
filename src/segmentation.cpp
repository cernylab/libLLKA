/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_segmentation.h>

#include <segmentation.h>

#include <cassert>
#include <map>
#include <string>

void LLKA_CC LLKA_destroyStructureSegments(const LLKA_StructureSegments *segs)
{
    if (segs->nModels == 0)
        return;

    for (size_t modelIdx = 0; modelIdx < segs->nModels; modelIdx++) {
        auto &model = segs->models[modelIdx];

        for (size_t chainIdx = 0; chainIdx < model.nChains; chainIdx++) {
            auto &chain = model.chains[chainIdx];

            for (size_t residueIdx = 0; residueIdx < chain.nResidues; residueIdx++) {
                auto &residue = chain.residues[residueIdx];

                delete [] residue.atoms;
            }

            delete [] chain.residues;
            delete [] chain.atoms;
        }

        delete [] model.chains;
        delete [] model.atoms;
    }

    delete [] segs->models;
}

LLKA_StructureSegments LLKA_CC LLKA_structureSegments(LLKA_Structure *stru)
{
    if (stru->nAtoms == 0)
        return {};

    // Count number of chains, residues and atoms in all models of the structure
    auto counts = LLKAInternal::countAtomsForSegmentation(stru->atoms, stru->nAtoms);

    LLKA_StructureSegments segments{
        .models = new LLKA_Model[counts.size()],
        .nModels = counts.size()
    };

    // Initialize the mapping arrays
    for (auto modelIt = counts.cbegin(); modelIt != counts.cend(); modelIt++) {
        auto &model = modelIt->second;
        size_t llkaModelIdx = std::distance(counts.cbegin(), modelIt);
        auto &llkaModel = segments.models[llkaModelIdx];

        llkaModel.chains = new LLKA_Chain[model.chains.size()];
        llkaModel.nChains = model.chains.size();

        llkaModel.atoms = new LLKA_Atom *[model.nAtoms];
        llkaModel.nAtoms = 0; // Set during the mapping

        for (auto chainIt = model.chains.cbegin(); chainIt != model.chains.cend(); chainIt++) {
            auto &chain = chainIt->second;
            size_t llkaChainIdx = std::distance(model.chains.cbegin(), chainIt);
            auto &llkaChain = llkaModel.chains[llkaChainIdx];

            llkaChain.residues = new LLKA_Residue[chain.residues.size()];
            llkaChain.nResidues = chain.residues.size();

            llkaChain.atoms = new LLKA_Atom *[chain.nAtoms];
            llkaChain.nAtoms = 0; // Set during the mapping

            for (auto residueIt = chain.residues.cbegin(); residueIt != chain.residues.cend(); residueIt++) {
                auto &nAtoms = residueIt->second;
                size_t llkaResidueIdx = std::distance(chain.residues.cbegin(), residueIt);
                auto &llkaResidue = llkaChain.residues[llkaResidueIdx];

                llkaResidue.atoms = new LLKA_Atom *[nAtoms];
                llkaResidue.nAtoms = 0; // Set during the mapping
            }
        }
    }

    // Do the actual mapping by walking the structure atom by atom.
    for (size_t atomIdx = 0; atomIdx < stru->nAtoms; atomIdx++) {
        auto atom = &stru->atoms[atomIdx];

        auto modelIt = counts.find(atom->pdbx_PDB_model_num);

        assert(modelIt != counts.cend());
        auto chainIt = modelIt->second.chains.find(atom->label_asym_id);
        assert(chainIt != modelIt->second.chains.cend());
        auto residueIt = chainIt->second.residues.find(atom->label_seq_id);
        assert(residueIt != chainIt->second.residues.cend());

        size_t llkaModelIdx = std::distance(counts.begin(), modelIt);
        size_t llkaChainIdx = std::distance(modelIt->second.chains.begin(), chainIt);
        size_t llkaResidueIdx = std::distance(chainIt->second.residues.begin(), residueIt);

        auto &model = segments.models[llkaModelIdx];
        auto &chain = model.chains[llkaChainIdx];
        auto &residue = chain.residues[llkaResidueIdx];

        model.atoms[model.nAtoms++] = atom;

        chain.atoms[chain.nAtoms++] = atom;

        residue.atoms[residue.nAtoms++] = atom;
    }

    return segments;
}
