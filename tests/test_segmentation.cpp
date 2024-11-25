/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "effedup.hpp"
#include "testing_structures.h"

#include <llka_segmentation.h>

#include <array>

static
auto testSegmentation()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1DK1_B_26_28_ATOMS, REAL_1DK1_B_26_28_ATOMS_LEN);
    LLKA_StructureSegments segs = LLKA_structureSegments(&stru);

    EFF_expect(segs.nModels, 1UL, "Wrong number of models");
    auto &model = segs.models[0];
    EFF_expect(model.nAtoms, 113UL, "Wrong number of atoms in model");
    EFF_expect(model.atoms[0]->id, 541U, "Wrong atom id of the first atom in model");

    EFF_expect(model.nChains, 1UL, "Wrong number of chains in model");
    auto &chain = model.chains[0];
    EFF_expect(chain.nResidues, 3UL, "Wrong number of residues in chain");
    EFF_expect(chain.nAtoms, 113UL, "Wrong number of atoms in chain");
    EFF_expect(chain.atoms[0]->id, 541U, "Wrong atom id of the first atom in chain");

    auto &residueA = chain.residues[0];
    EFF_expect(residueA.nAtoms, 23UL, "Wrong number of atoms in residue");
    auto &residueAfirstAtom = residueA.atoms[0];
    EFF_expect(residueAfirstAtom->id, 541U, "Wrong atom id of first atom in residue");
    EFF_expect(residueAfirstAtom->label_seq_id, 26, "Wrong label_seq_id of first atom in residue");
    auto &residueAlastAtom = residueA.atoms[22];
    EFF_expect(residueAlastAtom->id, 563U, "Wrong atom id of last atom in residue");
    EFF_expect(residueAlastAtom->label_seq_id, 26, "Wrong label_seq_id of last atom in residue");

    auto &residueB = chain.residues[2];
    EFF_expect(residueB.nAtoms, 46UL, "Wrong number of atoms in residue");
    auto &residueBfirstAtom = residueB.atoms[0];
    EFF_expect(residueBfirstAtom->id, 608U, "Wrong atom id of first atom in residue");
    EFF_expect(residueBfirstAtom->label_seq_id, 28, "Wrong label_seq_id of first atom in residue");
    EFF_expect(residueBfirstAtom->label_alt_id, 'A', "Wrong label_alt_id of first atom in residue");
    auto &residueBlastAtom = residueB.atoms[45];
    EFF_expect(residueBlastAtom->id, 653U, "Wrong atom id of last atom in residue");
    EFF_expect(residueBlastAtom->label_seq_id, 28, "Wrong label_seq_id of last atom in residue");
    EFF_expect(residueBlastAtom->label_alt_id, 'B', "Wrong label_alt_id of last atom in residue");

    LLKA_destroyStructureSegments(&segs);
    LLKA_destroyStructure(&stru);
}

static
auto testSegmentationMultipleChains()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1CJG_A_B_1_2_M1, REAL_1CJG_A_B_1_2_M1_LEN);
    LLKA_StructureSegments segs = LLKA_structureSegments(&stru);

    EFF_expect(segs.nModels, 1UL, "Wrong number of models");

    auto &model = segs.models[0];
    EFF_expect(model.nAtoms, 82UL, "Wrong number of atoms in model");
    EFF_expect(model.atoms[0]->id, 1397U, "Wrong atom id of the first atom in model");

    EFF_expect(model.nChains, 2UL, "Wrong number of chains in model");

    // Chain A
    auto &chainA = model.chains[0];
    EFF_expect(chainA.nResidues, 2UL, "Wrong number of residues in chain A");
    EFF_expect(chainA.nAtoms, 41UL, "Wrong number of atoms in chain A");

    auto &firstResidueA = chainA.residues[0];
    EFF_expect(firstResidueA.nAtoms, 19UL, "Wrong number of atoms in first residue of chain A");
    EFF_expect(firstResidueA.atoms[0]->id, 1397U, "Wrong atom id on first residue of chain A");
    EFF_expect(firstResidueA.atoms[0]->label_asym_id, "A", "Wrong label_asym_id on first residue of chain A");

    auto &secondResidueA = chainA.residues[1];
    EFF_expect(secondResidueA.nAtoms, 22UL, "Wrong number of atoms in second residue of chain A");
    EFF_expect(secondResidueA.atoms[0]->id, 1416U, "Wrong atom id on second residue of chain A");
    EFF_expect(secondResidueA.atoms[0]->label_asym_id, "A", "Wrong atom id on second residue of chain A");

    // Chain B
    auto &chainB = model.chains[1];
    EFF_expect(chainB.nResidues, 2UL, "Wrong number of residues in chain B");
    EFF_expect(chainB.nAtoms, 41UL, "Wrong number of atoms in chain B");

    auto &firstResidueB = chainB.residues[0];
    EFF_expect(firstResidueB.nAtoms, 19UL, "Wrong number of atoms in first residue of chain B");
    EFF_expect(firstResidueB.atoms[0]->id, 2363U, "Wrong atom id on first residue of chain B");
    EFF_expect(firstResidueB.atoms[0]->label_asym_id, "B", "Wrong label_asym_id on first residue of chain B");

    auto &secondResidueB = chainB.residues[1];
    EFF_expect(secondResidueB.nAtoms, 22UL, "Wrong number of atoms in second residue of chain B");
    EFF_expect(secondResidueB.atoms[0]->id, 2382U, "Wrong atom id on second residue of chain B");
    EFF_expect(secondResidueB.atoms[0]->label_asym_id, "B", "Wrong label_asym_id on second residue of chain B");

    LLKA_destroyStructureSegments(&segs);
    LLKA_destroyStructure(&stru);
}

static
auto testSegmentationMultipleModels()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1CJG_A_1_6_M6_M1_M2_M3, REAL_1CJG_A_1_6_M6_M1_M2_M3_LEN);
    LLKA_StructureSegments segs = LLKA_structureSegments(&stru);

    EFF_expect(segs.nModels, 4UL, "Wrong number of models");

    std::array<int32_t, 4> MODEL_NUMBERS = { 1, 2, 3, 6 };
    std::array<int32_t, 2> RESIDUE_SEQ_IDS = { 1, 6 };

    for (size_t modelIdx = 0; modelIdx < segs.nModels; modelIdx++) {
        auto &model = segs.models[modelIdx];
        EFF_expect(model.nAtoms, 38UL, "Wrong number of atoms in model");
        EFF_expect(model.atoms[0]->label_atom_id, "N", "Wrong atom id of the first atom in model");
        EFF_expect(model.atoms[0]->pdbx_PDB_model_num, MODEL_NUMBERS[modelIdx], "Wrong pdbx_PDB_model_num on first atom in model");

        EFF_expect(model.nChains, 1UL, "Wrong number of chains in model");
        auto &chain = model.chains[0];
        EFF_expect(chain.nResidues, 2UL, "Wrong number of residues in chain");
        EFF_expect(chain.nAtoms, 38UL, "Wrong number of atoms in chain");
        EFF_expect(chain.atoms[0]->label_asym_id, "A", "Wrong label_asym_id on first atom in chain");

        auto &firstResidue = chain.residues[0];
        EFF_expect(firstResidue.nAtoms, 19UL, "Wrong number of atoms in first residue");
        EFF_expect(firstResidue.atoms[0]->label_atom_id, "N", "Wrong label_atom_id on first atom in residue");
        EFF_expect(firstResidue.atoms[18]->label_atom_id, "HE3", "Wrong label_atom_id on last atom in residue");
        EFF_expect(firstResidue.atoms[0]->label_seq_id, RESIDUE_SEQ_IDS[0], "Wrong label_seq_id on first atom in residue");
        EFF_expect(firstResidue.atoms[18]->label_seq_id, RESIDUE_SEQ_IDS[0], "Wrong label_seq_id on last atom in residue");

        auto &secondResidue = chain.residues[1];
        EFF_expect(secondResidue.nAtoms, 19UL, "Wrong number of atoms in second residue");
        EFF_expect(secondResidue.atoms[0]->label_atom_id, "N", "Wrong label_atom_id on first atom in residue");
        EFF_expect(secondResidue.atoms[18]->label_atom_id, "HD23", "Wrong label_atom_id on last atom in residue");
        EFF_expect(secondResidue.atoms[0]->label_seq_id, RESIDUE_SEQ_IDS[1], "Wrong label_seq_id on first atom in residue");
        EFF_expect(secondResidue.atoms[18]->label_seq_id, RESIDUE_SEQ_IDS[1], "Wrong label_seq_id on last atom in residue");
    }

    LLKA_destroyStructureSegments(&segs);
    LLKA_destroyStructure(&stru);
}

auto main(int, char **) -> int
{
    testSegmentation();
    testSegmentationMultipleChains();
    testSegmentationMultipleModels();
}
