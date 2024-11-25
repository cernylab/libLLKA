/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_cpp.h>

#include "effedup.hpp"
#include "testing_structures.h"

#include <array>

#ifdef LLKA_PLATFORM_EMSCRIPTEN
    #define _a(v) (&v)
#else
    #define _a(v) (v)
#endif // LLKA_PLATFORM_EMSCRIPTEN

static
auto testSegmentationMultipleChains()
{
    auto stru = LLKA::makeStructure(REAL_1CJG_A_B_1_2_M1, REAL_1CJG_A_B_1_2_M1_LEN);
    auto segs = LLKA::StructureSegments(stru);

    EFF_expect(segs.models.size(), 1UL, "Wrong number of models");

    auto &model = segs.models.at(1);
#ifndef LLKA_PLATFORM_EMSCRIPTEN
    EFF_expect(model.atoms.size(), 82UL, "Wrong number of atoms in model");
    EFF_expect(_a(model.atoms.front())->id, 1397U, "Wrong atom id of the first atom in model");
#endif // LLKA_PLATFORM_EMSCRIPTEN

    EFF_expect(model.chains.size(), 2UL, "Wrong number of chains in model");

    // Chain A
    auto &chainA = model.chains.at("A");
    EFF_expect(chainA.residues.size(), 2UL, "Wrong number of residues in chain A");
#ifndef LLKA_PLATFORM_EMSCRIPTEN
    EFF_expect(_a(chainA.atoms.size()), 41UL, "Wrong number of atoms in chain A");
#endif // LLKA_PLATFORM_EMSCRIPTEN

    auto &firstResidueA = chainA.residues.cbegin()->second;
    EFF_expect(firstResidueA.atoms.size(), 19UL, "Wrong number of atoms in first residue of chain A");
    EFF_expect(_a(firstResidueA.atoms.front())->id, 1397U, "Wrong atom id on first residue of chain A");
    EFF_expect(_a(firstResidueA.atoms.front())->label_asym_id, "A", "Wrong label_asym_id on first residue of chain A");

    auto &secondResidueA = (++chainA.residues.cbegin())->second;
    EFF_expect(secondResidueA.atoms.size(), 22UL, "Wrong number of atoms in second residue of chain A");
    EFF_expect(_a(secondResidueA.atoms.front())->id, 1416U, "Wrong atom id on second residue of chain A");
    EFF_expect(_a(secondResidueA.atoms.front())->label_asym_id, "A", "Wrong atom id on second residue of chain A");

    // Chain B
    auto &chainB = model.chains.at("B");
    EFF_expect(chainB.residues.size(), 2UL, "Wrong number of residues in chain B");
#ifndef LLKA_PLATFORM_EMSCRIPTEN
    EFF_expect(chainB.atoms.size(), 41UL, "Wrong number of atoms in chain B");
#endif // LLKA_PLATFORM_EMSCRIPTEN

    auto &firstResidueB = chainB.residues.cbegin()->second;
    EFF_expect(firstResidueB.atoms.size(), 19UL, "Wrong number of atoms in first residue of chain B");
    EFF_expect(_a(firstResidueB.atoms.front())->id, 2363U, "Wrong atom id on first residue of chain B");
    EFF_expect(_a(firstResidueB.atoms.front())->label_asym_id, "B", "Wrong label_asym_id on first residue of chain B");

    auto &secondResidueB = (++chainB.residues.cbegin())->second;
    EFF_expect(secondResidueB.atoms.size(), 22UL, "Wrong number of atoms in second residue of chain B");
    EFF_expect(_a(secondResidueB.atoms.front())->id, 2382U, "Wrong atom id on second residue of chain B");
    EFF_expect(_a(secondResidueB.atoms.front())->label_asym_id, "B", "Wrong label_asym_id on second residue of chain B");
}

static
auto testSegmentationMultipleModels()
{
    std::array<int32_t, 2> RESIDUE_SEQ_IDS = { 1, 6 };

    auto stru = LLKA::makeStructure(REAL_1CJG_A_1_6_M6_M1_M2_M3, REAL_1CJG_A_1_6_M6_M1_M2_M3_LEN);
    auto segs = LLKA::StructureSegments(stru);

    EFF_expect(segs.models.size(), 4UL, "Wrong number of models");

    for (const auto &modelIt : segs.models) {
        const auto &model = modelIt.second;

#ifndef LLKA_PLATFORM_EMSCRIPTEN
        EFF_expect(model.atoms.size(), 38UL, "Wrong number of atoms in model");
        EFF_expect(_a(model.atoms.front())->label_atom_id, "N", "Wrong atom id of the first atom in model");
        EFF_expect(_a(model.atoms.front())->pdbx_PDB_model_num, modelIt.first, "Wrong pdbx_PDB_model_num on first atom in model");
#endif // LLKA_PLATFORM_EMSCRIPTEN

        EFF_expect(model.chains.size(), 1UL, "Wrong number of chains in model");
        auto &chain = model.chains.at("A");
        EFF_expect(chain.residues.size(), 2UL, "Wrong number of residues in chain");
#ifndef LLKA_PLATFORM_EMSCRIPTEN
        EFF_expect(chain.atoms.size(), 38UL, "Wrong number of atoms in chain");
        EFF_expect(_a(chain.atoms.front())->label_asym_id, "A", "Wrong label_asym_id on first atom in chain");
#endif // LLKA_PLATFORM_EMSCRIPTEN

        auto &firstResidue = chain.residues.cbegin()->second;
        EFF_expect(firstResidue.atoms.size(), 19UL, "Wrong number of atoms in first residue");
        EFF_expect(_a(firstResidue.atoms.front())->label_atom_id, "N", "Wrong label_atom_id on first atom in residue");
        EFF_expect(_a(firstResidue.atoms.back())->label_atom_id, "HE3", "Wrong label_atom_id on last atom in residue");
        EFF_expect(_a(firstResidue.atoms.front())->label_seq_id, RESIDUE_SEQ_IDS[0], "Wrong label_seq_id on first atom in residue");
        EFF_expect(_a(firstResidue.atoms.back())->label_seq_id, RESIDUE_SEQ_IDS[0], "Wrong label_seq_id on last atom in residue");

        auto &secondResidue = (++chain.residues.cbegin())->second;
        EFF_expect(secondResidue.atoms.size(), 19UL, "Wrong number of atoms in second residue");
        EFF_expect(_a(secondResidue.atoms.front())->label_atom_id, "N", "Wrong label_atom_id on first atom in residue");
        EFF_expect(_a(secondResidue.atoms.back())->label_atom_id, "HD23", "Wrong label_atom_id on last atom in residue");
        EFF_expect(_a(secondResidue.atoms.front())->label_seq_id, RESIDUE_SEQ_IDS[1], "Wrong label_seq_id on first atom in residue");
        EFF_expect(_a(secondResidue.atoms.back())->label_seq_id, RESIDUE_SEQ_IDS[1], "Wrong label_seq_id on last atom in residue");
    }
}

static
auto testSugarPucker()
{
    auto stru = LLKA::makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);

    // Test LLKA_Structure
    auto nucl = LLKA::extractNucleotide(stru, 1, "A", 3);
    auto res = LLKA::sugarPucker(nucl);
    EFF_expect(res.isSuccess(), true, "LLKA::sugarPucker() failed with return code " + LLKA::errorToString(res.failure()));
    EFF_expect(res.success(), LLKA_O4_ENDO, "Wrong sugar pucker");

    // Test LLKA_StructureView
    auto nuclView = LLKA::extractNucleotideView(stru, 1, "A", 4);
    auto res2 = LLKA::sugarPucker(nuclView);
    EFF_expect(res2.isSuccess(), true, "LLKA::sugarPucker() failed with return code " + LLKA::errorToString(res2.failure()));
    EFF_expect(res2.success(), LLKA_C2_ENDO, "Wrong sugar pucker");
}

auto main() -> int
{
    testSegmentationMultipleChains();
    testSegmentationMultipleModels();

    testSugarPucker();
}
