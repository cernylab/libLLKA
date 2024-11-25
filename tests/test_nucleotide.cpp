/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_nucleotide.h>

#include "effedup.hpp"
#include "testing_structures.h"

#include <array>

static
auto testExtractRibose(const LLKA_Structure *stru)
{
    auto nucl = LLKA_extractNucleotide(stru, 1, "A", 3);

    LLKA_Structure riboseStru;
    auto tRet = LLKA_extractRibose(&nucl, &riboseStru);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractRibose returned unexpected value");
    EFF_expect(riboseStru.nAtoms, 5UL, "wrong ribose core atoms");

    constexpr std::array<const char *, 5> names{ "C4'", "O4'", "C1'", "C2'", "C3'" };
    for (size_t idx = 0; idx < 5; idx++)
        EFF_expect(riboseStru.atoms[idx].label_atom_id, names[idx], "wrong atom name");

    LLKA_destroyStructure(&riboseStru);
    LLKA_destroyStructure(&nucl);
}

static
auto testExtractNucleotide(const LLKA_Structure *stru)
{
    auto firstNucl = LLKA_extractNucleotide(stru, 1, "A", 3);
    EFF_expect(firstNucl.nAtoms, 19UL, "wrong number of atoms in nucleotide");
    EFF_expect(firstNucl.atoms[0].label_seq_id, firstNucl.atoms[firstNucl.nAtoms-1].label_seq_id, "sequence id of first and last atom in nucleotide do not match");

    auto secondNucl = LLKA_extractNucleotide(stru, 1, "A", 4);
    EFF_expect(secondNucl.nAtoms, 22UL, "wrong number of atoms in nucleotide");
    EFF_expect(secondNucl.atoms[0].label_seq_id, secondNucl.atoms[secondNucl.nAtoms-1].label_seq_id, "sequence id of first and last atom in nucleotide do not match");

    EFF_expect(firstNucl.atoms[0].label_seq_id + 1 == secondNucl.atoms[0].label_seq_id, true , "First atoms of two different nucleotides have the same sequence id");

    LLKA_destroyStructure(&secondNucl);
    LLKA_destroyStructure(&firstNucl);
}

static
auto testSugarPucker(const LLKA_Structure *stru)
{
    LLKA_SugarPucker pucker;

    // Test LLKA_Structure
    auto nucl = LLKA_extractNucleotide(stru, 1, "A", 3);
    auto tRet = LLKA_sugarPucker(&nucl, &pucker);
    EFF_expect(tRet, LLKA_OK, "Unexpected return code from LLKA_sugarPucker()");
    EFF_expect(pucker, LLKA_O4_ENDO, "Wrong sugar pucker");

    LLKA_destroyStructure(&nucl);

    // Test LLKA_StructureView
    auto nuclView = LLKA_extractNucleotideView(stru, 1, "A", 4);
    tRet = LLKA_sugarPuckerView(&nuclView, &pucker);
    EFF_expect(tRet, LLKA_OK, "Unexpected return code from LLKA_sugarPuckerView()");
    EFF_expect(pucker, LLKA_C2_ENDO, "Wrong sugar pucker");

    LLKA_destroyStructureView(&nuclView);
}

auto main(int, char *[]) -> int
{
    auto stru = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);

    testExtractNucleotide(&stru);
    testExtractRibose(&stru);
    testSugarPucker(&stru);

    LLKA_destroyStructure(&stru);

    return 0;
}
