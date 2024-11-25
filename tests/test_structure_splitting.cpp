/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_structure.h>

#include "effedup.hpp"
#include "testing_structures.h"

static
auto testSplitAltIds()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1DK1_B_26_28_ATOMS, REAL_1DK1_B_26_28_ATOMS_LEN);

    LLKA_AlternatePositions alts;
    LLKA_Structures splittedStrus = LLKA_splitByAltIds(&stru, &alts);

    EFF_expect(splittedStrus.nStrus, 2UL, "wrong number of alt-id-splitted structures")
    EFF_expect(alts.nPositions, 2UL, "wrong number of alterate position identifiers")

    EFF_expect(splittedStrus.strus[0].nAtoms, 68UL, "wrong number of atoms in splitted structure")
    EFF_expect(splittedStrus.strus[1].nAtoms, 68UL, "wrong number of atoms in splitted structure")

    EFF_expect(alts.positions[0], 'A', "wrong alternate position id")
    EFF_expect(alts.positions[1], 'B', "wrong alternate position id")

    EFF_expect(splittedStrus.strus[0].atoms[0].label_alt_id, LLKA_NO_ALTID, "wrong alternate position id")
    EFF_expect(splittedStrus.strus[0].atoms[67].label_alt_id, 'A', "wrong alternate position id")
    EFF_expect(splittedStrus.strus[1].atoms[0].label_alt_id, LLKA_NO_ALTID, "wrong alternate position id")
    EFF_expect(splittedStrus.strus[1].atoms[67].label_alt_id, 'B', "wrong alternate position id")

    LLKA_destroyStructures(&splittedStrus);
    LLKA_destroyAlternatePositions(&alts);
    LLKA_destroyStructure(&stru);
}

static
auto testSplitDinucleotides()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1DK1_B_26_28_ATOMS, REAL_1DK1_B_26_28_ATOMS_LEN);

    LLKA_Structures steps;
    auto tRet = LLKA_splitStructureToDinucleotideSteps(&stru, &steps);
    EFF_expect(tRet, LLKA_OK, "unexpected return value from LLKA_splitStructureToDinucleotideSteps()")
    EFF_expect(steps.nStrus, 3UL, "wrong number of steps")

    EFF_expect(steps.strus[0].atoms[0].label_seq_id, 26, "wrong sequence id")
    EFF_expect(steps.strus[0].atoms[steps.strus[0].nAtoms - 1].label_seq_id, 27, "wrong sequence id")

    EFF_expect(steps.strus[2].atoms[0].label_seq_id, 27, "wrong sequence id")
    EFF_expect(steps.strus[2].atoms[0].label_alt_id, 'B', "wrong alternate position id")
    EFF_expect(steps.strus[2].atoms[steps.strus[2].nAtoms - 1].label_seq_id, 28, "wrong sequence id")
    EFF_expect(steps.strus[2].atoms[steps.strus[2].nAtoms - 1].label_alt_id, 'B', "wrong alternate position id")

    LLKA_destroyStructures(&steps);
    LLKA_destroyStructure(&stru);
}

auto main(int, char **) -> int
{
    testSplitAltIds();
    testSplitDinucleotides();

    return EXIT_SUCCESS;
}
