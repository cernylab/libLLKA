/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "effedup.hpp"
#include "testing_structures.h"

#include <llka_measurements.h>
#include <llka_structure.h>

static
auto testAngles(LLKA_Structure &stru)
{
    auto atomA = LLKA_findAtom(&stru, "O5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    auto atomB = LLKA_findAtom(&stru, "C5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    auto atomC = LLKA_findAtom(&stru, "C4'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    assert(atomA != nullptr); assert(atomB != nullptr); assert(atomC != nullptr);
    auto ang = LLKA_measureAngle(atomA, atomB, atomC);

    EFF_cmpFlt(ang, 1.885722272657443, "unexpected angle between O5' - C5' - C4' atoms");

    atomA = LLKA_findAtom(&stru, "C5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomB = LLKA_findAtom(&stru, "C4'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomC = LLKA_findAtom(&stru, "C3'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    assert(atomA != nullptr); assert(atomB != nullptr); assert(atomC != nullptr);
    ang = LLKA_measureAngle(atomA, atomB, atomC);

    EFF_cmpFlt(ang, 2.0342231990866746, "unexpected angle between C5' - C4' - C3' atoms");

    atomA = LLKA_findAtom(&stru, "O2", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomB = LLKA_findAtom(&stru, "C2", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomC = LLKA_findAtom(&stru, "C5", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    assert(atomA != nullptr); assert(atomB != nullptr); assert(atomC != nullptr);
    ang = LLKA_measureAngle(atomA, atomB, atomC);

    EFF_cmpFlt(ang, 3.1036907363579922, "unexpected angle between O2 - C2 - C5 atoms");
}

static
auto testDihedrals(LLKA_Structure &stru)
{
    auto atomA = LLKA_findAtom(&stru, "O5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    auto atomB = LLKA_findAtom(&stru, "C5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    auto atomC = LLKA_findAtom(&stru, "C4'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    auto atomD = LLKA_findAtom(&stru, "O4'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    assert(atomA != nullptr); assert(atomB != nullptr); assert(atomC != nullptr); assert(atomD != nullptr);
    auto dih = LLKA_measureDihedral(atomA, atomB, atomC, atomD);

    EFF_cmpFlt(dih, 0.93724348718692196, "unexpected dihedral angle between O5' - C5' - C4' - O4' atoms");

    atomA = LLKA_findAtom(&stru, "O5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomB = LLKA_findAtom(&stru, "C5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomC = LLKA_findAtom(&stru, "C4'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomD = LLKA_findAtom(&stru, "C3'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    assert(atomA != nullptr); assert(atomB != nullptr); assert(atomC != nullptr); assert(atomD != nullptr);
    dih = LLKA_measureDihedral(atomA, atomB, atomC, atomD);

    EFF_cmpFlt(dih, 3.0409263486359723, "unexpected dihedral angle between O5' - C5' - C4' - C3' atoms");

    atomA = LLKA_findAtom(&stru, "C4", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomB = LLKA_findAtom(&stru, "N3", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomC = LLKA_findAtom(&stru, "C2", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomD = LLKA_findAtom(&stru, "N1", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    assert(atomA != nullptr); assert(atomB != nullptr); assert(atomC != nullptr); assert(atomD != nullptr);
    dih = LLKA_measureDihedral(atomA, atomB, atomC, atomD);

    EFF_cmpFlt(dih, 0.0075379644566811376, "unexpected dihedral angle between C4 - N3 - C2 - N1 atoms");
}

static
auto testDistances(LLKA_Structure &stru)
{
    auto atomA = LLKA_findAtom(&stru, "O5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    auto atomB = LLKA_findAtom(&stru, "C5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    assert(atomA != nullptr); assert(atomB != nullptr);
    auto dist = LLKA_measureDistance(atomA, atomB);

    EFF_cmpFlt(dist, 1.4378821231241452, "unexpected distance between O5' -> C5' atoms");

    atomA = LLKA_findAtom(&stru, "C1'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomB = LLKA_findAtom(&stru, "C2", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    assert(atomA != nullptr); assert(atomB != nullptr);
    dist = LLKA_measureDistance(atomA, atomB);

    EFF_cmpFlt(dist, 2.5047702489450008, "unexpected distance between C1' -> C2 atoms");

    atomA = LLKA_findAtom(&stru, "O5'", nullptr, nullptr, 1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    atomB = LLKA_findAtom(&stru, "O5'", nullptr, nullptr, 2, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    assert(atomA != nullptr); assert(atomB != nullptr);
    dist = LLKA_measureDistance(atomA, atomB);

    EFF_cmpFlt(dist, 7.2086283022500197, "unexpected distance between O5'_1 -> O5'_2 atoms");
}

auto main(int , char **) -> int
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1BNA_A_1_2_ATOMS, REAL_1BNA_A_1_2_ATOMS_LEN);

    testDistances(stru);
    testAngles(stru);
    testDihedrals(stru);

    LLKA_destroyStructure(&stru);
}
