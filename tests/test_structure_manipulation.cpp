// vim: set sw=4 ts=4 sts=4 expandtab :

#include "effedup.hpp"
#include "testing_structures.h"

#include <llka_structure.h>

static
auto testAppendAtom()
{
    const LLKA_Atom atoms[] = {
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 1, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 2, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 3, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 4, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID }
    };

    const LLKA_Atom toAppend = { "C", "C4'", "1", "DC", "A", "C4'", "DC", "A", { 0.551, -3.732, 0.531 }, 5, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID };

    LLKA_Structure stru = LLKA_makeStructure(atoms, 4);
    LLKA_appendAtom(&toAppend, &stru);

    EFF_expect(stru.nAtoms, 5UL, "append atom, count");
    EFF_expect(LLKA_compareAtoms(&stru.atoms[4], &toAppend, LLKA_FALSE), LLKA_TRUE, "append atom, compare appended");

    LLKA_Point c{ 10.0, 11.0, 12.0 };
    LLKA_appendAtomFromParams(
        6,
        "C", "C4'", "1", "DA", "A",
        nullptr, nullptr, nullptr,
        6, LLKA_NO_ALTID, 6,
        LLKA_NO_INSCODE,
        1,
        &c,
        &stru
    );

    EFF_expect(stru.nAtoms, 6UL, "append atom by params, count");
    EFF_expect(std::string{stru.atoms[5].type_symbol}, std::string{"C"}, "append atom by params, type_symbol");
    EFF_expect(std::string{stru.atoms[5].label_comp_id}, std::string{"DA"}, "append atom by params, label_comp_id");
    EFF_expect(stru.atoms[5].coords.x, 10.0, "append atom by params, X coordinate");

    LLKA_destroyStructure(&stru);
}

static
auto testAppendStructure()
{
    const LLKA_Atom atomsFirst[] = {
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 1, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 2, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 3, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 4, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID }
    };

    const LLKA_Atom atomsSecond[] = {
        { "N", "N9", "1", "DA", "A", "N9", "DA", "A", { 1.059, -3.973, 1.922 }, 5, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "N", "N9", "1", "DA", "A", "N9", "DA", "A", { 1.059, -3.973, 1.922 }, 6, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "N", "N9", "1", "DA", "A", "N9", "DA", "A", { 1.059, -3.973, 1.922 }, 7, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "N", "N9", "1", "DA", "A", "N9", "DA", "A", { 1.059, -3.973, 1.922 }, 8, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "N", "N9", "1", "DA", "A", "N9", "DA", "A", { 1.059, -3.973, 1.922 }, 9, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID }
    };

    LLKA_Structure firstStru = LLKA_makeStructure(atomsFirst, 4);
    LLKA_Structure secondStru = LLKA_makeStructure(atomsSecond, 5);

    LLKA_appendStructure(&secondStru, &firstStru);

    EFF_expect(firstStru.nAtoms, 9UL, "wrong number of atoms");
    EFF_expect(firstStru.atoms[0].auth_atom_id, "C5'", "wrong atom name");
    EFF_expect(firstStru.atoms[8].auth_atom_id, "N9", "wrong atom name");
    EFF_expect(firstStru.atoms[8].id, 9U, "wrong atom id");

    LLKA_destroyStructure(&secondStru);
    LLKA_destroyStructure(&firstStru);
}

static
auto testComparison()
{
    LLKA_Point coordsA = { 10, 11, 12 };
    auto atomA = LLKA_makeAtom(1, "C", "C5'", "1", "DA", "A", nullptr, nullptr, nullptr, 1, LLKA_NO_ALTID, 1, LLKA_NO_INSCODE, 1, &coordsA);

    LLKA_Point coordsAA = { 10, 11, 12 };
    auto atomAA = LLKA_makeAtom(2, "C", "C5'", "1", "DA", "A", nullptr, nullptr, nullptr, 1, LLKA_NO_ALTID, 1, LLKA_NO_INSCODE, 1, &coordsAA);

    LLKA_Point coordsB = { 10, 11, 12 };
    auto atomB = LLKA_makeAtom(3, "C", "C4'", "1", "DA", "A", nullptr, nullptr, nullptr, 2, LLKA_NO_ALTID, 2, LLKA_NO_INSCODE, 1, &coordsB);
    auto atomB2 = LLKA_makeAtom(4, "C", "C4'", "1", "DA", "A", nullptr, nullptr, nullptr, 2, LLKA_NO_ALTID, 2, LLKA_NO_INSCODE, 2, &coordsB);

    EFF_expect(LLKA_compareAtoms(&atomA, &atomAA, LLKA_FALSE), LLKA_FALSE, "compare identical atoms");
    EFF_expect(LLKA_compareAtoms(&atomA, &atomAA, LLKA_TRUE), LLKA_TRUE, "compare idential atoms (ignore ID)");
    EFF_expect(LLKA_compareAtoms(&atomA, &atomB, LLKA_TRUE), LLKA_FALSE, "compare different atoms");
    EFF_expect(LLKA_compareAtoms(&atomB, &atomB2, LLKA_FALSE), LLKA_FALSE, "compare atoms from different models (ignore ID)");

    LLKA_destroyAtom(&atomA);
    LLKA_destroyAtom(&atomAA);
    LLKA_destroyAtom(&atomB);
    LLKA_destroyAtom(&atomB2);
}

static
auto testDuplication()
{
    LLKA_Point coordsA = { 10, 11, 12 };
    auto atomA = LLKA_makeAtom(1, "C", "C5'", "1", "DA", "A", nullptr, nullptr, nullptr, 1, LLKA_NO_ALTID, 1, LLKA_NO_INSCODE, 1, &coordsA);

    LLKA_Atom atomB;
    LLKA_duplicateAtom(&atomA, &atomB);

    EFF_expect(LLKA_compareAtoms(&atomA, &atomB, LLKA_FALSE), LLKA_TRUE, "duplicate atom");

    LLKA_destroyAtom(&atomA);
    LLKA_destroyAtom(&atomB);
}

static
auto testFindAtoms()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1DK1_B_27_28_ATOMS, REAL_1DK1_B_27_28_ATOMS_LEN);
    LLKA_Atom *atoms = stru.atoms;

    auto atom = LLKA_findAtom(&stru, "OP2", "A", "B", 27, 'B', LLKA_NO_INSCODE, 1);
    EFF_expect(atom, &atoms[5] , "expected atom not found");

    atom = LLKA_findAtom(&stru, "OP2", "A", "B", 27, 'A', LLKA_NO_INSCODE, 1);
    EFF_expect(atom, &atoms[4] , "expected atom not found");

    auto atom2 = LLKA_findAtom(&stru, "OP2", "A", "B", 27, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    EFF_expect(atom, atom2, "expected atom not found");

    atom = LLKA_findAtom(&stru, "O5'", nullptr, nullptr, -1, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    EFF_expect(atom, &atoms[6], "expected atom not found");

    atom = LLKA_findAtom(&stru, "O5'", nullptr, nullptr, 29, LLKA_NO_ALTID, LLKA_NO_INSCODE, 1);
    EFF_expect(atom, nullptr, "did not expect to find an atom");

    LLKA_destroyStructure(&stru);
}

static
auto testRemoval()
{
    const LLKA_Atom atoms[] = {
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 1, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 2, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 3, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID },
        { "C", "C5'", "1", "DC", "A", "C5'", "DC", "A", { 1.059, -3.973, 1.922 }, 4, 1, 1, 1, LLKA_NO_INSCODE, LLKA_NO_ALTID }
    };

    // Test removal from the front
    {
        LLKA_Structure stru = LLKA_makeStructure(atoms, 4);
        LLKA_removeAtomById(1, &stru);

        EFF_expect(stru.nAtoms, 3UL, "remove atom from front, count");
        EFF_expect(stru.atoms[0].id, 2U, "remove atom from the front, first atom");
        EFF_expect(stru.atoms[1].id, 3U, "remove atom from the front, second atom");
        EFF_expect(stru.atoms[2].id, 4U, "remove atom from the front, third atom");

        LLKA_destroyStructure(&stru);
    }

    // Test removal from the end
    {
        LLKA_Structure stru = LLKA_makeStructure(atoms, 4);
        LLKA_removeAtomById(4, &stru);

        EFF_expect(stru.nAtoms, 3UL, "remove atom from the end, count");
        EFF_expect(stru.atoms[0].id, 1U, "remove atom from the end, first atom");
        EFF_expect(stru.atoms[1].id, 2U, "remove atom from the end, second atom");
        EFF_expect(stru.atoms[2].id, 3U, "remove atom from the end, third atom");

        LLKA_destroyStructure(&stru);
    }

    // Test removal from the middle
    {
        LLKA_Structure stru = LLKA_makeStructure(atoms, 4);
        LLKA_removeAtomById(2, &stru);

        EFF_expect(stru.nAtoms, 3UL, "remove atom from the middle, count");
        EFF_expect(stru.atoms[0].id, 1U, "remove atom from the middle, first atom");
        EFF_expect(stru.atoms[1].id, 3U, "remove atom from the middle, second atom");
        EFF_expect(stru.atoms[2].id, 4U, "remove atom from the middle, third atom");

        LLKA_destroyStructure(&stru);
    }
}

auto main(int, char **) -> int
{
    testComparison();
    testDuplication();
    testAppendAtom();
    testAppendStructure();
    testRemoval();
    testFindAtoms();

    return EXIT_SUCCESS;
}
