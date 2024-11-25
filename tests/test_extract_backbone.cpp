/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "effedup.hpp"
#include "testing_structures.h"

#include <llka_ntc.h>

#include <array>
#include <tuple>

static const std::array<std::tuple<std::string, int32_t>, 12> EXPECTED_ORDER{{
    {"C5'", 1},
    {"C4'", 1},
    {"O4'", 1},
    {"C3'", 1},
    {"O3'", 1},
    {"P", 2},
    {"O5'", 2},
    {"C5'", 2},
    {"C4'", 2},
    {"O4'", 2},
    {"C3'", 2},
    {"O3'", 2},
}};

static
auto testPlainBackbone()
{
    LLKA_Structure stru = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA_Structure backbone;

    EFF_expect(LLKA_extractBackbone(&stru, &backbone), LLKA_OK, "LLKA_extractBackbone() returned unexpected value")
    EFF_expect(backbone.nAtoms, EXPECTED_ORDER.size(), "number of backbone atoms is wrong")

    for (size_t idx = 0; idx < backbone.nAtoms; idx++) {
        const auto &atom = backbone.atoms[idx];
        const auto &[name, seqId] = EXPECTED_ORDER[idx];
        EFF_expect(std::string{atom.label_atom_id}, name, "label_atom_id is wrong");
        EFF_expect(atom.label_seq_id, seqId, "label_seq_id is wrong");
    }

    LLKA_destroyStructure(&backbone);
    LLKA_destroyStructure(&stru);
}

static
auto testPlainBackboneView()
{
    LLKA_Structure stru = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA_StructureView backbone;

    EFF_expect(LLKA_extractBackboneView(&stru, &backbone), LLKA_OK, "LLKA_extractBackbone() returned unexpected value")
    EFF_expect(backbone.nAtoms, EXPECTED_ORDER.size(), "number of backbone atoms is wrong")

    for (size_t idx = 0; idx < backbone.nAtoms; idx++) {
        const auto atom = backbone.atoms[idx];
        const auto &[name, seqId] = EXPECTED_ORDER[idx];
        EFF_expect(std::string{atom->label_atom_id}, name, "label_atom_id is wrong");
        EFF_expect(atom->label_seq_id, seqId, "label_seq_id is wrong");
    }

    LLKA_destroyStructureView(&backbone);
    LLKA_destroyStructure(&stru);
}

static
auto testPlainBackboneAltPos()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1DK1_B_27_28_ATOMS, REAL_1DK1_B_27_28_ATOMS_LEN);
    LLKA_Structure backbone;

    EFF_expect(
        LLKA_extractBackbone(&stru, &backbone),
        LLKA_E_MULTIPLE_ALT_IDS,
        "LLKA_extractBackbone() was supposed to detect multiple alternate positions but it did not"
    )

    LLKA_destroyStructure(&stru);

    stru = LLKA_makeStructure(REAL_1DK1_B_27_28_ALT_A_ONLY_ATOMS, REAL_1DK1_B_27_28_ALT_A_ONLY_ATOMS_LEN);

    EFF_expect(
        LLKA_extractBackbone(&stru, &backbone),
        LLKA_OK,
        "LLKA_extractBackbone() returned unexpected value"
    )

    LLKA_destroyStructure(&backbone);
    LLKA_destroyStructure(&stru);
}

static
auto testExtendedBackbone()
{
    std::array<std::tuple<std::string, int32_t>, 18> EXTENDED_EXPECTED_ORDER{{
        {"C5'", 1},
        {"C4'", 1},
        {"O4'", 1},
        {"C3'", 1},
        {"O3'", 1},
        {"C1'", 1},
        {"XXX", 1}, // Either N1 or N1, depepens on the base
        {"XXX", 1}, // Either C4 or C3, depends on the base
        {"P", 2},
        {"O5'", 2},
        {"C5'", 2},
        {"C4'", 2},
        {"O4'", 2},
        {"C3'", 2},
        {"O3'", 2},
        {"C1'", 2},
        {"XXX", 2}, // Either N1 or N1, depepens on the base
        {"XXX", 2}  // Either C4 or C3, depends on the base
    }};

    LLKA_Structure stru = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA_Structure backbone;
    LLKA_StepInfo info;

    auto tRet = LLKA_structureIsStep(&stru, &info);
    EFF_expect(tRet, LLKA_OK, "LLKA_structureIsStep() returned unexpected value")

    if (info.firstBaseKind == LLKA_PURINE) {
        std::get<0>(EXTENDED_EXPECTED_ORDER[6]) = "N9";
        std::get<0>(EXTENDED_EXPECTED_ORDER[7]) = "C4";
    } else {
        std::get<0>(EXTENDED_EXPECTED_ORDER[6]) = "N1";
        std::get<0>(EXTENDED_EXPECTED_ORDER[7]) = "C2";
    }
    if (info.secondBaseKind == LLKA_PURINE) {
        std::get<0>(EXTENDED_EXPECTED_ORDER[16]) = "N9";
        std::get<0>(EXTENDED_EXPECTED_ORDER[17]) = "C4";
    } else {
        std::get<0>(EXTENDED_EXPECTED_ORDER[16]) = "N1";
        std::get<0>(EXTENDED_EXPECTED_ORDER[17]) = "C2";
    }

    EFF_expect(LLKA_extractExtendedBackbone(&stru, &backbone), LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value");
    EFF_expect(backbone.nAtoms, EXTENDED_EXPECTED_ORDER.size(), "number of backbone atoms is wrong");

    for (size_t idx = 0; idx < backbone.nAtoms; idx++) {
        const auto &atom = backbone.atoms[idx];
        const auto &[name, seqId] = EXTENDED_EXPECTED_ORDER[idx];
        EFF_expect(std::string{atom.label_atom_id}, name, "label_atom_id is wrong");
        EFF_expect(atom.label_seq_id, seqId, "label_seq_id is wrong");
    }

    LLKA_destroyStructure(&backbone);
    LLKA_destroyStructure(&stru);
}

static
auto testExtendedBackboneView()
{
    std::array<std::tuple<std::string, int32_t>, 18> EXTENDED_EXPECTED_ORDER{{
        {"C5'", 1},
        {"C4'", 1},
        {"O4'", 1},
        {"C3'", 1},
        {"O3'", 1},
        {"C1'", 1},
        {"XXX", 1}, // Either N1 or N1, depepens on the base
        {"XXX", 1}, // Either C4 or C3, depends on the base
        {"P", 2},
        {"O5'", 2},
        {"C5'", 2},
        {"C4'", 2},
        {"O4'", 2},
        {"C3'", 2},
        {"O3'", 2},
        {"C1'", 2},
        {"XXX", 2}, // Either N1 or N1, depepens on the base
        {"XXX", 2}  // Either C4 or C3, depends on the base
    }};

    LLKA_Structure stru = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA_StructureView backbone;
    LLKA_StepInfo info;

    auto tRet = LLKA_structureIsStep(&stru, &info);
    EFF_expect(tRet, LLKA_OK, "LLKA_structureIsStep() returned unexpected value")

    if (info.firstBaseKind == LLKA_PURINE) {
        std::get<0>(EXTENDED_EXPECTED_ORDER[6]) = "N9";
        std::get<0>(EXTENDED_EXPECTED_ORDER[7]) = "C4";
    } else {
        std::get<0>(EXTENDED_EXPECTED_ORDER[6]) = "N1";
        std::get<0>(EXTENDED_EXPECTED_ORDER[7]) = "C2";
    }
    if (info.secondBaseKind == LLKA_PURINE) {
        std::get<0>(EXTENDED_EXPECTED_ORDER[16]) = "N9";
        std::get<0>(EXTENDED_EXPECTED_ORDER[17]) = "C4";
    } else {
        std::get<0>(EXTENDED_EXPECTED_ORDER[16]) = "N1";
        std::get<0>(EXTENDED_EXPECTED_ORDER[17]) = "C2";
    }

    EFF_expect(LLKA_extractExtendedBackboneView(&stru, &backbone), LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value");
    EFF_expect(backbone.nAtoms, EXTENDED_EXPECTED_ORDER.size(), "number of backbone atoms is wrong");

    for (size_t idx = 0; idx < backbone.nAtoms; idx++) {
        const auto &atom = backbone.atoms[idx];
        const auto &[name, seqId] = EXTENDED_EXPECTED_ORDER[idx];
        EFF_expect(std::string{atom->label_atom_id}, name, "label_atom_id is wrong");
        EFF_expect(atom->label_seq_id, seqId, "label_seq_id is wrong");
    }

    LLKA_destroyStructureView(&backbone);
    LLKA_destroyStructure(&stru);
}

static
auto testExtendedBackbonePseudouridine()
{
    std::array<std::tuple<std::string, int32_t>, 18> EXTENDED_EXPECTED_ORDER{{
        {"C5'", 12},
        {"C4'", 12},
        {"O4'", 12},
        {"C3'", 12},
        {"O3'", 12},
        {"C1'", 12},
        {"XXX", 12}, // Either N1 or N1, depepens on the base
        {"XXX", 12}, // Either C4 or C3, depends on the base
        {"P", 13},
        {"O5'", 13},
        {"C5'", 13},
        {"C4'", 13},
        {"O4'", 13},
        {"C3'", 13},
        {"O3'", 13},
        {"C1'", 13},
        {"C5", 13},
        {"C4", 13}
    }};

    LLKA_Structure stru = LLKA_makeStructure(REAL_1BZT_A_12_13_ATOMS, REAL_1BZT_A_12_13_ATOMS_LEN);
    LLKA_Structure backbone;
    LLKA_StepInfo info;

    auto tRet = LLKA_structureIsStep(&stru, &info);
    EFF_expect(tRet, LLKA_OK, "LLKA_structureIsStep() returned unexpected value")
    EFF_expect(info.secondBaseKind, LLKA_NON_STANDARD_BASE, "Wrong kind of second base")

    if (info.firstBaseKind == LLKA_PURINE) {
        std::get<0>(EXTENDED_EXPECTED_ORDER[6]) = "N9";
        std::get<0>(EXTENDED_EXPECTED_ORDER[7]) = "C4";
    } else {
        std::get<0>(EXTENDED_EXPECTED_ORDER[6]) = "N1";
        std::get<0>(EXTENDED_EXPECTED_ORDER[7]) = "C2";
    }

    EFF_expect(LLKA_extractExtendedBackbone(&stru, &backbone), LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value");
    EFF_expect(backbone.nAtoms, EXTENDED_EXPECTED_ORDER.size(), "number of backbone atoms is wrong");

    for (size_t idx = 0; idx < backbone.nAtoms; idx++) {
        const auto &atom = backbone.atoms[idx];
        const auto &[name, seqId] = EXTENDED_EXPECTED_ORDER[idx];
        EFF_expect(std::string{atom.label_atom_id}, name, "label_atom_id is wrong");
        EFF_expect(atom.label_seq_id, seqId, "label_seq_id is wrong");
    }

    LLKA_destroyStructure(&backbone);
    LLKA_destroyStructure(&stru);
}

auto main() -> int
{
    testPlainBackbone();
    testPlainBackboneView();
    testPlainBackboneAltPos();
    testExtendedBackbone();
    testExtendedBackboneView();
    testExtendedBackbonePseudouridine();
    // TODO: We should add a test with a residue with non-standard backbone
    return EXIT_SUCCESS;
}
