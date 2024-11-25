/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _NUCLEOTIDE_HPP
#define _NUCLEOTIDE_HPP

#include <llka_nucleotide.h>

#include <extract.hpp>
#include <residues.h>

#include <structure.hpp>
#include <util/geometry.h>
#include <util/templates.hpp>

#include <array>
#include <cstring>
#include <string>

template <typename StructureType> requires LLKAInternal::LLKAStructureType<StructureType>
inline
auto _extractNucleotideAtoms(const StructureType *stru, int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id)
{
    const auto ate = LLKAInternal::AtomToExtract{pdbx_PDB_model_num, label_asym_id, label_seq_id, {}};
    const auto allAtoms = LLKAInternal::getAllMatchingAtoms(stru, ate);

    std::vector<LLKA_Atom *> filteredAtoms{};

    // Find first atom whose compId resembles a known residue
    LLKABones::ANString compId{};
    for (const auto &at : allAtoms) {
        if (LLKAInternal::isKnownResidue(at->label_comp_id)) {
            compId = LLKABones::ANString{at->label_comp_id};
            break;
        }
    }

    // Do we even have anything that might be a nucleotide?
    if (compId.empty())
        return filteredAtoms;

    // Filter out all atoms that do not belong to a residue
    // NOTE: Do we need to care for microheterogentiy here?
    filteredAtoms.reserve(allAtoms.size());
    for (const auto &at : allAtoms) {
        if (std::strcmp(at->label_comp_id, compId.c_str()) == 0)
            filteredAtoms.push_back(at);
    }

    return filteredAtoms;
}

namespace LLKAInternal {

inline const std::map<std::string, LLKA_SugarPucker> NAME_TO_SUGAR_PUCKER_MAPPING{{
    // C1'
    { "C1end", LLKA_C1_ENDO }, { "C1endo", LLKA_C1_ENDO }, { "C1'end", LLKA_C1_ENDO }, { "C1'endo", LLKA_C1_ENDO }, { "C1' endo", LLKA_C1_ENDO },
    { "C1exo", LLKA_C1_EXO }, { "C1'exo", LLKA_C1_EXO }, { "C1' exo", LLKA_C1_EXO },
    // C2'
    { "C2end", LLKA_C2_ENDO }, { "C2endo", LLKA_C2_ENDO }, { "C2'end", LLKA_C2_ENDO }, { "C2'endo", LLKA_C2_ENDO }, { "C2' endo", LLKA_C2_ENDO },
    { "C2exo", LLKA_C2_EXO }, { "C2'exo", LLKA_C2_EXO }, { "C2' exo", LLKA_C2_EXO },
    // C3'
    { "C3end", LLKA_C3_ENDO }, { "C3endo", LLKA_C3_ENDO }, { "C3'end", LLKA_C3_ENDO }, { "C3'endo", LLKA_C3_ENDO }, { "C3' endo", LLKA_C3_ENDO },
    { "C3exo", LLKA_C3_EXO }, { "C3'exo", LLKA_C3_EXO }, { "C3' exo", LLKA_C3_EXO },
    // C4'
    { "C4end", LLKA_C4_ENDO }, { "C4endo", LLKA_C4_ENDO }, { "C4'end", LLKA_C4_ENDO }, { "C4'endo", LLKA_C4_ENDO }, { "C4' endo", LLKA_C4_ENDO },
    { "C4exo", LLKA_C4_EXO }, { "C4'exo", LLKA_C4_EXO }, { "C4' exo", LLKA_C4_EXO },
    // O4'
    { "O4end", LLKA_O4_ENDO }, { "O4endo", LLKA_O4_ENDO }, { "O4'end", LLKA_O4_ENDO }, { "O4'endo", LLKA_O4_ENDO }, { "O4' endo", LLKA_O4_ENDO },
    { "O1end", LLKA_O4_ENDO }, { "O1endo", LLKA_O4_ENDO }, { "O1'end", LLKA_O4_ENDO }, { "O1'endo", LLKA_O4_ENDO }, { "O1' endo", LLKA_O4_ENDO },
    { "O4exo", LLKA_O4_EXO }, { "O4'exo", LLKA_O4_EXO }, { "O4' exo", LLKA_O4_EXO },
    { "O1exo", LLKA_O4_EXO }, { "O1'exo", LLKA_O4_EXO }, { "O1' exo", LLKA_O4_EXO }
}};


// Items in this array must be ordered in the same way as items in LLKA_SugarPucker enum!
inline constinit std::array<char[7], 10> SUGAR_PUCKER_NAMES_VERY_TERSE{
    "C3end",
    "C4exo",
    "O4end",
    "C1exo",
    "C2end",
    "C3exo",
    "C4end",
    "O4exo",
    "C1end",
    "C2exo"
};

// Items in this array must be ordered in the same way as items in LLKA_SugarPucker enum!
inline constinit std::array<char[7], 10> SUGAR_PUCKER_NAMES_TERSE{
    "C3endo",
    "C4exo",
    "O4endo",
    "C1exo",
    "C2endo",
    "C3exo",
    "C4endo",
    "O4exo",
    "C1endo",
    "C2exo"
};

// Items in this array must be ordered in the same way as items in LLKA_SugarPucker enum!
inline constinit std::array<char[9], 10> SUGAR_PUCKER_NAMES_FANCY{
    "C3' endo",
    "C4' exo",
    "O4' endo",
    "C1' exo",
    "C2' endo",
    "C3' exo",
    "C4' endo",
    "O4' exo",
    "C1' endo",
    "C2' exo"
};

inline const std::array<std::string, 5> RIBOSE_CORE_ATOMS{
    "C4'", "O4'", "C1'", "C2'", "C3'"
};

inline
LLKA_SAD_CONSTEXPR_FUNC
auto calcRibosePseudorotation(double tau0, double tau1, double tau2, double tau3, double tau4)
{
    if (LLKA_WITHIN_INCLUSIVE(-5.0e-5, tau2, 5.0e-5))
        tau2 = sign(tau2) * 5.0e-5;

    // PERF: Simplify the sine thing!
    const auto tanP = (tau4 + tau1 - tau3 - tau0) / (2.0 * tau2 * (std::sin(D2R(36.0)) + std::sin(D2R(72.0))));
    auto P = std::atan(tanP);

    if (tau2 < 0.0)
        P += M_PI;
    else if (tanP < 0.0)
        P += TWO_PI;

    const auto tMax = std::abs(tau2 / std::cos(P));

    return std::make_tuple(P, tMax);
}

inline
auto extractNucleotide(const LLKA_Structure *stru, int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id)
{
    auto filteredAtoms = _extractNucleotideAtoms(stru, pdbx_PDB_model_num, label_asym_id, label_seq_id);

    if (filteredAtoms.size() == 0)
        return LLKA_Structure{ .atoms = nullptr, .nAtoms = 0 };

    return LLKA_makeStructureFromPtrs(filteredAtoms.data(), filteredAtoms.size());
}

template <typename StructureType> requires LLKAStructureType<StructureType>
inline
auto extractNucleotideView(const StructureType *stru, int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id)
{
    auto filteredAtoms = _extractNucleotideAtoms(stru, pdbx_PDB_model_num, label_asym_id, label_seq_id);

    if (filteredAtoms.size() == 0)
        return LLKA_StructureView{ .atoms = nullptr, .nAtoms = 0, .capacity = 0 };

    auto view = makeStructureView(filteredAtoms.size());
    for (size_t idx = 0; idx < filteredAtoms.size(); idx++)
        view.atoms[idx] = filteredAtoms[idx];

    view.nAtoms = filteredAtoms.size();

    return view;
}

template <typename StructureTypeSrc, typename StructureTypeDst> requires LLKAStructureType<StructureTypeSrc> && LLKAStructureType<StructureTypeDst>
inline
auto extractRibose(const StructureTypeSrc *stru, StructureTypeDst *riboseStru)
{
    std::vector<LLKAInternal::AtomToExtract> toExtract{RIBOSE_CORE_ATOMS.size()};

    const auto &atom = getAtom(*stru, 0);
    const auto modelNum = atom.pdbx_PDB_model_num;
    const std::string asymId = atom.label_asym_id;
    const auto seqId = atom.label_seq_id;
    const auto compId = atom.label_comp_id;

    for (size_t idx = 0; idx < RIBOSE_CORE_ATOMS.size(); idx++) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = seqId;
        toExtract[idx].compId = compId;
        toExtract[idx].name = RIBOSE_CORE_ATOMS[idx];
    }

    if constexpr (std::is_same_v<StructureTypeDst, LLKA_StructureView>) {
        *riboseStru = makeStructureView(toExtract.size());
    }

    auto tRet = extractAtoms(stru, toExtract, riboseStru);
    if constexpr (std::is_same_v<StructureTypeDst, LLKA_StructureView>) {
        if (tRet != LLKA_OK)
            LLKA_destroyStructureView(riboseStru);
    }

    return tRet;
}

inline
bool isNucleotideCompound(const std::string &compId)
{
    return LLKAInternal::isKnownResidue(compId);
}

inline
bool isNucleotideCompound(const char *compId)
{
    return isNucleotideCompound(std::string{compId});
}

template <typename StructureType> requires LLKAStructureType<StructureType>
inline
auto measureNuAngles(const StructureType &riboseStru)
{
    LLKA_NuAngles angles;

    // Rely on the ordering of atoms.
    // The order is reliable if riboseStru was extracted with LLKA_extractRibose()
    angles.nu_0 = dihedralAngle<double>(getAtom(riboseStru, 0).coords, getAtom(riboseStru, 1).coords, getAtom(riboseStru, 2).coords, getAtom(riboseStru, 3).coords);
    angles.nu_1 = dihedralAngle<double>(getAtom(riboseStru, 1).coords, getAtom(riboseStru, 2).coords, getAtom(riboseStru, 3).coords, getAtom(riboseStru, 4).coords);
    angles.nu_2 = dihedralAngle<double>(getAtom(riboseStru, 2).coords, getAtom(riboseStru, 3).coords, getAtom(riboseStru, 4).coords, getAtom(riboseStru, 0).coords);
    angles.nu_3 = dihedralAngle<double>(getAtom(riboseStru, 3).coords, getAtom(riboseStru, 4).coords, getAtom(riboseStru, 0).coords, getAtom(riboseStru, 1).coords);
    angles.nu_4 = dihedralAngle<double>(getAtom(riboseStru, 4).coords, getAtom(riboseStru, 0).coords, getAtom(riboseStru, 1).coords, getAtom(riboseStru, 2).coords);

    return angles;
}

inline
auto pseudorotationToSugarPucker(double P)
{
    // TODO: Covert to -PI <-> PI range
    P = angleAsFull(P);

    if (LLKA_WITHIN_LINCL_REXCL(0, P, D2R(36.0)))
        return LLKA_C3_ENDO;
    else if (LLKA_WITHIN_LINCL_REXCL(D2R(36.0), P, D2R(72.0)))
        return LLKA_C4_EXO;
    else if (LLKA_WITHIN_LINCL_REXCL(D2R(72.0), P, D2R(108.0)))
        return LLKA_O4_ENDO;
    else if (LLKA_WITHIN_LINCL_REXCL(D2R(108.0), P, D2R(144.0)))
        return LLKA_C1_EXO;
    else if (LLKA_WITHIN_LINCL_REXCL(D2R(144.0), P, D2R(180.0)))
        return LLKA_C2_ENDO;
    else if (LLKA_WITHIN_LINCL_REXCL(D2R(180.0), P, D2R(216.0)))
        return LLKA_C3_EXO;
    else if (LLKA_WITHIN_LINCL_REXCL(D2R(216.0), P, D2R(252.0)))
        return LLKA_C4_ENDO;
    else if (LLKA_WITHIN_LINCL_REXCL(D2R(252.0), P, D2R(288.0)))
        return LLKA_O4_EXO;
    else if (LLKA_WITHIN_LINCL_REXCL(D2R(288.0), P, D2R(324.0)))
        return LLKA_C1_ENDO;
    else if (LLKA_WITHIN_LECXL_RINCL(D2R(324.0), P, D2R(360.0)))
        return LLKA_C2_EXO;

    // This should never really happen but to be safe
    return LLKA_C3_ENDO;
}

template <typename StructureType> requires LLKAStructureType<StructureType>
inline
auto riboseMetrics(const StructureType &stru, LLKA_RiboseMetrics &metrics)
{
    metrics.nus = measureNuAngles(stru);
    auto [P, tMax] = calcRibosePseudorotation(
        metrics.nus.nu_0,
        metrics.nus.nu_1,
        metrics.nus.nu_2,
        metrics.nus.nu_3,
        metrics.nus.nu_4
    );

    metrics.P = P;
    metrics.tMax = tMax;
    metrics.pucker = pseudorotationToSugarPucker(metrics.P);
}

} // namespace LLKAInternal

#endif // _NUCLEOTIDE_HPP
