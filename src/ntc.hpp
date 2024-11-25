/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _NTC_HPP
#define _NTC_HPP

#include <llka_ntc.h>

#include <extract.hpp>
#include "residues.h"
#include "structure.hpp"
#include "util/geometry.h"

#include <array>
#include <cassert>
#include <memory>
#include <string>
#include <tuple>

namespace LLKAInternal {

inline const std::array<std::pair<LLKA_DinucleotideTorsion, double LLKA_StepMetrics::*>, 9> DINU_TORSIONS{{
    { LLKA_TOR_DELTA_1, &LLKA_StepMetrics::delta_1 },
    { LLKA_TOR_EPSILON_1, &LLKA_StepMetrics::epsilon_1 },
    { LLKA_TOR_ZETA_1, &LLKA_StepMetrics::zeta_1 },
    { LLKA_TOR_ALPHA_2, &LLKA_StepMetrics::alpha_2 },
    { LLKA_TOR_BETA_2, &LLKA_StepMetrics::beta_2 },
    { LLKA_TOR_GAMMA_2, &LLKA_StepMetrics::gamma_2 },
    { LLKA_TOR_DELTA_2, &LLKA_StepMetrics::delta_2 },
    { LLKA_TOR_CHI_1, &LLKA_StepMetrics::chi_1 },
    { LLKA_TOR_CHI_2, &LLKA_StepMetrics::chi_2 }
}};

inline
auto dinucleotideTorsion_unchecked(LLKA_DinucleotideTorsion torsion, const LLKA_Structure *stru, LLKA_StructureView &view)
{
    // This is supposed to be unsafe so we do sanity checks only through asserts
    assert(size_t(torsion) < 9); // Backbone torsions plus two base torsions
    assert(view.capacity >= 4);

    const auto firstSeqId = stru->atoms[0].label_seq_id;
    const std::string firstCompId = stru->atoms[0].label_comp_id;
    const auto secondSeqId = stru->atoms[stru->nAtoms - 1].label_seq_id;
    const std::string secondCompId = stru->atoms[stru->nAtoms - 1].label_comp_id;
    // Something to think about: Is the expectation that the second residue must
    // always have a "higher" sequence id that the first one really reasonable?
    // Would not a check that the sequence ids differ be better?
    assert(firstSeqId < secondSeqId);

    const auto &boneFirst = LLKAInternal::findBone(stru->atoms[0].label_comp_id);
    const auto &boneSecond = LLKAInternal::findBone(stru->atoms[stru->nAtoms - 1].label_comp_id);

    size_t atomIdx = 0;

    if (torsion < LLKA_TOR_CHI_1) {
        // We need just the backbone atoms to calculate backbone torsions
        assert(stru->nAtoms >= LLKABones::NUM_PLAIN_ATOMS);

        LLKAInternal::AtomToExtract ate{stru->atoms[0].pdbx_PDB_model_num, stru->atoms[0].label_asym_id, -1, firstCompId, ""};

        const auto &quads = LLKABones::BACKBONE_QUADS(boneFirst, boneSecond);
        const auto &atomsToGet = quads[torsion];
        for (const auto &[name, resNo] : atomsToGet) {
            ate.seqId = resNo == 1 ? firstSeqId : secondSeqId;
            ate.compId = resNo == 1 ? firstCompId : secondCompId;
            ate.name = name;
            auto atom = getMatchingAtom(stru, ate);
            if (!atom)
                return LLKA_E_INVALID_ARGUMENT;

            view.atoms[atomIdx++] = atom;
        }
    } else {
        // We must have the full "metrics substructure" including the four base-dependent atoms
        assert(stru->nAtoms >= LLKABones::NUM_TORSIONS_ATOMS);

        LLKA_BaseKind baseKind;
        const auto baseDetectionAtom = torsion == LLKA_TOR_CHI_1 ? &stru->atoms[0] : &stru->atoms[stru->nAtoms - 1];
        const auto seqId = torsion == LLKA_TOR_CHI_1 ? firstSeqId : secondSeqId;
        std::string compId = torsion == LLKA_TOR_CHI_1 ? firstCompId : secondCompId;

        if (LLKA_baseKind(baseDetectionAtom->label_comp_id, &baseKind) != LLKA_OK)
            return LLKA_E_INVALID_ARGUMENT;

        const auto &baseQuad = (torsion == LLKA_TOR_CHI_1 ? boneFirst : boneSecond).baseQuad;

        LLKAInternal::AtomToExtract ate{stru->atoms[0].pdbx_PDB_model_num, stru->atoms[0].label_asym_id, seqId, std::move(compId)};
        for (const auto &name : baseQuad) {
            ate.name = name;
            auto atom = LLKAInternal::getMatchingAtom(stru, ate);
            if (!atom)
                return LLKA_E_INVALID_ARGUMENT;

            view.atoms[atomIdx++] = atom;
        }
    }

    assert(atomIdx == 4);

    view.nAtoms = 4;

    return LLKA_OK;
}

inline
auto crossResidueMetric(LLKA_CrossResidueMetric metric, const LLKA_Structure *stru, LLKA_StructureView &view)
{
    assert(metric <= LLKA_XR_TOR_MU);
    assert(view.capacity >= 4);

    // We must have the full "metrics substructure" including the four base-dependent atoms
    if (stru->nAtoms < LLKABones::NUM_TORSIONS_ATOMS)
        return LLKA_E_INVALID_ARGUMENT;

    const auto modelNum = stru->atoms[0].pdbx_PDB_model_num;
    std::string asymId = stru->atoms[0].label_asym_id;
    const auto seqIdFirst = stru->atoms[0].label_seq_id;
    const std::string compIdFirst = stru->atoms[0].label_comp_id;
    const auto seqIdSecond = stru->atoms[stru->nAtoms - 1].label_seq_id;
    const std::string compIdSecond = stru->atoms[stru->nAtoms - 1].label_comp_id;

    const auto &boneFirst = LLKAInternal::findBone(stru->atoms[0].label_comp_id);
    const auto &boneSecond = LLKAInternal::findBone(stru->atoms[stru->nAtoms - 1].label_comp_id);

    LLKAInternal::AtomToExtract ateFirst{modelNum, asymId, seqIdFirst, ""};
    LLKAInternal::AtomToExtract ateSecond{modelNum, asymId, seqIdSecond, ""};

    if (metric == LLKA_XR_DIST_CC) {
        std::unique_ptr<LLKA_Atom[]> atoms{new LLKA_Atom[2]};
        ateFirst.compId = compIdFirst;
        ateFirst.name = "C1'";
        auto atom = getMatchingAtom(stru, ateFirst);
        if (!atom)
            return LLKA_E_INVALID_ARGUMENT;

        view.atoms[0] = atom;

        ateSecond.compId = compIdSecond;
        ateSecond.name = "C1'";
        atom = LLKAInternal::getMatchingAtom(stru, ateSecond);
        if (!atom)
            return LLKA_E_INVALID_ARGUMENT;

        view.atoms[1] = atom;

        view.nAtoms = 2;

        return LLKA_OK;
    } else if (metric == LLKA_XR_DIST_NN) {
        const auto &firstResidueBaseDepAtoms = boneFirst.base;
        const auto &secondResidueBaseDepAtoms = boneSecond.base;

        std::unique_ptr<LLKA_Atom[]> atoms{new LLKA_Atom[2]};
        ateFirst.compId = compIdFirst;
        ateFirst.name = firstResidueBaseDepAtoms[0];
        auto atom = LLKAInternal::getMatchingAtom(stru, ateFirst);
        if (!atom)
            return LLKA_E_INVALID_ARGUMENT;

        view.atoms[0] = atom;

        ateSecond.compId = compIdSecond;
        ateSecond.name = secondResidueBaseDepAtoms[0];
        atom = getMatchingAtom(stru, ateSecond);
        if (!atom)
            return LLKA_E_INVALID_ARGUMENT;

        view.atoms[1] = atom;

        view.nAtoms = 2;

        return LLKA_OK;
    } else {
        const auto &firstResidueBaseDepAtoms = boneFirst.base;
        const auto &secondResidueBaseDepAtoms = boneSecond.base;

        // PERF: Look into having a static array
        std::vector<AtomToExtract> toExtract{
            { modelNum, asymId, seqIdFirst, compIdFirst, firstResidueBaseDepAtoms[0] },
            { modelNum, asymId, seqIdFirst, compIdFirst, "C1'" },
            { modelNum, asymId, seqIdSecond, compIdSecond, "C1'" },
            { modelNum, asymId, seqIdSecond, compIdSecond, secondResidueBaseDepAtoms[0] },
        };

        return extractAtoms(stru, toExtract, &view);
    }
}

inline
auto calculateStepMetrics_unchecked(const LLKA_Structure *stru, LLKA_StepMetrics *metrics) -> LLKA_RetCode
{
    LLKA_RetCode tRet;
    LLKA_StructureView view = makeStructureView(4);

    for (const auto &[tor, clsPtr] : DINU_TORSIONS) {
        tRet = dinucleotideTorsion_unchecked(tor, stru, view);
        if (tRet != LLKA_OK) {
            LLKA_destroyStructureView(&view);
            return tRet;
        }

        auto v = dihedralAngle<double>(view);

        if (std::isnan(v)) {
            LLKA_destroyStructureView(&view);
            return LLKA_E_BAD_DATA;
        }

        metrics->*clsPtr = v;
    }

    tRet = crossResidueMetric(LLKA_XR_DIST_CC, stru, view);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&view);
        return tRet;
    }
    metrics->CC = spatialDistance<double>(*view.atoms[0], *view.atoms[1]);
    if (std::isnan(metrics->CC)) {
        LLKA_destroyStructureView(&view);
        return LLKA_E_BAD_DATA;
    }

    tRet = crossResidueMetric(LLKA_XR_DIST_NN, stru, view);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&view);
        return tRet;
    }
    metrics->NN = spatialDistance<double>(*view.atoms[0], *view.atoms[1]);
    if (std::isnan(metrics->NN)) {
        LLKA_destroyStructureView(&view);
        return LLKA_E_BAD_DATA;
    }

    tRet = crossResidueMetric(LLKA_XR_TOR_MU, stru, view);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&view);
        return tRet;
    }
    metrics->mu = dihedralAngle<double>(view);
    if (std::isnan(metrics->mu)) {
        LLKA_destroyStructureView(&view);
        return LLKA_E_BAD_DATA;
    }

    LLKA_destroyStructureView(&view);
    return LLKA_OK;
}

} // namespace LLKA

#endif // _NTC_HPP
