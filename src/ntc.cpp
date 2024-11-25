/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "ntc.hpp"

#include "ntc_constants.h"
#include "residues.h"
#include "structure.hpp"
#include "util/elementaries.h"
#include "util/templates.hpp"

#include <cstring>
#include <limits>
#include <map>
#include <vector>

namespace LLKAInternal {

#define MK_BKBN_ATOM_TAG_FIRST(name) (name + "f")
#define MK_BKBN_ATOM_TAG_SECOND(name) (name + "s")

template <typename T>
static
auto atomPropertyMatches(const std::vector<const LLKA_Atom *> &atoms, T LLKA_Atom::* prop)
{
    if (atoms.empty())
        return false;

    const T &checkAgainst = atoms.front()->*prop;
    for (size_t idx = 1; idx < atoms.size(); idx++) {
        if (checkAgainst != atoms[idx]->*prop)
            return false;
    }

    return true;
}

static
auto hasSingleAltId(const std::vector<const LLKA_Atom *> &atoms)
{
    char altIdCheck = LLKA_NO_ALTID;
    for (const auto &atom : atoms) {
        if (atom->label_alt_id != LLKA_NO_ALTID) {
            if (altIdCheck == LLKA_NO_ALTID)
                altIdCheck = atom->label_alt_id;
            else if (altIdCheck != atom->label_alt_id)
                return false;
        }
    }

    return true;
}

static
auto crossResidueMetricAtoms(const char *firstBase, const char *secondBase, LLKA_CrossResidueMetric metric, LLKA_AtomNameQuad &quad)
{
    const auto &boneFirst = findBone(firstBase);
    const auto &boneSecond = findBone(secondBase);

    const auto &baseQuadFirst = boneFirst.baseQuad;
    const auto &baseQuadSecond = boneSecond.baseQuad;

    assert(baseQuadFirst[0].size() < 8); assert(baseQuadFirst[1].size() < 8); assert(baseQuadFirst[2].size() < 8); assert(baseQuadFirst[3].size() < 8);
    assert(baseQuadSecond[0].size() < 8); assert(baseQuadSecond[1].size() < 8); assert(baseQuadSecond[2].size() < 8); assert(baseQuadSecond[3].size() < 8);

    if (metric == LLKA_XR_DIST_CC) {
        std::strncpy(quad.a, baseQuadFirst[3].c_str(), sizeof(quad.a));
        std::strncpy(quad.b, baseQuadSecond[3].c_str(), sizeof(quad.b));
    } else if (metric == LLKA_XR_DIST_NN) {
        std::strncpy(quad.a, baseQuadFirst[2].c_str(), sizeof(quad.a));
        std::strncpy(quad.b, baseQuadSecond[2].c_str(), sizeof(quad.b));
    } else {
        std::strncpy(quad.a, baseQuadFirst[3].c_str(), sizeof(quad.a));
        std::strncpy(quad.b, baseQuadFirst[2].c_str(), sizeof(quad.b));
        std::strncpy(quad.c, baseQuadSecond[2].c_str(), sizeof(quad.c));
        std::strncpy(quad.d, baseQuadSecond[3].c_str(), sizeof(quad.d));
    }

    return LLKA_OK;
}

static
auto dinucleotideTorsionAtoms(const char *firstBase, const char *secondBase, LLKA_DinucleotideTorsion torsion, LLKA_AtomNameQuad &quad)
{
    const auto &boneFirst = LLKAInternal::findBone(firstBase);
    const auto &boneSecond = LLKAInternal::findBone(secondBase);

    const LLKABones::ANString *a;
    const LLKABones::ANString *b;
    const LLKABones::ANString *c;
    const LLKABones::ANString *d;

    if (torsion < LLKA_TOR_CHI_1) {
        const auto &_quads = LLKABones::BACKBONE_QUADS(boneFirst, boneSecond);
        const auto &_quad = _quads[torsion];

        a = &std::get<0>(_quad[0]);
        b = &std::get<0>(_quad[1]);
        c = &std::get<0>(_quad[2]);
        d = &std::get<0>(_quad[3]);
    } else {
        const auto &_quad = (torsion == LLKA_TOR_CHI_1 ? boneFirst : boneSecond).baseQuad;

        a = &_quad[0];
        b = &_quad[1];
        c = &_quad[2];
        d = &_quad[3];
    }

    // Exported quads have a limit of 7 characters for atom name
    assert(a->length() < 8); assert(b->length() < 8); assert(c->length() < 8); assert(d->length() < 8);

    std::strncpy(quad.a, a->c_str(), sizeof(quad.a));
    std::strncpy(quad.b, b->c_str(), sizeof(quad.b));
    std::strncpy(quad.c, c->c_str(), sizeof(quad.c));
    std::strncpy(quad.d, d->c_str(), sizeof(quad.d));

    return LLKA_OK;
}

template <LLKAStructureType T>
static
LLKA_RetCode LLKA_CC extractBackbone(const LLKA_Structure *stru, T *backbone)
{
    LLKA_StepInfo info;
    LLKA_RetCode tRet;

    tRet = LLKA_structureIsStep(stru, &info);
    if (tRet != LLKA_OK)
        return tRet;

    const auto &boneFirst = LLKAInternal::findBone(stru->atoms[0].label_comp_id);
    const auto &boneSecond = LLKAInternal::findBone(stru->atoms[stru->nAtoms - 1].label_comp_id);

    std::vector<LLKAInternal::AtomToExtract> toExtract{LLKABones::NUM_PLAIN_ATOMS};

    auto modelNum = stru->atoms[0].pdbx_PDB_model_num;
    std::string asymId = stru->atoms[0].label_asym_id;
    std::string firstCompId = stru->atoms[0].label_comp_id;
    std::string secondCompId = stru->atoms[stru->nAtoms - 1].label_comp_id;

    size_t idx = 0;
    for (const auto &atomName : LLKABones::PLAIN_BACKBONE(boneFirst.firstResidue)) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.firstSeqId;
        toExtract[idx].compId = firstCompId;
        toExtract[idx].name = atomName;
        idx++;
    }
    for (const auto &atomName : LLKABones::PLAIN_BACKBONE(boneSecond.secondResidue)) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.secondSeqId;
        toExtract[idx].compId = secondCompId;
        toExtract[idx].name = atomName;
        idx++;
    }
    assert(idx == LLKABones::NUM_PLAIN_ATOMS);

    if constexpr (std::is_same_v<T, LLKA_StructureView>)
        initStructureView(*backbone, toExtract.size());

    return LLKAInternal::extractAtoms(stru, toExtract, backbone);
}

template <LLKAStructureType T>
static
LLKA_RetCode LLKA_CC extractExtendedBackbone(const LLKA_Structure *stru, T *backbone)
{
    LLKA_StepInfo info;
    LLKA_RetCode tRet;

    tRet = LLKA_structureIsStep(stru, &info);
    if (tRet != LLKA_OK)
        return tRet;

    const auto &boneFirst = LLKAInternal::findBone(stru->atoms[0].label_comp_id);
    const auto &boneSecond = LLKAInternal::findBone(stru->atoms[stru->nAtoms - 1].label_comp_id);

    std::vector<LLKAInternal::AtomToExtract> toExtract{LLKABones::NUM_EXTENDED_ATOMS};

    auto modelNum = stru->atoms[0].pdbx_PDB_model_num;
    std::string asymId = stru->atoms[0].label_asym_id;
    std::string firstCompId = stru->atoms[0].label_comp_id;
    std::string secondCompId = stru->atoms[stru->nAtoms - 1].label_comp_id;

    size_t idx = 0;

    // First residue
    for (const auto &atomName : LLKABones::EXTENDED_BACKBONE(boneFirst.firstResidue)) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.firstSeqId;
        toExtract[idx].compId = firstCompId;
        toExtract[idx].name = atomName;
        idx++;
    }

    for (const auto &baseAtom : boneFirst.base) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.firstSeqId;
        toExtract[idx].compId = firstCompId;
        toExtract[idx].name = baseAtom;
        idx++;
    }

    // Second residue
    for (const auto &atomName : LLKABones::EXTENDED_BACKBONE(boneSecond.secondResidue)) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.secondSeqId;
        toExtract[idx].compId = secondCompId;
        toExtract[idx].name = atomName;
        idx++;
    }

    for (const auto &baseAtom : boneSecond.base) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.secondSeqId;
        toExtract[idx].compId = secondCompId;
        toExtract[idx].name = baseAtom;
        idx++;
    }
    assert(idx == LLKABones::NUM_EXTENDED_ATOMS);

    if constexpr (std::is_same_v<T, LLKA_StructureView>)
        initStructureView(*backbone, toExtract.size());

    return LLKAInternal::extractAtoms(stru, toExtract, backbone);
}

} // namespace LLKAInternal

size_t LLKA_CC LLKA_backboneAtomIndex(const LLKA_BackboneAtom *bkbnAtom, const LLKA_Structure *backbone)
{
    if (backbone->nAtoms < LLKABones::NUM_PLAIN_ATOMS)
        return LLKA_INVALID_BKBN_ATOM_INDEX;

    // We expect backbone to be well-ordered
    const auto firstSeqId = backbone->atoms[0].label_seq_id;
    const auto secondSeqId = backbone->atoms[backbone->nAtoms - 1].label_seq_id;
    if (firstSeqId >= secondSeqId)
        return LLKA_INVALID_BKBN_ATOM_INDEX;

    const auto seqId = bkbnAtom->residue == LLKA_FIRST_RESIDUE ? firstSeqId : secondSeqId;
    for (size_t idx = 0; idx < backbone->nAtoms; idx++) {
        const auto &atom = backbone->atoms[idx];
        if (atom.label_seq_id != seqId)
            continue;
        if (!std::strcmp(atom.label_atom_id, bkbnAtom->name))
            return idx;
    }

    return LLKA_INVALID_BKBN_ATOM_INDEX;
}

LLKA_RetCode LLKA_CC LLKA_calculateStepMetrics(const LLKA_Structure *stru, LLKA_StepMetrics *metrics)
{
    LLKA_StepInfo info;

    auto tRet = LLKA_structureIsStep(stru, &info);
    if (tRet != LLKA_OK)
        return tRet;

    return LLKAInternal::calculateStepMetrics_unchecked(stru, metrics);
}

LLKA_RetCode LLKA_CC LLKA_calculateStepMetricsDifferenceAgainstReference(const LLKA_Structure *stru, LLKA_NtC ntc, LLKA_StepMetrics *metrics)
{
    if (ntc == LLKA_INVALID_NTC)
        return LLKA_E_INVALID_ARGUMENT;

    auto tRet = LLKA_calculateStepMetrics(stru, metrics);
    if (tRet != LLKA_OK)
        return tRet;

    const auto &ntcRef = LLKAInternal::NTC_AVERAGES[ntc];
    metrics->delta_1 = LLKAInternal::angleDifference(metrics->delta_1, ntcRef.delta_1);
    metrics->epsilon_1 = LLKAInternal::angleDifference(metrics->epsilon_1, ntcRef.epsilon_1);
    metrics->zeta_1 = LLKAInternal::angleDifference(metrics->zeta_1, ntcRef.zeta_1);
    metrics->alpha_2 = LLKAInternal::angleDifference(metrics->alpha_2, ntcRef.alpha_2);
    metrics->beta_2 = LLKAInternal::angleDifference(metrics->beta_2, ntcRef.beta_2);
    metrics->gamma_2 = LLKAInternal::angleDifference(metrics->gamma_2, ntcRef.gamma_2);
    metrics->delta_2 = LLKAInternal::angleDifference(metrics->delta_2, ntcRef.delta_2);
    metrics->chi_1 = LLKAInternal::angleDifference(metrics->chi_1, ntcRef.chi_1);
    metrics->chi_2 = LLKAInternal::angleDifference(metrics->chi_2, ntcRef.chi_2);
    metrics->CC -= ntcRef.CC;
    metrics->NN -= ntcRef.NN;
    metrics->mu = LLKAInternal::angleDifference(metrics->mu, ntcRef.mu);

    return LLKA_OK;
}

const char * LLKA_CC LLKA_CANAToName(LLKA_CANA cana)
{
    if (cana == LLKA_INVALID_CANA)
        return "NAN";
    const auto idx = static_cast<size_t>(cana);
    return LLKAInternal::CANA_NAMES[idx];
}

LLKA_RetCode LLKA_CC LLKA_crossResidueMetric(LLKA_CrossResidueMetric metric, const LLKA_Structure *stru, LLKA_Structure *metricStru)
{
    assert(metric <= LLKA_XR_TOR_MU);

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
        auto atom = LLKAInternal::getMatchingAtom(stru, ateFirst);
        if (!atom)
            return LLKA_E_INVALID_ARGUMENT;
        LLKA_duplicateAtom(atom, &atoms[0]);

        ateSecond.compId = compIdSecond;
        ateSecond.name = "C1'";
        atom = LLKAInternal::getMatchingAtom(stru, ateSecond);
        if (!atom) {
            LLKA_destroyAtom(&atoms[0]);
            return LLKA_E_INVALID_ARGUMENT;
        }

        LLKA_duplicateAtom(atom, &atoms[1]);

        metricStru->atoms = atoms.release();
        metricStru->nAtoms = 2;

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
        LLKA_duplicateAtom(atom, &atoms[0]);

        ateSecond.compId = compIdSecond;
        ateSecond.name = secondResidueBaseDepAtoms[0];
        atom = LLKAInternal::getMatchingAtom(stru, ateSecond);
        if (!atom) {
            LLKA_destroyAtom(&atoms[0]);
            return LLKA_E_INVALID_ARGUMENT;
        }

        LLKA_duplicateAtom(atom, &atoms[1]);

        metricStru->atoms = atoms.release();
        metricStru->nAtoms = 2;

        return LLKA_OK;
    } else {
        const auto &firstResidueBaseDepAtoms = boneFirst.base;
        const auto &secondResidueBaseDepAtoms = boneSecond.base;

        // PERF: Look into having a static array
        std::vector<LLKAInternal::AtomToExtract> toExtract{
            { modelNum, asymId, seqIdFirst, compIdFirst, firstResidueBaseDepAtoms[0] },
            { modelNum, asymId, seqIdFirst, compIdFirst, "C1'" },
            { modelNum, asymId, seqIdSecond, compIdSecond, "C1'" },
            { modelNum, asymId, seqIdSecond, compIdSecond, secondResidueBaseDepAtoms[0] },
        };

        return LLKAInternal::extractAtoms(stru, toExtract, metricStru);
    }
}

LLKA_RetCode LLKA_CC LLKA_crossResidueMetricAtomsFromBases(const char *firstBase, const char *secondBase, LLKA_CrossResidueMetric metric, LLKA_AtomNameQuad *quad)
{
    if (!LLKAInternal::isKnownResidue(firstBase) || !LLKAInternal::isKnownResidue(secondBase))
        return LLKA_E_INVALID_ARGUMENT;

    return LLKAInternal::crossResidueMetricAtoms(firstBase, secondBase,  metric, *quad);
}

LLKA_RetCode LLKA_CC LLKA_crossResidueMetricAtomsFromStructure(const LLKA_Structure *stru, LLKA_CrossResidueMetric metric, LLKA_AtomNameQuad *quad)
{
    LLKA_StepInfo info;

    auto tRet = LLKA_structureIsStep(stru, &info);
    if (tRet != LLKA_OK)
        return tRet;

    return LLKAInternal::crossResidueMetricAtoms(stru->atoms[0].label_comp_id, stru->atoms[stru->nAtoms - 1].label_comp_id, metric, *quad);
}

const char * LLKA_CC LLKA_crossResidueMetricName(LLKA_CrossResidueMetric metric, LLKA_Bool greek)
{
    switch (metric) {
    case LLKA_XR_DIST_CC:
        return greek ? "CC" : "CC";
    case LLKA_XR_DIST_NN:
        return greek ? "NN" : "NN";
    case LLKA_XR_TOR_MU:
        return greek ? "\xCE\xBC" : "mu";
    default:
        assert(false);
        return "Unknown cross-residue metric";
    }
}

LLKA_RetCode LLKA_CC LLKA_dinucleotideTorsion(LLKA_DinucleotideTorsion torsion, const LLKA_Structure *stru, LLKA_Structure *torsionStru)
{
    assert(size_t(torsion) < 9); // Backbone torsions plus two base torsions

    if (torsion < LLKA_TOR_CHI_1 && stru->nAtoms < LLKABones::NUM_PLAIN_ATOMS)
        return LLKA_E_INVALID_ARGUMENT; // We need the full backbone to calculate backbone torsions
    else if (stru->nAtoms < LLKABones::NUM_TORSIONS_ATOMS)
        return LLKA_E_INVALID_ARGUMENT; // We need the full "metrics" structure to calculate chi torsions

    const auto firstSeqId = stru->atoms[0].label_seq_id;
    const auto secondSeqId = stru->atoms[stru->nAtoms - 1].label_seq_id;
    if (firstSeqId >= secondSeqId)
        return LLKA_E_INVALID_ARGUMENT;

    auto view = LLKAInternal::makeStructureView(4);

    auto tRet = LLKAInternal::dinucleotideTorsion_unchecked(torsion, stru, view);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&view);
        return tRet;
    }

    LLKAInternal::structureViewToStructure(view, *torsionStru);
    LLKA_destroyStructureView(&view);

    return LLKA_OK;
}

LLKA_RetCode LLKA_CC LLKA_dinucleotideTorsionAtomsFromBases(const char *firstBase, const char *secondBase, LLKA_DinucleotideTorsion torsion, LLKA_AtomNameQuad *quad)
{
    if (!LLKAInternal::isKnownResidue(firstBase) || !LLKAInternal::isKnownResidue(secondBase))
        return LLKA_E_INVALID_ARGUMENT;

    return LLKAInternal::dinucleotideTorsionAtoms(firstBase, secondBase, torsion, *quad);
}

LLKA_RetCode LLKA_CC LLKA_dinucleotideTorsionAtomsFromStructure(const LLKA_Structure *stru, LLKA_DinucleotideTorsion torsion, LLKA_AtomNameQuad *quad)
{
    LLKA_StepInfo info;

    auto tRet = LLKA_structureIsStep(stru, &info);
    if (tRet != LLKA_OK)
        return tRet;

    return LLKAInternal::dinucleotideTorsionAtoms(stru->atoms[0].label_comp_id, stru->atoms[stru->nAtoms - 1].label_comp_id, torsion, *quad);
}

const char * LLKA_CC LLKA_dinucleotideTorsionName(LLKA_DinucleotideTorsion torsion, LLKA_Bool greek)
{
    switch (torsion) {
    case LLKA_TOR_DELTA_1:
        return greek ? "\xCE\xB4_1" : "delta_1";
    case LLKA_TOR_EPSILON_1:
        return greek ? "\xCE\xB5_1" : "epsilon_1";
    case LLKA_TOR_ZETA_1:
        return greek ? "\xCE\xB6_1" : "zeta_1";
    case LLKA_TOR_ALPHA_2:
        return greek ? "\xCE\xB1_2" : "alpha_2";
    case LLKA_TOR_BETA_2:
        return greek ? "\xCE\xB2_2" : "beta_2";
    case LLKA_TOR_GAMMA_2:
        return greek ? "\xCE\xB3_2" : "gamma_2";
    case LLKA_TOR_DELTA_2:
        return greek ? "\xCE\xB4_2" : "delta_2";
    case LLKA_TOR_CHI_1:
        return greek ? "\xCF\x87_1" : "chi_1";
    case LLKA_TOR_CHI_2:
        return greek ? "\xCF\x87_2" : "chi_2";
    default:
        assert(false);
        return "Unknown torsion";
    }
}

LLKA_RetCode LLKA_CC LLKA_extractBackbone(const LLKA_Structure *stru, LLKA_Structure *backbone)
{
    return LLKAInternal::extractBackbone(stru, backbone);
}

LLKA_RetCode LLKA_CC LLKA_extractBackboneView(const LLKA_Structure *stru, LLKA_StructureView *backbone)
{
    return LLKAInternal::extractBackbone(stru, backbone);
}

LLKA_RetCode LLKA_CC LLKA_extractExtendedBackbone(const LLKA_Structure *stru, LLKA_Structure *backbone)
{
    return LLKAInternal::extractExtendedBackbone(stru, backbone);
}

LLKA_RetCode LLKA_CC LLKA_extractExtendedBackboneView(const LLKA_Structure *stru, LLKA_StructureView *backbone)
{
    return LLKAInternal::extractExtendedBackbone(stru, backbone);
}

LLKA_RetCode LLKA_CC LLKA_extractMetricsStructure(const LLKA_Structure *stru, LLKA_Structure *metricsStru)
{
    LLKA_StepInfo info;
    LLKA_RetCode tRet;

    tRet = LLKA_structureIsStep(stru, &info);
    if (tRet != LLKA_OK)
        return tRet;

    const auto &boneFirst = LLKAInternal::findBone(stru->atoms[0].label_comp_id);
    const auto &boneSecond = LLKAInternal::findBone(stru->atoms[stru->nAtoms - 1].label_comp_id);

    std::vector<LLKAInternal::AtomToExtract> toExtract{LLKABones::NUM_TORSIONS_ATOMS};

    auto modelNum = stru->atoms[0].pdbx_PDB_model_num;
    std::string asymId = stru->atoms[0].label_asym_id;
    std::string firstCompId = stru->atoms[0].label_comp_id;
    std::string secondCompId = stru->atoms[stru->nAtoms - 1].label_comp_id;

    size_t idx = 0;

    // First residue
    for (const auto &atomName : LLKABones::TORSIONS_BACKBONE(boneFirst.firstResidue)) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.firstSeqId;
        toExtract[idx].compId = firstCompId;
        toExtract[idx].name = atomName;
        idx++;
    }

    for (const auto &baseAtom : boneFirst.base) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.firstSeqId;
        toExtract[idx].compId = firstCompId;
        toExtract[idx].name = baseAtom;
        idx++;
    }

    // Second residue
    for (const auto &atomName : LLKABones::TORSIONS_BACKBONE(boneSecond.secondResidue)) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.secondSeqId;
        toExtract[idx].compId = secondCompId;
        toExtract[idx].name = atomName;
        idx++;
    }

    for (const auto &baseAtom : boneSecond.base) {
        toExtract[idx].modelNum = modelNum;
        toExtract[idx].asymId = asymId;
        toExtract[idx].seqId = info.secondSeqId;
        toExtract[idx].compId = secondCompId;
        toExtract[idx].name = baseAtom;
        idx++;
    }
    assert(idx == LLKABones::NUM_TORSIONS_ATOMS);

    return LLKAInternal::extractAtoms(stru, toExtract, metricsStru);
}

LLKA_CANA LLKA_CC LLKA_nameToCANA(const char *name)
{
    auto it = std::find_if(LLKAInternal::CANA_NAMES.cbegin(), LLKAInternal::CANA_NAMES.cend(), [name](const auto &item) { return std::strcmp(item, name) == 0; });
    if (it == LLKAInternal::CANA_NAMES.cend())
        return LLKA_INVALID_CANA;
    return static_cast<LLKA_CANA>(std::distance(LLKAInternal::CANA_NAMES.cbegin(), it));
}

LLKA_NtC LLKA_CC LLKA_nameToNtC(const char *name)
{
    auto it = LLKAInternal::NAME_TO_NTC_INDEX_MAPPING.find(name);
    if (it == LLKAInternal::NAME_TO_NTC_INDEX_MAPPING.cend())
        return LLKA_INVALID_NTC;
    return static_cast<LLKA_NtC>(it->second);
}

const char * LLKA_CC LLKA_NtCToName(LLKA_NtC ntc)
{
    if (ntc == LLKA_INVALID_NTC)
        return "NANT";
    const auto idx = std::underlying_type_t<LLKA_NtC>(ntc);
    return LLKAInternal::NTC_AVERAGES[idx].name.c_str();
}

LLKA_Structure LLKA_CC LLKA_NtCStructure(LLKA_NtC ntc)
{
    if (ntc == LLKA_INVALID_NTC)
        return { .atoms = nullptr, .nAtoms = 0 };

    const auto &ref = LLKAInternal::NTC_REFERENCES[ntc];

    LLKA_Structure stru;
    LLKA_duplicateStructure(&ref, &stru);

    return stru;
}

LLKA_RetCode LLKA_CC LLKA_structureIsStep(const LLKA_Structure *stru, LLKA_StepInfo *info)
{
    static const auto atomComparator = [](const std::string &label_atom_id, const LLKA_Atom *const &atom) -> bool {
        return std::strcmp(atom->label_atom_id, label_atom_id.c_str()) == 0;
    };

    if (stru->nAtoms == 0)
        return LLKA_E_MISSING_ATOMS;

    const auto &boneFirst = LLKAInternal::findBone(stru->atoms[0].label_comp_id);
    const auto &boneSecond = LLKAInternal::findBone(stru->atoms[stru->nAtoms - 1].label_comp_id);

    std::vector<const LLKA_Atom *> atomsFirstResidue{};
    std::vector<const LLKA_Atom *> atomsSecondResidue{};

    // Gather all atoms of interest
    for (size_t idx = 0; idx < stru->nAtoms; idx++) {
        const auto atom = &stru->atoms[idx];

        for (const auto &atomName : LLKABones::EXTENDED_BACKBONE(boneFirst.firstResidue)) {
            if (!std::strcmp(atomName.c_str(), atom->label_atom_id)) {
                atomsFirstResidue.push_back(atom);
                break;
            }
        }

        for (const auto &atomName : LLKABones::EXTENDED_BACKBONE(boneSecond.secondResidue)) {
            if (!std::strcmp(atomName.c_str(), atom->label_atom_id)) {
                atomsSecondResidue.push_back(atom);
                break;
            }
        }
    }

    int32_t firstSeqId = std::numeric_limits<int32_t>::max();
    for (const auto &atom : atomsFirstResidue) {
        auto seqId = atom->label_seq_id;
        if (seqId < firstSeqId)
            firstSeqId = seqId;
    }

    LLKAInternal::filter(atomsFirstResidue, [firstSeqId](const LLKA_Atom *atom) { return atom->label_seq_id != firstSeqId; });
    LLKAInternal::filter(atomsSecondResidue, [firstSeqId](const LLKA_Atom *atom) { return atom->label_seq_id == firstSeqId; });

    // Do altId check. Since we allow the first and second step to have different altIds,
    // we need to check each nucleotide individually.
    if (!LLKAInternal::hasSingleAltId(atomsFirstResidue))
        return LLKA_E_MULTIPLE_ALT_IDS;
    if (!LLKAInternal::hasSingleAltId(atomsSecondResidue))
        return LLKA_E_MULTIPLE_ALT_IDS;

    // Check that the structure contains all atoms necessary to form the step backbone
    for (const auto &atomName : LLKABones::PLAIN_BACKBONE(boneFirst.firstResidue)) {
        if (!LLKAInternal::contains(atomsFirstResidue, atomName, atomComparator))
            return LLKA_E_MISSING_ATOMS;
    }
    for (const auto &atomName : LLKABones::PLAIN_BACKBONE(boneSecond.secondResidue)) {
        if (!LLKAInternal::contains(atomsSecondResidue, atomName, atomComparator))
            return LLKA_E_MISSING_ATOMS;
    }

    if (atomsFirstResidue.size() != LLKABones::NUM_FIRST_EXTENDED_BACKONE_ATOMS || atomsSecondResidue.size() != LLKABones::NUM_SECOND_EXTENDED_BACKONE_ATOMS)
        return LLKA_E_MISMATCHING_SIZES;

    // Check that we have just two residues
    if (!LLKAInternal::atomPropertyMatches(atomsSecondResidue, &LLKA_Atom::label_seq_id))
        return LLKA_E_MISMATCHING_DATA;

    int32_t secondSeqId = atomsSecondResidue.front()->label_seq_id;

    LLKA_BaseKind firstBaseKind;
    LLKA_BaseKind secondBaseKind;

    if (LLKA_baseKind(atomsFirstResidue.front()->label_comp_id, &firstBaseKind) != LLKA_OK)
        return LLKA_E_INVALID_ARGUMENT;
    if (LLKA_baseKind(atomsSecondResidue.front()->label_comp_id, &secondBaseKind) != LLKA_OK)
        return LLKA_E_INVALID_ARGUMENT;

    // Check that we have base atoms
    for (const auto &res : { std::make_pair(firstSeqId, boneFirst), std::make_pair(secondSeqId, boneSecond) }) {
        const auto seqId = res.first;
        const auto &baseAtoms = res.second.base;
        for (const auto &atomName : baseAtoms) {
            size_t idx = 0;
            for (; idx < stru->nAtoms; idx++) {
                const auto atom = &stru->atoms[idx];
                if (atom->label_seq_id != seqId)
                    continue;
                if (atomName == atom->label_atom_id)
                    break;
            }
            if (idx == stru->nAtoms)
                return LLKA_E_MISSING_ATOMS; // Atom not found
        }
    }

    if (info != nullptr) {
        info->firstSeqId = firstSeqId;
        info->firstSeqIdAuth = atomsFirstResidue.front()->auth_seq_id;
        info->secondSeqId = secondSeqId;
        info->secondSeqIdAuth = atomsSecondResidue.front()->auth_seq_id;
        info->firstBaseKind = firstBaseKind;
        info->secondBaseKind = secondBaseKind;
    }

    return LLKA_OK;
}
