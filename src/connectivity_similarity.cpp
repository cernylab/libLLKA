// vim: set sw=4 ts=4 sts=4 expandtab :

#include <llka_connectivity_similarity.h>

#include <llka_superposition.h>

#include "ntc_constants.h"
#include "similarity.h"
#include "util/geometry.h"
#include "util/templates.hpp"

#include <cassert>
#include <memory>

namespace LLKAInternal {

static constexpr LLKA_BackboneAtom BKBN_C5_PRIME_FIRST = { LLKA_FIRST_RESIDUE, "C5'" };
static constexpr LLKA_BackboneAtom BKBN_C5_PRIME_SECOND = { LLKA_SECOND_RESIDUE, "C5'" };
static constexpr LLKA_BackboneAtom BKBN_O3_PRIME_FIRST = { LLKA_FIRST_RESIDUE, "O3'" };
static constexpr LLKA_BackboneAtom BKBN_O3_PRIME_SECOND = { LLKA_SECOND_RESIDUE, "O3'" };

static
auto averagesToMetrics(const NtCAverages &avgs) -> LLKA_StepMetrics
{
    return {
        avgs.delta_1,
        avgs.epsilon_1,
        avgs.zeta_1,
        avgs.alpha_2,
        avgs.beta_2,
        avgs.gamma_2,
        avgs.delta_2,
        avgs.chi_1,
        avgs.chi_2,
        avgs.CC,
        avgs.NN,
        avgs.mu
    };
}

template <LLKAStructureType TA, LLKAStructureType TB>
static
auto measureStepConnectivity(const TA *bkbnPositionFirst, LLKA_Structure *bkbnFirst, const TB *bkbnPositionSecond, LLKA_Structure *bkbnSecond, LLKA_Connectivity *result)
{
    LLKA_RetCode tRet;

    double rmsd;
    if constexpr (std::is_same_v<TA, LLKA_Structure>)
        tRet = LLKA_superposeStructures(bkbnFirst, bkbnPositionFirst, &rmsd);
    else
        tRet = LLKA_superposeStructuresView(bkbnFirst, bkbnPositionFirst, &rmsd);

    if (tRet != LLKA_OK)
        return tRet;

    if constexpr (std::is_same_v<TB, LLKA_Structure>)
        tRet = LLKA_superposeStructures(bkbnSecond, bkbnPositionSecond, &rmsd);
    else
        tRet = LLKA_superposeStructuresView(bkbnSecond, bkbnPositionSecond, &rmsd);

    if (tRet != LLKA_OK)
        return tRet;

    // We need *second residue* C5' and O3' from the first step
    const auto C5PrimeFirstIdx = LLKA_backboneAtomIndex(&LLKAInternal::BKBN_C5_PRIME_SECOND, bkbnFirst);
    assert(C5PrimeFirstIdx != LLKA_INVALID_BKBN_ATOM_INDEX);
    const auto O3PrimeFirstIdx = LLKA_backboneAtomIndex(&LLKAInternal::BKBN_O3_PRIME_SECOND, bkbnFirst);
    assert(O3PrimeFirstIdx != LLKA_INVALID_BKBN_ATOM_INDEX);

    // We need *first residue* C5' and O3' from the second step
    const auto C5PrimeSecondIdx = LLKA_backboneAtomIndex(&LLKAInternal::BKBN_C5_PRIME_FIRST, bkbnSecond);
    assert(C5PrimeSecondIdx != LLKA_INVALID_BKBN_ATOM_INDEX);
    const auto O3PrimeSecondIdx = LLKA_backboneAtomIndex(&LLKAInternal::BKBN_O3_PRIME_FIRST, bkbnSecond);
    assert(O3PrimeSecondIdx != LLKA_INVALID_BKBN_ATOM_INDEX);

    result->C5PrimeDistance = LLKAInternal::spatialDistance<double>(bkbnFirst->atoms[C5PrimeFirstIdx], bkbnSecond->atoms[C5PrimeSecondIdx]);
    result->O3PrimeDistance = LLKAInternal::spatialDistance<double>(bkbnFirst->atoms[O3PrimeFirstIdx], bkbnSecond->atoms[O3PrimeSecondIdx]);

    return LLKA_OK;
}

template <LLKAStructureType T>
static
auto measureStepSimilarity(const LLKA_StepMetrics &stepMetrics, const T *bkbnStepStru, const LLKA_Structure *rmsdRefStru, const LLKA_StepMetrics &refMetrics, LLKA_Similarity *result)
{
    LLKA_RetCode tRet;
    LLKA_Structure bkbnRmsdRefStru;

    tRet = LLKA_extractExtendedBackbone(rmsdRefStru, &bkbnRmsdRefStru);
    if (tRet != LLKA_OK)
        return tRet;

    if constexpr (std::is_same_v<T, LLKA_Structure>)
        tRet = LLKA_superposeStructures(&bkbnRmsdRefStru, bkbnStepStru, &result->rmsd);
    else
        tRet = LLKA_superposeStructuresView(&bkbnRmsdRefStru, bkbnStepStru, &result->rmsd);

    if (tRet != LLKA_OK) {
        LLKA_destroyStructure(&bkbnRmsdRefStru);
        return tRet;
    }
    LLKA_destroyStructure(&bkbnRmsdRefStru);

    result->euclideanDistance = 0;
    // Dinucleotide torsions
    for (
        const auto &clsPtr :
        { &LLKA_StepMetrics::delta_1, &LLKA_StepMetrics::epsilon_1, &LLKA_StepMetrics::zeta_1, &LLKA_StepMetrics::alpha_2, &LLKA_StepMetrics::beta_2,
          &LLKA_StepMetrics::gamma_2, &LLKA_StepMetrics::delta_2, &LLKA_StepMetrics::chi_1, &LLKA_StepMetrics::chi_2 }
    ) {
        auto angDiff = LLKAInternal::R2D(LLKAInternal::angleDifference(stepMetrics.*clsPtr, refMetrics.*clsPtr));
        result->euclideanDistance += angDiff * angDiff;
    }
    // Cross-residue torsion
    auto aux = LLKAInternal::R2D(LLKAInternal::angleDifference(stepMetrics.mu, refMetrics.mu));
    result->euclideanDistance += aux * aux;

    // Cross-residue distances
    aux = XR_DISTANCE_MULTIPLIER * std::abs(stepMetrics.CC - refMetrics.CC);
    result->euclideanDistance += aux * aux;

    aux = XR_DISTANCE_MULTIPLIER * std::abs(stepMetrics.NN - refMetrics.NN);
    result->euclideanDistance += aux * aux;

    result->euclideanDistance = std::sqrt(result->euclideanDistance);

    return LLKA_OK;
}

} // namespace LLKAInternal

LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityNtCs(const LLKA_Structure *positionFirst, LLKA_NtC ntcFirst, const LLKA_Structure *positionSecond, LLKA_NtC ntcSecond, LLKA_Connectivity *result)
{
    return LLKA_measureStepConnectivityStructures(positionFirst, &LLKAInternal::NTC_REFERENCES[ntcFirst], positionSecond, &LLKAInternal::NTC_REFERENCES[ntcSecond], result);
}

LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityNtCsMultipleFirst(const LLKA_Structure *positionFirst, const LLKA_NtC ntcsFirst[], const LLKA_Structure *positionSecond, LLKA_NtC ntcSecond, LLKA_Connectivities *results)
{
    size_t ntcCount = 0;
    while (ntcsFirst[ntcCount] != LLKA_INVALID_NTC) ntcCount++;

    if (ntcCount != results->nConns)
        return LLKA_E_MISMATCHING_SIZES;

    LLKA_RetCode tRet;
    LLKA_StructureView bkbnPositionFirst;
    LLKA_StructureView bkbnPositionSecond;
    LLKA_Structure bkbnSecond;

    tRet = LLKA_extractExtendedBackboneView(positionFirst, &bkbnPositionFirst);
    if (tRet != LLKA_OK)
        return tRet;

    tRet = LLKA_extractExtendedBackboneView(positionSecond, &bkbnPositionSecond);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&bkbnPositionFirst);
        return tRet;
    }

    tRet = LLKA_extractExtendedBackbone(&LLKAInternal::NTC_REFERENCES[ntcSecond], &bkbnSecond);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&bkbnPositionFirst);
        LLKA_destroyStructureView(&bkbnPositionSecond);
        return tRet;
    }

    for (size_t idx = 0; idx < ntcCount; idx++) {
        LLKA_Structure bkbnFirst;

        tRet = LLKA_extractExtendedBackbone(&LLKAInternal::NTC_REFERENCES[ntcsFirst[idx]], &bkbnFirst);
        if (tRet != LLKA_OK)
            break;

        tRet = LLKAInternal::measureStepConnectivity(&bkbnPositionFirst, &bkbnFirst, &bkbnPositionSecond, &bkbnSecond, &results->conns[idx]);
        LLKA_destroyStructure(&bkbnFirst);

        if (tRet != LLKA_OK)
            break;
    }

    LLKA_destroyStructureView(&bkbnPositionSecond);
    LLKA_destroyStructureView(&bkbnPositionFirst);
    LLKA_destroyStructure(&bkbnSecond);

    return tRet;
}

LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityNtCsMultipleSecond(const LLKA_Structure *positionFirst, LLKA_NtC ntcFirst, const LLKA_Structure *positionSecond, const LLKA_NtC ntcsSecond[], LLKA_Connectivities *results)
{
    size_t ntcCount = 0;
    while (ntcsSecond[ntcCount] != LLKA_INVALID_NTC) ntcCount++;

    if (ntcCount != results->nConns)
        return LLKA_E_MISMATCHING_SIZES;

    LLKA_RetCode tRet;
    LLKA_Structure bkbnPositionFirst;
    LLKA_Structure bkbnFirst;
    LLKA_Structure bkbnPositionSecond;

    tRet = LLKA_extractExtendedBackbone(positionFirst, &bkbnPositionFirst);
    if (tRet != LLKA_OK)
        return tRet;

    tRet = LLKA_extractExtendedBackbone(&LLKAInternal::NTC_REFERENCES[ntcFirst], &bkbnFirst);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructure(&bkbnPositionFirst);
        return tRet;
    }

    tRet = LLKA_extractExtendedBackbone(positionSecond, &bkbnPositionSecond);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructure(&bkbnFirst);
        LLKA_destroyStructure(&bkbnPositionFirst);
        return tRet;
    }

    for (size_t idx = 0; idx < ntcCount; idx++) {
        LLKA_Structure bkbnSecond;

        tRet = LLKA_extractExtendedBackbone(&LLKAInternal::NTC_REFERENCES[ntcsSecond[idx]], &bkbnSecond);
        if (tRet != LLKA_OK)
            break;

        tRet = LLKAInternal::measureStepConnectivity(&bkbnPositionFirst, &bkbnFirst, &bkbnPositionSecond, &bkbnSecond, &results->conns[idx]);
        LLKA_destroyStructure(&bkbnSecond);

        if (tRet != LLKA_OK)
            break;
    }

    LLKA_destroyStructure(&bkbnPositionSecond);
    LLKA_destroyStructure(&bkbnFirst);
    LLKA_destroyStructure(&bkbnPositionFirst);

    return tRet;
}

LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityStructures(const LLKA_Structure *positionFirst, const LLKA_Structure *dinuFirst, const LLKA_Structure *positionSecond, const LLKA_Structure *dinuSecond, LLKA_Connectivity *result)
{
    LLKA_RetCode tRet;
    LLKA_StructureView bkbnPositionFirst;
    LLKA_Structure bkbnFirst;
    LLKA_StructureView bkbnPositionSecond;
    LLKA_Structure bkbnSecond;

    tRet = LLKA_extractExtendedBackboneView(positionFirst, &bkbnPositionFirst);
    if (tRet != LLKA_OK)
        return tRet;

    tRet = LLKA_extractExtendedBackbone(dinuFirst, &bkbnFirst);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&bkbnPositionFirst);
        return tRet;
    }

    tRet = LLKA_extractExtendedBackboneView(positionSecond, &bkbnPositionSecond);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructure(&bkbnFirst);
        LLKA_destroyStructureView(&bkbnPositionFirst);
        return tRet;
    }

    tRet = LLKA_extractExtendedBackbone(dinuSecond, &bkbnSecond);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&bkbnPositionSecond);
        LLKA_destroyStructure(&bkbnFirst);
        LLKA_destroyStructureView(&bkbnPositionFirst);
        return tRet;
    }

    tRet = LLKAInternal::measureStepConnectivity(&bkbnPositionFirst, &bkbnFirst, &bkbnPositionSecond, &bkbnSecond, result);

    LLKA_destroyStructure(&bkbnSecond);
    LLKA_destroyStructureView(&bkbnPositionSecond);
    LLKA_destroyStructure(&bkbnFirst);
    LLKA_destroyStructureView(&bkbnPositionFirst);

    return tRet;
}

// TODO: We need functions that take multiple references to avoid the function call overhead

LLKA_API LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityStructuresMultiple(const LLKA_Structure *positionFirst, const LLKA_Structure *dinuFirst, const LLKA_Structure *positionSecond, const LLKA_Structures *dinusSecond, LLKA_Connectivities *results)
{
    if (dinusSecond->nStrus != results->nConns)
        return LLKA_E_MISMATCHING_SIZES;

    const auto nStrus = dinusSecond->nStrus;

    LLKA_RetCode tRet;
    LLKA_StructureView bkbnPositionFirst;
    LLKA_Structure bkbnFirst;
    LLKA_StructureView bkbnPositionSecond;

    tRet = LLKA_extractExtendedBackboneView(positionFirst, &bkbnPositionFirst);
    if (tRet != LLKA_OK)
        return tRet;

    tRet = LLKA_extractExtendedBackbone(dinuFirst, &bkbnFirst);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&bkbnPositionFirst);
        return tRet;
    }

    tRet = LLKA_extractExtendedBackboneView(positionSecond, &bkbnPositionSecond);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructure(&bkbnFirst);
        LLKA_destroyStructureView(&bkbnPositionFirst);
        return tRet;
    }

    for (size_t idx = 0; idx < nStrus; idx++) {
        LLKA_Structure bkbnSecond;

        tRet = LLKA_extractExtendedBackbone(&dinusSecond->strus[idx], &bkbnSecond);
        if (tRet != LLKA_OK)
            break;

        tRet = LLKAInternal::measureStepConnectivity(&bkbnPositionFirst, &bkbnFirst, &bkbnPositionSecond, &bkbnSecond, &results->conns[idx]);
        LLKA_destroyStructure(&bkbnSecond);

        if (tRet != LLKA_OK)
            break;
    }

    LLKA_destroyStructureView(&bkbnPositionSecond);
    LLKA_destroyStructure(&bkbnFirst);
    LLKA_destroyStructureView(&bkbnPositionFirst);

    return tRet;
}

LLKA_RetCode LLKA_CC LLKA_measureStepSimilarityNtC(const LLKA_Structure *stepStru, LLKA_NtC ntc, LLKA_Similarity *result)
{
    LLKA_StepMetrics stepMetrics;
    LLKA_StructureView bkbnStepStru;

    auto tRet = LLKA_calculateStepMetrics(stepStru, &stepMetrics);
    if (tRet != LLKA_OK)
        return tRet;

    tRet = LLKA_extractExtendedBackboneView(stepStru, &bkbnStepStru);
    if (tRet != LLKA_OK)
        return tRet;

    const auto &refStru = LLKAInternal::NTC_REFERENCES[ntc];
    const auto &refMetrics = averagesToMetrics(LLKAInternal::NTC_AVERAGES[ntc]);

    tRet = LLKAInternal::measureStepSimilarity(stepMetrics, &bkbnStepStru, &refStru, refMetrics, result);
    LLKA_destroyStructureView(&bkbnStepStru);

    return tRet;
}

LLKA_RetCode LLKA_CC LLKA_measureStepSimilarityNtCMultiple(const LLKA_Structure *stepStru, const LLKA_NtC ntcs[], LLKA_Similarities *results)
{
    LLKA_StepMetrics stepMetrics;
    LLKA_StructureView bkbnStepStru;

    auto tRet = LLKA_calculateStepMetrics(stepStru, &stepMetrics);
    if (tRet != LLKA_OK)
        return tRet;

    tRet = LLKA_extractExtendedBackboneView(stepStru, &bkbnStepStru);
    if (tRet != LLKA_OK)
        return tRet;

    size_t idx = 0;
    while (ntcs[idx] != LLKA_INVALID_NTC) {
        if (idx >= results->nSimilars)
            return LLKA_E_MISMATCHING_SIZES;

        const auto ntc = ntcs[idx];
        const auto &refStru = LLKAInternal::NTC_REFERENCES[ntc];
        const auto &refMetrics = averagesToMetrics(LLKAInternal::NTC_AVERAGES[ntc]);

        tRet = LLKAInternal::measureStepSimilarity(stepMetrics, &bkbnStepStru, &refStru, refMetrics, &results->similars[idx]);
        if (tRet != LLKA_OK)
            goto out;

        idx++;
    }

out:
    LLKA_destroyStructureView(&bkbnStepStru);

    return tRet;
}

LLKA_RetCode LLKA_CC LLKA_measureStepSimilarityStructure(const LLKA_Structure *stepStru, const LLKA_Structure *refStru, LLKA_Similarity *result)
{
    LLKA_StepMetrics stepMetrics;
    LLKA_StepMetrics refMetrics;
    LLKA_StructureView bkbnStepStru;

    auto tRet = LLKA_calculateStepMetrics(stepStru, &stepMetrics);
    if (tRet != LLKA_OK)
        return tRet;

    tRet = LLKA_calculateStepMetrics(refStru, &refMetrics);
    if (tRet != LLKA_OK)
        return tRet;

    tRet = LLKA_extractExtendedBackboneView(stepStru, &bkbnStepStru);
    if (tRet != LLKA_OK)
        return tRet;

    tRet = LLKAInternal::measureStepSimilarity(stepMetrics, &bkbnStepStru, refStru, refMetrics, result);
    LLKA_destroyStructureView(&bkbnStepStru);

    return tRet;
}
