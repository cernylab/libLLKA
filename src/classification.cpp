/* vim: set sw=4 ts=4 sts=4 expandtab : */


#include <llka_classification.h>
#include <llka_nucleotide.h>
#include <llka_superposition.h>

#include "util/arch.hpp"
#include "util/elementaries.h"
#include "util/geometry.h"
#include "util/printers.hpp"


#include "tracing/llka_tracer.h"


#include <array>
#include <cassert>
#include <map>
#include <memory>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>


#include "nucleotide.hpp"
#include "ntc.hpp"
#include "ntc_references.h"
#include "similarity.h"


#define CLASSIFICATION_VIOLATION_STR(x) case x: return #x

struct LLKA_ClassificationContext {
    // TODO: This is a very interim solution.
    // This will be reworked once we get on with SIMD

    ~LLKA_ClassificationContext()
    {
        for (auto &gs : goldenSteps)
            LLKAInternal::destroyString(gs.name);

        for (auto &p : ntcExtendedBackbonePoints)
            delete [] p.points;
    }

    std::vector<LLKA_GoldenStep> goldenSteps{};
    std::vector<LLKA_ClassificationCluster> clusters{};
    std::vector<LLKA_Confal> confals{};
    std::vector<double> confalPercentiles{};
    LLKA_ClassificationLimits limits;

    double maxCloseEnoughRmsd;

    std::array<LLKA_Points, 96> ntcExtendedBackbonePoints{};
};

namespace LLKAInternal {

inline constinit double MINIMUM_ALLOWED_DELTA = D2R(55.0);
inline constinit double MAXIMUM_ALLOWED_DELTA = D2R(185.0);

inline constinit std::array<double LLKA_StepMetrics::*, 9> ALL_TORSIONS_STEP_METRIC_CLSPTRS{
    &LLKA_StepMetrics::delta_1, &LLKA_StepMetrics::epsilon_1, &LLKA_StepMetrics::zeta_1, &LLKA_StepMetrics::alpha_2, &LLKA_StepMetrics::beta_2,
    &LLKA_StepMetrics::gamma_2, &LLKA_StepMetrics::delta_2, &LLKA_StepMetrics::chi_1, &LLKA_StepMetrics::chi_2
};
inline constinit std::array<LLKA_ClassificationMetric LLKA_ClassificationCluster::*, 9> ALL_TORSIONS_CLASSIFICATION_METRIC_CLSPTRS{
    &LLKA_ClassificationCluster::delta_1, &LLKA_ClassificationCluster::epsilon_1, &LLKA_ClassificationCluster::zeta_1, &LLKA_ClassificationCluster::alpha_2, &LLKA_ClassificationCluster::beta_2,
    &LLKA_ClassificationCluster::gamma_2, &LLKA_ClassificationCluster::delta_2, &LLKA_ClassificationCluster::chi_1, &LLKA_ClassificationCluster::chi_2
};
inline constinit std::array<double LLKA_Confal::*, 9> ALL_TORSIONS_CONFAL_CLSPTRS{
    &LLKA_Confal::delta_1, &LLKA_Confal::epsilon_1, &LLKA_Confal::zeta_1, &LLKA_Confal::alpha_2, &LLKA_Confal::beta_2,
    &LLKA_Confal::gamma_2, &LLKA_Confal::delta_2, &LLKA_Confal::chi_1, &LLKA_Confal::chi_2
};
inline constinit std::array<double LLKA_ConfalScore::*, 9> ALL_TORSIONS_CONFAL_SCORE_CLSPTRS{
    &LLKA_ConfalScore::delta_1, &LLKA_ConfalScore::epsilon_1, &LLKA_ConfalScore::zeta_1, &LLKA_ConfalScore::alpha_2, &LLKA_ConfalScore::beta_2,
    &LLKA_ConfalScore::gamma_2, &LLKA_ConfalScore::delta_2, &LLKA_ConfalScore::chi_1, &LLKA_ConfalScore::chi_2
};
static_assert(ALL_TORSIONS_STEP_METRIC_CLSPTRS.size() == ALL_TORSIONS_CLASSIFICATION_METRIC_CLSPTRS.size());
static_assert(ALL_TORSIONS_STEP_METRIC_CLSPTRS.size() == ALL_TORSIONS_CONFAL_CLSPTRS.size());
static_assert(ALL_TORSIONS_STEP_METRIC_CLSPTRS.size() == ALL_TORSIONS_CONFAL_SCORE_CLSPTRS.size());

inline constinit std::array<double LLKA_StepMetrics::*, 10> TORSION_STEP_METRIC_CLSPTRS{
    &LLKA_StepMetrics::delta_1, &LLKA_StepMetrics::epsilon_1, &LLKA_StepMetrics::zeta_1, &LLKA_StepMetrics::alpha_2, &LLKA_StepMetrics::beta_2,
    &LLKA_StepMetrics::gamma_2, &LLKA_StepMetrics::delta_2, &LLKA_StepMetrics::chi_1, &LLKA_StepMetrics::chi_2,
    &LLKA_StepMetrics::mu
};
inline constinit std::array<LLKA_ClassificationMetric LLKA_ClassificationCluster::*, 10> TORSION_CLASSIFICATION_METRIC_CLSPTRS{
    &LLKA_ClassificationCluster::delta_1, &LLKA_ClassificationCluster::epsilon_1, &LLKA_ClassificationCluster::zeta_1, &LLKA_ClassificationCluster::alpha_2, &LLKA_ClassificationCluster::beta_2,
    &LLKA_ClassificationCluster::gamma_2, &LLKA_ClassificationCluster::delta_2, &LLKA_ClassificationCluster::chi_1, &LLKA_ClassificationCluster::chi_2,
    &LLKA_ClassificationCluster::mu
};
static_assert(TORSION_STEP_METRIC_CLSPTRS.size() == TORSION_CLASSIFICATION_METRIC_CLSPTRS.size());

inline constinit std::array<LLKA_ClassificationMetric LLKA_NuAnglesMetrics::*, 5> NU_ANGLES_METRICS_CLSPTRS{
    &LLKA_NuAnglesMetrics::nu_0, &LLKA_NuAnglesMetrics::nu_1, &LLKA_NuAnglesMetrics::nu_2, &LLKA_NuAnglesMetrics::nu_3, &LLKA_NuAnglesMetrics::nu_4
};
inline constinit std::array<double LLKA_NuAngles::*, 5> NU_ANGLES_CLSPTRS{
    &LLKA_NuAngles::nu_0, &LLKA_NuAngles::nu_1, &LLKA_NuAngles::nu_2, &LLKA_NuAngles::nu_3, &LLKA_NuAngles::nu_4
};
static_assert(NU_ANGLES_METRICS_CLSPTRS.size() == NU_ANGLES_CLSPTRS.size());

class NearestNeighbor {
public:
    NearestNeighbor() :
        euclideanDistance{-1},
        goldenStepIdx{Arch::INVALID_SIZE_T}
    {}

    auto isValid() { return euclideanDistance >= 0.0 && goldenStepIdx != Arch::INVALID_SIZE_T; }

    LLKA_StepMetrics metricsDifference;
    double euclideanDistance;
    size_t goldenStepIdx;
};

using NearestNeighborsView = std::span<const LLKAInternal::NearestNeighbor, Arch::ArraySize::MAX>;

static
auto calcConfalScore(const LLKA_StepMetrics &differencesFromNtCAverages, const LLKA_Confal &confal, bool noViolations)
{
    constexpr auto calcScore = [](auto diff, auto confal) {
        auto diffSq = diff * diff;
        auto confalSq = confal * confal; // PERF: Precalculate in init

        // PERF: Look into replacing division by multiplication
        return 100.0 * std::exp(-diffSq / (2.0 * confalSq));
    };

    double invTotal = 0;
    LLKA_ConfalScore score;

    for (size_t idx = 0; idx < ALL_TORSIONS_STEP_METRIC_CLSPTRS.size(); idx++) {
        auto metricClsPtr = ALL_TORSIONS_STEP_METRIC_CLSPTRS[idx];
        auto confalClsPtr = ALL_TORSIONS_CONFAL_CLSPTRS[idx];
        auto confalScoreClsPtr = ALL_TORSIONS_CONFAL_SCORE_CLSPTRS[idx];


        auto s = calcScore(R2D(differencesFromNtCAverages.*metricClsPtr), confal.*confalClsPtr);
        invTotal += 1.0 / s;

        score.*confalScoreClsPtr = s;
    }

    auto calcAndAccum = [&invTotal, &calcScore](auto diff, auto confal) {
        auto s = calcScore(diff, confal);
        invTotal += 1.0 / s;

        return s;
    };

    score.CC = calcAndAccum(differencesFromNtCAverages.CC, confal.CC);
    score.NN = calcAndAccum(differencesFromNtCAverages.NN, confal.NN);
    score.mu = calcAndAccum(R2D(differencesFromNtCAverages.mu), confal.mu);

    score.total = ((12.0 / invTotal) + 0.5) * decltype(score.total)(noViolations);

    return score;
}

static
auto calcDistancesFromNtCAverages(const LLKA_StepMetrics &stepMetrics, const LLKA_ClassificationCluster &cluster)
{
    LLKA_StepMetrics dists{};
    double totalDiffSq = 0;

    for (size_t idx = 0; idx < ALL_TORSIONS_STEP_METRIC_CLSPTRS.size(); idx++) {
        auto stepClsPtr = ALL_TORSIONS_STEP_METRIC_CLSPTRS[idx];
        auto clustClsPtr = ALL_TORSIONS_CLASSIFICATION_METRIC_CLSPTRS[idx];

        auto angDiff = angleDifference(stepMetrics.*stepClsPtr, (cluster.*clustClsPtr).meanValue);

        dists.*stepClsPtr = angDiff;
        totalDiffSq += R2D(angDiff) * R2D(angDiff);
    }
    dists.CC = stepMetrics.CC - cluster.CC.meanValue;
    totalDiffSq += XR_DISTANCE_MULTIPLIER * XR_DISTANCE_MULTIPLIER * (dists.CC * dists.CC);

    dists.NN = stepMetrics.NN - cluster.NN.meanValue;
    totalDiffSq += XR_DISTANCE_MULTIPLIER * XR_DISTANCE_MULTIPLIER * (dists.NN * dists.NN);

    dists.mu = angleDifference(stepMetrics.mu, cluster.mu.meanValue);
    totalDiffSq += R2D(dists.mu) * R2D(dists.mu);

    return std::make_tuple(dists, std::sqrt(totalDiffSq));
}

static
auto calcRmsdToClosestNtC(const LLKA_Structure *stru, const LLKA_Points &ntcExtBkbnPts)
{
    // We expect that no structure will ever have extended backbone more that 50 atoms large
    LLKA_Point scratchArray[50];
    LLKA_Points scratchPts {
        { scratchArray },
        50
    };

    LLKA_Structure extBkbn{};

    // PERF: This code right here indicates that we should split atoms and coordinates
    // to two arrays
    LLKA_extractExtendedBackbone(stru, &extBkbn);
    for (size_t idx = 0; idx < extBkbn.nAtoms; idx++)
        scratchPts.points[idx] = extBkbn.atoms[idx].coords;

    LLKA_Points pts;
    pts.raw = scratchPts.raw;
    pts.nPoints = extBkbn.nAtoms;

    double rmsd = 0;
    auto tRet = LLKA_superposePoints(&pts, &ntcExtBkbnPts, &rmsd);
    assert(tRet == LLKA_OK);
#ifdef NDEBUG
    (void)tRet;
#endif // NDEBUG

    LLKA_destroyStructure(&extBkbn);

    return rmsd;
}

static
auto checkNtCTolerances(
    const LLKA_StepMetrics &stepMetrics,
    const LLKA_ClassificationCluster &bestieCluster,
    const NearestNeighborsView &nearestNeighbors,
    const double pseudorotation_1,
    const double pseudorotation_2,
    const LLKA_ClassificationContext *ctx
) -> std::tuple<int32_t, int16_t, int16_t>
{
    assert(nearestNeighbors.size() > 0);

    // PERF: This will really benefit of SIMD once we get around to doing it

    LLKA_StepMetrics sumsOfSinesOfTorsions = {};
    LLKA_StepMetrics sumsOfCosinesOfTorsions = {};
    LLKA_StepMetrics averageNeighborTorsions = {};

    int32_t violations = 0;
    int16_t violatingTorsionsAverage = 0;
    int16_t violatingTorsionsNearest = 0;

    // We rely on this to be 9
    assert(ALL_TORSIONS_STEP_METRIC_CLSPTRS.size() == 9);

    for (const auto &neigh : nearestNeighbors) {
        const auto &gs = ctx->goldenSteps[neigh.goldenStepIdx];

        for (const auto &clsPtr : ALL_TORSIONS_STEP_METRIC_CLSPTRS) {
            sumsOfSinesOfTorsions.*clsPtr += std::sin(angleAsFull(gs.metrics.*clsPtr));
            sumsOfCosinesOfTorsions.*clsPtr += std::cos(angleAsFull(gs.metrics.*clsPtr));
        }
    }

    for (const auto &clsPtr : ALL_TORSIONS_STEP_METRIC_CLSPTRS) {
        auto ang = angleAsFull(std::atan2(sumsOfSinesOfTorsions.*clsPtr, sumsOfCosinesOfTorsions.*clsPtr));
        averageNeighborTorsions.*clsPtr = ang;
    }

    for (size_t idx = 0; idx < ALL_TORSIONS_STEP_METRIC_CLSPTRS.size(); idx++) {
        auto clsPtr = ALL_TORSIONS_STEP_METRIC_CLSPTRS[idx];
        auto angDiff = angleDifference(stepMetrics.*clsPtr, averageNeighborTorsions.*clsPtr);
        if (std::abs(angDiff) > ctx->limits.averageNeighborsTorsionCutoff) {
            violations |= LLKA_CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT;
            // This relies in BACKBONE_TORSIONS_STEP_METRIC_CLSPTRS and the errors enum to be in the same order
            violatingTorsionsAverage |= (1 << idx);
        }
    }

    const auto &bestieGoldenStep = ctx->goldenSteps[nearestNeighbors[0].goldenStepIdx];
    for (size_t idx = 0; idx < ALL_TORSIONS_STEP_METRIC_CLSPTRS.size(); idx++) {
        auto clsPtr = ALL_TORSIONS_STEP_METRIC_CLSPTRS[idx];
        auto angDiff = angleDifference(stepMetrics.*clsPtr, bestieGoldenStep.metrics.*clsPtr);
        if (std::abs(angDiff) > ctx->limits.nearestNeighborTorsionsCutoff) {
            violations |= LLKA_CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT;
            // This relies in BACKBONE_TORSIONS_STEP_METRIC_CLSPTRS and the errors enum to be in the same order
            violatingTorsionsNearest |= (1 << idx);
        }
    }

    if (violations)
        return { violations, violatingTorsionsAverage, violatingTorsionsNearest };

    LLKA_StepMetrics differencesFromClusterAverage{};

    for (size_t idx = 0; idx < ALL_TORSIONS_STEP_METRIC_CLSPTRS.size(); idx++) {
        const auto stepMetricsClsPtr = ALL_TORSIONS_STEP_METRIC_CLSPTRS[idx];
        const auto clusterMetricsClsPtr = ALL_TORSIONS_CLASSIFICATION_METRIC_CLSPTRS[idx];

        // Mean values in cluster metrics are in 0 -> 2PI range.
        differencesFromClusterAverage.*stepMetricsClsPtr = angleDifference(angleAsFull(stepMetrics.*stepMetricsClsPtr), (bestieCluster.*clusterMetricsClsPtr).meanValue);
    }
    differencesFromClusterAverage.mu = angleDifference(angleAsFull(stepMetrics.mu), bestieCluster.mu.meanValue);

    if (stepMetrics.CC < bestieCluster.CC.minValue)
        violations |= LLKA_CLASSIFICATION_E_CC_TOO_LOW;
    else if (stepMetrics.CC > bestieCluster.CC.maxValue)
        violations |= LLKA_CLASSIFICATION_E_CC_TOO_HIGH;

    if (stepMetrics.NN < bestieCluster.NN.minValue)
        violations |= LLKA_CLASSIFICATION_E_NN_TOO_LOW;
    else if (stepMetrics.NN > bestieCluster.NN.maxValue)
        violations |= LLKA_CLASSIFICATION_E_NN_TOO_HIGH;

    if (angleDifference(stepMetrics.mu, bestieCluster.mu.minValue) < 0.0)
        violations |= LLKA_CLASSIFICATION_E_MU_TOO_LOW;
    else if (angleDifference(stepMetrics.mu, bestieCluster.mu.maxValue) > 0.0)
        violations |= LLKA_CLASSIFICATION_E_MU_TOO_HIGH;

    // NOTE: Original implementation does some additional "golden step mode" checks here.

    double totalDistance = 0;
    for (size_t idx = 0; idx < 7; idx++) { // Use only first 7 backbone torsions
        const auto stepMetricsClsPtr = ALL_TORSIONS_STEP_METRIC_CLSPTRS[idx];
        const auto clusterMetricsClsPtr = ALL_TORSIONS_CLASSIFICATION_METRIC_CLSPTRS[idx];

        // Simple way to calculate angle difference for 0 -> 2PI ranges
        totalDistance += angleAsFull(stepMetrics.*stepMetricsClsPtr) - (bestieCluster.*clusterMetricsClsPtr).meanValue;
    }
    if (std::abs(totalDistance) > ctx->limits.totalDistanceCutoff)
        violations |= LLKA_CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH;

    auto pseudo1Diff = std::abs(angleDifference(pseudorotation_1, bestieCluster.ribosePseudorotation_1));
    auto pseudo2Diff = std::abs(angleDifference(pseudorotation_2, bestieCluster.ribosePseudorotation_2));
    if (pseudo1Diff > ctx->limits.pseudorotationCutoff)
        violations |= LLKA_CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT;
    if (pseudo2Diff > ctx->limits.pseudorotationCutoff)
        violations |= LLKA_CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT;

    if (violations & LLKA_CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT || violations & LLKA_CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT) {
        ECHMET_TRACE(
            LLKATracing, PSEUDOROTATION_TOO_DIFFERENT,
            pseudorotation_1, bestieCluster.ribosePseudorotation_1, pseudo1Diff,
            pseudorotation_2, bestieCluster.ribosePseudorotation_2, pseudo2Diff,
            ctx->limits.pseudorotationCutoff
        );
    }

    return { violations, violatingTorsionsAverage, violatingTorsionsNearest };
}

static
auto determineBestieClusterIdx(const NearestNeighborsView &nearestNeighbors, const LLKA_ClassificationContext *ctx) -> std::tuple<size_t, double>
{
    std::map<size_t, double> clusterVotes{}; // Key is clusterIdx, value the number of votes for each cluster

    if (nearestNeighbors.size() == 0)
        return { 0, -1 };

    for (const auto &nn : nearestNeighbors) {
        assert(nn.goldenStepIdx != Arch::INVALID_SIZE_T);

        const auto clusterIdx = ctx->goldenSteps[nn.goldenStepIdx].clusterIdx;

        // PERF: This is calculated in degrees instead of radians to match the scoring
        // behavior of the reference implementation
        double totalDiff = 0;
        for (const auto &clsPtr : ALL_TORSIONS_STEP_METRIC_CLSPTRS) {
            auto angDiffDeg = R2D(nn.metricsDifference.*clsPtr);
            totalDiff += angDiffDeg * angDiffDeg;
        }
        totalDiff += XR_DISTANCE_MULTIPLIER * XR_DISTANCE_MULTIPLIER * nn.metricsDifference.CC * nn.metricsDifference.CC;
        totalDiff += XR_DISTANCE_MULTIPLIER * XR_DISTANCE_MULTIPLIER * nn.metricsDifference.NN * nn.metricsDifference.NN;
        totalDiff += R2D(nn.metricsDifference.mu) * R2D(nn.metricsDifference.mu);

        auto score = 1.0 / totalDiff;
        clusterVotes[clusterIdx] += score;
    }

    // TODO: Rework
    // Clusters sorted by the highest score
    std::vector<std::pair<size_t, double>> sorted{}; // first: clusterNumber, second: score
    for (const auto &[clusterNumber, votes] : clusterVotes)
        sorted.emplace_back(clusterNumber, votes);
    std::sort(sorted.begin(), sorted.end(), [](const auto &lhs, const auto &rhs) { return rhs.second < lhs.second; });

    return { sorted[0].first, sorted[0].second };
}

/*
 * Inserts \p NearestNeighbor into a vector of NearestNeighbors
 * in a way that the vector is always sorted by nearestDistance.
 *
 * We expect the vector to be allocated in advance with a known size.
 * Therefore, in may contain invalid (empty) values. Since we always keep
 * the vector sorted, we only need to look at the last valid element.
 *
 * NOTE: This function is defined out-of-order to avoid forward declarations.
 */
static
auto insertNearestNeighbor(const NearestNeighbor &nn, std::vector<NearestNeighbor> &vec, const size_t nValidNeighbors) -> size_t
{
    assert(vec.size() >= nValidNeighbors);

    if (nValidNeighbors == 0) {
        vec[0] = nn;
        return 1UL;
    }

    const size_t lastValidIdx = nValidNeighbors - 1;

    size_t idx = 0;
    for (; idx <= lastValidIdx; idx++) {
        if (nn.euclideanDistance < vec[idx].euclideanDistance)
            break;
    }

    if (idx > lastValidIdx) {
        if (idx < vec.size()) {
            // Append
            vec[idx] = nn;
            return nValidNeighbors + 1;
        } else {
            // Vector is full and this item has greater nearestNeighborDistance than all elements it contains so just disregard it
            return nValidNeighbors;
        }
    } else {
        // Insert-sort
        const bool vecIsFull = nValidNeighbors == vec.size();
        // Make room to insert our item by moving things up
        for (size_t jdx = (vecIsFull ? lastValidIdx : lastValidIdx + 1); jdx > idx; jdx--)
            vec[jdx] = vec[jdx - 1];
        vec[idx] = nn;

        return nValidNeighbors + !vecIsFull;
    }
}

static
auto findClosestNtC(const LLKA_StepMetrics &stepMetrics, const LLKA_ClassificationContext *ctx)
{
    bool rejectDelta = false;
    if (
        !LLKA_WITHIN_EXCLUSIVE(MINIMUM_ALLOWED_DELTA, angleAsFull(stepMetrics.delta_1), MAXIMUM_ALLOWED_DELTA) ||
        !LLKA_WITHIN_EXCLUSIVE(MINIMUM_ALLOWED_DELTA, angleAsFull(stepMetrics.delta_2), MAXIMUM_ALLOWED_DELTA)
    )
        rejectDelta = true;

    // PERF: This is sadly slow. Rewrite this to use SIMD once we have the basic implementation working and tested

    std::vector<NearestNeighbor> nearestNeighbors(ctx->limits.numberOfUsedNearestNeighbors);
    size_t nValidNearestNeighbors = 0;
    double shortestEuclideanDistance = std::numeric_limits<double>::max();
    size_t closestGoldenStepIdx = Arch::INVALID_SIZE_T;
    int32_t lastClusterNumber = -1;
    bool rejectCluster = rejectDelta; // Do not calculate any metrics because we will reject the step
                                      // anyway if the deltas are outside tolerances
    NearestNeighbor emergencyNearestNeighbor{}; // Used as the nearest neighbor of there are no candidates that meet all
                                                // of the matching criteria. This emergency neighbor will have the lowest euclidean distance
                                                // to the measured step.
    for (size_t gsIdx = 0; gsIdx < ctx->goldenSteps.size(); gsIdx++) {
        const auto &gs = ctx->goldenSteps[gsIdx];
        const auto &gsMetrics = gs.metrics;

        NearestNeighbor nearestNeighbor{};

        double totalDiffSquared = 0;
        for (const auto &clsPtr : ALL_TORSIONS_STEP_METRIC_CLSPTRS) {
            auto angDiff = angleDifference(stepMetrics.*clsPtr, gsMetrics.*clsPtr);
            nearestNeighbor.metricsDifference.*clsPtr = angDiff;
            totalDiffSquared += angDiff * angDiff;
        }

        // Original implementation tries to detect if the classified step is actually a golden step.
        // It does that by calculating the metrics for the first 9 torsions and compares against zero.
        // The original implementation compares euclidan distance, the notManhattan distance and estimated st.dev.
        // It does not account for floating point arithmetics rounding which makes the original implementation
        // suspicious.
        // Here, we apply the following logic:
        // - We calculate the squared angular difference for the first 9 torsions.
        // - We assume the standard mmCIF coordinate precision of 3 decimal places.
        // - We set the tolerance to (0.0005 * 9)^2 = 0.00002025.
        // - If we fall within the tolerance, we consider the classified step to be
        //   the same as the golden step and skip it.
        if (compareWithTolerance(totalDiffSquared, 0.0, 0.00002025))
            continue;

        // If this fails, we must have screwed up at the sanity check during context initialization
        assert(gs.clusterIdx != Arch::INVALID_SIZE_T);
        const auto &cluster = ctx->clusters[gs.clusterIdx];
        assert(cluster.number == gs.clusterNumber);

        const auto ccDiff = stepMetrics.CC - gsMetrics.CC;
        const auto nnDiff = stepMetrics.NN - gsMetrics.NN;
        const auto muDiff = angleDifference(stepMetrics.mu, gsMetrics.mu);

        // PERF: This is sort of excessive but we will not keep this code anyway
        nearestNeighbor.metricsDifference.CC = ccDiff;
        nearestNeighbor.metricsDifference.NN = nnDiff;
        nearestNeighbor.metricsDifference.mu = muDiff;

        ECHMET_TRACE(LLKATracing, CLASSIFICATION_METRICS_DIFFERENCES, stepMetrics, gs, nearestNeighbor.metricsDifference);

        const auto cc = D2R(XR_DISTANCE_MULTIPLIER) * ccDiff;
        const auto nn = D2R(XR_DISTANCE_MULTIPLIER) * nnDiff;

        totalDiffSquared += (cc * cc) + (nn * nn) + (muDiff * muDiff);

        const auto totalEuclideanDistance = std::sqrt(totalDiffSquared);

        nearestNeighbor.euclideanDistance = totalEuclideanDistance;
        nearestNeighbor.goldenStepIdx = gsIdx;

        if (totalEuclideanDistance < shortestEuclideanDistance) {
            closestGoldenStepIdx = gsIdx;
            shortestEuclideanDistance = totalEuclideanDistance;
            emergencyNearestNeighbor = nearestNeighbor;
        }

        // Quick cluster rejection relies on golden steps being grouped by clusterNumber.
        // LLKA_initializeClassificationContext() takes care of this grouping.
        if (rejectCluster && lastClusterNumber == gs.clusterNumber)
            goto reject_golden_step;

        rejectCluster = false;
        lastClusterNumber = gs.clusterNumber;

        for (size_t idx = 0; idx < TORSION_CLASSIFICATION_METRIC_CLSPTRS.size(); idx++) {
            // Spread around mean angle difference cannot be easily measured in the -PI <-> PI representation.
            const auto &clsfMetric = cluster.*TORSION_CLASSIFICATION_METRIC_CLSPTRS[idx];
            auto actual = angleAsFull(stepMetrics.*TORSION_STEP_METRIC_CLSPTRS[idx]);

            auto low = clsfMetric.minValue;
            auto high = clsfMetric.maxValue;

            // Branchless switch to "inverted" arc.
            // Whether this is actually faster than an if-else variant is debatable...
            using T = decltype(low);
            T inverted = (low > high);
            T notInverted = !(low > high);
            T shift = inverted * high;
            low -= shift;
            actual -= shift;
            actual += TWO_PI * (actual < 0.0);
            high = inverted * TWO_PI + notInverted * high;

            if (!LLKA_WITHIN_EXCLUSIVE(low, actual, high)) {
                ECHMET_TRACE(LLKATracing, GOLDEN_STEP_REJECTED_TOLERANCE_EXCEEDED, actual, low, high, idx, gs.clusterIdx, ctx->clusters);
                rejectCluster = true;
                goto reject_golden_step;
            }
        }

        {
            const auto low = cluster.CC.minValue;
            const auto high = cluster.CC.maxValue;

            if (!LLKA_WITHIN_EXCLUSIVE(low, stepMetrics.CC, high)) {
                ECHMET_TRACE(LLKATracing, GOLDEN_STEP_REJECTED_TOLERANCE_EXCEEDED, stepMetrics.CC, low, high, 10, gs.clusterIdx, ctx->clusters);
                rejectCluster = true;
                goto reject_golden_step;
            }
        }
        {
            const auto low = cluster.NN.minValue;
            const auto high = cluster.NN.maxValue;

            if (!LLKA_WITHIN_EXCLUSIVE(low, stepMetrics.NN, high)) {
                ECHMET_TRACE(LLKATracing, GOLDEN_STEP_REJECTED_TOLERANCE_EXCEEDED, stepMetrics.NN, low, high, 11, gs.clusterIdx, ctx->clusters);
                rejectCluster = true;
                goto reject_golden_step;
            }
        }

        // NOTE: We could do the insertion without copying but it probably would not make much difference
        nValidNearestNeighbors = insertNearestNeighbor(nearestNeighbor, nearestNeighbors, nValidNearestNeighbors);
reject_golden_step:;
        // Jump right to the end of the loop if we reject the golden step
    }

    if (closestGoldenStepIdx == Arch::INVALID_SIZE_T)
        throw LLKA_CLASSIFICATION_E_WRONG_METRICS;

    ECHMET_TRACE(LLKATracing, ALL_NEAREST_NEIGHBORS, nearestNeighbors, nValidNearestNeighbors, ctx);

    if (nValidNearestNeighbors == 0) {
        assert(emergencyNearestNeighbor.goldenStepIdx != Arch::INVALID_SIZE_T);

        nearestNeighbors[0] = std::move(emergencyNearestNeighbor);
        nValidNearestNeighbors = 1;
    }

    return std::make_tuple(std::move(nearestNeighbors), nValidNearestNeighbors, closestGoldenStepIdx, rejectDelta);
}

static
auto invalidateClassifiedStep(LLKA_ClassifiedStep &classifiedStep)
{
    classifiedStep = {};
    classifiedStep.assignedNtC = LLKA_INVALID_NTC;
    classifiedStep.assignedCANA = LLKA_INVALID_CANA;
    classifiedStep.closestNtC = LLKA_INVALID_NTC;
    classifiedStep.closestCANA = LLKA_INVALID_CANA;
    classifiedStep.closestGoldenStep = "";
}

static
auto precalculateNtCExtendedBackbones()
{
    std::array<LLKA_Points, 96> precalculated{};

    for (size_t idx = 0; idx < 96; idx++) {
        LLKA_Structure ntcStru{
            const_cast<LLKA_Atom *>(NTC_RAW_REFS[idx].atoms),
            NTC_RAW_REFS[idx].nAtoms,
        };
        LLKA_Structure extBkbnStru{};

        auto tRet = LLKA_extractExtendedBackbone(&ntcStru, &extBkbnStru);
        assert(tRet == LLKA_OK);
#ifdef NDEBUG
        (void)tRet;
#endif // NDEBUG

        precalculated[idx].points = new LLKA_Point[extBkbnStru.nAtoms];
        for (size_t pdx = 0; pdx < extBkbnStru.nAtoms; pdx++)
            precalculated[idx].points[pdx] = extBkbnStru.atoms[pdx].coords;
        precalculated[idx].nPoints = extBkbnStru.nAtoms;

        LLKA_destroyStructure(&extBkbnStru);
    }

    return precalculated;
}

// NOTE: This function is defined out-of-order to avoid forward declarations */
static
LLKA_RetCode classifyStep(const LLKA_Structure &stru, const LLKA_ClassificationContext *ctx, LLKA_ClassifiedStep &classifiedStep)
{
    invalidateClassifiedStep(classifiedStep);

    // PERF: We should internalize these functions to avoid doing unnecessary checks over and over again

    // Make sure that we have something resembling a step
    LLKA_StepInfo info;
    auto tRet = LLKA_structureIsStep(&stru, &info);
    if (tRet != LLKA_OK)
        return tRet;

    // These two must be the same for all atoms in the structure (unless there is something really odd going on)
    const auto modelNum = stru.atoms[0].pdbx_PDB_model_num;
    const auto asymId = stru.atoms[0].label_asym_id;

    // Split the step up to individual nucleotides
    LLKA_StructureView nuclFirst = LLKAInternal::extractNucleotideView(&stru, modelNum, asymId, info.firstSeqId);
    LLKA_StructureView nuclSecond = LLKAInternal::extractNucleotideView(&stru, modelNum, asymId, info.secondSeqId);

    LLKA_StructureView riboseViewFirst = {};
    LLKA_StructureView riboseViewSecond = {};

    tRet = LLKAInternal::extractRibose(&nuclFirst, &riboseViewFirst);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&nuclSecond);
        LLKA_destroyStructureView(&nuclFirst);

        return tRet;
    }

    tRet = LLKAInternal::extractRibose(&nuclSecond, &riboseViewSecond);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&riboseViewFirst);
        LLKA_destroyStructureView(&nuclSecond);
        LLKA_destroyStructureView(&nuclFirst);

        return tRet;
    }

    LLKA_StepMetrics stepMetrics{};

    tRet = LLKAInternal::calculateStepMetrics_unchecked(&stru, &stepMetrics);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&riboseViewSecond);
        LLKA_destroyStructureView(&riboseViewFirst);
        LLKA_destroyStructureView(&nuclSecond);
        LLKA_destroyStructureView(&nuclFirst);

        return tRet;
    }

    try {
        // Measure ribose geometry. We will use this regardless of how the NtC assignment turns out
        LLKA_RiboseMetrics riboseMetrics;

        LLKAInternal::riboseMetrics(riboseViewFirst, riboseMetrics);
        classifiedStep.nuAngles_1 = riboseMetrics.nus;
        classifiedStep.ribosePseudorotation_1 = riboseMetrics.P;
        classifiedStep.tau_1 = riboseMetrics.tMax;
        classifiedStep.sugarPucker_1 = riboseMetrics.pucker;

        LLKAInternal::riboseMetrics(riboseViewSecond, riboseMetrics);
        classifiedStep.nuAngles_2 = riboseMetrics.nus;
        classifiedStep.ribosePseudorotation_2 = riboseMetrics.P;
        classifiedStep.tau_2 = riboseMetrics.tMax;
        classifiedStep.sugarPucker_2 = riboseMetrics.pucker;

        // Look for best matching NtC and golden step
        auto [ nearestNeighbors, nValidNearestNeighbors, closestGoldenStepIdx, rejectDelta ] = findClosestNtC(stepMetrics, ctx);
        auto [ bestieClusterIdx, bestieVotes ] = determineBestieClusterIdx(std::views::counted(nearestNeighbors.cbegin(), nValidNearestNeighbors), ctx);

        if (nValidNearestNeighbors == 0)
            ECHMET_TRACE(LLKATracing, DETAILS_STEPS_WITH_NO_NEIGHBORS, stru, stepMetrics);

        const bool notEnoughNearestNeigbors = nValidNearestNeighbors < ctx->limits.minimumNearestNeighbors;

        // NOTE:
        // When the original implementation rejects delta torions, it calculates differences
        // between the actual and cluster mean nu angles (that is why it cannot abort early)
        // and writes them out.
        // It is unclear whether these angles then get written to the Cif file or of it is
        // just an internal sanity check.
        // We should consider removing the early abort with odd delta and add a tracepoint
        // that would calculate and log the nu angles.
        if (rejectDelta)
            classifiedStep.violations |= LLKA_CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED;

        if (notEnoughNearestNeigbors) {
            bestieClusterIdx = ctx->goldenSteps[closestGoldenStepIdx].clusterIdx;
            classifiedStep.violations |= LLKA_CLASSIFICATION_E_NOT_ENOUGH_NEAREST_NEIGHBORS;
        }

        // There is no cluster that matches well enough to properly classify the step.
        // Use the cluster from the most closely matching golden step instead.
        bool notEnoughVotes = bestieVotes < ctx->limits.minimumClusterVotes && bestieVotes > 0;
        if (notEnoughVotes) {
            bestieClusterIdx = ctx->goldenSteps[closestGoldenStepIdx].clusterIdx;
            classifiedStep.violations |= LLKA_CLASSIFICATION_E_BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES;
        }

        const auto &bestieCluster = ctx->clusters[bestieClusterIdx];

        for (size_t idx = 0; idx < LLKAInternal::NU_ANGLES_CLSPTRS.size(); idx++) {
            const auto clsPtrDiff = LLKAInternal::NU_ANGLES_CLSPTRS[idx];
            const auto clsPtrMetrics = LLKAInternal::NU_ANGLES_METRICS_CLSPTRS[idx];

            classifiedStep.nuAngleDifferences_1.*clsPtrDiff = angleDifference(classifiedStep.nuAngles_1.*clsPtrDiff, (bestieCluster.nusFirst.*clsPtrMetrics).meanValue);
            classifiedStep.nuAngleDifferences_2.*clsPtrDiff = angleDifference(classifiedStep.nuAngles_2.*clsPtrDiff, (bestieCluster.nusSecond.*clsPtrMetrics).meanValue);
        }

        ECHMET_TRACE(LLKATracing, BESTIE_CLUSTER_INFO, bestieCluster, notEnoughNearestNeigbors, notEnoughVotes);
        ECHMET_TRACE(LLKATracing, CLOSEST_GOLDEN_STEP_INFO, ctx->goldenSteps[closestGoldenStepIdx], ctx->clusters);

        classifiedStep.closestNtC = bestieCluster.NtC;
        classifiedStep.closestCANA = bestieCluster.CANA;
        classifiedStep.confalScore = {};
        classifiedStep.closestGoldenStep = ctx->goldenSteps[closestGoldenStepIdx].name;

        auto [ distancesFromNtCAverages, euclideanDistanceIdeal ] = calcDistancesFromNtCAverages(stepMetrics, bestieCluster);

        ECHMET_TRACE(LLKATracing, DIFFERENCES_FROM_NTC_AVERAGES, distancesFromNtCAverages, classifiedStep.closestNtC);

        classifiedStep.metrics = stepMetrics;
        classifiedStep.differencesFromNtCAverages = distancesFromNtCAverages;
        classifiedStep.euclideanDistanceNtCIdeal = euclideanDistanceIdeal;

        classifiedStep.rmsdToClosestNtC = calcRmsdToClosestNtC(&stru, ctx->ntcExtendedBackbonePoints[classifiedStep.closestNtC]);

        if (nValidNearestNeighbors > 0) {
            auto [ violations, violatingTorsionsAverage, violatingTorsionsNearest ] = checkNtCTolerances(
                stepMetrics,
                bestieCluster,
                std::views::counted(nearestNeighbors.cbegin(), nValidNearestNeighbors),
                classifiedStep.ribosePseudorotation_1,
                classifiedStep.ribosePseudorotation_2,
                ctx
            );

            classifiedStep.violations |= violations;
            classifiedStep.violatingTorsionsAverage = violatingTorsionsAverage;
            classifiedStep.violatingTorsionsNearest = violatingTorsionsNearest;

            classifiedStep.confalScore = calcConfalScore(
                classifiedStep.differencesFromNtCAverages,
                ctx->confals[bestieClusterIdx],
                classifiedStep.violations == LLKA_CLASSIFICATION_OK
            );

            if (classifiedStep.violations == LLKA_CLASSIFICATION_OK) {
                classifiedStep.assignedNtC = bestieCluster.NtC;
                classifiedStep.assignedCANA = bestieCluster.CANA;
            } else {
                if (classifiedStep.rmsdToClosestNtC <= ctx->maxCloseEnoughRmsd)
                    classifiedStep.violations |= LLKA_CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH;
            }
        }

        LLKA_destroyStructureView(&riboseViewSecond);
        LLKA_destroyStructureView(&riboseViewFirst);

        LLKA_destroyStructureView(&nuclSecond);
        LLKA_destroyStructureView(&nuclFirst);

        return LLKA_OK;
    } catch (const decltype(LLKA_CLASSIFICATION_OK) violations) {
        LLKAInternal::invalidateClassifiedStep(classifiedStep);
        classifiedStep.violations = violations;

        LLKA_destroyStructureView(&riboseViewSecond);
        LLKA_destroyStructureView(&riboseViewFirst);

        LLKA_destroyStructureView(&nuclSecond);
        LLKA_destroyStructureView(&nuclFirst);

        return LLKA_OK;
    }
}

} // namespace LLKAInternal

LLKA_AverageConfal LLKA_CC LLKA_averageConfal(const LLKA_ClassifiedStep *classifiedSteps, size_t nClassifiedSteps, const LLKA_ClassificationContext *ctx)
{
    double totalConfal = 0.0;

    if (nClassifiedSteps== 0) {
        return { .score = 0.0, .percentile = 0.0 };
    }

    for (size_t idx = 0; idx < nClassifiedSteps; idx++) {
        const auto &step = classifiedSteps[idx];
        totalConfal += step.confalScore.total;
    }

    double averageConfal = totalConfal / nClassifiedSteps;
    int percentileIndex = std::floor(averageConfal);
    assert(percentileIndex >= 0 && percentileIndex<= 100);

    return {
        .score = averageConfal,
        .percentile = ctx->confalPercentiles.at(percentileIndex)
    };
}

LLKA_API LLKA_AverageConfal LLKA_CC LLKA_averageConfalAttempted(const LLKA_ClassifiedSteps *attemptedClassifiedSteps, const LLKA_ClassificationContext *ctx)
{
    double totalConfal = 0.0;
    size_t usedSteps = 0;

    for (size_t idx = 0; idx < attemptedClassifiedSteps->nAttemptedSteps; idx++) {
        const auto &step = attemptedClassifiedSteps->attemptedSteps[idx];
        if (step.status == LLKA_OK) {
            totalConfal += step.step.confalScore.total;
            usedSteps++;
        }
    }

    if (usedSteps == 0) {
        return { .score = 0.0, .percentile = 0.0 };
    }

    double averageConfal = totalConfal / usedSteps;
    int percentileIndex = std::floor(averageConfal);
    assert(percentileIndex >= 0 && percentileIndex<= 100);

    return {
        .score = averageConfal,
        .percentile = ctx->confalPercentiles.at(percentileIndex)
    };
}

LLKA_RetCode LLKA_CC LLKA_classificationClusterForNtC(LLKA_NtC ntc, const LLKA_ClassificationContext *ctx, LLKA_ClassificationCluster *cluster)
{
    for (const auto &c : ctx->clusters) {
        if (c.NtC == ntc) {
            *cluster = c;
            return LLKA_OK;
        }
    }

    return LLKA_E_INVALID_ARGUMENT;
}

const char * LLKA_CC LLKA_classificationViolationToName(int32_t violation)
{
    switch (violation) {
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_OK);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_SCORE_TOO_LOW);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_NOT_ENOUGH_NEAREST_NEIGHBORS);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_CC_TOO_LOW);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_CC_TOO_HIGH);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_NN_TOO_LOW);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_NN_TOO_HIGH);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_MU_TOO_LOW);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_MU_TOO_HIGH);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED);
    CLASSIFICATION_VIOLATION_STR(LLKA_CLASSIFICATION_E_WRONG_METRICS);
    }

    return "Unknown classification violation";
}

LLKA_RetCode LLKA_CC LLKA_classifyStep(const LLKA_Structure *stru, const LLKA_ClassificationContext *ctx, LLKA_ClassifiedStep *classifiedStep)
{
    return LLKAInternal::classifyStep(*stru, ctx, *classifiedStep);
}

LLKA_RetCode LLKA_CC LLKA_classifyStepsMultiple(const LLKA_Structures *strus, const LLKA_ClassificationContext *ctx, LLKA_ClassifiedSteps *classifiedSteps)
{
    if (strus->nStrus == 0)
        return LLKA_E_NOTHING_TO_CLASSIFY;

    classifiedSteps->attemptedSteps = new LLKA_AttemptedClassifiedStep[strus->nStrus];
    classifiedSteps->nAttemptedSteps = strus->nStrus;

    for (size_t idx = 0; idx < strus->nStrus; idx++) {
        ECHMET_TRACE(LLKATracing, BEGIN_STEP_CLASSIFICATION_MULTIPLE, idx);

        const auto &stru = strus->strus[idx];
        auto &attempt = classifiedSteps->attemptedSteps[idx];

        attempt.status = LLKAInternal::classifyStep(stru, ctx, attempt.step);
    }

    return LLKA_OK;
}

LLKA_RetCode LLKA_CC LLKA_confalForNtC(LLKA_NtC ntc, const LLKA_ClassificationContext *ctx, LLKA_Confal *confal)
{
    int32_t number = -1;

    for (const auto &c : ctx->clusters) {
        if (c.NtC == ntc) {
            number = c.number;
            break;
        }
    }

    if (number == -1)
        return LLKA_E_INVALID_ARGUMENT;

    for (const auto &c : ctx->confals) {
        if (c.clusterNumber == number) {
            *confal = c;
            return LLKA_OK;
        }
    }

    return LLKA_E_INVALID_ARGUMENT;
}

double LLKA_CC LLKA_confalPercentile(double confalScore, const LLKA_ClassificationContext *ctx)
{
    if (confalScore < 0 || confalScore > 100)
        return -1;

    return ctx->confalPercentiles.at(std::floor(confalScore));
}

void LLKA_CC LLKA_destroyClassificationContext(LLKA_ClassificationContext *ctx)
{
    delete ctx;
}

void LLKA_CC LLKA_destroyClassifiedSteps(LLKA_ClassifiedSteps *classifiedSteps)
{
    if (classifiedSteps->nAttemptedSteps > 0)
        delete [] classifiedSteps->attemptedSteps;
}

LLKA_RetCode LLKA_CC LLKA_initializeClassificationContext(
    const LLKA_ClassificationCluster *clusters, size_t nClusters,
    const LLKA_GoldenStep *goldenSteps, size_t nGoldenSteps,
    const LLKA_Confal *confals, size_t nConfals,
    const LLKA_ClusterNuAngles *clusterNuAngles, size_t nClusterNuAngles,
    const LLKA_ConfalPercentile *confalPercentiles, size_t nConfalPercentiles,
    const LLKA_ClassificationLimits *limits,
    double maxCloseEnoughRmsd,
    LLKA_ClassificationContext **ctx
) {
    if (nGoldenSteps == 0 || nClusters == 0 || nConfals == 0 || nClusterNuAngles == 0)
        return LLKA_E_INVALID_ARGUMENT;
    if (maxCloseEnoughRmsd <= 0.0)
        return LLKA_E_INVALID_ARGUMENT;

    if (nClusters != nConfals)
        return LLKA_E_MISMATCHING_SIZES;
    if (nClusters != nClusterNuAngles)
        return LLKA_E_MISMATCHING_SIZES;

    if (nConfalPercentiles != 101)
        return LLKA_E_BAD_DATA;

    auto _ctx = std::make_unique<LLKA_ClassificationContext>();

    _ctx->clusters.resize(nClusters);

    // Track which cluster numbers correspond to indices in the clusters vector.
    // We can use this to tell golden steps where to get their respective clusters.
    std::map<int32_t, size_t> clusterNumberToIndexMapping{};

    for (size_t idx = 0; idx < nClusters; idx++) {
        const auto &c = clusters[idx];

        if (LLKAInternal::containsKey(clusterNumberToIndexMapping, c.number))
            return LLKA_E_BAD_CLASSIFICATION_CLUSTERS;
        clusterNumberToIndexMapping[c.number] = idx;

        _ctx->clusters[idx] = c;

        auto &cluster = _ctx->clusters[idx];
        for (const auto clsPtr : LLKAInternal::ALL_TORSIONS_CLASSIFICATION_METRIC_CLSPTRS) {
            auto &metric = cluster.*clsPtr;

            if (metric.meanValue < 0.0 || metric.deviation < 0.0)
                return LLKA_E_BAD_CLASSIFICATION_CLUSTERS;

            metric.deviation *= LLKAInternal::BACKBONE_TORSIONS_DEVIATION_MULTIPLIER;
            metric.minValue = LLKAInternal::clampAngle(LLKAInternal::angleAsFull(metric.meanValue - metric.deviation));
            metric.maxValue = LLKAInternal::clampAngle(LLKAInternal::angleAsFull(metric.meanValue + metric.deviation));
        }

        cluster.CC.deviation *= LLKAInternal::XR_DISTANCE_DEVIATION_MULTIPLIER;
        cluster.CC.minValue = cluster.CC.meanValue - cluster.CC.deviation;
        cluster.CC.maxValue = cluster.CC.meanValue + cluster.CC.deviation;

        cluster.NN.deviation *= LLKAInternal::XR_DISTANCE_DEVIATION_MULTIPLIER;
        cluster.NN.minValue = cluster.NN.meanValue - cluster.NN.deviation;
        cluster.NN.maxValue = cluster.NN.meanValue + cluster.NN.deviation;

        cluster.mu.deviation *= LLKAInternal::MU_TORSION_DEVIATION_MULTIPLIER;
        cluster.mu.minValue = LLKAInternal::clampAngle(LLKAInternal::angleAsFull(cluster.mu.meanValue - cluster.mu.deviation));
        cluster.mu.maxValue = LLKAInternal::clampAngle(LLKAInternal::angleAsFull(cluster.mu.meanValue + cluster.mu.deviation));
    }

    _ctx->goldenSteps.resize(nGoldenSteps);
    for (size_t idx = 0; idx < nGoldenSteps; idx++) {
        const auto &gs = goldenSteps[idx];

        if (!LLKAInternal::containsKey(clusterNumberToIndexMapping, gs.clusterNumber))
            return LLKA_E_BAD_GOLDEN_STEPS;

        _ctx->goldenSteps[idx] = gs;
        _ctx->goldenSteps[idx].name = LLKAInternal::duplicateString(gs.name);
        _ctx->goldenSteps[idx].clusterIdx = clusterNumberToIndexMapping[gs.clusterNumber];
    }

    // Sort the golden steps by cluster number. This will allow us to quicky reject an entire cluster
    // in the classification process.
    std::sort(
        _ctx->goldenSteps.begin(),
        _ctx->goldenSteps.end(),
        [](const auto &lhs, const auto &rhs) {
            return lhs.clusterNumber < rhs.clusterNumber;
        }
    );

    _ctx->confals.resize(nConfals);
    for (size_t idx = 0; idx < nConfals; idx++) {
        const auto &confal = confals[idx];

        if (!LLKAInternal::containsKey(clusterNumberToIndexMapping, confal.clusterNumber))
            return LLKA_E_BAD_CONFALS;
        const auto clusterIdx = clusterNumberToIndexMapping[confal.clusterNumber];

        _ctx->confals[clusterIdx] = confal;
    }

    for (size_t idx = 0; idx < nClusterNuAngles; idx++) {
        const auto &nus = clusterNuAngles[idx];

        if (!LLKAInternal::containsKey(clusterNumberToIndexMapping, nus.clusterNumber))
            return LLKA_E_BAD_AVERAGE_NU_ANGLES;
        const auto clusterIdx = clusterNumberToIndexMapping[nus.clusterNumber];
        auto &cluster = _ctx->clusters[clusterIdx];

        cluster.nusFirst = nus.firstNucleotide;
        cluster.nusSecond = nus.secondNucleotide;

        // Initialize Nu angles data
        for (const auto &clsPtr : LLKAInternal::NU_ANGLES_METRICS_CLSPTRS) {
            auto &first = cluster.nusFirst.*clsPtr;
            first.meanValue = LLKAInternal::angleAsFull(first.meanValue);
            first.minValue = LLKAInternal::clampAngle(LLKAInternal::angleAsFull(first.meanValue - first.deviation));
            first.maxValue = LLKAInternal::clampAngle(LLKAInternal::angleAsFull(first.meanValue + first.deviation));

            auto &second = cluster.nusSecond.*clsPtr;
            second.meanValue = LLKAInternal::angleAsFull(second.meanValue);
            second.minValue = LLKAInternal::clampAngle(LLKAInternal::angleAsFull(second.meanValue - second.deviation));
            second.maxValue = LLKAInternal::clampAngle(LLKAInternal::angleAsFull(second.meanValue + second.deviation));
        }
    }

    _ctx->confalPercentiles.resize(nConfalPercentiles);
    for (size_t idx = 0; idx < nConfalPercentiles; idx++)
        _ctx->confalPercentiles[idx] = confalPercentiles[idx].value;

    if (
        limits->averageNeighborsTorsionCutoff <= 0.0 || limits->nearestNeighborTorsionsCutoff <= 0.0 || limits->totalDistanceCutoff <= 0.0 || limits->pseudorotationCutoff <= 0.0 ||
        limits->minimumClusterVotes <= 0.0 ||
        limits->minimumNearestNeighbors < 1 || limits->numberOfUsedNearestNeighbors < limits->minimumNearestNeighbors
    )
        return LLKA_E_BAD_CLASSIFICATION_LIMITS;

    _ctx->limits = *limits;
    _ctx->maxCloseEnoughRmsd = maxCloseEnoughRmsd;

    _ctx->ntcExtendedBackbonePoints = LLKAInternal::precalculateNtCExtendedBackbones();

    *ctx = _ctx.release();

    return LLKA_OK;
}

#ifndef ECHMET_TRACER_DISABLE_TRACING

#include <iomanip>
#include <sstream>

namespace ECHMET {

ECHMET_MAKE_TRACEPOINT_NOINLINE(LLKATracing, CLASSIFICATION_METRICS_DIFFERENCES, "Classification metrics differences")
ECHMET_BEGIN_MAKE_LOGGER(LLKATracing, CLASSIFICATION_METRICS_DIFFERENCES, const LLKA_StepMetrics &stepMetrics, const LLKA_GoldenStep &gs, const LLKA_StepMetrics &differences)
{
    const auto &gsMetrics = gs.metrics;

    std::ostringstream oss{};

    oss
        << "Measuring against golden step " << gs.name << "\n"
        << "delta_1:   " << stepMetrics.delta_1 << ", " << gsMetrics.delta_1 << ", " << differences.delta_1 << "\n"
        << "epsilon_1: " << stepMetrics.epsilon_1 << ", " << gsMetrics.epsilon_1  << ", "  << differences.epsilon_1 << "\n"
        << "zeta_1:    " << stepMetrics.zeta_1 << ", " << gsMetrics.zeta_1 << ", " << differences.zeta_1 << "\n"
        << "alpha_2:   " << stepMetrics.alpha_2 << ", " << gsMetrics.alpha_2 << ", " << differences.alpha_2 << "\n"
        << "beta_2:    " << stepMetrics.beta_2 << ", " << gsMetrics.beta_2 << ", "  << differences.beta_2 << "\n"
        << "gamma_2:   " << stepMetrics.gamma_2 << ", " << gsMetrics.gamma_2 << ", " << differences.gamma_2 << "\n"
        << "delta_2:   " << stepMetrics.delta_2 << ", " << gsMetrics.delta_2 << ", " << differences.delta_2 << "\n"
        << "chi_1:     " << stepMetrics.chi_1 << ", " << gsMetrics.chi_1 << ", " << differences.chi_1 << "\n"
        << "chi_2:     " << stepMetrics.chi_2 << ", " << gsMetrics.chi_2 << ", " << differences.chi_2 << "\n"
        << "CC:        " << stepMetrics.CC << ", " << gsMetrics.CC << ", " << differences.CC << "\n"
        << "NN:        " << stepMetrics.NN << ", " << gsMetrics.NN << ", " << differences.NN << "\n"
        << "mu:        " << stepMetrics.mu << ", " << gsMetrics.mu << ", "  << differences.mu << "\n";

	return oss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT_NOINLINE(LLKATracing, ALL_NEAREST_NEIGHBORS, "All nearest neighbors of a step that will be voted for")
ECHMET_BEGIN_MAKE_LOGGER(LLKATracing, ALL_NEAREST_NEIGHBORS, const std::vector<LLKAInternal::NearestNeighbor> &nearestNeighbors, size_t nValidNearestNeighbors, const LLKA_ClassificationContext *ctx)
{
    std::ostringstream oss{};

    oss << "--- Nearest neighbors to vote for ---\n";
    oss << "No.\tEucl. dist.\tGolden step\n";

    for (size_t idx = 0; idx < nValidNearestNeighbors; idx++) {
        const auto &nn = nearestNeighbors[idx];

        oss << idx << "\t" << std::setprecision(7) << nn.euclideanDistance << "\t" << ctx->goldenSteps[nn.goldenStepIdx].name << "\n";
    }

    return oss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT_NOINLINE(LLKATracing, BEGIN_STEP_CLASSIFICATION_MULTIPLE, "Report the beginning of a step classification from a function that classifies multiple steps in sequence")
ECHMET_BEGIN_MAKE_LOGGER(LLKATracing, BEGIN_STEP_CLASSIFICATION_MULTIPLE, size_t idx)
{
    return "Attempting to classify step " + std::to_string(idx) + "\n";
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT_NOINLINE(LLKATracing, DIFFERENCES_FROM_NTC_AVERAGES, "Differences between step metrics values and averages for the closest NtC class")
ECHMET_BEGIN_MAKE_LOGGER(LLKATracing, DIFFERENCES_FROM_NTC_AVERAGES, LLKA_StepMetrics &diffs, const LLKA_NtC ntc)
{
    std::ostringstream oss{};

    oss
        << "Differences between NtC averages (closest NtC " << LLKA_NtCToName(ntc) << ")\n"
        << "delta_1:   " << diffs.delta_1 << "\n"
        << "epsilon_1: " << diffs.epsilon_1 << "\n"
        << "zeta_1:    " << diffs.zeta_1 << "\n"
        << "alpha_2:   " << diffs.alpha_2 << "\n"
        << "beta_2:    " << diffs.beta_2 << "\n"
        << "gamma_2:   " << diffs.gamma_2 << "\n"
        << "delta_2:   " << diffs.delta_2 << "\n"
        << "chi_1:     " << diffs.chi_1 << "\n"
        << "chi_2:     " << diffs.chi_2 << "\n"
        << "CC:        " << diffs.CC << "\n"
        << "NN:        " << diffs.NN << "\n"
        << "mu:        " << diffs.mu << "\n";

    return oss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT_NOINLINE(LLKATracing, GOLDEN_STEP_REJECTED_TOLERANCE_EXCEEDED, "Golden step was rejected as neighbor because some metrics exceeds tolerance")
ECHMET_BEGIN_MAKE_LOGGER(LLKATracing, GOLDEN_STEP_REJECTED_TOLERANCE_EXCEEDED, double actual, double low, double high, size_t metricsIdx, size_t clusterIdx, const std::vector<LLKA_ClassificationCluster> &clusters)
{
    std::ostringstream oss{};

    oss
        << "Rejecting cluster " << clusterIdx << " (" << LLKA_NtCToName(clusters[clusterIdx].NtC) << ") " << " because metrics " << metricsIdx << " exceeds tolerance (low, actual, high): [ "
        << low << "; " << actual << "; " << high << " ]\n";

    return oss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT_NOINLINE(LLKATracing, DETAILS_STEPS_WITH_NO_NEIGHBORS, "Print detailed information about step that does not have any neighbors")
ECHMET_BEGIN_MAKE_LOGGER(LLKATracing, DETAILS_STEPS_WITH_NO_NEIGHBORS, const LLKA_Structure &step, const LLKA_StepMetrics &metrics)
{
    std::ostringstream oss{};

    oss
        << "--- This step does not have any neighbors ---\n" << step << "\n"
        << "delta_1:   " << metrics.delta_1 << "\n"
        << "epsilon_1: " << metrics.epsilon_1 << "\n"
        << "zeta_1:    " << metrics.zeta_1 << "\n"
        << "alpha_2:   " << metrics.alpha_2 << "\n"
        << "beta_2:    " << metrics.beta_2 << "\n"
        << "gamma_2:   " << metrics.gamma_2 << "\n"
        << "delta_2:   " << metrics.delta_2 << "\n"
        << "chi_1:     " << metrics.chi_1 << "\n"
        << "chi_2:     " << metrics.chi_2 << "\n"
        << "CC:        " << metrics.CC << "\n"
        << "NN:        " << metrics.NN << "\n"
        << "mu:        " << metrics.mu << "\n";

    return oss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT_NOINLINE(LLKATracing, BESTIE_CLUSTER_INFO, "Information about the chosen classification cluster")
ECHMET_BEGIN_MAKE_LOGGER(LLKATracing, BESTIE_CLUSTER_INFO, const LLKA_ClassificationCluster &cluster, const bool notEnoughNN, const bool notEnoughVotes)
{
    std::ostringstream oss{};

    oss
        << "Selected best matching cluster: " << cluster.number << ", " << LLKA_NtCToName(cluster.NtC) << ", " << LLKA_CANAToName(cluster.CANA) << "\n"
        << "Selected from enough neighbors: " << (notEnoughNN ? "No" : "Yes") << "\n"
        << "Got enough votes: " << (notEnoughVotes ? "No" : "Yes") << "\n";

    return oss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT_NOINLINE(LLKATracing, PSEUDOROTATION_TOO_DIFFERENT, "Details about pseudorotations exceeding tolerance")
ECHMET_BEGIN_MAKE_LOGGER(LLKATracing, PSEUDOROTATION_TOO_DIFFERENT, double actual1, double reference1, double difference1, double actual2, double reference2, double difference2, double tolerance)
{
    std::ostringstream oss{};

    oss
        << "Pseudorotation 1 (actual, ref, diff): " << actual1 << ", " << reference1 << ", " << difference1 << (difference1 > tolerance ? " (exceeds)" : "") << "\n"
        << "Pseudorotation 2 (actual, ref, diff): " << actual2 << ", " << reference2 << ", " << difference2 << (difference2 > tolerance ? " (exceeds)" : "") << "\n"
        << "Maximum difference: " << tolerance << "\n";

    return oss.str();
}
ECHMET_END_MAKE_LOGGER

ECHMET_MAKE_TRACEPOINT_NOINLINE(LLKATracing, CLOSEST_GOLDEN_STEP_INFO, "Details about the closest golden step")
ECHMET_BEGIN_MAKE_LOGGER(LLKATracing, CLOSEST_GOLDEN_STEP_INFO, const LLKA_GoldenStep &gs, const std::vector<LLKA_ClassificationCluster> &clusters)
{
    std::ostringstream oss{};

    oss
        << " --- Closest golden step ---\n"
        << "Name: " << gs.name << "\n"
        << "Cluster NtC: " << LLKA_NtCToName(clusters[gs.clusterIdx].NtC) << "\n"
        << "Cluster CANA: " << LLKA_CANAToName(clusters[gs.clusterIdx].CANA) << "\n";

    return oss.str();
}
ECHMET_END_MAKE_LOGGER

} // namespace ECHMET

#endif // ECHMET_TRACER_DISABLE_TRACING
