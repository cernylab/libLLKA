/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "testing_structures.h"

#include <llka_classification.h>
#include <llka_resource_loaders.h>
#include <llka_tracing.h>

#include "effedup.hpp"

#include "../src/util/elementaries.h"

namespace EffedUp {

template <>
auto stringify<decltype(LLKA_CLASSIFICATION_OK)>(const decltype(LLKA_CLASSIFICATION_OK) &v) -> std::string
{
    std::vector<std::string> errs{};

    for (int bit = 0; bit < 32; bit++) {
        const auto singleViolation = (1 << bit);
        if (v & singleViolation)
            errs.push_back(LLKA_classificationViolationToName(singleViolation));
    }

    std::string out{};
    for (const auto &e : errs)
        out += e + ", ";

    return out;
}

} // namespace EffedUp

static
auto testClassifyNotClassifiable(const LLKA_ClassificationContext *ctx)
{
    LLKA_Structure stru = LLKA_makeStructure(
        REAL_1DK1_B_27_28_ATOMS,
        REAL_1DK1_B_27_28_ATOMS_LEN
    );

    LLKA_ClassifiedStep classifiedStep{};
    auto tRet = LLKA_classifyStep(&stru, ctx, &classifiedStep);
    EFF_expect(tRet, LLKA_E_MULTIPLE_ALT_IDS, "unexpected return value");

    LLKA_destroyStructure(&stru);
}

static
auto testClassifySimple(const LLKA_ClassificationContext *ctx)
{
    using VT = decltype(LLKA_CLASSIFICATION_OK);

    LLKA_Structure stru = LLKA_makeStructure(
        REAL_1BNA_A_1_2_ATOMS,
        REAL_1BNA_A_1_2_ATOMS_LEN
    );

    LLKA_ClassifiedStep classifiedStep{};
    auto tRet = LLKA_classifyStep(&stru, ctx, &classifiedStep);
    EFF_expect(tRet, LLKA_OK, "unable to classify step");

    EFF_expect(VT(classifiedStep.violations), LLKA_CLASSIFICATION_OK, "detected classification violations where there were not supposed to be any");
    EFF_expect(classifiedStep.assignedNtC, LLKA_BB04, "wrong NtC class");
    EFF_expect(classifiedStep.closestGoldenStep, "1mjo_F_DG_10_DA_11", "wrong golden step name");
    EFF_cmpFlt(classifiedStep.rmsdToClosestNtC, 0.366977877615573, "RMSD to closest NtC exceeds tolerance");

    EFF_cmpFlt(classifiedStep.confalScore.CC, 80.5998065953948, "wrong CC confal score");
    EFF_cmpFlt(classifiedStep.confalScore.NN, 99.9952502902733, "wrong NN confal score");
    EFF_cmpFlt(classifiedStep.confalScore.total, 32.9484362098594, "wrong total confal score");

    LLKA_destroyStructure(&stru);
}

static
auto testClassifyThorough(const LLKA_ClassificationContext *ctx)
{
    using VT = decltype(LLKA_CLASSIFICATION_OK);

    LLKA_Structure stru = LLKA_makeStructure(
        REAL_3VOK_U_1_2_ATOMS,
        REAL_3VOK_U_1_2_ATOMS_LEN
    );

    LLKA_ClassifiedStep classifiedStep{};
    auto tRet = LLKA_classifyStep(&stru, ctx, &classifiedStep);
    EFF_expect(tRet, LLKA_OK, "unable to classify step");

    EFF_expect(VT(classifiedStep.violations), LLKA_CLASSIFICATION_OK, "detected classification violations where there were not supposed to be any");
    EFF_expect(classifiedStep.assignedNtC, LLKA_BBS1, "wrong NtC class");
    EFF_expect(classifiedStep.assignedCANA, LLKA_SYN, "wrong CANA class");
    EFF_expect(classifiedStep.assignedNtC, classifiedStep.closestNtC, "assigned and closest NtCs do not match");
    EFF_expect(classifiedStep.assignedCANA, classifiedStep.closestCANA, "assigned and closest CANAs do not match");
    EFF_expect(classifiedStep.closestGoldenStep, "1lv5_F_DA_1_DC_2", "wrong golden step");

    EFF_cmpFlt(classifiedStep.confalScore.CC, 35.6080290838014, "wrong CC confal score");
    EFF_cmpFlt(classifiedStep.confalScore.NN, 43.1082983683617, "wrong NN confal score");
    EFF_cmpFlt(classifiedStep.confalScore.total, 72.3729761097013, "wrong total confal score");

    EFF_expect(classifiedStep.sugarPucker_1, LLKA_C2_ENDO, "wrong sugar pucker on first nucleotide");
    EFF_expect(classifiedStep.sugarPucker_2, LLKA_C1_EXO, "wrong sugar pucker on second nucleotide");
    EFF_cmpFlt(classifiedStep.ribosePseudorotation_1, 2.84035778887486, "wrong ribose pseudorotation on first ribose");
    EFF_cmpFlt(classifiedStep.tau_2, 0.7674388643588, "wrong tau_2");
    EFF_cmpFlt(classifiedStep.nuAngles_2.nu_2, -0.608524908531341, "wrong nu2 angle on second nucleotide");

    auto avg = LLKA_averageConfal(&classifiedStep, 1, ctx);
    EFF_cmpFlt(avg.score, 72.372976109701284, "wrong average confal score");
    EFF_cmpFlt(avg.percentile, 94.599999999999994, "wrong confal percentile");

    LLKA_destroyStructure(&stru);
}

static
auto testClassifyViolations(const LLKA_ClassificationContext *ctx)
{
    LLKA_Structure stru = LLKA_makeStructure(
        REAL_1DK1_B_5_6_ATOMS,
        REAL_1DK1_B_5_6_ATOMS_LEN
    );

    LLKA_ClassifiedStep classifiedStep{};
    auto tRet = LLKA_classifyStep(&stru, ctx, &classifiedStep);
    EFF_expect(tRet, LLKA_OK, "unable to classify step");

    EFF_expect(
        classifiedStep.violations,
        LLKA_CLASSIFICATION_E_NOT_ENOUGH_NEAREST_NEIGHBORS |
        LLKA_CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT |
        LLKA_CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT |
        LLKA_CLASSIFICATION_E_BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES,
        "wrong violations"
    );
    EFF_expect(classifiedStep.assignedNtC, LLKA_INVALID_NTC, "no NtC should be assigned");
    EFF_expect(classifiedStep.closestNtC, LLKA_IC03, "wrong NtC class");

    LLKA_destroyStructure(&stru);
}

static
auto testGetCluster(const LLKA_ClassificationContext *ctx)
{
    LLKA_ClassificationCluster cluster;

    auto tRet = LLKA_classificationClusterForNtC(LLKA_AA00, ctx, &cluster);
    EFF_expect(tRet, LLKA_OK, "LLKA_classificationClusterForNtC() returned unexpected value");
    EFF_expect(cluster.NtC, LLKA_AA00, "wrong cluster NtC")

    tRet = LLKA_classificationClusterForNtC(LLKA_INVALID_NTC, ctx, &cluster);
    EFF_expect(tRet, LLKA_E_INVALID_ARGUMENT, "LLKA_classificationClusterForNtC() returned unexpected value");
}

static
auto testSugarPuckerNaming()
{
    auto pucker = LLKA_nameToSugarPucker("C1end");
    EFF_expect(pucker, LLKA_C1_ENDO, "wrong sugar pucker");

    pucker = LLKA_nameToSugarPucker("C3'exo");
    EFF_expect(pucker, LLKA_C3_EXO, "wrong sugar pucker");

    pucker = LLKA_nameToSugarPucker("C4' exo");
    EFF_expect(pucker, LLKA_C4_EXO, "wrong sugar pucker");

    pucker = LLKA_nameToSugarPucker("O4' endo");
    EFF_expect(pucker, LLKA_O4_ENDO, "wrong sugar pucker");

    pucker = LLKA_nameToSugarPucker("O1'exo");
    EFF_expect(pucker, LLKA_O4_EXO, "wrong sugar pucker");

    pucker = LLKA_nameToSugarPucker("C3 exo");
    EFF_expect(pucker, LLKA_INVALID_SUGAR_PUCKER, "got valid suger pucker for an invalid name");

    std::string name = LLKA_sugarPuckerToName(LLKA_C1_ENDO, LLKA_SPN_TERSE);
    EFF_expect(name, std::string{"C1endo"}, "wrong sugar pucker name");

    name = LLKA_sugarPuckerToName(LLKA_C1_EXO, LLKA_SPN_FANCY);
    EFF_expect(name, std::string{"C1' exo"}, "wrong sugar pucker name");

    name = LLKA_sugarPuckerToName(LLKA_INVALID_SUGAR_PUCKER, LLKA_SPN_VERY_TERSE);
    EFF_expect(name, std::string{""}, "wrong sugar pucker name");
}

static
auto initializeClassificationContext()
{
    LLKA_Resource goldenSteps = {};
    goldenSteps.type = LLKA_RES_GOLDEN_STEPS;
    auto tRet = LLKA_loadResourceFile(LLKA_PathLiteral("./golden_steps.csv"), &goldenSteps);
    EFF_expect(tRet, LLKA_OK, "could not load golden steps definitions");

    LLKA_Resource clusters = {};
    clusters.type = LLKA_RES_CLUSTERS;
    tRet = LLKA_loadResourceFile(LLKA_PathLiteral("./clusters.csv"), &clusters);
    EFF_expect(tRet, LLKA_OK, "could not load clusters definitions");

    LLKA_Resource confals = {};
    confals.type= LLKA_RES_CONFALS;
    tRet = LLKA_loadResourceFile(LLKA_PathLiteral("./confals.csv"), &confals);
    EFF_expect(tRet, LLKA_OK, "could not load confals definitions");

    LLKA_Resource nuAngles = {};
    nuAngles.type = LLKA_RES_AVERAGE_NU_ANGLES;
    tRet = LLKA_loadResourceFile(LLKA_PathLiteral("./nu_angles.csv"), &nuAngles);
    EFF_expect(tRet, LLKA_OK, "could not load average nu angles definitions");

    LLKA_Resource confalPercentiles = {};
    confalPercentiles.type = LLKA_RES_CONFAL_PERCENTILES;
    tRet = LLKA_loadResourceFile(LLKA_PathLiteral("./confal_percentiles.csv"), &confalPercentiles);
    EFF_expect(tRet, LLKA_OK, "could not load confal percentiles definitions");

    LLKA_ClassificationLimits limits = {};
    limits.averageNeighborsTorsionCutoff = LLKAInternal::D2R(28.0);
    limits.nearestNeighborTorsionsCutoff = LLKAInternal::D2R(28.0);
    limits.totalDistanceCutoff = LLKAInternal::D2R(60.0);
    limits.pseudorotationCutoff = LLKAInternal::D2R(72.0);
    limits.minimumClusterVotes = 0.001111;
    limits.minimumNearestNeighbors = 7;
    limits.numberOfUsedNearestNeighbors = 11;
    const double MAX_CLOSE_ENOUGH_RMSD = 0.5;

    LLKA_ClassificationContext *ctx;
    tRet = LLKA_initializeClassificationContext(
        clusters.data.clusters, clusters.count,
        goldenSteps.data.goldenSteps, goldenSteps.count,
        confals.data.confals, confals.count,
        nuAngles.data.clusterNuAngles, nuAngles.count,
        confalPercentiles.data.confalPercentiles, confalPercentiles.count,
        &limits,
        MAX_CLOSE_ENOUGH_RMSD,
        &ctx
    );
    EFF_expect(tRet, LLKA_OK, "failed to initialize classification context");

    LLKA_destroyResource(&nuAngles);
    LLKA_destroyResource(&confals);
    LLKA_destroyResource(&clusters);
    LLKA_destroyResource(&goldenSteps);
    LLKA_destroyResource(&confalPercentiles);

    return ctx;
}

auto main(int, char **) -> int
{
    testSugarPuckerNaming();

    auto ctx = initializeClassificationContext();

    testClassifySimple(ctx);
    testClassifyThorough(ctx);
    testClassifyViolations(ctx);
    testClassifyNotClassifiable(ctx);
    testGetCluster(ctx);

    LLKA_destroyClassificationContext(ctx);
}
