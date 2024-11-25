/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_cpp.h>

#include "effedup.hpp"
#include "testing_structures.h"

#include "../src/util/elementaries.h"

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED

static
auto classification_initContext()
{
    auto resCl = LLKA::loadClusters(std::filesystem::path{"./clusters.csv"});
    EFF_expect(resCl.isSuccess(), true, "failed to load clusters.csv: " + LLKA::errorToString(resCl.failure()));
    auto resGs = LLKA::loadGoldenSteps(std::filesystem::path{"./golden_steps.csv"});
    EFF_expect(resGs.isSuccess(), true, "failed to load golden_steps.csv: " + LLKA::errorToString(resGs.failure()));
    auto resCf = LLKA::loadConfals(std::filesystem::path{"./confals.csv"});
    EFF_expect(resCf.isSuccess(), true, "failed to load confals.csv: " + LLKA::errorToString(resCf.failure()));
    auto resNu = LLKA::loadClusterNuAngles(std::filesystem::path{"./nu_angles.csv"});
    EFF_expect(resNu.isSuccess(), true, "failed to load nu_angles.csv: " + LLKA::errorToString(resNu.failure()));
    auto resCp = LLKA::loadConfalPercentiles(std::filesystem::path{"./confal_percentiles.csv"});
    EFF_expect(resCp.isSuccess(), true, "failed to load confal_percentiles.csv: " + LLKA::errorToString(resCp.failure()));

    LLKA_ClassificationLimits limits = {};
    limits.averageNeighborsTorsionCutoff = LLKAInternal::D2R(28.0);
    limits.nearestNeighborTorsionsCutoff = LLKAInternal::D2R(28.0);
    limits.totalDistanceCutoff = LLKAInternal::D2R(60.0);
    limits.pseudorotationCutoff = LLKAInternal::D2R(72.0);
    limits.minimumClusterVotes = 0.001111;
    limits.minimumNearestNeighbors = 7;
    limits.numberOfUsedNearestNeighbors = 11;
    const double MAX_CLOSE_ENOUGH_RMSD = 0.5;

    const auto &cl = resCl.success();
    const auto &gs = resGs.success();
    const auto &cf = resCf.success();
    const auto &nu = resNu.success();
    const auto &cp = resCp.success();

    auto resCtx = LLKA::initializeClassificationContext(
        cl, gs, cf, nu, cp, limits, MAX_CLOSE_ENOUGH_RMSD
    );

    EFF_expect(resCtx.isSuccess(), true, "failed to initialize classification context " + LLKA::errorToString(resCtx.failure()));

    return resCtx.success();
}

static
auto classification_testClassifySimple(const LLKA::ClassificationContext &ctx)
{
    using VT = decltype(LLKA_CLASSIFICATION_OK);

    auto stru = LLKA::makeStructure(
        REAL_1BNA_A_1_2_ATOMS,
        REAL_1BNA_A_1_2_ATOMS_LEN
    );

    auto res = LLKA::classifyStep(stru, ctx);
    EFF_expect(res.isSuccess(), true, "unable to classify step " + LLKA::errorToString(res.failure()));
    const auto &classifiedStep = res.success();

    EFF_expect(VT(classifiedStep.violations), LLKA_CLASSIFICATION_OK, "detected classification violations where there were not supposed to be any");
    EFF_expect(classifiedStep.assignedNtC, LLKA_BB04, "wrong NtC class");
    EFF_expect(classifiedStep.closestGoldenStep, "1mjo_F_DG_10_DA_11", "wrong golden step name");
    EFF_cmpFlt(classifiedStep.rmsdToClosestNtC, 0.366977877615573, "RMSD to closest NtC exceeds tolerance");

    EFF_cmpFlt(classifiedStep.confalScore.CC, 80.5998065953948, "wrong CC confal score");
    EFF_cmpFlt(classifiedStep.confalScore.NN, 99.9952502902733, "wrong NN confal score");
    EFF_cmpFlt(classifiedStep.confalScore.total, 32.9484362098594, "wrong total confal score");
}

static
auto classification_testClassifyThorough(const LLKA::ClassificationContext &ctx)
{
    using VT = decltype(LLKA_CLASSIFICATION_OK);

    auto stru = LLKA::makeStructure(
        REAL_3VOK_U_1_2_ATOMS,
        REAL_3VOK_U_1_2_ATOMS_LEN
    );

    auto res = LLKA::classifyStep(stru, ctx);
    EFF_expect(res.isSuccess(), true, "unable to classify step " + LLKA::errorToString(res.failure()));
    const auto &classifiedStep = res.success();

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

    auto avg = LLKA::averageConfal({ classifiedStep }, ctx);
    EFF_cmpFlt(avg.score, 72.372976109701284, "wrong average confal score");
    EFF_cmpFlt(avg.percentile, 94.599999999999994, "wrong confal percentile");
}

static
auto classification_testClassifyViolations(const LLKA::ClassificationContext &ctx)
{
    auto stru = LLKA::makeStructure(
        REAL_1DK1_B_5_6_ATOMS,
        REAL_1DK1_B_5_6_ATOMS_LEN
    );

    auto res = LLKA::classifyStep(stru, ctx);
    EFF_expect(res.isSuccess(), true, "unable to classify step " + LLKA::errorToString(res.failure()));
    const auto &classifiedStep = res.success();

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
}

static
auto minicif_to_string()
{
    auto res = LLKA::cifToStructure(std::filesystem::path{"./1BNA.cif"}, LLKA_MINICIF_GET_CIFDATA);
    EFF_expect(res.isSuccess(), true, "unexpected return value " + LLKA::errorToString(res.failure().tRet) + ", " + res.failure().error)
    const auto &cifData = res.success().cifData;

    auto res2 = LLKA::cifDataToString(cifData, true);
    EFF_expect(res2.isSuccess(), true, "unexpected return value " + LLKA::errorToString(res.failure().tRet) + ", " + res.failure().error);
}

static
auto minicif_structure()
{
    auto res = LLKA::cifToStructure(std::filesystem::path{"./1BNA.cif"});
    EFF_expect(res.isSuccess(), true, "unexpected return value " + LLKA::errorToString(res.failure().tRet) + ", " + res.failure().error)
    const auto &importedStru = res.success();
    EFF_expect(importedStru.structure.size(), 566UL, "wrong number of atoms in structure");
    EFF_expect(importedStru.structure[165].label_seq_id, 9, "unexpected label_seq_id");

    auto res2 = LLKA::cifToStructure(std::filesystem::path{"./1BNA_broken.cif"});
    EFF_expect(res2.isSuccess(), false, "CIF parser succeeded when it should not have")
    const auto &fail = res2.failure();
    EFF_expect(fail.tRet, LLKA_E_BAD_DATA, "wrong error code")
    EFF_expect(fail.error, "Malformed loop on line 1089, unexpected token kind 8", "unexpected error message")
}

static
auto resource_loaders_load()
{
    auto res = LLKA::loadGoldenSteps(std::filesystem::path{"./golden_steps.csv"});
    EFF_expect(res.isSuccess(), true, "failed to load golden_steps.csv: " + LLKA::errorToString(res.failure()));
    const auto &resource = res.success();

    // NOTE: This may need adjusting if the golden steps definitions change
    EFF_expect(resource.empty(), false, "golden_steps are not supposed to be empty");
    EFF_expect(resource.front().clusterNumber, 901, "wrong cluster number");
}

#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

static
auto splitting_testSplitAltIds()
{
    auto stru = LLKA::makeStructure(REAL_1DK1_B_26_28_ATOMS, REAL_1DK1_B_26_28_ATOMS_LEN);

    auto splitted = LLKA::splitByAltIds(stru);
    EFF_expect(splitted.size(), 2UL, "wrong number of splitted structures");
    EFF_expect(splitted[0].altId, 'A', "wrong alternate position id");
    EFF_expect(splitted[1].altId, 'B', "wrong alternate position id");

    EFF_expect(splitted[0].structure.front().label_alt_id, LLKA_NO_ALTID, "wrong alternate position id");
    EFF_expect(splitted[0].structure.back().label_alt_id, 'A', "wrong alternate position id");
    EFF_expect(splitted[1].structure.back().label_alt_id, 'B', "wrong alternate position id");
}

static
auto splitting_testSplitDinucleotides()
{
    auto stru = LLKA::makeStructure(REAL_1DK1_B_26_28_ATOMS, REAL_1DK1_B_26_28_ATOMS_LEN);

    auto res = LLKA::splitStructureToDinucleotideSteps(stru);
    EFF_expect(res.isSuccess(), true, "unexpected return value " + LLKA::errorToString(res.failure()))
    const auto &steps = res.success();

    EFF_expect(steps.size(), 3UL, "wrong number of steps")
    EFF_expect(steps[0][0].label_seq_id, 26, "wrong sequence id");
    EFF_expect(steps[2].back().label_seq_id, 28, "wrong sequence id");
}

auto main() -> int
{
#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
    minicif_structure();
    minicif_to_string();

    resource_loaders_load();

    auto ctx = classification_initContext();
    classification_testClassifySimple(ctx);
    classification_testClassifyThorough(ctx);
    classification_testClassifyViolations(ctx);
#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

    splitting_testSplitAltIds();
    splitting_testSplitDinucleotides();
}

