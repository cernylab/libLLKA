/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_connectivity_similarity.h>
#include <llka_measurements.h>

#include "effedup.hpp"
#include "testing_structures.h"

#include <util/geometry.h>
#include <util/printers.hpp>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <type_traits>

static
auto prnSimilarity(const std::string &desc, const LLKA_Similarity &similar)
{
    std::cout
        << desc << "\n"
        << "Euclidean distance: " << similar.euclideanDistance << "\n"
        << "RMSD: " << similar.rmsd << "\n";
}

static
auto testShouldFail(const LLKA_Structure &stru)
{
    LLKA_RetCode tRet;
    LLKA_Structure dinucleotide;
    LLKA_Structure torsion;

    tRet = LLKA_extractExtendedBackbone(&stru, &dinucleotide);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value")

    tRet = LLKA_dinucleotideTorsion(LLKA_TOR_CHI_1, &dinucleotide, &torsion);
    EFF_expect(tRet, LLKA_E_INVALID_ARGUMENT, "LLKA_NtCTorsion returned unexpected value")

    LLKA_destroyStructure(&dinucleotide);
}

static
auto testAllMetrics(const LLKA_Structure &stru)
{
    LLKA_RetCode tRet;
    LLKA_Structure allMetricsStru;
    LLKA_Structure metricStru;
    double EXPECTED_TORSIONS[] = {
        2.11367168175336,
        3.03091618163782,
       -1.54423502253339,
       -0.98817404111552,
       -3.12713306928336,
        0.91145018867645,
        1.72629995189555,
       -2.13194310890551,
       -2.22249417305504
    };

    tRet = LLKA_extractMetricsStructure(&stru, &allMetricsStru);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractMetricsStructure() returned unexpected value")

    // Test torsions
    std::underlying_type_t<LLKA_DinucleotideTorsion> tor = LLKA_TOR_DELTA_1;
    for (; tor <= LLKA_TOR_CHI_2; tor++) {
        const auto tTor = static_cast<LLKA_DinucleotideTorsion>(tor);
        tRet = LLKA_dinucleotideTorsion(tTor, &allMetricsStru, &metricStru);
        EFF_expect(tRet, LLKA_OK, "LLKA_dinucleotideTorsion() returned unexpected value")

        auto ang = LLKA_measureDihedral(&metricStru.atoms[0], &metricStru.atoms[1], &metricStru.atoms[2], &metricStru.atoms[3]);
        std::cout << LLKA_dinucleotideTorsionName(tTor, LLKA_TRUE) << " angle = " << ang << " (" << ang * 180.0 / M_PI << ")\n";
        std::cout << LLKA_dinucleotideTorsionName(tTor, LLKA_FALSE) << " atoms\n" << metricStru << "\n";

        EFF_cmpFlt(ang, EXPECTED_TORSIONS[tor], LLKA_dinucleotideTorsionName(tTor, LLKA_FALSE))

        LLKA_destroyStructure(&metricStru);
    }

    // Test crossresidue torsions and distances
    tRet = LLKA_crossResidueMetric(LLKA_XR_DIST_CC, &allMetricsStru, &metricStru);
    EFF_expect(tRet, LLKA_OK, "LLKA_crossResidueMetric() returned unexpected value")
    auto dist = LLKA_measureDistance(&metricStru.atoms[0], &metricStru.atoms[1]);
    std::cout << LLKA_crossResidueMetricName(LLKA_XR_DIST_CC, LLKA_TRUE) << " distance " << dist << "\n";
    EFF_cmpFlt(dist, 4.78625636588764, LLKA_crossResidueMetricName(LLKA_XR_DIST_CC, LLKA_FALSE))
    LLKA_destroyStructure(&metricStru);

    tRet = LLKA_crossResidueMetric(LLKA_XR_DIST_NN, &allMetricsStru, &metricStru);
    EFF_expect(tRet, LLKA_OK, "LLKA_crossResidueMetric() returned unexpected value")
    dist = LLKA_measureDistance(&metricStru.atoms[0], &metricStru.atoms[1]);
    std::cout << LLKA_crossResidueMetricName(LLKA_XR_DIST_NN, LLKA_TRUE) << " distance " << dist << "\n";
    EFF_cmpFlt(dist, 4.22325786567669, LLKA_crossResidueMetricName(LLKA_XR_DIST_NN, LLKA_FALSE))
    LLKA_destroyStructure(&metricStru);

    tRet = LLKA_crossResidueMetric(LLKA_XR_TOR_MU, &allMetricsStru, &metricStru);
    EFF_expect(tRet, LLKA_OK, "LLKA_crossResidueMetric() returned unexpected value")
    auto ang = LLKA_measureDihedral(&metricStru.atoms[0], &metricStru.atoms[1], &metricStru.atoms[2], &metricStru.atoms[3]);
    EFF_cmpFlt(ang, 0.45071211146595, LLKA_crossResidueMetricName(LLKA_XR_TOR_MU, LLKA_FALSE))
    std::cout << LLKA_crossResidueMetricName(LLKA_XR_TOR_MU, LLKA_TRUE) << " torsion angle " << ang << " (" << (180.0 * ang / M_PI) << ")\n";
    LLKA_destroyStructure(&metricStru);

    LLKA_destroyStructure(&allMetricsStru);
}

static
auto testAllMetrics2()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1BNA_B_16_17_ATOMS, REAL_1BNA_B_16_17_ATOMS_LEN);
    LLKA_StepMetrics metrics;

    EFF_expect(LLKA_calculateStepMetrics(&stru, &metrics), LLKA_OK, "LLKA_calculateStepMetrics() returned unexpected value")
    std::cout << "\n" << metrics << "\n";

    EFF_cmpFlt(metrics.delta_1,   2.37223786208374, LLKA_dinucleotideTorsionName(LLKA_TOR_DELTA_1, LLKA_FALSE))
    EFF_cmpFlt(metrics.epsilon_1, 3.03929251506922, LLKA_dinucleotideTorsionName(LLKA_TOR_EPSILON_1, LLKA_FALSE))
    EFF_cmpFlt(metrics.zeta_1,   -1.71720034870345, LLKA_dinucleotideTorsionName(LLKA_TOR_ZETA_1, LLKA_FALSE))
    EFF_cmpFlt(metrics.alpha_2,  -0.98818512153652, LLKA_dinucleotideTorsionName(LLKA_TOR_ALPHA_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.beta_2,   -2.95903498082587, LLKA_dinucleotideTorsionName(LLKA_TOR_BETA_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.gamma_2,   0.93968862228876, LLKA_dinucleotideTorsionName(LLKA_TOR_GAMMA_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.delta_2,   2.55870427626711, LLKA_dinucleotideTorsionName(LLKA_TOR_DELTA_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.chi_1,    -2.00305820661614, LLKA_dinucleotideTorsionName(LLKA_TOR_CHI_1, LLKA_FALSE))
    EFF_cmpFlt(metrics.chi_2,    -1.85775524738115, LLKA_dinucleotideTorsionName(LLKA_TOR_CHI_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.CC,        5.04421341737243, LLKA_crossResidueMetricName(LLKA_XR_DIST_CC, LLKA_FALSE))
    EFF_cmpFlt(metrics.NN,        4.45444407305783, LLKA_crossResidueMetricName(LLKA_XR_DIST_NN, LLKA_FALSE))
    EFF_cmpFlt(metrics.mu,        0.58047587191140, LLKA_crossResidueMetricName(LLKA_XR_TOR_MU, LLKA_FALSE))

    LLKA_destroyStructure(&stru);
}

static
auto testAllMetricsDifferenceAgainstReference()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1BNA_B_16_17_ATOMS, REAL_1BNA_B_16_17_ATOMS_LEN);
    LLKA_StepMetrics metrics;

    EFF_expect(LLKA_calculateStepMetricsDifferenceAgainstReference(&stru, LLKA_BB00, &metrics), LLKA_OK, "LLKA_calculateStepMetricsDifferenceAgainstReference returned unexpected value")
    std::cout << "\n" << metrics << "\n";

    EFF_cmpFlt(metrics.delta_1,   -0.033174913014841, LLKA_dinucleotideTorsionName(LLKA_TOR_DELTA_1, LLKA_FALSE))
    EFF_cmpFlt(metrics.epsilon_1, -0.156579878257593, LLKA_dinucleotideTorsionName(LLKA_TOR_EPSILON_1, LLKA_FALSE))
    EFF_cmpFlt(metrics.zeta_1,     0.059893895677181, LLKA_dinucleotideTorsionName(LLKA_TOR_ZETA_1, LLKA_FALSE))
    EFF_cmpFlt(metrics.alpha_2,   -0.005913818514107, LLKA_dinucleotideTorsionName(LLKA_TOR_ALPHA_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.beta_2,     0.189713522697103, LLKA_dinucleotideTorsionName(LLKA_TOR_BETA_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.gamma_2,    0.167729494131665, LLKA_dinucleotideTorsionName(LLKA_TOR_GAMMA_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.delta_2,    0.147706447562145, LLKA_dinucleotideTorsionName(LLKA_TOR_DELTA_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.chi_1,     -0.127701925348230, LLKA_dinucleotideTorsionName(LLKA_TOR_CHI_1, LLKA_FALSE))
    EFF_cmpFlt(metrics.chi_2,     -0.081882733476915, LLKA_dinucleotideTorsionName(LLKA_TOR_CHI_2, LLKA_FALSE))
    EFF_cmpFlt(metrics.CC,         0.097213417372426, LLKA_crossResidueMetricName(LLKA_XR_DIST_CC, LLKA_FALSE))
    EFF_cmpFlt(metrics.NN,         0.081444073057826, LLKA_crossResidueMetricName(LLKA_XR_DIST_NN, LLKA_FALSE))
    EFF_cmpFlt(metrics.mu,         0.134718780952045, LLKA_crossResidueMetricName(LLKA_XR_TOR_MU, LLKA_FALSE))

    LLKA_destroyStructure(&stru);
}

    static
auto testSimilarity()
{
    LLKA_RetCode tRet;
    LLKA_Structure stru = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    LLKA_Similarity similar;

    tRet = LLKA_measureStepSimilarityNtC(&stru, LLKA_AB01, &similar);
    EFF_expect(tRet, LLKA_OK, "LLKA_measureStepSimilarity() returned unexpected value")

    prnSimilarity("1BNA 3_4 (AB01) vs. reference AB01", similar);
    EFF_cmpFlt(similar.rmsd, 0.176588293714412, "wrong RMSD")
    EFF_cmpFlt(similar.euclideanDistance, 25.0592315110376, "wrong euclidean distance")

    LLKA_destroyStructure(&stru);
}

static
auto testSimilarityMultiple()
{
    const LLKA_NtC ntcs[] = { LLKA_AB01, LLKA_AB03, LLKA_BB00, LLKA_INVALID_NTC };
    LLKA_Structure stru = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);

    LLKA_Similarities results{
        new LLKA_Similarity[3],
        3
    };

    auto tRet = LLKA_measureStepSimilarityNtCMultiple(&stru, ntcs, &results);
    EFF_expect(tRet, LLKA_OK, "measureStepSimilarityNtCMultiple() returned unexpected value")

    prnSimilarity("1BNA 3_4 (AB01) vs. reference AB01", results.similars[0]);
    EFF_cmpFlt(results.similars[0].rmsd, 0.176588293714412, "wrong RMSD")
    EFF_cmpFlt(results.similars[0].euclideanDistance, 25.0592315110376, "wrong euclidean distance")

    prnSimilarity("1BNA 3_4 (AB01) vs. reference AB03", results.similars[1]);
    EFF_cmpFlt(results.similars[1].rmsd, 0.24827518856867, "wrong RMSD")
    EFF_cmpFlt(results.similars[1].euclideanDistance, 46.3910103838844, "wrong euclidean distance")

    prnSimilarity("1BNA 3_4 (AB01) vs. reference BB00", results.similars[2]);
    EFF_cmpFlt(results.similars[2].rmsd, 0.348261210522226, "wrong RMSD")
    EFF_cmpFlt(results.similars[2].euclideanDistance, 56.5297526013971, "wrong euclidean distance")

    delete[] results.similars;
    LLKA_destroyStructure(&stru);
}

auto main(int, char**) -> int
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1BNA_A_6_7_ATOMS, REAL_1BNA_A_6_7_ATOMS_LEN);

    std::cout << std::setprecision(15);

    testShouldFail(stru);
    testAllMetrics(stru);
    testAllMetrics2();
    testAllMetricsDifferenceAgainstReference();
    testSimilarity();
    testSimilarityMultiple();

    LLKA_destroyStructure(&stru);

    return EXIT_SUCCESS;
}
