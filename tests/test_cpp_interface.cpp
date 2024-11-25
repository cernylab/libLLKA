// vim: set sw=4 ts=4 sts=4 expandtab :

#include <llka_cpp.h>

#include "effedup.hpp"
#include "testing_structures.h"

#include <util/geometry.h>
#include <util/geometry_cpp.h>
#include <util/printers.hpp>

#include <array>
#include <cmath>
#include <utility>

static
auto NtC_testAllMetrics()
{
    double EXPECTED_TORSIONS[] = {
        2.11367168175336,
        3.03091618163782,
       -1.54423502253339,
       -0.988174041115517,
       -3.127133069283356,
        0.911450188676453,
        1.72629995189555,
       -2.13194310890551,
       -2.22249417305504
    };
    auto stru = LLKA::makeStructure(REAL_1BNA_A_6_7_ATOMS, REAL_1BNA_A_6_7_ATOMS_LEN);

    auto tRet = LLKA::extractMetricsStructure(stru);
    EFF_expect(tRet.isSuccess(), true, "LLKA_extractMetricsStructure() returned unexpected value " + LLKA::errorToString(tRet.failure()));
    const auto &allMetricsStru = tRet.success();

    // Test torsions
    std::underlying_type_t<LLKA_DinucleotideTorsion> tor = LLKA_TOR_DELTA_1;
    for (; tor <= LLKA_TOR_CHI_2; tor++) {
        const auto tTor = static_cast<LLKA_DinucleotideTorsion>(tor);
        auto tRet2 = LLKA::dinucleotideTorsion(tTor, allMetricsStru);
        EFF_expect(tRet2.isSuccess(), true, "LLKA_dinucleotideTorsion() returned unexpected value " + LLKA::errorToString(tRet2.failure()));
        const auto &metricStru = tRet2.success();

        auto ang = LLKA::measureDihedral<double>(metricStru[0], metricStru[1], metricStru[2], metricStru[3]);
        std::cout << LLKA::dinucleotideTorsionName(tTor, LLKA_TRUE) << " angle = " << ang << " (" << ang * 180.0 / M_PI << ")\n";
        std::cout << LLKA::dinucleotideTorsionName(tTor, LLKA_FALSE) << " atoms\n" << metricStru << "\n";

        EFF_cmpFlt(ang, EXPECTED_TORSIONS[tor], LLKA_dinucleotideTorsionName(tTor, LLKA_FALSE))
    }

    // Test crossresidue torsions and distances
    auto tRet3 = LLKA::crossResidueMetric(LLKA_XR_DIST_CC, allMetricsStru);
    EFF_expect(tRet3.isSuccess(), true, "LLKA_crossResidueMetric() returned unexpected value " + LLKA::errorToString(tRet3.failure()))
    auto metricStru = tRet3.success();
    auto dist = LLKA::measureDistance<double>(metricStru[0], metricStru[1]);
    std::cout << LLKA::crossResidueMetricName(LLKA_XR_DIST_CC, LLKA_TRUE) << " distance " << dist << "\n";
    EFF_cmpFlt(dist, 4.78625636588764, LLKA_crossResidueMetricName(LLKA_XR_DIST_CC, LLKA_FALSE))

    tRet3 = LLKA::crossResidueMetric(LLKA_XR_DIST_NN, allMetricsStru);
    EFF_expect(tRet3.isSuccess(), true, "LLKA_crossResidueMetric() returned unexpected value " + LLKA::errorToString(tRet3.failure()))
    metricStru = tRet3.success();
    dist = LLKA::measureDistance<double>(metricStru[0], metricStru[1]);
    std::cout << LLKA::crossResidueMetricName(LLKA_XR_DIST_NN, LLKA_TRUE) << " distance " << dist << "\n";
    EFF_cmpFlt(dist, 4.22325786567669, LLKA_crossResidueMetricName(LLKA_XR_DIST_NN, LLKA_FALSE))

    tRet3 = LLKA::crossResidueMetric(LLKA_XR_TOR_MU, allMetricsStru);
    EFF_expect(tRet3.isSuccess(), true, "LLKA_crossResidueMetric() returned unexpected value " + LLKA::errorToString(tRet3.failure()))
    metricStru = tRet3.success();
    auto ang = LLKA::measureDihedral<double>(metricStru[0], metricStru[1], metricStru[2], metricStru[3]);
    EFF_cmpFlt(ang, 0.45071211146595, LLKA_crossResidueMetricName(LLKA_XR_TOR_MU, LLKA_FALSE))
    std::cout << LLKA::crossResidueMetricName(LLKA_XR_TOR_MU, LLKA_TRUE) << " torsion angle " << ang << " (" << (180.0 * ang / M_PI) << ")\n";
}

static
auto NtC_testExtractExtendedBackbone()
{
    std::array<std::tuple<std::string, int32_t>, 18> EXTENDED_EXPECTED_ORDER{{
        {"C5'", 1},
        {"C4'", 1},
        {"O4'", 1},
        {"C3'", 1},
        {"O3'", 1},
        {"C1'", 1},
        {"XXX", 1}, // Either N1 or N1, depends on the base
        {"XXX", 1}, // Either C4 or C3, depends on the base
        {"P",   2},
        {"O5'", 2},
        {"C5'", 2},
        {"C4'", 2},
        {"O4'", 2},
        {"C3'", 2},
        {"O3'", 2},
        {"C1'", 2},
        {"XXX", 2}, // Either N1 or N1, depends on the base
        {"XXX", 2}  // Either C4 or C3, depends on the base
    }};

    auto stru = LLKA::makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);

    auto siRet = LLKA::structureIsStep(stru);
    EFF_expect(siRet.isSuccess(), true, "structureIsStep returned false")
    const auto &info = siRet.success();

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

    auto tRet = LLKA::extractExtendedBackbone(stru);
    EFF_expect(tRet.isSuccess(), true, "extractExtendedBackbone() returned unexpected value " + LLKA::errorToString(tRet.failure()));
    const auto &backbone = tRet.success();

    EFF_expect(backbone.size(), EXTENDED_EXPECTED_ORDER.size(), "number of backbone atoms is wrong");

    for (size_t idx = 0; idx < backbone.size(); idx++) {
        const auto &atom = backbone[idx];
        const auto &[name, seqId] = EXTENDED_EXPECTED_ORDER[idx];
        EFF_expect(std::string{atom.label_atom_id}, name, "label_atom_id is wronh");
        EFF_expect(atom.label_seq_id, seqId, "label_seq_id is wrong");
    }
}

static
auto NtC_testExtractPlainBackbone()
{
    static const std::array<std::tuple<std::string, int32_t>, 12> EXPECTED_ORDER{{
        {"C5'", 1},
        {"C4'", 1},
        {"O4'", 1},
        {"C3'", 1},
        {"O3'", 1},
        {"P",   2},
        {"O5'", 2},
        {"C5'", 2},
        {"C4'", 2},
        {"O4'", 2},
        {"C3'", 2},
        {"O3'", 2},
    }};

    auto stru = LLKA::makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);

    auto tRet = LLKA::extractBackbone(stru);
    EFF_expect(tRet.isSuccess(), true,  "extractBackbone() returned unexpected value " + LLKA::errorToString(tRet.failure()))
    const auto &backbone = tRet.success();

    EFF_expect(backbone.size(), EXPECTED_ORDER.size(), "number of backbone atoms is wrong");

    for (size_t idx = 0; idx < backbone.size(); idx++) {
        const auto &atom = backbone[idx];
        const auto &[name, seqId] = EXPECTED_ORDER[idx];
        EFF_expect(std::string{atom.label_atom_id}, name, "atom_id");
        EFF_expect(atom.label_seq_id, seqId, "seq_id");
    }
}

static
auto structure_testComparison()
{
    LLKA::Atom atomA = LLKA::makeAtom(
        "C", "C5", "1", "DA", "A",
        {}, {}, {},
        { 10, 11, 12 },
        1,
        1, 1,
        1,
        LLKA_NO_INSCODE,
        LLKA_NO_ALTID
    );

    LLKA::Atom atomAA = LLKA::makeAtom(
        "C", "C5", "1", "DA", "A",
        {}, {}, {},
        { 10, 11, 12 },
        2,
        1, 1,
        1,
        LLKA_NO_INSCODE,
        LLKA_NO_ALTID
    );

    LLKA::Atom atomB = LLKA::makeAtom(
        "C", "C4'", "1", "DA", "A",
        {}, {}, {},
        { 10, 11, 12 },
        2,
        1, 1,
        1,
        LLKA_NO_INSCODE,
        LLKA_NO_ALTID
    );

    EFF_expect(LLKA::atomsEqual(atomA, atomAA), false, "check if atoms are truly equal (does not ignore ID)")
    EFF_expect(LLKA::compareAtoms(atomA, atomAA, false), false , "compare identical atoms")
    EFF_expect(LLKA::compareAtoms(atomA, atomAA, true), true, "compare idential atoms (ignore ID)")
}


#define CHECK_POINT(pt, ref) \
    do { \
        EFF_cmpFlt(pt.x, ref.x, "Point X coordinate is wrong"); \
        EFF_cmpFlt(pt.y, ref.y, "Point Y coordinate is wrong"); \
        EFF_cmpFlt(pt.z, ref.z, "Point Z coordinate is wrong"); \
    } while (false)

static
auto superposition_testShifted()
{
    LLKA::Points a{
        { 1, 1, 1 },
        { 2, 2, 2 },
        { 3, 3, 3 },
        { 4, 4, 4 }
    };
    LLKA::Points b{
        { 10 + 1, 10 + 1, 10 + 1 },
        { 10 + 2, 10 + 2, 10 + 2 },
        { 10 + 3, 10 + 3, 10 + 3 },
        { 10 + 4, 10 + 4, 10 + 4 }
    };

    auto tRet = LLKA::superpose(a, b);
    EFF_expect(tRet.isSuccess(), true, "LLKA_superposePoints() returned unexpected value " + LLKA::errorToString(tRet.failure()))
    auto rmsd = tRet.success();

    CHECK_POINT(a[0], b[0]);
    CHECK_POINT(a[1], b[1]);
    CHECK_POINT(a[2], b[2]);
    CHECK_POINT(a[3], b[3]);
    EFF_cmpFlt(rmsd, 0.0, "RMSD is wrong")
}

static
auto superposition_testStructure()
{
    LLKA::Structure ref_AB01 = LLKA::makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA::Structure real_AB01 = LLKA::makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);

    auto tRet = LLKA::extractBackbone(ref_AB01);
    EFF_expect(tRet.isSuccess(), true, "LLKA_extractBackbone() returned unexpected value " + LLKA::errorToString(tRet.failure()))
    auto ref_AB01_bkbn = tRet.success();

    tRet = LLKA::extractBackbone(real_AB01);
    EFF_expect(tRet.isSuccess(), true, "LLKA_extractBackbone() returned unexpected value " + LLKA::errorToString(tRet.failure()));
    auto real_AB01_bkbn = tRet.success();

    auto tRet2 = LLKA::superpose(real_AB01_bkbn, ref_AB01_bkbn);
    EFF_expect(tRet2.isSuccess(), true, "LLKA_superposePoints() returned unexpected value " + LLKA::errorToString(tRet2.failure()))
    auto rmsd = tRet2.success();

    LLKA::Points expected{
        { 0.950860168366356, -4.10490242639581, 1.88352814333396 },
        { 0.50054990344796, -3.69993907500445, 0.514678133530478 },
        { -0.347244679657512, -4.71592059085125, -0.00580232173625572 },
        { -0.318861255763144, -2.43716579580966, 0.473326017357792 },
        { 0.506320258152502, -1.32662194302102, 0.113179074950623 },
        { -0.0811260469759572, 0.154118620506023, 0.0681432878857751},
        { -0.647318283320744, 0.253840078930681, -1.40665838767353 },
        { 0.237737393186242, 0.145726060576478, -2.52768907793944 },
        { -0.603024580420967, 0.282861734520597, -3.76182161161158 },
        { -1.56591988716846, -0.783314033669751, -3.80479671565305 },
        { -1.41697603269085, 1.56517005180424, -3.86385794226367 },
        { -1.71799695715546, 1.88014731841389, -5.20722860018113 }
    };

    for (size_t idx = 0; idx < 12; idx++)
        CHECK_POINT(real_AB01_bkbn[idx].coords, expected[idx]);

    EFF_cmpFlt(rmsd, 0.147657613575976, "RMSD is wrong");
}

static
auto connectivity_testConnectivity()
{
    // Test based on results for 1BNA molecule
    const LLKA::Structure realFirst = LLKA::makeStructure(REAL_1BNA_A_2_3_ATOMS, REAL_1BNA_A_2_3_ATOMS_LEN);
    const LLKA::Structure realSecond = LLKA::makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    const LLKA::Structure refFirst = LLKA::makeStructure(REF_BA01_ATOMS, REF_BA01_ATOMS_LEN);
    const LLKA::Structure refSecond = LLKA::makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);

    auto tRet = LLKA::measureStepConnectivity(realFirst, refFirst, realSecond, refSecond);
    EFF_expect(tRet.isSuccess(), true, "measureStepConnectivity() returned unexpected value " + LLKA::errorToString(tRet.failure()));
    const auto &conn = tRet.success();

    EFF_cmpFlt(conn.C5PrimeDistance, 0.0778497140143741, "wrong C5' distance")
    EFF_cmpFlt(conn.O3PrimeDistance, 0.110550489844076, "wrong O3' distance")
}

static
auto connectivity_testConnectivityMultiple()
{
    // Test based on results for 1BNA molecule
    const LLKA::Structure realFirst = LLKA::makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    const LLKA::Structure realSecond = LLKA::makeStructure(REAL_1BNA_A_4_5_ATOMS, REAL_1BNA_A_4_5_ATOMS_LEN);

    const LLKA::Structure refFirst = LLKA::makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);

    const LLKA::Structures refsSecond{
        LLKA::makeStructure(REF_BA01_ATOMS, REF_BA01_ATOMS_LEN),
        LLKA::makeStructure(REF_BB00_ATOMS, REF_BB00_ATOMS_LEN),
        LLKA::makeStructure(REF_BB04_ATOMS, REF_BB04_ATOMS_LEN)
    };

    auto tRet = LLKA::measureStepConnectivity(realFirst, refFirst, realSecond, refsSecond);
    EFF_expect(tRet.isSuccess(), true, "measureStepConnectivity() returned unexpected value " + LLKA::errorToString(tRet.failure()))
    const auto &results = tRet.success();

    EFF_cmpFlt(results[0].C5PrimeDistance, 0.222853768522876, "wrong C5' distance");
    EFF_cmpFlt(results[0].O3PrimeDistance, 0.332800460520272, "wrong O3' distance")

    EFF_cmpFlt(results[1].C5PrimeDistance, 0.398108806100796, "wrong C5' distance");
    EFF_cmpFlt(results[1].O3PrimeDistance, 0.29940718941947, "wrong O3' distance")

    EFF_cmpFlt(results[2].C5PrimeDistance, 0.379975441318778, "wrong C5' distance");
    EFF_cmpFlt(results[2].O3PrimeDistance, 0.372926953223463, "wrong O3' distance")
}


static
auto nucleotide_testExtractRibose()
{
    auto stru = LLKA::makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);

    auto nucl = LLKA::extractNucleotide(stru, 1, "A", 3);

    auto tRet = LLKA::extractRibose(nucl);
    EFF_expect(tRet.isSuccess(), true, "LLKA_extractRibose returned unexpected value " + LLKA::errorToString(tRet.failure()));
    const auto &riboseStru = tRet.success();
    EFF_expect(riboseStru.size(), 5UL, "wrong ribose core atoms");

    const std::array<std::string, 5> names{ "C4'", "O4'", "C1'", "C2'", "C3'" };
    for (size_t idx = 0; idx < 5; idx++)
        EFF_expect(riboseStru[idx].label_atom_id, names[idx], "wrong atom name");
}

static
auto nucleotide_testExtractNucleotide()
{
    auto stru = LLKA::makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);

    auto firstNucl = LLKA::extractNucleotide(stru, 1, "A", 3);
    EFF_expect(firstNucl.size(), 19UL, "wrong number of atoms in nucleotide");
    EFF_expect(firstNucl.front().label_seq_id, firstNucl.back().label_seq_id, "sequence id of first and last atom in nucleotide do not match");

    auto secondNucl = LLKA::extractNucleotide(stru, 1, "A", 4);
    EFF_expect(secondNucl.size(), 22UL, "wrong number of atoms in nucleotide");
    EFF_expect(secondNucl.front().label_seq_id, secondNucl.back().label_seq_id, "sequence id of first and last atom in nucleotide do not match");

    EFF_expect(firstNucl.front().label_seq_id + 1 == secondNucl.front().label_seq_id, true , "First atoms of two different nucleotides have the same sequence id");
}

static
auto similarity_testSimilarity()
{
    LLKA::Structure stru = LLKA::makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);

    auto tRet = LLKA::measureStepSimilarity(stru, LLKA_AB01);
    EFF_expect(tRet.isSuccess(), true, "LLKA::measureStepSimilarity() returned unexpected value " + LLKA::errorToString(tRet.failure()));
    const auto &similar = tRet.success();

    EFF_cmpFlt(similar.rmsd, 0.176588293714412, "wrong RMSD")
    EFF_cmpFlt(similar.euclideanDistance, 25.0592315110376, "wrong euclidean distance")
}

static
auto similarity_testSimilarityMultiple()
{
    const std::vector<LLKA_NtC> ntcs = { LLKA_AB01, LLKA_AB03, LLKA_BB00 };
    LLKA::Structure stru = LLKA::makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);

    auto tRet = LLKA::measureStepSimilarity(stru, ntcs);
    EFF_expect(tRet.isSuccess(), true, "LLKA::measureStepSimilarity() returned unexpected value " + LLKA::errorToString(tRet.failure()))
    const auto &results = tRet.success();

    EFF_cmpFlt(results[0].rmsd, 0.176588293714412, "wrong RMSD")
    EFF_cmpFlt(results[0].euclideanDistance, 25.0592315110376, "wrong euclidean distance")

    EFF_cmpFlt(results[1].rmsd, 0.24827518856867, "wrong RMSD")
    EFF_cmpFlt(results[1].euclideanDistance, 46.3910103838844, "wrong euclidean distance")

    EFF_cmpFlt(results[2].rmsd, 0.348261210522226, "wrong RMSD")
    EFF_cmpFlt(results[2].euclideanDistance, 56.5297526013971, "wrong euclidean distance")
}

auto main() -> int
{
    structure_testComparison();

    superposition_testShifted();
    superposition_testStructure();

    NtC_testExtractPlainBackbone();
    NtC_testExtractExtendedBackbone();
    NtC_testAllMetrics();

    connectivity_testConnectivity();
    connectivity_testConnectivityMultiple();

    nucleotide_testExtractNucleotide();
    nucleotide_testExtractRibose();

    similarity_testSimilarity();
    similarity_testSimilarityMultiple();

    return EXIT_SUCCESS;
}
