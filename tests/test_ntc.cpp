// vim: set sw=4 ts=4 sts=4 expandtab :

#include "effedup.hpp"
#include "testing_structures.h"

#include <llka_ntc.h>

static
auto testCANA()
{
    auto cana = LLKA_nameToCANA("AAw");
    EFF_expect(cana, LLKA_AAw, "wrong conversion from name to CANA item");

    cana = LLKA_nameToCANA("this is wrong");
    EFF_expect(cana, LLKA_INVALID_CANA, "wrong conversion from invalid CANA name to CANA item");

    std::string name = LLKA_CANAToName(LLKA_miB);
    EFF_expect(name, std::string{"miB"}, "wrong conversion from CANA item to name");

    name = LLKA_CANAToName(LLKA_INVALID_CANA);
    EFF_expect(name, std::string{"NAN"}, "wrong conversion from invalid CANA item to name");
}

static
auto testNtC()
{
    auto ntc = LLKA_nameToNtC("BA17");
    EFF_expect(ntc, LLKA_BA17, "wrong conversion from name to NtC item");

    ntc = LLKA_nameToNtC("this is wrong");
    EFF_expect(ntc, LLKA_INVALID_NTC, "wrong conversion from invalid NtC name to NtC item");

    std::string name = LLKA_NtCToName(LLKA_BB1S);
    EFF_expect(name, std::string{"BB1S"}, "wrong conversion from NtC item to name");

    name = LLKA_NtCToName(LLKA_INVALID_NTC);
    EFF_expect(name, std::string{"NANT"}, "wrong conversion from invalid NtC item to name");
}

static
auto testDifferenceAgainstReference()
{
    LLKA_Structure stru = LLKA_makeStructure(REAL_1BNA_A_1_2_ATOMS, REAL_1BNA_A_1_2_ATOMS_LEN);
    LLKA_StepMetrics metrics;

    auto tRet = LLKA_calculateStepMetricsDifferenceAgainstReference(&stru, LLKA_BB04, &metrics);
    EFF_expect(tRet, LLKA_OK, "LLKA_calculateStepMetricsDifferenceAgainstReference() returned unexpected value")
    EFF_cmpFlt(metrics.delta_1,   0.291304749269346, "delta_1 metric is wrong");
    EFF_cmpFlt(metrics.epsilon_1, 0.303056944684855, "epsilon_1 metric is wrong");
    EFF_cmpFlt(metrics.zeta_1,    0.032871663591650, "zeta_1 metric is wrong");
    EFF_cmpFlt(metrics.alpha_2,  -0.355335537090472, "alpha_2 metric is wrong");
    EFF_cmpFlt(metrics.beta_2,    0.299443860928775, "beta_2 metric is wrong");
    EFF_cmpFlt(metrics.gamma_2,  -0.104152469812837, "gamma_2 metric is wrong");
    EFF_cmpFlt(metrics.delta_2,  -0.206873055232490, "delta_2 metric is wrong");
    EFF_cmpFlt(metrics.chi_1,    -0.131009933521641, "chi_1 metric is wrong");
    EFF_cmpFlt(metrics.chi_2,    -0.053816152662312, "chi_2 metric is wrong");
    EFF_cmpFlt(metrics.CC,        0.228744666804201, "CC metric is wrong");
    EFF_cmpFlt(metrics.NN,        0.003706138611939, "NN metric is wrong");
    EFF_cmpFlt(metrics.mu,        0.064645999883097, "mu metric is wrong");

    LLKA_destroyStructure(&stru);
}

#define CHECK_QUAD(_a, _b, _c, _d, tor, quad) \
    EFF_expect(quad.a, std::string{_a}, std::string{"Atom a on torsion "} + LLKA_dinucleotideTorsionName(tor, false) + " is wrong") \
    EFF_expect(quad.b, std::string{_b}, std::string{"Atom b on torsion "} + LLKA_dinucleotideTorsionName(tor, false) + " is wrong") \
    EFF_expect(quad.c, std::string{_c}, std::string{"Atom c on torsion "} + LLKA_dinucleotideTorsionName(tor, false) + " is wrong") \
    EFF_expect(quad.d, std::string{_d}, std::string{"Atom d on torsion "} + LLKA_dinucleotideTorsionName(tor, false) + " is wrong")

#define CHECK_DINU_TOR_ATOMS_BASES(fb, sb, tor, quad, a, b, c, d) \
    tRet = LLKA_dinucleotideTorsionAtomsFromBases(fb, sb, tor, &quad); \
    EFF_expect(tRet, LLKA_OK, "LLKA_dinucleotideTorsionAtomsFromStructure() returned unexpected value") \
    CHECK_QUAD(a, b, c, d, tor, quad)

#define CHECK_XR_DIST_ATOMS_BASES(fb, sb, xr, quad, _a, _b) \
    tRet = LLKA_crossResidueMetricAtomsFromBases(fb, sb, xr, &quad); \
    EFF_expect(tRet, LLKA_OK, "LLKA_crossResidueMetricAtoms() returned unexpected value") \
    EFF_expect(quad.a, std::string{_a}, std::string{"Atom a on distance "} + LLKA_crossResidueMetricName(xr, false) + " is wrong") \
    EFF_expect(quad.b, std::string{_b}, std::string{"Atom b on distance "} + LLKA_crossResidueMetricName(xr, false) + " is wrong")

#define CHECK_XR_TOR_ATOMS_BASES(fb, sb, xr, quad, _a, _b, _c, _d) \
    tRet = LLKA_crossResidueMetricAtomsFromBases(fb, sb, xr, &quad); \
    EFF_expect(tRet, LLKA_OK, "LLKA_crossResidueMetricAtomsFromStructure() returned unexpected value") \
    EFF_expect(quad.a, std::string{_a}, std::string{"Atom a on torsion "} + LLKA_crossResidueMetricName(xr, false) + " is wrong") \
    EFF_expect(quad.b, std::string{_b}, std::string{"Atom b on torsion "} + LLKA_crossResidueMetricName(xr, false) + " is wrong") \
    EFF_expect(quad.c, std::string{_c}, std::string{"Atom c on torsion "} + LLKA_crossResidueMetricName(xr, false) + " is wrong") \
    EFF_expect(quad.d, std::string{_d}, std::string{"Atom d on torsion "} + LLKA_crossResidueMetricName(xr, false) + " is wrong")

static
auto testTorsionNamesBases()
{
    LLKA_AtomNameQuad quad;
    LLKA_RetCode tRet;

    CHECK_DINU_TOR_ATOMS_BASES("DC", "DG", LLKA_TOR_DELTA_1,   quad, "C5'", "C4'", "C3'", "O3'");
    CHECK_DINU_TOR_ATOMS_BASES("DC", "DG", LLKA_TOR_EPSILON_1, quad, "C4'", "C3'", "O3'", "P");
    CHECK_DINU_TOR_ATOMS_BASES("DC", "DG", LLKA_TOR_ZETA_1,    quad, "C3'", "O3'", "P", "O5'");
    CHECK_DINU_TOR_ATOMS_BASES("DC", "DG", LLKA_TOR_ALPHA_2,   quad, "O3'", "P", "O5'", "C5'");
    CHECK_DINU_TOR_ATOMS_BASES("DC", "DG", LLKA_TOR_BETA_2,    quad, "P", "O5'", "C5'", "C4'");
    CHECK_DINU_TOR_ATOMS_BASES("DC", "DG", LLKA_TOR_GAMMA_2,   quad, "O5'", "C5'", "C4'", "C3'");
    CHECK_DINU_TOR_ATOMS_BASES("DC", "DG", LLKA_TOR_DELTA_2,   quad, "C5'", "C4'", "C3'", "O3'");
    CHECK_DINU_TOR_ATOMS_BASES("DC", "DG", LLKA_TOR_CHI_1,     quad, "O4'", "C1'", "N1", "C2");
    CHECK_DINU_TOR_ATOMS_BASES("DC", "DG", LLKA_TOR_CHI_2,     quad, "O4'", "C1'", "N9", "C4");
    CHECK_XR_DIST_ATOMS_BASES("DC", "DG",  LLKA_XR_DIST_CC,    quad, "C2", "C4");
    CHECK_XR_DIST_ATOMS_BASES("DC", "DG",  LLKA_XR_DIST_NN,    quad, "N1", "N9");
    CHECK_XR_TOR_ATOMS_BASES("DC", "DG",   LLKA_XR_TOR_MU,     quad, "C2", "N1", "N9", "C4");

    CHECK_DINU_TOR_ATOMS_BASES("A", "PSU", LLKA_TOR_DELTA_1,   quad, "C5'", "C4'", "C3'", "O3'");
    CHECK_DINU_TOR_ATOMS_BASES("A", "PSU", LLKA_TOR_EPSILON_1, quad, "C4'", "C3'", "O3'", "P");
    CHECK_DINU_TOR_ATOMS_BASES("A", "PSU", LLKA_TOR_ZETA_1,    quad, "C3'", "O3'", "P", "O5'");
    CHECK_DINU_TOR_ATOMS_BASES("A", "PSU", LLKA_TOR_ALPHA_2,   quad, "O3'", "P", "O5'", "C5'");
    CHECK_DINU_TOR_ATOMS_BASES("A", "PSU", LLKA_TOR_BETA_2,    quad, "P", "O5'", "C5'", "C4'");
    CHECK_DINU_TOR_ATOMS_BASES("A", "PSU", LLKA_TOR_GAMMA_2,   quad, "O5'", "C5'", "C4'", "C3'");
    CHECK_DINU_TOR_ATOMS_BASES("A", "PSU", LLKA_TOR_DELTA_2,   quad, "C5'", "C4'", "C3'", "O3'");
    CHECK_DINU_TOR_ATOMS_BASES("A", "PSU", LLKA_TOR_CHI_1,     quad, "O4'", "C1'", "N9", "C4");
    CHECK_DINU_TOR_ATOMS_BASES("A", "PSU", LLKA_TOR_CHI_2,     quad, "O4'", "C1'", "C5", "C4");
    CHECK_XR_DIST_ATOMS_BASES("A", "PSU",  LLKA_XR_DIST_CC,    quad, "C4", "C4");
    CHECK_XR_DIST_ATOMS_BASES("A", "PSU",  LLKA_XR_DIST_NN,    quad, "N9", "C5");
    CHECK_XR_TOR_ATOMS_BASES("A", "PSU",   LLKA_XR_TOR_MU,     quad, "C4", "N9", "C5", "C4");

    tRet = LLKA_dinucleotideTorsionAtomsFromBases("Organic", "A", LLKA_TOR_DELTA_1, &quad);
    EFF_expect(tRet, LLKA_E_INVALID_ARGUMENT, std::string{"Expected LLKA_E_INVALID_ARGUMENT, got "} + LLKA_errorToString(tRet))

    tRet = LLKA_dinucleotideTorsionAtomsFromBases("A", "Unicorn", LLKA_TOR_DELTA_1, &quad);
    EFF_expect(tRet, LLKA_E_INVALID_ARGUMENT, std::string{"Expected LLKA_E_INVALID_ARGUMENT, got "} + LLKA_errorToString(tRet))
}

#define CHECK_DINU_TOR_ATOMS_STRU(stru, tor, quad, a, b, c, d) \
    tRet = LLKA_dinucleotideTorsionAtomsFromStructure(&stru, tor, &quad); \
    EFF_expect(tRet, LLKA_OK, "LLKA_dinucleotideTorsionAtomsFromStructure() returned unexpected value") \
    CHECK_QUAD(a, b, c, d, tor, quad)

#define CHECK_XR_DIST_ATOMS_STRU(stru, xr, quad, _a, _b) \
    tRet = LLKA_crossResidueMetricAtomsFromStructure(&stru, xr, &quad); \
    EFF_expect(tRet, LLKA_OK, "LLKA_crossResidueMetricAtoms() returned unexpected value") \
    EFF_expect(quad.a, std::string{_a}, std::string{"Atom a on distance "} + LLKA_crossResidueMetricName(xr, false) + " is wrong") \
    EFF_expect(quad.b, std::string{_b}, std::string{"Atom b on distance "} + LLKA_crossResidueMetricName(xr, false) + " is wrong")

#define CHECK_XR_TOR_ATOMS_STRU(stru, xr, quad, _a, _b, _c, _d) \
    tRet = LLKA_crossResidueMetricAtomsFromStructure(&stru, xr, &quad); \
    EFF_expect(tRet, LLKA_OK, "LLKA_crossResidueMetricAtomsFromStructure() returned unexpected value") \
    EFF_expect(quad.a, std::string{_a}, std::string{"Atom a on torsion "} + LLKA_crossResidueMetricName(xr, false) + " is wrong") \
    EFF_expect(quad.b, std::string{_b}, std::string{"Atom b on torsion "} + LLKA_crossResidueMetricName(xr, false) + " is wrong") \
    EFF_expect(quad.c, std::string{_c}, std::string{"Atom c on torsion "} + LLKA_crossResidueMetricName(xr, false) + " is wrong") \
    EFF_expect(quad.d, std::string{_d}, std::string{"Atom d on torsion "} + LLKA_crossResidueMetricName(xr, false) + " is wrong")

static
auto testTorsionNamesStructure()
{
    LLKA_AtomNameQuad quad;
    LLKA_RetCode tRet;

    LLKA_Structure stru = LLKA_makeStructure(REAL_1BNA_A_1_2_ATOMS, REAL_1BNA_A_1_2_ATOMS_LEN);

    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_DELTA_1,   quad, "C5'", "C4'", "C3'", "O3'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_EPSILON_1, quad, "C4'", "C3'", "O3'", "P");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_ZETA_1,    quad, "C3'", "O3'", "P", "O5'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_ALPHA_2,   quad, "O3'", "P", "O5'", "C5'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_BETA_2,    quad, "P", "O5'", "C5'", "C4'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_GAMMA_2,   quad, "O5'", "C5'", "C4'", "C3'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_DELTA_2,   quad, "C5'", "C4'", "C3'", "O3'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_CHI_1,     quad, "O4'", "C1'", "N1", "C2");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_CHI_2,     quad, "O4'", "C1'", "N9", "C4");
    CHECK_XR_DIST_ATOMS_STRU(stru,  LLKA_XR_DIST_CC,    quad, "C2", "C4");
    CHECK_XR_DIST_ATOMS_STRU(stru,  LLKA_XR_DIST_NN,    quad, "N1", "N9");
    CHECK_XR_TOR_ATOMS_STRU(stru,   LLKA_XR_TOR_MU,     quad, "C2", "N1", "N9", "C4");

    LLKA_destroyStructure(&stru);

    stru = LLKA_makeStructure(REAL_1BZT_A_12_13_ATOMS, REAL_1BZT_A_12_13_ATOMS_LEN);

    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_DELTA_1,   quad, "C5'", "C4'", "C3'", "O3'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_EPSILON_1, quad, "C4'", "C3'", "O3'", "P");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_ZETA_1,    quad, "C3'", "O3'", "P", "O5'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_ALPHA_2,   quad, "O3'", "P", "O5'", "C5'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_BETA_2,    quad, "P", "O5'", "C5'", "C4'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_GAMMA_2,   quad, "O5'", "C5'", "C4'", "C3'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_DELTA_2,   quad, "C5'", "C4'", "C3'", "O3'");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_CHI_1,     quad, "O4'", "C1'", "N9", "C4");
    CHECK_DINU_TOR_ATOMS_STRU(stru, LLKA_TOR_CHI_2,     quad, "O4'", "C1'", "C5", "C4");
    CHECK_XR_DIST_ATOMS_STRU(stru,  LLKA_XR_DIST_CC,    quad, "C4", "C4");
    CHECK_XR_DIST_ATOMS_STRU(stru,  LLKA_XR_DIST_NN,    quad, "N9", "C5");
    CHECK_XR_TOR_ATOMS_STRU(stru,   LLKA_XR_TOR_MU,     quad, "C4", "N9", "C5", "C4");

    LLKA_destroyStructure(&stru);
}

auto main(int, char **) -> int
{
    testCANA();
    testNtC();

    testDifferenceAgainstReference();

    testTorsionNamesStructure();
    testTorsionNamesBases();
}
