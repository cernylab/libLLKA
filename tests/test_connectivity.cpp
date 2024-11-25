/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_connectivity_similarity.h>
#include <llka_superposition.h>
#include <util/printers.hpp>

#include "effedup.hpp"
#include "testing_structures.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>

static
auto prnConnectivity(const std::string &desc, const LLKA_Connectivity &conn)
{
    std::cout
        << std::setprecision(15)
        << desc << "\n"
        << "C5' distance: " << conn.C5PrimeDistance
        << ", O3' distance: " << conn.O3PrimeDistance
        << "\n";
}

static
auto testConnectivity()
{
    // Test based on results for 1BNA molecule
    const LLKA_Structure realFirst = LLKA_makeStructure(REAL_1BNA_A_2_3_ATOMS, REAL_1BNA_A_2_3_ATOMS_LEN);
    const LLKA_Structure realSecond = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    const LLKA_Structure refFirst = LLKA_makeStructure(REF_BA01_ATOMS, REF_BA01_ATOMS_LEN);
    const LLKA_Structure refSecond = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);

    LLKA_RetCode tRet;
    LLKA_Connectivity conn;

    tRet = LLKA_measureStepConnectivityStructures(&realFirst, &refFirst, &realSecond, &refSecond, &conn);
    EFF_expect(tRet, LLKA_OK, "measureStepConnectivityStructures() returned unexpected value");
    prnConnectivity("1BNA 2_3 -> 3_4, BA01 (tested), AB01 (actual)" , conn);
    EFF_cmpFlt(conn.C5PrimeDistance, 0.0778497140143741, "wrong C5' distance")
    EFF_cmpFlt(conn.O3PrimeDistance, 0.110550489844076, "wrong O3' distance")

    LLKA_destroyStructure(&refSecond);
    LLKA_destroyStructure(&refFirst);
    LLKA_destroyStructure(&realSecond);
    LLKA_destroyStructure(&realFirst);
}

static
auto testConnectivityMultiple()
{
    // Test based on results for 1BNA molecule
    const LLKA_Structure realFirst = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    const LLKA_Structure realSecond = LLKA_makeStructure(REAL_1BNA_A_4_5_ATOMS, REAL_1BNA_A_4_5_ATOMS_LEN);

    const LLKA_Structure refFirst = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);

    LLKA_Structures refsSecond{
        new LLKA_Structure[3],
        3
    };
    LLKA_Connectivities results{
        new LLKA_Connectivity[3],
        3
    };

    refsSecond.strus[0] = LLKA_makeStructure(REF_BA01_ATOMS, REF_BA01_ATOMS_LEN);
    refsSecond.strus[1] = LLKA_makeStructure(REF_BB00_ATOMS, REF_BB00_ATOMS_LEN);
    refsSecond.strus[2] = LLKA_makeStructure(REF_BB04_ATOMS, REF_BB04_ATOMS_LEN);

    auto tRet = LLKA_measureStepConnectivityStructuresMultiple(&realFirst, &refFirst, &realSecond, &refsSecond, &results);
    EFF_expect(tRet, LLKA_OK, "measureStepConnectivityMultiple() returned unexpected value")

    prnConnectivity("1BNA 3_4 -> 4_5, AB01 (actual) -> BA01 (tested)", results.conns[0]);
    EFF_cmpFlt(results.conns[0].C5PrimeDistance, 0.222853768522876, "wrong C5' distance");
    EFF_cmpFlt(results.conns[0].O3PrimeDistance, 0.332800460520272, "wrong O3' distance")

    prnConnectivity("1BNA 3_4 -> 4_5, AB01 (actual) -> BB00 (tested)", results.conns[1]);
    EFF_cmpFlt(results.conns[1].C5PrimeDistance, 0.398108806100796, "wrong C5' distance");
    EFF_cmpFlt(results.conns[1].O3PrimeDistance, 0.29940718941947, "wrong O3' distance")

    prnConnectivity("1BNA 3_4 -> 4_5, BB04 (actual) -> BB04 (tested)", results.conns[2]);
    EFF_cmpFlt(results.conns[2].C5PrimeDistance, 0.379975441318778, "wrong C5' distance");
    EFF_cmpFlt(results.conns[2].O3PrimeDistance, 0.372926953223463, "wrong O3' distance")

    delete [] results.conns;

    // We have allocated LLKA_Structures manually so we cannot call LLKA_destroyStructures() safely
    // Let's free the resources ourselves
    LLKA_destroyStructure(&refsSecond.strus[2]);
    LLKA_destroyStructure(&refsSecond.strus[1]);
    LLKA_destroyStructure(&refsSecond.strus[0]);
    delete [] refsSecond.strus;

    LLKA_destroyStructure(&refFirst);
    LLKA_destroyStructure(&realSecond);
    LLKA_destroyStructure(&realFirst);
}

static
auto testConnectivityNtCs()
{
    // Test based on results for 1BNA molecule
    const LLKA_Structure realFirst = LLKA_makeStructure(REAL_1BNA_A_2_3_ATOMS, REAL_1BNA_A_2_3_ATOMS_LEN);
    const LLKA_Structure realSecond = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);

    LLKA_RetCode tRet;
    LLKA_Connectivity conn;

    tRet = LLKA_measureStepConnectivityNtCs(&realFirst, LLKA_BA01, &realSecond, LLKA_AB01, &conn);
    EFF_expect(tRet, LLKA_OK, "measureStepConnectivityNtCs() returned unexpected value");
    prnConnectivity("1BNA 2_3 -> 3_4, BA01 (tested), AB01 (actual)" , conn);
    EFF_cmpFlt(conn.C5PrimeDistance, 0.0778497140143741, "wrong C5' distance")
    EFF_cmpFlt(conn.O3PrimeDistance, 0.110550489844076, "wrong O3' distance")

    LLKA_destroyStructure(&realSecond);
    LLKA_destroyStructure(&realFirst);
}

static
auto testConnectivityMultipleNtCsFirst()
{
    const LLKA_Structure realFirst = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    const LLKA_Structure realSecond = LLKA_makeStructure(REAL_1BNA_A_4_5_ATOMS, REAL_1BNA_A_4_5_ATOMS_LEN);
    const LLKA_NtC ntcsFirst[] = { LLKA_AB02, LLKA_BBS1, LLKA_BB07, LLKA_INVALID_NTC };

    LLKA_Connectivities results{
        new LLKA_Connectivity[3],
        3
    };

    auto tRet = LLKA_measureStepConnectivityNtCsMultipleFirst(&realFirst, ntcsFirst, &realSecond, LLKA_BB04, &results);
    EFF_expect(tRet, LLKA_OK, "measureStepConnectivityNtCsMultipleFirst() returned unexpected value");

    prnConnectivity("1BNA 3_4 -> 4_5, AB02 (tested) -> BB04 (actual)", results.conns[0]);
    EFF_cmpFlt(results.conns[0].C5PrimeDistance, 0.172379441824832, "wrong C5' distance");
    EFF_cmpFlt(results.conns[0].O3PrimeDistance, 0.218599268902842, "wrong O3' distance")

    prnConnectivity("1BNA 3_4 -> 4_5, BBS1 (tested) -> BB04 (actual)", results.conns[1]);
    EFF_cmpFlt(results.conns[1].C5PrimeDistance, 0.22807852040331, "wrong C5' distance");
    EFF_cmpFlt(results.conns[1].O3PrimeDistance, 0.118886962245301, "wrong O3' distance")

    prnConnectivity("1BNA 3_4 -> 4_5, BB07 (tested) -> BB04 (actual)", results.conns[2]);
    EFF_cmpFlt(results.conns[2].C5PrimeDistance, 0.289601415171216, "wrong C5' distance");
    EFF_cmpFlt(results.conns[2].O3PrimeDistance, 0.106908891049413, "wrong O3' distance")

    delete[] results.conns;

    LLKA_destroyStructure(&realSecond);
    LLKA_destroyStructure(&realFirst);
}

static
auto testConnectivityMultipleNtCsSecond()
{
    // Test based on results for 1BNA molecule
    const LLKA_Structure realFirst = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    const LLKA_Structure realSecond = LLKA_makeStructure(REAL_1BNA_A_4_5_ATOMS, REAL_1BNA_A_4_5_ATOMS_LEN);
    const LLKA_NtC ntcsSecond[] = { LLKA_BA01, LLKA_BB00, LLKA_BB04, LLKA_INVALID_NTC };

    LLKA_Connectivities results{
        new LLKA_Connectivity[3],
        3
    };

    auto tRet = LLKA_measureStepConnectivityNtCsMultipleSecond(&realFirst, LLKA_AB01, &realSecond, ntcsSecond, &results);
    EFF_expect(tRet, LLKA_OK, "measureStepConnectivityMultipleSecond() returned unexpected value")

    prnConnectivity("1BNA 3_4 -> 4_5, AB01 (actual) -> BA01 (tested)", results.conns[0]);
    EFF_cmpFlt(results.conns[0].C5PrimeDistance, 0.222853768522876, "wrong C5' distance");
    EFF_cmpFlt(results.conns[0].O3PrimeDistance, 0.332800460520272, "wrong O3' distance")

    prnConnectivity("1BNA 3_4 -> 4_5, AB01 (actual) -> BB00 (tested)", results.conns[1]);
    EFF_cmpFlt(results.conns[1].C5PrimeDistance, 0.398108806100796, "wrong C5' distance");
    EFF_cmpFlt(results.conns[1].O3PrimeDistance, 0.29940718941947, "wrong O3' distance")

    EFF_cmpFlt(results.conns[2].C5PrimeDistance, 0.379975441318778, "wrong C5' distance");
    EFF_cmpFlt(results.conns[2].O3PrimeDistance, 0.372926953223463, "wrong O3' distance")

    delete[] results.conns;

    LLKA_destroyStructure(&realSecond);
    LLKA_destroyStructure(&realFirst);
}

auto main(int, char**) -> int
{
    testConnectivity();
    testConnectivityMultiple();
    testConnectivityNtCs();
    testConnectivityMultipleNtCsFirst();
    testConnectivityMultipleNtCsSecond();

    return EXIT_SUCCESS;
}
