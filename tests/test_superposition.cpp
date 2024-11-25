/* vim: set sw=4 ts=4 sts=4 expandtab : */


#include <llka_ntc.h>
#include <llka_superposition.h>

#include "effedup.hpp"
#include "testing_structures.h"

#include <sstream>

#define CHECK_POINT(pt, ref) \
    do { \
        EFF_cmpFlt(pt.x, ref.x, "Point X coordinate is wrong"); \
        EFF_cmpFlt(pt.y, ref.y, "Point Y coordinate is wrong"); \
        EFF_cmpFlt(pt.z, ref.z, "Point Z coordinate is wrong"); \
    } while (false)

static
auto strMatrix(const LLKA_Matrix &matrix)
{
    std::ostringstream oss{};

    for (size_t r = 0; r < matrix.nRows; r++) {
        for (size_t c = 0; c < matrix.nCols; c++) {
            oss << LLKA_matrixGet(r, c, &matrix) << "; ";
        }
        oss << "\n";
    }

    return oss.str();
}

static
auto prnPt(const LLKA_Point &pt)
{
    std::cout << std::setprecision(15) << "{" << pt.x << ", " << pt.y << ", " << pt.z << "}\n";
}

static
auto prnPoints(const LLKA_Points &pts)
{
    for (size_t idx = 0; idx < pts.nPoints; idx++)
        prnPt(pts.points[idx]);
}

static
auto prnStru(const LLKA_Structure &stru)
{
    for (size_t idx = 0; idx < stru.nAtoms; idx++)
        prnPt(stru.atoms[idx].coords);
}

static
auto prnStru(const LLKA_StructureView &stru)
{
    for (size_t idx = 0; idx < stru.nAtoms; idx++)
        prnPt(stru.atoms[idx]->coords);
}

template <typename T>
static
auto prnSuperposition(const std::string &desc, const LLKA_Structure &superposed, const T &reference, double rmsd)
{
    std::cout << "Superposition " << desc << "\n";
    std::cout << "Superposed:\n";
    prnStru(superposed);
    std::cout << "Reference:\n";
    prnStru(reference);
    std::cout << "RMSD = " << rmsd << "\n\n";
}

static
auto prnSuperposition(const std::string &desc, const LLKA_Points &superposed, const LLKA_Points &reference, double rmsd)
{
    std::cout << "Superposition " << desc << "\n";
    std::cout << "Superposed:\n";
    prnPoints(superposed);
    std::cout << "Reference:\n";
    prnPoints(reference);
    std::cout << "RMSD = " << rmsd << "\n\n";
}

static
auto testShifted()
{
    LLKA_Points a{
        {new LLKA_Point[4]},
        4
    };
    LLKA_Points b{
        {new LLKA_Point[4]},
        4
    };

    a.points[0] = { 1, 1, 1 };
    a.points[1] = { 2, 2, 2 };
    a.points[2] = { 3, 3, 3 };
    a.points[3] = { 4, 4, 4 };

    b.points[0] = { 10 + 1, 10 + 1, 10 + 1 };
    b.points[1] = { 10 + 2, 10 + 2, 10 + 2 };
    b.points[2] = { 10 + 3, 10 + 3, 10 + 3 };
    b.points[3] = { 10 + 4, 10 + 4, 10 + 4 };

    double rmsd;
    auto tRet = LLKA_superposePoints(&a, &b, &rmsd);
    EFF_expect(tRet, LLKA_OK, "LLKA_superposePoints() returned unexpected value")

    prnSuperposition("Shifted square", a, b, rmsd);

    CHECK_POINT(a.points[0], b.points[0]);
    CHECK_POINT(a.points[1], b.points[1]);
    CHECK_POINT(a.points[2], b.points[2]);
    CHECK_POINT(a.points[3], b.points[3]);
    EFF_cmpFlt(rmsd, 0.0, "RMSD is wrong")

    delete[] a.points;
    delete[] b.points;
}

auto testShifted2()
{
    LLKA_Points a{
        {new LLKA_Point[4]},
        4
    };
    LLKA_Points b{
        {new LLKA_Point[4]},
        4
    };

    a.points[0] = { 0, 0, 0 };
    a.points[1] = { 1, 0, 0 };
    a.points[2] = { 1, 1, 0 };
    a.points[3] = { 0, 1, 0 };

    b.points[0] = { 10 + 0, 10 + 0, 10 + 0 };
    b.points[1] = { 10 + 0, 10 + 1, 10 + 0 };
    b.points[2] = { 10 + 0, 10 + 1, 10 + 1 };
    b.points[3] = { 10 + 0, 10 + 0, 10 + 1 };

    double rmsd;
    auto tRet = LLKA_superposePoints(&a, &b, &rmsd);
    EFF_expect(tRet, LLKA_OK, "LLKA_superposePoints() returned unexpected value")

    prnSuperposition("Shifted square - Z", a, b, rmsd);

    CHECK_POINT(a.points[0], b.points[0]);
    CHECK_POINT(a.points[1], b.points[1]);
    CHECK_POINT(a.points[2], b.points[2]);
    CHECK_POINT(a.points[3], b.points[3]);
    EFF_cmpFlt(rmsd, 0.0, "RMSD is wrong");

    delete[] a.points;
    delete[] b.points;
}

auto testOneMismatching()
{
    LLKA_Points a{
        {new LLKA_Point[4]},
        4
    };
    LLKA_Points b{
        {new LLKA_Point[4]},
        4
    };

    a.points[0] = { 1, 1, 1 };
    a.points[1] = { 2, 2, 2 };
    a.points[2] = { 3, 3, 3 };
    a.points[3] = { 4, 4, 4 };

    b.points[0] = { 10 + 1, 10 + 1, 10 + 1 };
    b.points[1] = { 10 + 2, 10 + 2, 10 + 2 };
    b.points[2] = { 10 + 3, 10 + 3, 10 + 3 };
    b.points[3] = { 10 + 4, 10 + 4, 10 + 7 };

    double rmsd;
    auto tRet = LLKA_superposePoints(&a, &b, &rmsd);
    EFF_expect(tRet, LLKA_OK, "LLKA_superposePoints() returned unexpected value")

    prnSuperposition("One mismatching", a, b, rmsd);

    LLKA_Point expected0 = { 11.4030913638093, 11.4030913638093, 11.1658735912377 };
    LLKA_Point expected1 = { 12.1343637879364, 12.1343637879364, 12.5552911970792 };
    LLKA_Point expected2 = { 12.8656362120636, 12.8656362120636, 13.9447088029208 };
    LLKA_Point expected3 = { 13.5969086361907, 13.5969086361907, 15.3341264087623 };

    CHECK_POINT(a.points[0], expected0);
    CHECK_POINT(a.points[1], expected1);
    CHECK_POINT(a.points[2], expected2);
    CHECK_POINT(a.points[3], expected3);
    EFF_cmpFlt(rmsd, 1.0869242161333, "RMSD is wrong")

    delete[] a.points;
    delete[] b.points;
}

static
auto testRealPlainBackbone()
{
    LLKA_Structure ref_AB01 = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA_Structure real_AB01 = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    LLKA_Structure ref_AB01_bkbn;
    LLKA_Structure real_AB01_bkbn;

    auto tRet = LLKA_extractBackbone(&ref_AB01, &ref_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractBackbone() returned unexpected value")

    tRet = LLKA_extractBackbone(&real_AB01, &real_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractBackbone() returned unexpected value")

    double rmsd;
    tRet = LLKA_superposeStructures(&real_AB01_bkbn, &ref_AB01_bkbn, &rmsd);
    EFF_expect(tRet, LLKA_OK, "LLKA_superposePoints() returned unexpected value")

    prnSuperposition("Superposition of plain backbone on AB01 on 1BNA_A_3_4", real_AB01_bkbn, ref_AB01_bkbn, rmsd);

    LLKA_Points expected{
        {new LLKA_Point[12]},
        12
    };
    expected.points[0] = { 0.950860168366356, -4.10490242639581, 1.88352814333396 };
    expected.points[1] = { 0.50054990344796, -3.69993907500445, 0.514678133530478 };
    expected.points[2] = { -0.347244679657512, -4.71592059085125, -0.00580232173625572 };
    expected.points[3] = { -0.318861255763144, -2.43716579580966, 0.473326017357792 };
    expected.points[4] = { 0.506320258152502, -1.32662194302102, 0.113179074950623 };
    expected.points[5] = { -0.0811260469759572, 0.154118620506023, 0.0681432878857751};
    expected.points[6] = { -0.647318283320744, 0.253840078930681, -1.40665838767353 };
    expected.points[7] = { 0.237737393186242, 0.145726060576478, -2.52768907793944 };
    expected.points[8] = { -0.603024580420967, 0.282861734520597, -3.76182161161158 };
    expected.points[9] = { -1.56591988716846, -0.783314033669751, -3.80479671565305 };
    expected.points[10] = { -1.41697603269085, 1.56517005180424, -3.86385794226367 };
    expected.points[11] = { -1.71799695715546, 1.88014731841389, -5.20722860018113 };

    for (size_t idx = 0; idx < 12; idx++)
        CHECK_POINT(real_AB01_bkbn.atoms[idx].coords, expected.points[idx]);

    EFF_cmpFlt(rmsd, 0.147657613575976, "RMSD is wrong");

    delete[] expected.points;

    LLKA_destroyStructure(&ref_AB01_bkbn);
    LLKA_destroyStructure(&real_AB01_bkbn);
    LLKA_destroyStructure(&ref_AB01);
    LLKA_destroyStructure(&real_AB01);
}

static
auto testRealExtendedBackbone()
{
    LLKA_Structure ref_AB01 = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA_Structure real_AB01 = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    LLKA_Structure ref_AB01_bkbn;
    LLKA_Structure real_AB01_bkbn;

    auto tRet = LLKA_extractExtendedBackbone(&ref_AB01, &ref_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value")

    tRet = LLKA_extractExtendedBackbone(&real_AB01, &real_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value")

    double rmsd;
    tRet = LLKA_superposeStructures(&real_AB01_bkbn, &ref_AB01_bkbn, &rmsd);
    EFF_expect(tRet, LLKA_OK, "LLKA_superposePoints() returned unexpected value")

    prnSuperposition("Superposition of extended backbone of AB01 on 1BNA_A_3_4", real_AB01_bkbn, ref_AB01_bkbn, rmsd);

    LLKA_Points expected{
        {new LLKA_Point[18]},
        18
    };
    expected.points[0] = { 0.987110766185577, -4.09071771489897, 1.9236822102414 };
    expected.points[1] = { 0.498993067687067, -3.70426391350659, 0.562462326533447 };
    expected.points[2] = { -0.388755973663699, -4.71109589650857, 0.0933514680094707 };
    expected.points[3] = { -0.292926693772989, -2.42405654051827, 0.52155788233221 };
    expected.points[4] = { 0.54395097846502, -1.34065484814116, 0.109472290902603 };
    expected.points[5] = { -1.17780200834469, -4.12249011848516, -0.933417051780487 };
    expected.points[6] = { -2.43228362842333, -4.89306776930891, -0.85803171543564 };
    expected.points[7] = { -3.01959843323915, -5.3413321995708, -1.99405078244382 };
    expected.points[8] = { -0.011455206873104, 0.151981739114691, 0.052253520613716 };
    expected.points[9] = { -0.625836124161192, 0.233795179516065, -1.40424936281698 };
    expected.points[10] = { 0.217366956113348, 0.0817923356419141, -2.55212814599476 };
    expected.points[11] = { -0.662209557696309, 0.212387424607011, -3.75962973153933 };
    expected.points[12] = { -1.64960487803313, -0.831943050574988, -3.7466166618369 };
    expected.points[13] = { -1.45033183726805, 1.51072162058409, -3.86177956266596 };
    expected.points[14] = { -1.79041850416729, 1.80432120889078, -5.2006631799918 };
    expected.points[15] = { -2.93672182938542, -0.21171102176145, -3.72940666545558 };
    expected.points[16] = { -3.83674604788897, -1.07392594379514, -2.92836744398792 };
    expected.points[17] = { -4.64473104553373, -2.04674049128451, -3.42943939468367 };

    for (size_t idx = 0; idx < 14; idx++)
        CHECK_POINT(real_AB01_bkbn.atoms[idx].coords, expected.points[idx]);

    EFF_cmpFlt(rmsd, 0.176588293714412, "RMSD is wrong");

    delete[] expected.points;

    LLKA_destroyStructure(&ref_AB01_bkbn);
    LLKA_destroyStructure(&real_AB01_bkbn);
    LLKA_destroyStructure(&ref_AB01);
    LLKA_destroyStructure(&real_AB01);
}

static
auto testRealExtendedBackboneView()
{
    LLKA_Structure ref_AB01 = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA_Structure real_AB01 = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    LLKA_StructureView ref_AB01_bkbn;
    LLKA_Structure real_AB01_bkbn;

    auto tRet = LLKA_extractExtendedBackboneView(&ref_AB01, &ref_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value")

    tRet = LLKA_extractExtendedBackbone(&real_AB01, &real_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value")

    double rmsd;
    tRet = LLKA_superposeStructuresView(&real_AB01_bkbn, &ref_AB01_bkbn, &rmsd);
    EFF_expect(tRet, LLKA_OK, "LLKA_superposePoints() returned unexpected value")

    prnSuperposition("Superposition of extended backbone of AB01 on 1BNA_A_3_4", real_AB01_bkbn, ref_AB01_bkbn, rmsd);

    LLKA_Points expected{
        {new LLKA_Point[18]},
        18
    };
    expected.points[0] = { 0.987110766185577, -4.09071771489897, 1.9236822102414 };
    expected.points[1] = { 0.498993067687067, -3.70426391350659, 0.562462326533447 };
    expected.points[2] = { -0.388755973663699, -4.71109589650857, 0.0933514680094707 };
    expected.points[3] = { -0.292926693772989, -2.42405654051827, 0.52155788233221 };
    expected.points[4] = { 0.54395097846502, -1.34065484814116, 0.109472290902603 };
    expected.points[5] = { -1.17780200834469, -4.12249011848516, -0.933417051780487 };
    expected.points[6] = { -2.43228362842333, -4.89306776930891, -0.85803171543564 };
    expected.points[7] = { -3.01959843323915, -5.3413321995708, -1.99405078244382 };
    expected.points[8] = { -0.011455206873104, 0.151981739114691, 0.052253520613716 };
    expected.points[9] = { -0.625836124161192, 0.233795179516065, -1.40424936281698 };
    expected.points[10] = { 0.217366956113348, 0.0817923356419141, -2.55212814599476 };
    expected.points[11] = { -0.662209557696309, 0.212387424607011, -3.75962973153933 };
    expected.points[12] = { -1.64960487803313, -0.831943050574988, -3.7466166618369 };
    expected.points[13] = { -1.45033183726805, 1.51072162058409, -3.86177956266596 };
    expected.points[14] = { -1.79041850416729, 1.80432120889078, -5.2006631799918 };
    expected.points[15] = { -2.93672182938542, -0.21171102176145, -3.72940666545558 };
    expected.points[16] = { -3.83674604788897, -1.07392594379514, -2.92836744398792 };
    expected.points[17] = { -4.64473104553373, -2.04674049128451, -3.42943939468367 };

    for (size_t idx = 0; idx < 14; idx++)
        CHECK_POINT(real_AB01_bkbn.atoms[idx].coords, expected.points[idx]);

    EFF_cmpFlt(rmsd, 0.176588293714412, "RMSD is wrong");

    delete[] expected.points;

    LLKA_destroyStructureView(&ref_AB01_bkbn);
    LLKA_destroyStructure(&real_AB01_bkbn);
    LLKA_destroyStructure(&ref_AB01);
    LLKA_destroyStructure(&real_AB01);
}

static
auto testTransformationMatrixPoints()
{
    LLKA_Points a{
        {new LLKA_Point[4]},
        4
    };
    LLKA_Points b{
        {new LLKA_Point[4]},
        4
    };

    a.points[0] = { 1, 1, 1 };
    a.points[1] = { 2, 2, 2 };
    a.points[2] = { 3, 3, 3 };
    a.points[3] = { 4, 4, 4 };

    b.points[0] = { 10 + 1, 10 + 1, 10 + 1 };
    b.points[1] = { 10 + 2, 10 + 2, 10 + 2 };
    b.points[2] = { 10 + 3, 10 + 3, 10 + 3 };
    b.points[3] = { 10 + 4, 10 + 4, 10 + 7 };

    LLKA_Matrix transformation{};
    auto tRet = LLKA_superpositionMatrixPoints(&a, &b, &transformation);
    EFF_expect(tRet, LLKA_OK, "LLKA_superpositionMatrixPoints() returned unexpected value")

    tRet = LLKA_applyTransformationPoints(&a, &transformation);
    EFF_expect(tRet, LLKA_OK, "LLKA_applyTransformationPoints() returned unexpected value")

    prnSuperposition("Transformation matrix - one mismatching", a, b, 0.0);

    LLKA_Point expected0 = { 11.4030913638093, 11.4030913638093, 11.1658735912377 };
    LLKA_Point expected1 = { 12.1343637879364, 12.1343637879364, 12.5552911970792 };
    LLKA_Point expected2 = { 12.8656362120636, 12.8656362120636, 13.9447088029208 };
    LLKA_Point expected3 = { 13.5969086361907, 13.5969086361907, 15.3341264087623 };

    CHECK_POINT(a.points[0], expected0);
    CHECK_POINT(a.points[1], expected1);
    CHECK_POINT(a.points[2], expected2);
    CHECK_POINT(a.points[3], expected3);

    LLKA_destroyMatrix(&transformation);
    delete[] a.points;
    delete[] b.points;
}

static
auto testTransformationMatrixStructure()
{
    LLKA_Structure ref_AB01 = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA_Structure real_AB01 = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    LLKA_Structure ref_AB01_bkbn;
    LLKA_Structure real_AB01_bkbn;

    auto tRet = LLKA_extractExtendedBackbone(&ref_AB01, &ref_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value")

    tRet = LLKA_extractExtendedBackbone(&real_AB01, &real_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value")

    LLKA_Matrix transformation{};
    tRet = LLKA_superpositionMatrixStructures(&real_AB01_bkbn, &ref_AB01_bkbn, &transformation);
    EFF_expect(tRet, LLKA_OK, "LLKA_superpositionMatrixStructures() returned unexpected value")

    tRet = LLKA_applyTransformationStructure(&real_AB01_bkbn, &transformation);
    EFF_expect(tRet, LLKA_OK, "LLKA_applyTransformationStructure() returned unexpected value")

    prnSuperposition("Transformation matrix - superposition of extended backbone of AB01 on 1BNA_A_3_4", real_AB01_bkbn, ref_AB01_bkbn, 0.0);

    LLKA_Points expected{
        {new LLKA_Point[18]},
        18
    };
    expected.points[0] = { 0.987110766185577, -4.09071771489897, 1.9236822102414 };
    expected.points[1] = { 0.498993067687067, -3.70426391350659, 0.562462326533447 };
    expected.points[2] = { -0.388755973663699, -4.71109589650857, 0.0933514680094707 };
    expected.points[3] = { -0.292926693772989, -2.42405654051827, 0.52155788233221 };
    expected.points[4] = { 0.54395097846502, -1.34065484814116, 0.109472290902603 };
    expected.points[5] = { -1.17780200834469, -4.12249011848516, -0.933417051780487 };
    expected.points[6] = { -2.43228362842333, -4.89306776930891, -0.85803171543564 };
    expected.points[7] = { -3.01959843323915, -5.3413321995708, -1.99405078244382 };
    expected.points[8] = { -0.011455206873104, 0.151981739114691, 0.052253520613716 };
    expected.points[9] = { -0.625836124161192, 0.233795179516065, -1.40424936281698 };
    expected.points[10] = { 0.217366956113348, 0.0817923356419141, -2.55212814599476 };
    expected.points[11] = { -0.662209557696309, 0.212387424607011, -3.75962973153933 };
    expected.points[12] = { -1.64960487803313, -0.831943050574988, -3.7466166618369 };
    expected.points[13] = { -1.45033183726805, 1.51072162058409, -3.86177956266596 };
    expected.points[14] = { -1.79041850416729, 1.80432120889078, -5.2006631799918 };
    expected.points[15] = { -2.93672182938542, -0.21171102176145, -3.72940666545558 };
    expected.points[16] = { -3.83674604788897, -1.07392594379514, -2.92836744398792 };
    expected.points[17] = { -4.64473104553373, -2.04674049128451, -3.42943939468367 };

    for (size_t idx = 0; idx < 14; idx++)
        CHECK_POINT(real_AB01_bkbn.atoms[idx].coords, expected.points[idx]);

    delete[] expected.points;

    LLKA_destroyMatrix(&transformation);
    LLKA_destroyStructure(&ref_AB01_bkbn);
    LLKA_destroyStructure(&real_AB01_bkbn);
    LLKA_destroyStructure(&ref_AB01);
    LLKA_destroyStructure(&real_AB01);

}

static
auto testTransformationMatrixStructureView()
{
    LLKA_Structure ref_AB01 = LLKA_makeStructure(REF_AB01_ATOMS, REF_AB01_ATOMS_LEN);
    LLKA_Structure real_AB01 = LLKA_makeStructure(REAL_1BNA_A_3_4_ATOMS, REAL_1BNA_A_3_4_ATOMS_LEN);
    LLKA_StructureView ref_AB01_bkbn;
    LLKA_StructureView real_AB01_bkbn;

    auto tRet = LLKA_extractExtendedBackboneView(&ref_AB01, &ref_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value")

    tRet = LLKA_extractExtendedBackboneView(&real_AB01, &real_AB01_bkbn);
    EFF_expect(tRet, LLKA_OK, "LLKA_extractExtendedBackbone() returned unexpected value")

    LLKA_Matrix transformation{};
    tRet = LLKA_superpositionMatrixStructureViews(&real_AB01_bkbn, &ref_AB01_bkbn, &transformation);
    EFF_expect(tRet, LLKA_OK, "LLKA_superpositionMatrixStructures() returned unexpected value")

    static const double expected[4][4] = {
        { 0.6781837391735, -0.278187893249037, 0.680204610371256, -23.0397139352514 },
        { 0.362607977150515, -0.678396158822793, -0.638978956305519, 15.7850194094942 },
        { 0.639204404574515, 0.679992755655633, -0.359204094392586, -22.3197349277867 },
        { 0, 0, 0, 1 }
    };

    EFF_expect(transformation.nRows, size_t(4), "Transormation matrix has invalid number of rows");
    EFF_expect(transformation.nCols, size_t(4), "Transormation matrix has invalid number of columns");

    for (size_t r = 0; r < transformation.nRows; r++) {
        for (size_t c = 0; c < transformation.nCols; c++) {
            EFF_cmpFlt(
                LLKA_matrixGet(r, c, &transformation),
                expected[r][c],
                "Wrong value in transformation matrix at [" + std::to_string(r) + "; " + std::to_string(c) + "]\nFull matrix:\n" + strMatrix(transformation)
            )
        }
    }

    LLKA_destroyMatrix(&transformation);
    LLKA_destroyStructureView(&ref_AB01_bkbn);
    LLKA_destroyStructureView(&real_AB01_bkbn);
    LLKA_destroyStructure(&ref_AB01);
    LLKA_destroyStructure(&real_AB01);

}

auto main() -> int
{
    testShifted();
    testShifted2();
    testOneMismatching();
    testRealPlainBackbone();
    testRealExtendedBackbone();
    testRealExtendedBackboneView();
    testTransformationMatrixPoints();
    testTransformationMatrixStructure();
    testTransformationMatrixStructureView();

    return EXIT_SUCCESS;
}
