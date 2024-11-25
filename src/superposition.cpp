// vim: set sw=4 ts=4 sts=4 expandtab :

#include <llka_superposition.h>

#include "superposition.hpp"
#include "util/templates.hpp"

#include <memory>

namespace LLKAInternal {

inline
auto makePts(const LLKA_Structure *what, const LLKA_Structure *onto)
{
    // TODO: Use aligned allocations here to speed things up
    std::unique_ptr<double[]> whatPts{new double[what->nAtoms * 3]};
    std::unique_ptr<double[]> ontoPts{new double[onto->nAtoms * 3]};

    for (size_t idx = 0; idx < what->nAtoms; idx++) {
        size_t off = 3 * idx;
        whatPts[off] = what->atoms[idx].coords.x;
        whatPts[off + 1] = what->atoms[idx].coords.y;
        whatPts[off + 2] = what->atoms[idx].coords.z;

        ontoPts[off] = onto->atoms[idx].coords.x;
        ontoPts[off + 1] = onto->atoms[idx].coords.y;
        ontoPts[off + 2] = onto->atoms[idx].coords.z;
    }

    return std::make_tuple(std::move(whatPts), std::move(ontoPts));
}

inline
auto makePts(const LLKA_Structure *what, const LLKA_StructureView *onto)
{
    // TODO: Use aligned allocations here to speed things up
    std::unique_ptr<double[]> whatPts{new double[what->nAtoms * 3]};
    std::unique_ptr<double[]> ontoPts{new double[onto->nAtoms * 3]};

    for (size_t idx = 0; idx < what->nAtoms; idx++) {
        size_t off = 3 * idx;
        whatPts[off] = what->atoms[idx].coords.x;
        whatPts[off + 1] = what->atoms[idx].coords.y;
        whatPts[off + 2] = what->atoms[idx].coords.z;

        ontoPts[off] = onto->atoms[idx]->coords.x;
        ontoPts[off + 1] = onto->atoms[idx]->coords.y;
        ontoPts[off + 2] = onto->atoms[idx]->coords.z;
    }

    return std::make_tuple(std::move(whatPts), std::move(ontoPts));
}

inline
auto makePts(const LLKA_StructureView *what, const LLKA_StructureView *onto)
{
    // TODO: Use aligned allocations here to speed things up
    std::unique_ptr<double[]> whatPts{new double[what->nAtoms * 3]};
    std::unique_ptr<double[]> ontoPts{new double[onto->nAtoms * 3]};

    for (size_t idx = 0; idx < what->nAtoms; idx++) {
        size_t off = 3 * idx;
        whatPts[off] = what->atoms[idx]->coords.x;
        whatPts[off + 1] = what->atoms[idx]->coords.y;
        whatPts[off + 2] = what->atoms[idx]->coords.z;

        ontoPts[off] = onto->atoms[idx]->coords.x;
        ontoPts[off + 1] = onto->atoms[idx]->coords.y;
        ontoPts[off + 2] = onto->atoms[idx]->coords.z;
    }

    return std::make_tuple(std::move(whatPts), std::move(ontoPts));
}

template <LLKAStructureType T>
inline
auto superposeStructures(LLKA_Structure *what, const T *onto, double *rmsd)
{
    if (what->nAtoms != onto->nAtoms)
        return LLKA_E_MISMATCHING_SIZES;

    auto [whatPts, ontoPts] = makePts(what, onto);

    LLKAInternal::MappedPointsUnaligned whatMapped(whatPts.get(), 3, what->nAtoms);
    LLKAInternal::MappedPointsUnaligned ontoMapped(ontoPts.get(), 3, what->nAtoms);

    auto tRet = superpose(whatMapped, ontoMapped, rmsd);
    if (tRet != LLKA_OK)
        return tRet;

    for (size_t idx = 0; idx < what->nAtoms; idx++) {
        size_t off = 3 * idx;
        auto &atom = what->atoms[idx];
        atom.coords.x = whatPts[off];
        atom.coords.y = whatPts[off + 1];
        atom.coords.z = whatPts[off + 2];
    }

    return LLKA_OK;
}

template <LLKAStructureType T>
inline
auto superpositionMatrixStructures(const T *what, const T *onto, LLKA_Matrix *matrix)
{
    if (what->nAtoms != onto->nAtoms)
        return LLKA_E_MISMATCHING_SIZES;
    if (what->nAtoms == 0)
        return LLKA_E_INVALID_ARGUMENT;

    auto [whatPts, ontoPts] = LLKAInternal::makePts(what, onto);

    LLKAInternal::MappedPointsUnaligned whatMapped(whatPts.get(), 3, what->nAtoms);
    LLKAInternal::MappedPointsUnaligned ontoMapped(ontoPts.get(), 3, what->nAtoms);

    return LLKAInternal::superpositionMatrix(whatMapped, ontoMapped, matrix);
}

} // LLKAInternal

LLKA_RetCode LLKA_applyTransformationPoints(LLKA_Points *what, const LLKA_Matrix *matrix)
{
    if (matrix->nCols != 4 || matrix->nRows != 4)
        return LLKA_E_INVALID_ARGUMENT;

    // TODO: Use aligned allocations here to speed things up
    std::unique_ptr<double[]> whatPts{new double[what->nPoints * 4]};

    for (size_t idx = 0; idx < what->nPoints; idx++) {
        size_t off = 4 * idx;
        whatPts[off] = what->points[idx].x;
        whatPts[off + 1] = what->points[idx].y;
        whatPts[off + 2] = what->points[idx].z;
        whatPts[off + 3] = 1;
    }

    LLKAInternal::MappedQuadsUnaligned whatMapped(whatPts.get(), 4, what->nPoints);
    LLKAInternal::MappedMat4x4Unaligned matrixMapped(matrix->data);

    LLKAInternal::applySuperposition(whatMapped, matrixMapped);

    for (size_t idx = 0; idx < what->nPoints; idx++) {
        size_t off = 4 * idx;
        auto &point = what->points[idx];
        point.x = whatPts[off];
        point.y = whatPts[off + 1];
        point.z = whatPts[off + 2];
    }

    return LLKA_OK;
}

LLKA_RetCode LLKA_applyTransformationStructure(LLKA_Structure *what, const LLKA_Matrix *matrix)
{
    if (matrix->nCols != 4 || matrix->nRows != 4)
        return LLKA_E_INVALID_ARGUMENT;

    // TODO: Use aligned allocations here to speed things up
    std::unique_ptr<double[]> whatPts{new double[what->nAtoms * 4]};

    for (size_t idx = 0; idx < what->nAtoms; idx++) {
        size_t off = 4 * idx;
        whatPts[off] = what->atoms[idx].coords.x;
        whatPts[off + 1] = what->atoms[idx].coords.y;
        whatPts[off + 2] = what->atoms[idx].coords.z;
        whatPts[off + 3] = 1;
    }

    LLKAInternal::MappedQuadsUnaligned whatMapped(whatPts.get(), 4, what->nAtoms);
    LLKAInternal::MappedMat4x4Unaligned matrixMapped(matrix->data);

    LLKAInternal::applySuperposition(whatMapped, matrixMapped);

    for (size_t idx = 0; idx < what->nAtoms; idx++) {
        size_t off = 4 * idx;
        auto &atom = what->atoms[idx];
        atom.coords.x = whatPts[off];
        atom.coords.y = whatPts[off + 1];
        atom.coords.z = whatPts[off + 2];
    }

    return LLKA_OK;
}

LLKA_Point LLKA_CC LLKA_centroidPoints(const LLKA_Points *points)
{
    if (points->nPoints == 0)
        return { 0, 0, 0 };

    LLKAInternal::MappedPointsUnaligned pointsMapped(points->raw, 3, points->nPoints);
    const auto ctr = LLKAInternal::centroid(pointsMapped);
    return { ctr.x(), ctr.y(), ctr.z() };
}

LLKA_Point LLKA_CC LLKA_centroidStructure(const LLKA_Structure *stru)
{
    if (stru->nAtoms < 1)
        return { 0, 0, 0 };

    // TODO: Use aligned allocations here to speed things up
    std::unique_ptr<double[]> raw{new double[stru->nAtoms * 3]};

    for (size_t idx = 0; idx < stru->nAtoms; idx++) {
        const auto &coords = stru->atoms[idx].coords;
        const size_t off = 3 * idx;
        raw[off] = coords.x;
        raw[off + 1] = coords.y;
        raw[off + 2] = coords.z;
    }

    LLKAInternal::MappedPointsUnaligned pointsMapped(raw.get(), 3, stru->nAtoms);

    const auto ctr = LLKAInternal::centroid(pointsMapped);
    return { ctr.x(), ctr.y(), ctr.z() };
}

LLKA_RetCode LLKA_CC LLKA_rmsdPoints(const LLKA_Points *a, const LLKA_Points *b, double *rmsd)
{
    if (a->nPoints != b->nPoints)
        return LLKA_E_MISMATCHING_SIZES;

    LLKAInternal::MappedPointsUnaligned mA(a->raw, 3, a->nPoints);
    LLKAInternal::MappedPointsUnaligned mB(b->raw, 3, b->nPoints);

    return LLKAInternal::rmsd(mA, mB, rmsd);
}

LLKA_RetCode LLKA_CC LLKA_rmsdStructures(const LLKA_Structure *a, const LLKA_Structure *b, double *rmsd)
{
    if (a->nAtoms!= b->nAtoms)
        return LLKA_E_MISMATCHING_SIZES;

    auto [aPts, bPts] = LLKAInternal::makePts(a, b);

    LLKAInternal::MappedPointsUnaligned mA(aPts.get(), 3, a->nAtoms);
    LLKAInternal::MappedPointsUnaligned mB(bPts.get(), 3, b->nAtoms);

    return LLKAInternal::rmsd(mA, mB, rmsd);
}

LLKA_RetCode LLKA_CC LLKA_superposePoints(LLKA_Points *what, const LLKA_Points *onto, double *rmsd)
{
    LLKAInternal::MappedPointsUnaligned mWhat(what->raw, 3, what->nPoints);
    LLKAInternal::MappedPointsUnaligned mOnto(onto->raw, 3, onto->nPoints);

    return LLKAInternal::superpose(mWhat, mOnto, rmsd);
}

LLKA_RetCode LLKA_CC LLKA_superposeStructures(LLKA_Structure *what, const LLKA_Structure *onto, double *rmsd)
{
    return LLKAInternal::superposeStructures(what, onto, rmsd);
}

LLKA_RetCode LLKA_CC LLKA_superposeStructuresView(LLKA_Structure *what, const LLKA_StructureView *onto, double *rmsd)
{
    return LLKAInternal::superposeStructures(what, onto, rmsd);
}

LLKA_RetCode LLKA_CC LLKA_superpositionMatrixPoints(const LLKA_Points *what, const LLKA_Points *onto, LLKA_Matrix *matrix)
{
    if (what->nPoints != onto->nPoints)
        return LLKA_E_MISMATCHING_SIZES;
    if (what->nPoints == 0)
        return LLKA_E_INVALID_ARGUMENT;

    LLKAInternal::MappedPointsUnaligned mWhat(what->raw, 3, what->nPoints);
    LLKAInternal::MappedPointsUnaligned mOnto(onto->raw, 3, onto->nPoints);

    return LLKAInternal::superpositionMatrix(mWhat, mOnto, matrix);
}

LLKA_RetCode LLKA_CC LLKA_superpositionMatrixStructures(const LLKA_Structure *what, const LLKA_Structure *onto, LLKA_Matrix *matrix)
{
    return LLKAInternal::superpositionMatrixStructures(what, onto, matrix);
}

LLKA_RetCode LLKA_CC LLKA_superpositionMatrixStructureViews(const LLKA_StructureView *what, const LLKA_StructureView *onto, LLKA_Matrix *matrix)
{
    return LLKAInternal::superpositionMatrixStructures(what, onto, matrix);
}
