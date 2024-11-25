// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_SUPERPOSITION_HPP
#define _LLKA_SUPERPOSITION_HPP

#include <llka_main.h>

#include "kabsch.hpp"

#include <cmath>

namespace LLKAInternal {

using MappedPointsUnaligned = Eigen::Map<Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::ColMajor>, Eigen::Unaligned>;
using MappedQuadsUnaligned = Eigen::Map<Eigen::Matrix<double, 4, Eigen::Dynamic, Eigen::ColMajor>, Eigen::Unaligned>;
using MappedMat4x4Unaligned = Eigen::Map<Eigen::Matrix<double, 4, 4, Eigen::ColMajor>, Eigen::Unaligned>;
using MappedPoint = Eigen::Map<Eigen::Vector3d>;
using Mat3x3 = Eigen::Matrix<double, 3, 3>;

template <typename M>
static
auto centroid(const M &coords) -> Eigen::Vector3d
{
    const Eigen::Index nCols = coords.cols();
    Eigen::Vector3d centroid{0, 0, 0};

    for (Eigen::Index idx = 0; idx < nCols; idx++) {
        const auto &col = coords.col(idx);
        centroid.x() += col.x();
        centroid.y() += col.y();
        centroid.z() += col.z();
    }

    centroid.x() = centroid.x() / nCols;
    centroid.y() = centroid.y() / nCols;
    centroid.z() = centroid.z() / nCols;

    return centroid;
}

template <typename T>
static
auto centroidify(T &coords, const Eigen::Vector3d &centroid)
{
    for (Eigen::Index col = 0; col < coords.cols(); col++) {
        coords(0, col) -= centroid.x();
        coords(1, col) -= centroid.y();
        coords(2, col) -= centroid.z();
    }
}

template <typename MA, typename MB>
auto applySuperposition(MA &what, const MB &matrix) -> void
{
    what = matrix * what;
}

template <typename MA, typename MB>
auto rmsd(const MA &a, const MB &b, double *rmsd) -> LLKA_RetCode
{
    if (a.cols() != b.cols())
        return LLKA_E_MISMATCHING_SIZES;

    const Eigen::Index nPoints = a.cols();
    double sum = 0;
    for (Eigen::Index idx = 0; idx < nPoints; idx++) {
        const auto &colA = a.col(idx);
        const auto &colB = b.col(idx);

        Eigen::Vector3d diff = colA - colB;

        sum += diff.dot(diff);
    }

    *rmsd = std::sqrt(sum / nPoints);

    return LLKA_OK;
}

template <typename MA, typename MB>
auto superpose(MA &what, const MB &onto, double *_rmsd) -> LLKA_RetCode
{
    if (what.cols() != onto.cols())
        return LLKA_E_MISMATCHING_SIZES;

    Eigen::Matrix<double, 3, Eigen::Dynamic> cOnto = onto;
    auto centroidWhat = centroid(what);
    auto centroidOnto = centroid(onto);

    centroidify(what, centroidWhat);
    centroidify(cOnto, centroidOnto);

    Mat3x3 rot = kabsch(what, cOnto);
    what = rot * what;
    rmsd(what, cOnto, _rmsd);

    centroidify(what, -centroidOnto);

    return LLKA_OK;
}

template <typename MA, typename MB>
auto superpositionMatrix(const MA &what, const MB &onto, LLKA_Matrix *matrix) -> LLKA_RetCode
{
    using Mat4x4 = Eigen::Matrix<double, 4, 4>;

    if (what.cols() != onto.cols())
        return LLKA_E_MISMATCHING_SIZES;

    Eigen::Matrix<double, 3, Eigen::Dynamic> cWhat = what;
    Eigen::Matrix<double, 3, Eigen::Dynamic> cOnto = onto;
    auto centroidWhat = centroid(what);
    auto centroidOnto = centroid(onto);

    centroidify(cWhat, centroidWhat);
    centroidify(cOnto, centroidOnto);
    Mat3x3 rot = kabsch(cWhat, cOnto);

    // Translate to origin
    Mat4x4 translation = Mat4x4::Identity();
    translation(0, 3) = -centroidWhat(0);
    translation(1, 3) = -centroidWhat(1);
    translation(2, 3) = -centroidWhat(2);

    // Apply rotation
    Mat4x4 rotation = Mat4x4::Identity();
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++)
            rotation(r, c) = rot(r, c);
    }

    // Rotate and translate to origin
    Mat4x4 aux = rotation * translation;

    translation(0, 3) = centroidOnto(0);
    translation(1, 3) = centroidOnto(1);
    translation(2, 3) = centroidOnto(2);

    // Translate the product to the target position
    Mat4x4 transformation = translation * aux;

    LLKA_initMatrix(4, 4, matrix);
    std::memcpy(matrix->data, transformation.data(), sizeof(Mat4x4::Scalar) * 4 * 4);

    return LLKA_OK;
}

} // namespace LLKAInternal

#endif // _LLKA_SUPERPOSITION_HPP
