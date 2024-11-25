// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_KABSCH_HPP
#define _LLKA_KABSCH_HPP

#include <util/elementaries.h>

#include <Eigen/Dense>

namespace LLKAInternal {

template <typename MA, typename MB>
static
auto kabsch(const MA &a, const MB &b) -> Eigen::Matrix3d
{
    // NOTE:
    // This implementation of the Kabsch algorithm stores the coordinates
    // of each superposed point in columns instead of rows. This enables
    // seamless mapping to column-major ordered arrays preferred by Eigen
    Eigen::Matrix3d H = a * b.transpose();

    Eigen::JacobiSVD<Eigen::MatrixXd> svd{H, Eigen::ComputeThinU | Eigen::ComputeThinV};
    auto U = svd.matrixU();
    auto V = svd.matrixV();
    auto transpU = U.transpose();

    auto det = (V * transpU).determinant();
    auto sDet = sign(det);

    Eigen::Matrix3d F{
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, sDet }
    };

    return V * F * transpU;
}

} // namespace LLKAInternal

#endif // _LLKA_KABSCH_HPP
