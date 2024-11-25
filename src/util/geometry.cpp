// vim: set sw=4 ts=4 sts=4 expandtab :

#include "geometry.h"

#include <Eigen/Dense>

#include <cmath>

namespace LLKAInternal {

auto centroidify(const LLKA_Point &centroid, LLKA_Structure &stru) -> void
{
    const auto cX = centroid.x;
    const auto cY = centroid.y;
    const auto cZ = centroid.z;

    for (size_t idx = 0; idx < stru.nAtoms; idx++) {
        stru.atoms[idx].coords.x -= cX;
        stru.atoms[idx].coords.y -= cY;
        stru.atoms[idx].coords.z -= cZ;
    }
}

template <typename T> requires std::is_floating_point_v<T>
auto dihedralAngle(const LLKA_Point &A, const LLKA_Point &B, const LLKA_Point &C, const LLKA_Point &D) -> T
{
    Eigen::Vector3<T> ab{T(A.x - B.x), T(A.y - B.y), T(A.z - B.z)};
    Eigen::Vector3<T> cb{T(C.x - B.x), T(C.y - B.y), T(C.z - B.z)};
    Eigen::Vector3<T> dc{T(D.x - C.x), T(D.y - C.y), T(D.z - C.z)};

    Eigen::Vector3<T> abc = ab.cross(cb);
    Eigen::Vector3<T> bcd = (-cb).cross(dc);

    T norm = std::sqrt(abc.norm() * abc.norm() * bcd.norm() * bcd.norm());
    T dot = abc.dot(bcd) / norm;
    dot = ((-1.0 <= dot && dot <= 1.0) * dot) + ((-1.0 > dot) * -1.0) + ((dot > 1.0) * 1.0); // Clamp to <-1; 1>
    T ang = std::acos(dot);

    Eigen::Vector3<T> aux = abc.cross(bcd);
    T auxDot = cb.dot(aux);

    assert(!std::isnan(auxDot));
    assert(!std::isnan(ang));

    return sign(auxDot) * ang;
}
template float dihedralAngle<float>(const LLKA_Point &A, const LLKA_Point &B, const LLKA_Point &C, const LLKA_Point &D);
template double dihedralAngle<double>(const LLKA_Point &A, const LLKA_Point &B, const LLKA_Point &C, const LLKA_Point &D);
template long double dihedralAngle<long double>(const LLKA_Point &A, const LLKA_Point &B, const LLKA_Point &C, const LLKA_Point &D);

template <typename T> requires std::is_floating_point_v<T>
auto dihedralAngle(const LLKA_Structure &stru) -> T
{
    if (stru.nAtoms != 4)
        return NAN;

    LLKA_Points pts{
        { new LLKA_Point[4] },
        4
    };

    for (size_t idx = 0; idx < 4; idx++)
        pts.points[idx] = stru.atoms[idx].coords;

    auto angle = dihedralAngle<T>(pts);

    delete[] pts.points;

    return angle;
}
template float dihedralAngle<float>(const LLKA_Structure &stru);
template double dihedralAngle<double>(const LLKA_Structure &stru);
template long double dihedralAngle<long double>(const LLKA_Structure &stru);

} // namespace LLKAInternal
