// vim: set sw=4 ts=4 sts=4 expandtab :

#include <llka_cpp.h>

#include "geometry_cpp.h"
#include "geometry.h"

namespace LLKAInternal {

template <typename T> requires std::is_floating_point_v<T>
auto dihedralAngle(const LLKA::Structure &stru) -> T
{
    if (stru.size() != 4)
        return NAN;

    LLKA::Points pts{
        stru[0].coords, stru[1].coords, stru[2].coords, stru[3].coords
    };

    const LLKA_Points cPts{
        { pts.data() },
        pts.size()
    };

    auto angle = dihedralAngle<T>(cPts);

    return angle;
}
template float dihedralAngle<float>(const LLKA::Structure &stru);
template double dihedralAngle<double>(const LLKA::Structure &stru);
template long double dihedralAngle<long double>(const LLKA::Structure &stru);

template <typename T> requires std::is_floating_point_v<T>
auto spatialDistance(const LLKA::Atom &a, const LLKA::Atom &b) -> T
{
    return spatialDistance<T>(a.coords, b.coords);
}
template float spatialDistance<float>(const LLKA::Atom &a, const LLKA::Atom &b);
template double spatialDistance<double>(const LLKA::Atom &a, const LLKA::Atom &b);
template long double spatialDistance<long double>(const LLKA::Atom &a, const LLKA::Atom &b);

} // namespace LLKAInternal

