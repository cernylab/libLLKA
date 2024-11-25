// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_UTIL_GEOMETRY_CPP_H
#define _LLKA_UTIL_GEOMETRY_CPP_H

#include <llka_cpp.h>

namespace LLKAInternal {

template <typename T> requires std::is_floating_point_v<T>
LLKA_INTERNAL_API
auto dihedralAngle(const LLKA::Structure &stru) -> T;

template <typename T> requires std::is_floating_point_v<T>
LLKA_INTERNAL_API
auto spatialDistance(const LLKA::Atom &a, const LLKA::Atom &b) -> T;

} // namespace LLKAInternal

#endif // _LLKA_UTIL_GEOMETRY_CPP_H


