// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_UTIL_GEOMETRY_H
#define _LLKA_UTIL_GEOMETRY_H

#include <llka_structure.h>

#include "elementaries.h"

#include <cassert>

#ifdef LLKA_USE_SIMD_X86
    #include <xmmintrin.h>
#elif defined(LLKA_USE_SIMD_WASM)
    #include <wasm_simd128.h>
#endif // LLKA_USE_SIMD_*

namespace LLKAInternal {

template <typename T> requires std::is_floating_point_v<T>
inline constexpr
auto angleDifference(const T &a, const T &b)
{
    auto clampA = a - LLKA_TWO_PI * int(a * LLKA_INV_TWO_PI);
    auto clampB = b - LLKA_TWO_PI * int(b * LLKA_INV_TWO_PI);
    auto diff = clampA - clampB;

    if (diff > M_PI || diff < -M_PI)
        return diff - sign(diff) * TWO_PI;
    else
        return diff;
}

#ifdef LLKA_USE_SIMD_X86

typedef double LLKA_ALIGNED_BEF_16 VD[2] LLKA_ALIGNED_AFT_16;

static constexpr VD V_TWO_PI = { LLKA_TWO_PI, LLKA_TWO_PI };
static constexpr VD V_INV_TWO_PI = { LLKA_INV_TWO_PI, LLKA_INV_TWO_PI };

inline
auto angleDifference(const double a, const double b) -> double
{
    VD out = { 0, 0 };

    __m128d w = _mm_set_pd(b, a);
    __m128d w2 = _mm_load_pd(V_INV_TWO_PI);
    __m128d w3 = _mm_load_pd(V_TWO_PI);

    w2 = _mm_mul_pd(w, w2);
    __m128i iw = _mm_cvttpd_epi32(w2);
    w2 = _mm_cvtepi32_pd(iw);

    w2 = _mm_mul_pd(w3, w2);

    w = _mm_sub_pd(w, w2);

    _mm_store_pd(out, w);
    double diff = out[0] - out[1];

    if (diff > M_PI || diff < -M_PI) {
        return diff - sign(diff) * TWO_PI;
    } else {
        return diff;
    }
}

#elif defined(LLKA_USE_SIMD_WASM)

typedef double LLKA_ALIGNED_BEF_16 VD[2] LLKA_ALIGNED_AFT_16;

static constexpr VD V_TWO_PI = { LLKA_TWO_PI, LLKA_TWO_PI };
static constexpr VD V_INV_TWO_PI = { LLKA_INV_TWO_PI, LLKA_INV_TWO_PI };

inline
auto angleDifference(const double a, const double b) -> double
{
    VD out = { 0, 0 };

    v128_t w = wasm_f64x2_make(a, b);
    v128_t w2 = wasm_v128_load(V_INV_TWO_PI);
    v128_t w3 = wasm_v128_load(V_TWO_PI);

    w2 = wasm_f64x2_mul(w, w2);
    w2 = wasm_f64x2_trunc(w2);

    w2 = wasm_f64x2_mul(w3, w2);

    w = wasm_f64x2_sub(w, w2);

    wasm_v128_store(out, w);
    double diff = out[0] - out[1];

    if (diff > M_PI || diff < -M_PI) {
        return diff - sign(diff) * TWO_PI;
    } else {
        return diff;
    }
}

#else

inline constexpr
auto angleDifference(const double a, const double b)
{
    double clampA = a - TWO_PI * int(a * INV_TWO_PI);
    double clampB = b - TWO_PI * int(b * INV_TWO_PI);
    double diff = clampA - clampB;

    if (diff > M_PI || diff < -M_PI)
        return diff - sign(diff) * TWO_PI;
    else
        return diff;
}

#endif // LLKA_USE_SIMD_*

template <typename T> requires std::is_floating_point_v<T>
inline constexpr
auto angleAsFull(const T &ang)
{
    T neg = ang < T(0);

    return neg * (TWO_PI + ang) + (T(1) - neg) * ang;
}

template <typename T> requires std::is_floating_point_v<T>
inline constexpr
auto angleAsRange(const T &ang)
{
    T outward = ang > T{M_PI};

    return outward * (ang - T{LLKA_TWO_PI}) + (T{1} - outward) * ang;
}

auto centroidify(const LLKA_Point &centroid, LLKA_Structure &stru) -> void;

template <typename T> requires std::is_floating_point_v<T>
inline constexpr
auto clampAngle(const T &a)
{
    return a - LLKA_TWO_PI * int(a * LLKA_INV_TWO_PI);
}

inline constexpr
auto clampAngle(const double a)
{
    return a - TWO_PI * int(a * INV_TWO_PI);
}

template <typename T> requires std::is_floating_point_v<T>
LLKA_INTERNAL_API
auto dihedralAngle(const LLKA_Point &A, const LLKA_Point &B, const LLKA_Point &C, const LLKA_Point &D) -> T;

template <typename T> requires std::is_floating_point_v<T>
LLKA_INTERNAL_API
auto dihedralAngle(const LLKA_Structure &stru) -> T;

template <typename T> requires std::is_floating_point_v<T>
inline
auto dihedralAngle(const LLKA_Points &pts) -> T
{
    assert(pts.nPoints >= 4);

    return dihedralAngle<T>(pts.points[0], pts.points[1], pts.points[2], pts.points[3]);
}

template <typename T> requires std::is_floating_point_v<T>
inline
auto dihedralAngle(const LLKA_StructureView &view) -> T
{
    assert(view.nAtoms == 4);

    return dihedralAngle<T>(view.atoms[0]->coords, view.atoms[1]->coords, view.atoms[2]->coords, view.atoms[3]->coords);
}

template <typename T> requires std::is_floating_point_v<T>
inline
auto spatialDistance(const LLKA_Point &a, const LLKA_Point &b) -> T
{
    T dX = b.x - a.x;
    T dY = b.y - a.y;
    T dZ = b.z - a.z;

    return std::sqrt(dX*dX + dY*dY + dZ*dZ);
}

template <typename T> requires std::is_floating_point_v<T>
inline
auto spatialDistance(const LLKA_Atom &a, const LLKA_Atom &b) -> T
{
    return spatialDistance<T>(a.coords, b.coords);
}

template <typename T> requires std::is_floating_point_v<T>
inline
auto angle(const LLKA_Point &a, const LLKA_Point &b, const LLKA_Point &c) -> T
{
    T dot =
        (T(a.x - b.x) * T(c.x - b.x)) +
        (T(a.y - b.y) * T(c.y - b.y)) +
        (T(a.z - b.z) * T(c.z - b.z));

    T magBA = spatialDistance<T>(a, b);
    T magCB = spatialDistance<T>(b, c);

    T acos = dot / (magBA * magCB);
    acos = ((T(-1) <= acos && acos <= T(1)) * acos) + ((T(-1) > acos) * T(-1)) + ((acos > T(1)) * T(1)); // Clamp to <-1; 1>

    return std::acos(acos);
}

} // namespace LLKAInternal

#endif // _LLKA_UTIL_GEOMETRY_H
