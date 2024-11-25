/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_UTIL_H
#define _LLKA_UTIL_H

#include "llka_module.h"

#ifdef __cplusplus
    #include <cmath>
#else
    #include <math.h>
#endif /* __cplusplus */

LLKA_STATIC_INLINE
double LLKA_deg2rad(double x) {
    return M_PI * x / 180.0;
}

LLKA_STATIC_INLINE
float LLKA_deg2radf(float x) {
    return (float)(M_PI) * x / 180.0f;
}

LLKA_STATIC_INLINE
double LLKA_rad2deg(double x) {
    return 180.0 * x / M_PI;
}

LLKA_STATIC_INLINE
float LLKA_rad2degf(float x) {
    return 180.0f * x / (float)(M_PI);
}

LLKA_STATIC_INLINE
double LLKA_fullAngleFromRad(double x) {
    double neg = x < 0.0;

    return neg * ((2.0*M_PI) + x) + (1.0 - neg) * x;
}

LLKA_STATIC_INLINE
float LLKA_fullAngleFromRadf(float x) {
    float neg = x < 0.0f;

    return neg * ((2.0f * (float)M_PI) + x) + (1.0f - neg) * x;
}

LLKA_STATIC_INLINE
double LLKA_fullAngleFromDeg(double x) {
    double neg = x < 0.0;

    return neg * (360.0 + x) + (1.0 - neg) * x;
}

LLKA_STATIC_INLINE
float LLKA_fullAngleFromDegf(float x) {
    float neg = x < 0.0f;

    return neg * (360.0f + x) + (1.0f - neg) * x;
}

#endif /* _LLKA_UTIL_H */
