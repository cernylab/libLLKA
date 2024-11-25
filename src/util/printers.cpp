// vim: set sw=4 ts=4 sts=4 expandtab :

#include "printers.hpp"
#include "llka_ntc.h"

auto operator<<(std::ostream &os, const LLKA_Points &pts) -> std::ostream &
{
    for (size_t idx = 0; idx < pts.nPoints; idx++)
        os << pts.points[idx] << "\n";
    return os;
}

auto operator<<(std::ostream &os, const LLKA_Structure &stru) -> std::ostream &
{
    for (size_t idx = 0; idx < stru.nAtoms; idx++)
        os << stru.atoms[idx] << "\n";
    return os;
}

auto operator<<(std::ostream &os, const LLKA_StepMetrics &metrics) -> std::ostream &
{
    os
        << LLKA_dinucleotideTorsionName(LLKA_TOR_DELTA_1, LLKA_FALSE) << " " << metrics.delta_1 << "\n"
        << LLKA_dinucleotideTorsionName(LLKA_TOR_EPSILON_1, LLKA_FALSE) << " " << metrics.epsilon_1 << "\n"
        << LLKA_dinucleotideTorsionName(LLKA_TOR_ZETA_1, LLKA_FALSE) << " " << metrics.zeta_1 << "\n"
        << LLKA_dinucleotideTorsionName(LLKA_TOR_ALPHA_2, LLKA_FALSE) << " " << metrics.alpha_2 << "\n"
        << LLKA_dinucleotideTorsionName(LLKA_TOR_BETA_2, LLKA_FALSE) << " " << metrics.beta_2 << "\n"
        << LLKA_dinucleotideTorsionName(LLKA_TOR_GAMMA_2, LLKA_FALSE) << " " << metrics.gamma_2 << "\n"
        << LLKA_dinucleotideTorsionName(LLKA_TOR_DELTA_2, LLKA_FALSE) << " " << metrics.delta_2 << "\n"
        << LLKA_dinucleotideTorsionName(LLKA_TOR_CHI_1, LLKA_FALSE) << " " << metrics.chi_1 << "\n"
        << LLKA_dinucleotideTorsionName(LLKA_TOR_CHI_2, LLKA_FALSE) << " " << metrics.chi_2 << "\n"
        << LLKA_crossResidueMetricName(LLKA_XR_DIST_CC, LLKA_FALSE) << " " << metrics.CC << "\n"
        << LLKA_crossResidueMetricName(LLKA_XR_DIST_NN, LLKA_FALSE) << " " << metrics.NN << "\n"
        << LLKA_crossResidueMetricName(LLKA_XR_TOR_MU, LLKA_FALSE) << " " << metrics.mu << "\n";

    return os;
}
