// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_NTC_CONSTANT_H
#define _LLKA_NTC_CONSTANT_H

#include <llka_ntc.h>

#include <array>
#include <map>
#include <string>

namespace LLKAInternal {

class NtCAverages {
public:
    std::string name;
    double delta_1;
    double epsilon_1;
    double zeta_1;
    double alpha_2;
    double beta_2;
    double gamma_2;
    double delta_2;

    double chi_1;
    double chi_2;

    double CC;
    double NN;
    double mu;

    double nu0first;
    double nu1first;
    double nu2first;
    double nu3first;
    double nu4first;
    double nu0second;
    double nu1second;
    double nu2second;
    double nu3second;
    double nu4second;
};


// This array must follow the order of the items in the LLKA_CANA enum!
inline constexpr std::array<char[4], 14> CANA_NAMES {
    "AAA",
    "AAw",
    "AAu",
    "A-B",
    "B-A",
    "BBB",
    "BBw",
    "B12",
    "BB2",
    "miB",
    "ICL",
    "OPN",
    "SYN",
    "ZZZ"
};

extern const std::map<std::string, size_t> NAME_TO_NTC_INDEX_MAPPING;     /*! < Mapping of NtC names to indices into <tt>NTC_AVERAGES</tt> and <tt>NTC_REFERENCES</tt> */
extern const std::array<NtCAverages, 96> NTC_AVERAGES;                    /*! < Array of averaged torsion angles, cross-residue torsions and distances and nu angles of all NtCs */
extern const std::array<LLKA_Structure, 96> NTC_REFERENCES;               /*! < Array of structures used as reference NtC structures */

} // namespace LLKAInternal

#endif // _LLKA_NTC_CONSTANT_H
