// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_UTIL_PRINTERS_HPP
#define _LLKA_UTIL_PRINTERS_HPP

#include <llka_ntc.h>
#include <llka_structure.h>

#include <ostream>

inline
auto operator<<(std::ostream &os, const LLKA_Point &pt) -> std::ostream &
{
    os << "[ " << pt.x << "; " << pt.y << "; " << pt.z << " ]";
    return os;
}

LLKA_INTERNAL_API
auto operator<<(std::ostream &os, const LLKA_Points &pts) -> std::ostream &;

inline
auto operator<<(std::ostream &os, const LLKA_Atom &atom) -> std::ostream &
{
    os
        << atom.label_asym_id << " "
        << atom.label_comp_id << " "
        << atom.label_seq_id << "(" << atom.auth_seq_id << ") "
        << (atom.label_alt_id ? std::string{"["} + atom.label_alt_id + std::string{"] "} : "")
        << atom.label_atom_id
        << atom.coords;
    return os;
}

LLKA_INTERNAL_API
auto operator<<(std::ostream &os, const LLKA_Structure &stru) -> std::ostream &;

LLKA_INTERNAL_API
auto operator<<(std::ostream &os, const LLKA_StepMetrics &metrics) -> std::ostream &;

#endif // _LLKA_UTIL_PRINTERS_HPP
