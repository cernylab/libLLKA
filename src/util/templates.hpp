/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_UTIL_TEMPLATES_H
#define _LLKA_UTIL_TEMPLATES_H

#include <llka_structure.h>

#include <algorithm>
#include <cstddef>
#include <string>
#include <type_traits>

namespace LLKAInternal {

template <typename T>
concept LLKAStructureType = std::is_same_v<T, LLKA_Structure> || std::is_same_v<T, LLKA_StructureView>;

template <size_t N>
struct StringLiteral {
    constexpr StringLiteral(const char (&str)[N]) {
        std::copy_n(str, N, value);
    }
    char value[N];
};

template <typename T>
struct Stringifier {
    static auto call(const T &v) -> std::string
    {
        if constexpr (std::is_arithmetic_v<T>) {
            return std::to_string(v);
        } else {
            return v;
        }
    }
};

inline
LLKA_Atom & getAtom(const LLKA_Structure &stru, size_t idx)
{
    return stru.atoms[idx];
}

inline
const LLKA_Atom & getAtom(const LLKA_StructureView &stru, size_t idx)
{
    return *stru.atoms[idx];
}

inline
LLKA_Atom * getAtomPtr(const LLKA_Structure &stru, size_t idx)
{
    return &stru.atoms[idx];
}

inline
const LLKA_Atom * getAtomPtr(const LLKA_StructureView &stru, size_t idx)
{
    return stru.atoms[idx];
}

} // namespace LLKAInternal

#endif // _LLKA_UTIL_TEMPLATES_H
