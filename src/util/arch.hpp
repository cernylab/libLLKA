// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_UTIL_ARCH_HPP
#define _LLKA_UTIL_ARCH_HPP

#include <cstddef>

namespace LLKAInternal {

enum class _MemoryArchType {
    _32BIT,
    _64BIT_LP64,
    _64BIT_LLP64
};

// Determine memory architecture type

template <size_t sizeofPtr, size_t sizeofUL> struct _MemoryArchTypeGetter{};

template <> struct _MemoryArchTypeGetter<4, 4> {
    static constexpr _MemoryArchType TYPE = _MemoryArchType::_32BIT;
};

template <> struct _MemoryArchTypeGetter<8, 8> {
    static constexpr _MemoryArchType TYPE = _MemoryArchType::_64BIT_LP64;
};

template <> struct _MemoryArchTypeGetter<8, 4> {
    static constexpr _MemoryArchType TYPE = _MemoryArchType::_64BIT_LLP64;
};

template <_MemoryArchType> struct _ArraySize {};

template <> struct _ArraySize<_MemoryArchType::_32BIT> {
    static constexpr size_t MAX = -1UL;
};

template <> struct _ArraySize<_MemoryArchType::_64BIT_LP64> {
    static constexpr size_t MAX = -1UL;
};

template <> struct _ArraySize<_MemoryArchType::_64BIT_LLP64> {
    static constexpr size_t MAX = -1ULL; // Some compilers may give a spurious warning here.
};

struct Arch {
    using ArraySize = _ArraySize<_MemoryArchTypeGetter<sizeof(void *), sizeof(unsigned long)>::TYPE>;
    static constexpr size_t INVALID_SIZE_T = _ArraySize<_MemoryArchTypeGetter<sizeof(void *), sizeof(unsigned long)>::TYPE>::MAX;
};

} // namespace LLKAInternal

#endif // _LLKA_UTIL_ARCH_HPP
