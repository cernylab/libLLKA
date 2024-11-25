// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_UTIL_ELEMENTARIES_H
#define _LLKA_UTIL_ELEMENTARIES_H

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#define LLKA_TWO_PI (2*M_PI)
#define LLKA_INV_TWO_PI (1.0 / (2.0*M_PI))

#define LLKA_WITHIN_EXCLUSIVE(low, val, high) ((val > low) && (val < high))
#define LLKA_WITHIN_INCLUSIVE(low, val, high) ((val >= low) && (val <= high))
#define LLKA_WITHIN_LINCL_REXCL(low, val, high) ((val >= low) && (val < high))
#define LLKA_WITHIN_LECXL_RINCL(low, val, high) ((val > low) && (val <= high))

namespace LLKAInternal {

inline constexpr double TWO_PI{LLKA_TWO_PI};
inline constexpr double INV_TWO_PI{LLKA_INV_TWO_PI};

template <typename T>
inline constexpr
auto R2D(const T &x)
{
    return x * T{180}/T{M_PI};
}

template <typename T>
inline constexpr
auto D2R(const T &x)
{
    return x * T{M_PI}/T{180};
}

template <typename V, template <typename...> typename C>
inline
auto contains(const C<V> &container, const V &value) -> bool
{
    return std::find(container.cbegin(), container.cend(), value) != container.cend();
}

template <typename K, typename V, template <typename...> typename C>
inline
auto containsKey(const C<K, V> &container, const K &key) -> bool
{
    return container.find(key) != container.cend();
}

template <typename V, typename E, typename Cmp, template <typename...> typename C>
inline
auto contains(const C<E> &container, const V &value, const Cmp &comparator) -> bool
{
    const auto predicate = [&value, &comparator](const E &item) {
        return comparator(value, item);
    };

    return std::find_if(container.cbegin(), container.cend(), predicate) != container.cend();
}

template <typename K, typename V, typename E, typename Cmp, template <typename...> typename C>
inline
auto contains(const C<K, E> &container, const K &value, const Cmp &comparator) -> bool
{
    const auto predicate = [&value, &comparator](const E &item) {
        return comparator(value, item);
    };

    return std::find_if(container.cbegin(), container.cend(), predicate) != container.cend();
}

inline
auto dequote(const std::string &s)
{
    if (s.length() < 2)
        return s;

    char quoteChar = s[0];
    if (!(quoteChar == '"' || quoteChar == '\''))
        return s;

    if (!s.ends_with(quoteChar))
        return s;

    return s.substr(1, s.length() - 2);
}

template <typename E, typename Pred, template <typename...> typename C>
inline
auto filter(C<E> &container, const Pred &predicate)
{
    auto it = std::remove_if(container.begin(), container.end(), predicate);
    container.erase(it, container.end());
}

template <typename T>
inline
auto compareWithTolerance(const T &lhs, const T &rhs, const T &tolerance)
{
    assert(tolerance >= T{0});

    auto diff = std::abs(lhs - rhs);
    return diff <= tolerance;
}

inline
auto destroyString(const char *str)
{
    delete [] str;
}

inline
auto duplicateString(const char *str)
{
    const size_t len = std::strlen(str) + 1;
    char *dup = new char[len];

    std::copy_n(str, len, dup);

    return dup;
}

inline
auto duplicateString(const char *str, const size_t len)
{
    char *dup = new char[len + 1];

    std::copy_n(str, len, dup);
    dup[len] = '\0';

    return dup;
}

inline
auto duplicateString(const std::string &str)
{
    const size_t len = str.size();
    char *dup = new char[len + 1];

    std::copy_n(str.data(), len, dup);
    dup[len] = '\0';

    return dup;
}

template <typename T>
auto NaN() -> T;

template <>
inline
auto NaN<float>() -> float { return std::nanf(""); }

template <>
inline
auto NaN<double>() -> double { return std::nan(""); }

template <>
inline
auto NaN<long double>() -> long double { return std::nanl(""); }

template <typename T>
inline constexpr
auto sign(const T &v) -> T
{
    return (v > T(0)) - (v < T(0));
}

inline
auto sign(const double v) -> double
{
    union ID {
        int64_t i;
        double d;
    };
    static_assert(sizeof(ID::i) == sizeof(ID::d), "Integer and floating-point representations do not have the same size");

    ID asBits;
    ID sgn;

    int64_t signBitMask = 0x8000000000000000;
    asBits.d = v;
    sgn.i = (signBitMask & asBits.i) | (0x3FF0000000000000);

    return sgn.d;
}

inline
auto split(const std::string &str, const char delim, size_t sizeHint = 0)
{
    std::vector<std::string> tokens{};

    if (sizeHint)
        tokens.reserve(sizeHint);

    size_t pos = 0;
    while (true) {
        auto nextPos = str.find_first_of(delim, pos);
        tokens.push_back(str.substr(pos, nextPos - pos));
        if (nextPos != std::string::npos)
            pos = nextPos + 1;
        else
            break;
    }

    return tokens;
}

} // namespace LLKAInternal

#endif // _LLKA_UTIL_ELEMENTARIES_H
