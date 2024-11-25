/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_EFFEDUP_HPP
#define _LLKA_EFFEDUP_HPP

#include <llka_main.h>

#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace EffedUp {

static inline
auto cmpFlt(double actual, double expected) -> bool
{
    double TOL = 1.0e-11;
    auto diff = std::abs(expected - actual);
    return diff <= TOL;
}

template <typename T1, typename T2>
inline
auto compare(const T1 &lhs, const T2 &rhs) -> bool { return lhs == rhs; }

inline
auto compare(const char *lhs, const char *rhs) -> bool { return std::strcmp(lhs, rhs) == 0; }

inline
auto compare(char *lhs, char *rhs) -> bool { return std::strcmp(lhs, rhs) == 0; }

inline
auto compare(char *lhs, const char *rhs) -> bool { return std::strcmp(lhs, rhs) == 0; }

inline
auto compare(const char *lhs, char *rhs) -> bool { return std::strcmp(lhs, rhs) == 0; }

template <typename T>
inline
auto stringify(const T &v) -> std::string
{
    std::ostringstream oss{};
    oss << v;
    return oss.str();
}

template <>
inline auto stringify<LLKA_RetCode>(const LLKA_RetCode &tRet) -> std::string { return LLKA_errorToString(tRet); }
template <>
inline auto stringify<bool>(const bool &b) -> std::string { return b ? "true" : "false"; }
template <>
inline auto stringify<LLKA_Bool>(const LLKA_Bool &v) -> std::string
{
    if (v == LLKA_TRUE)
        return "true";
    else if (v == LLKA_FALSE)
        return "false";

    throw std::logic_error{"Unknown value of LLKA_Bool type"};
}

template <typename Result1, typename Result2>
static
void prnFail(const Result1 &actual, const Result2 &expected, const std::string &message)
{
    std::cout
        << "### Check has FAILED!\n"
        << "### " << message << "\n"
        << "### actual result   = " << stringify(actual) << "\n"
        << "### expected result = " << stringify(expected) << std::endl;
}

void prnFail(const std::string &message);

template <typename Flt>
static
auto prnFltFail(Flt actual, Flt expected, const std::string &msg)
{
    std::ostringstream oss;
    if (!msg.empty())
        oss << "### " << msg << "\n";
    oss << std::setprecision(15) << "### Difference between expected [" << expected << "] and actual [" << actual << "] values exceeds tolerance";
    std::cout << oss.str() << std::endl;
}

} // namespace EffedUp

#define EFF_cmpFlt(ac, ex, msg) \
    if (!EffedUp::cmpFlt((ac), ex)) { \
        std::cout << "### Aborting at " << __PRETTY_FUNCTION__ << ": " << __LINE__ << std::endl; \
        EffedUp::prnFltFail(ac, ex, msg); \
        std::abort(); \
    }

#define EFF_expect(ac, ex, msg) \
    if (!EffedUp::compare((ac), ex)) { \
        std::cout << "### Aborting at " << __PRETTY_FUNCTION__ << ": " << __LINE__ << std::endl; \
        EffedUp::prnFail(ac, ex, msg); \
        std::abort(); \
    }

#define EFF_fail(msg) \
    do { \
        std::cout << "### Aborting at " << __PRETTY_FUNCTION__ << ": " << __LINE__ << std::endl; \
        EffedUp::prnFail(msg); \
        std::abort(); \
    } while (false)

#endif // _LLKA_EFFEDUP_HPP
