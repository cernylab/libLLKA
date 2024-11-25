/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "llka_config.h"

/* Enforce calling convention */
#ifndef LLKA_CC
    #if defined LLKA_PLATFORM_WIN32
        #define LLKA_CC __stdcall
    #elif defined LLKA_PLATFORM_UNIX
        #ifdef LLKA_COMPILER_GCC_LIKE
            #ifdef __i386__
                #define LLKA_CC __attribute__((__cdecl__))
            #else
                #define LLKA_CC
            #endif /* __x86_64__ */
        #else
            #error "Unsupported or misdetected compiler"
        #endif /* LLKA_COMPILER_* */
    #elif defined LLKA_PLATFORM_EMSCRIPTEN
        #define LLKA_CC
    #else
        #error "Unsupported or misdetected target platform"
    #endif /* LLKA_PLATFORM_* */
#endif /* LLKA_CC */

/* Define force-inlined functions */
#ifndef LLKA_FORCE_INLINE
    #if defined LLKA_PLATFORM_WIN32
        #if defined LLKA_COMPILER_MINGW || defined LLKA_COMPILER_MSYS
            #define LLKA_FORCE_INLINE inline __attribute__((always_inline))
        #elif defined LLKA_COMPILER_MSVC
            #define LLKA_FORCE_INLINE __forceinline
        #else
            #error "Unsupported or misdetected compiler"
        #endif /* LLKA_COMPILER_* */
    #elif defined(LLKA_PLATFORM_UNIX) || defined(LLKA_PLATFORM_EMSCRIPTEN)
        #ifdef LLKA_COMPILER_GCC_LIKE
            #define LLKA_FORCE_INLINE inline __attribute__((always_inline))
        #else
            #error "Unsupported or misdetected compiler"
        #endif /* LLKA_COMPILER_* */
    #else
        #error "Unsupported or misdetected target platform"
    #endif /* LLKA_PLATFORM_* */
#endif /* LLKA_FORCE_INLINE */

/* Define non-aliased pointer */
#ifndef LLKA_RESTRICT_PTR
    #if defined LLKA_COMPILER_GCC_LIKE || defined LLKA_COMPILER_MINGW || defined LLKA_COMPILER_MSYS
        #define LLKA_RESTRICT_PTR __restrict__
    #elif defined LLKA_COMPILER_MSVC
        #define LLKA_RESTRICT_PTR __restrict
    #else
        #define LLKA_RESTRICT_PTR
    #endif /* LLKA_COMPILER_* */
#endif /* LLKA_RESTRICT_PTR */

/* Data alignment specifiers */
#ifndef LLKA_DATA_ALIGNMENT
    #define LLKA_DATA_ALIGNMENT
    #if defined LLKA_PLATFORM_WIN32
        #if defined LLKA_COMPILER_MINGW || defined LLKA_COMPILER_MSYS
            #define LLKA_ALIGNED_BEF_16
            #define LLKA_ALIGNED_BEF_32
            #define LLKA_ALIGNED_BEF_64
            #define LLKA_ALIGNED_AFT_16 __attribute__((aligned(16)))
            #define LLKA_ALIGNED_AFT_32 __attribute__((aligned(32)))
            #define LLKA_ALIGNED_AFT_64 __attribute__((aligned(64)))
        #elif defined LLKA_COMPILER_MSVC
            #define LLKA_ALIGNED_BEF_16 __declspec(align(16))
            #define LLKA_ALIGNED_BEF_32 __declspec(align(32))
            #define LLKA_ALIGNED_BEF_64 __declspec(align(64))
            #define LLKA_ALIGNED_AFT_16
            #define LLKA_ALIGNED_AFT_32
            #define LLKA_ALIGNED_AFT_64
        #else
            #error "Unsupported or misdetected compiler"
        #endif /* LLKA_COMPILER_* */
    #elif defined(LLKA_PLATFORM_UNIX) || defined(LLKA_PLATFORM_EMSCRIPTEN)
        #ifdef LLKA_COMPILER_GCC_LIKE
            #define LLKA_ALIGNED_BEF_16
            #define LLKA_ALIGNED_BEF_32
            #define LLKA_ALIGNED_BEF_64
            #define LLKA_ALIGNED_AFT_16 __attribute__((aligned(16)))
            #define LLKA_ALIGNED_AFT_32 __attribute__((aligned(32)))
            #define LLKA_ALIGNED_AFT_64 __attribute__((aligned(64)))
        #else
            #error "Unsupported or misdetected compiler"
        #endif /* LLKA_COMPILER_* */
    #else
        #error "Unsupported or misdetected target platform"
    #endif /* LLKA_PLATFORM_* */
#endif /* LLKA_DATA_ALIGNMENT */

/* Allow for redefinitions of LLKA_API as needed */
#ifdef LLKA_API
    #undef LLKA_API
#endif /* LLKA_API */

/* Export only symbols that are part of the public API */
#if defined LLKA_PLATFORM_WIN32
    #if defined LLKA_DLL_BUILD && !defined(LLKA_IMPORT_INTERNAL)
        #if defined LLKA_COMPILER_MINGW || defined LLKA_COMPILER_MSYS
            #define LLKA_API __attribute__ ((dllexport))
            #define LLKA_CPP_API __attribute__ ((dllexport))
            #define LLKA_INTERNAL_API __attribute__ ((dllexport))
        #elif defined LLKA_COMPILER_MSVC
            #define LLKA_API __declspec(dllexport)
            #define LLKA_CPP_API __declspec(dllexport)
            #define LLKA_INTERNAL_API __declspec(dllexport)
        #else
            #error "Unsupported or misdetected compiler"
        #endif /* LLKA_COMPILER_* */
    #else
        #if defined LLKA_COMPILER_MINGW || defined LLKA_COMPILER_MSYS
            #define LLKA_API __attribute__ ((dllimport))
            #define LLKA_CPP_API __attribute__ ((dllimport))
            #define LLKA_INTERNAL_API __attribute__ ((dllimport))
        #elif defined LLKA_COMPILER_MSVC
            #define LLKA_API __declspec(dllimport)
            #define LLKA_CPP_API __declspec(dllimport)
            #define LLKA_INTERNAL_API __declspec(dllimport)
        #else
            #error "Unsupported or misdetected compiler"
        #endif /* LLKA_COMPILER_* */
    #endif /* LLKA_DLL_BUILD */
#elif defined(LLKA_PLATFORM_UNIX) || defined(LLKA_PLATFORM_EMSCRIPTEN)
    #ifdef LLKA_COMPILER_GCC_LIKE
        #if defined LLKA_DLL_BUILD && !defined(LLKA_IMPORT_INTERNAL)
            #define LLKA_API __attribute__ ((visibility ("default")))
            #define LLKA_CPP_API __attribute__ ((visibility ("default")))
            #define LLKA_INTERNAL_API __attribute__ ((visibility ("default")))
        #else
            #define LLKA_API
            #define LLKA_CPP_API
            #define LLKA_INTERNAL_API
        #endif /* LLKA_DLL_BUILD */
    #else
        #error "Unsupported or misdetected compiler"
    #endif /* LLKA_COMPILER_* */
#else
    #error "Unsupported or misdetected target platform"
#endif /* LLKA_PLATFORM_* */

/* Do not expose "static inline" to C because that is, technically,
 * not supported by C compilers.
 * Besides handling the "static inline" itself we also need a way
 * to suppress "unused function" warning which some compilers
 * may issue.
 */

#if defined LLKA_PLATFORM_WIN32
    #ifdef __cplusplus
        #define LLKA_STATIC_INLINE static inline
    #else
        #define LLKA_STATIC_INLINE static
    #endif /* __cplusplus */
#elif defined(LLKA_PLATFORM_UNIX) || defined(LLKA_PLATFORM_EMSCRIPTEN)
    #ifdef LLKA_COMPILER_GCC_LIKE
        #ifdef __cplusplus
            #define LLKA_STATIC_INLINE static inline
        #else
            #define LLKA_STATIC_INLINE static __attribute__((unused))
        #endif /* __cplusplus */
    #else
        #error "Unsupported or misdetected compiler"
    #endif /* LLKA_COMPILER_* */
#else
    #error "Unsupported or misdetected target platform"
#endif /* LLKA_PLATFORM_* */
/*
 * EXPLANATION:
 * Okay, you have probably just finished examining this header file and now you
 * are asking yourself what the hell is this supposed to do and what was wrong
 * with the guy who wrote it this way? You can find the answer only to the
 * former here.
 *
 * All LLKA modules are intended to work on various platforms and be buildable
 * with different compilers. This is, obviously, a pain in the ass to take care
 * of properly. This header attempts to take care of four issues in a reasonably
 * reusable way, namely:
 *
 * - Target platform detection
 * - Target architecture detection (x86 and x86_64 only at this point)
 * - Compiler detection
 * - Symbol import/export mode.
 *
 * We need to know the target platform and architecture in order to enforce
 * correct calling convention. This is done in the first block of ifdefs
 * and it is pretty easy.
 *
 * Compiler detection is closely related to symbol import/export mode because
 * different compilers understand different directives that govern this.
 * Any symbol tagged with LLKA_API will be marked as global (=exported)
 * in the resulting binary as long as the following conditions are met:
 *
 *  - The binary is being built in DLL mode and LLKA_DLL_BUILD is defined
 *  - Export is not overridden by LLKA_IMPORT_INTERNAL
 *
 * If these conditions are not met, the symbol is marked for import (on Windows + MSVC)
 * or not marked at all (on UNIX or Windows + GCC).
 *
 * Notice that this file does not have include guards. This allows us to
 * switch the LLKA_IMPORT_INTERNAL flag on and off as needed by
 * defining or undefining it and reincluding this header.
 * This - although a bit messy - is useful when we are building a library A
 * that exports symbols through the LLKA_API tag but it also pulls in
 * other library B that uses the same mechanism to export symbols. If the symbols
 * exported by library B were not guarded by an LLKA_IMPORT_INTERNAL block,
 * library A would export symbols from both A and B which is not what we want.
 *
 * Clear? ...yeah, didn't think so.
 */

/* Additional platform-dependent definitions */
#if defined LLKA_PLATFORM_WIN32
    #define LLKA_PathChar wchar_t
    #define LLKA_PathLiteral(str) L##str
    #ifdef LLKA_COMPILER_MSVC
        #define __PRETTY_FUNCTION__ __FUNCSIG__
    #endif /* LLKA_COMPILER_MSVC */
#else
    #define LLKA_PathChar char
    #define LLKA_PathLiteral(str) str
#endif /* LLKA_PLATFORM_ */

/* Workarounds for troublesome compilers */
#ifdef LLKA_HAVE_SAD_LADY_COMPILER
    #define LLKA_SAD_CONSTEXPR_FUNC
    #define LLKA_SAD_CONSTEXPR_VAR const
#else
    #define LLKA_SAD_CONSTEXPR_FUNC constexpr
    #define LLKA_SAD_CONSTEXPR_VAR constexpr
#endif /* LLKA_HAVE_SAD_LADY_COMPILER */
