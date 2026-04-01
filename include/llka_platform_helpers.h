/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_PLATFORM_HELPERS_H
#define _LLKA_PLATFORM_HELPERS_H

#include <llka_module.h>

#ifdef LLKA_PLATFORM_WIN32
    #include <windows.h>

    #define LLKA_PHLP_STRCMP wcscmp
    #define LLKA_PHLP_STRLEN wcslen
    #define LLKA_PHLP_STRNCPY wcsncpy
    #define LLKA_PHLP_FOPEN(path, mode) _wfopen(path, L##mode)
    #define LLKA_PHLP_GETENV(variable_name) _wgetenv(L##variable_name)
    #define LLKA_PHLP_SNPRINTF swprintf
#else
    #define LLKA_PHLP_STRCMP strcmp
    #define LLKA_PHLP_STRLEN strlen
    #define LLKA_PHLP_STRNCPY strncpy
    #define LLKA_PHLP_FOPEN fopen
    #define LLKA_PHLP_GETENV getenv
    #define LLKA_PHLP_SNPRINTF snprintf
#endif

#endif /* _LLKA_PLATFORM_HELPERS_H */
