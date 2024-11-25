/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef MINICIF_P_H
#define MINICIF_P_H

#include <llka_minicif.h>

struct LLKA_CifDataPrivate {
    bool tainted;
};

struct LLKA_CifDataBlockPrivate {
    LLKA_CifData *root;
};

struct LLKA_CifDataCategoryPrivate {
    LLKA_CifDataCategory *prev;
    LLKA_CifDataCategory *next;
    LLKA_CifData *root;
    bool isLoop;
};

struct LLKA_CifDataItemPrivate {
    LLKA_CifDataItem *prev;
    LLKA_CifDataItem *next;
    LLKA_CifData *root;
};

#endif // MINICIF_P_H
