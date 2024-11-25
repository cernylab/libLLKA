/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_SEGMENTATION_H
#define _LLKA_SEGMENTATION_H

#include "llka_structure.h"

typedef struct LLKA_Residue {
    LLKA_Atom **atoms;
    size_t nAtoms;
} LLKA_Residue;
LLKA_IS_POD(LLKA_Residue)

typedef struct LLKA_Chain {
    LLKA_Residue *residues;
    size_t nResidues;

    LLKA_Atom **atoms;
    size_t nAtoms;
} LLKA_Chain;
LLKA_IS_POD(LLKA_Chain)

typedef struct LLKA_Model {
    LLKA_Chain *chains;
    size_t nChains;

    LLKA_Atom **atoms;
    size_t nAtoms;
} LLKA_Model;
LLKA_IS_POD(LLKA_Model)

typedef struct LLKA_StructureSegments {
    LLKA_Model *models;
    size_t nModels;
} LLKA_StructureSegments;
LLKA_IS_POD(LLKA_StructureSegments)

LLKA_BEGIN_API_FUNCTIONS

LLKA_API void LLKA_CC LLKA_destroyStructureSegments(const LLKA_StructureSegments *segs);
LLKA_API LLKA_StructureSegments LLKA_CC LLKA_structureSegments(LLKA_Structure *stru);

LLKA_END_API_FUNCTIONS

#endif /*_LLKA_SEGMENTATION_H */
