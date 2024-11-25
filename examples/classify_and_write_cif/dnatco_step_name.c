/* vim: set sw=4 ts=4 sts=4 expandtab : */

/****************************************************************

  This module contains utility functions that construct DNATCO
  step name. The logic is a bit lengthy so it was put into
  a separate file to keep the important example code brief
  and easy to follow.

****************************************************************/

#include <llka_structure.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static
void toLwrEntryId(char *id)
{
    id[1] += (id[1] >= 'A' && id[1] <= 'Z') * 32;
    id[2] += (id[2] >= 'A' && id[2] <= 'Z') * 32;
    id[3] += (id[3] >= 'A' && id[3] <= 'Z') * 32;
}

static
void stepAltIds(const LLKA_Structure *stru, char *altIdFirst, char *altIdSecond)
{
    size_t idx;
    int32_t seqIdFirst;

    if (stru->nAtoms == 0) {
        fprintf(stderr, "Step has no atoms, aborting.\n");
        abort(); /* Something must have gone very wrong if we ended up here */
    }

    seqIdFirst = stru->atoms[0].label_seq_id;
    *altIdFirst = stru->atoms[0].label_alt_id;

    /* Check if the first nucleotide in the step has any altId */
    for (idx = 1; idx < stru->nAtoms; idx++) {
        if (*altIdFirst != '\0' || stru->atoms[idx].label_seq_id != seqIdFirst)
            break;

        *altIdFirst = stru->atoms[idx].label_alt_id;
    }

    /* Scroll to the second nucleotide */
    while (idx < stru->nAtoms && stru->atoms[idx].label_seq_id == seqIdFirst)
        idx++;

    if (idx == stru->nAtoms) {
        fprintf(stderr, "Step seems to have only one nucleotide. That don't make no sense... aborting.\n");
        abort();
    }

    *altIdSecond = stru->atoms[idx].label_alt_id;
    idx++;

    /* Check if the second nucleotide has any any altId */
    for (; idx < stru->nAtoms; idx++) {
        if (*altIdSecond != '\0')
            break;
        *altIdSecond = stru->atoms[idx].label_alt_id;
    }
}

char * deriveDNATCOStepName(const char *entryId, int hasMultipleModels, const LLKA_Structure *step)
{
    char altIdFirst, altIdSecond;
    int strLen, strLenModels;
    char *str;
    char *strPtr;
    char lwrEntryId[5];
    const LLKA_Atom *firstAtom;
    const LLKA_Atom *lastAtom;

    /* Proper code should check that the entryId really is a valid PDB ID */
    strcpy(lwrEntryId, entryId);
    toLwrEntryId(lwrEntryId);

    stepAltIds(step, &altIdFirst, &altIdSecond);

    firstAtom = &step->atoms[0];
    lastAtom = &step->atoms[step->nAtoms - 1];

    strLenModels = hasMultipleModels ? snprintf(NULL, 0, "-m%d", firstAtom->pdbx_PDB_model_num) : 0;

    /* Calculate the maximum required length of the string */
    strLen = snprintf(
        NULL,
        0,
        "%s_%s_%s.x_%d.%s_%s.x_%d.%s",  /* Use placeholder characters to account for altIds should they be needed */
        lwrEntryId,                   /* PDB ID */
        firstAtom->auth_asym_id,      /* Chain */
        firstAtom->auth_comp_id,      /* First base name */
        /* First altId placeholder */
        firstAtom->auth_seq_id,       /* First author residue number */
        firstAtom->pdbx_PDB_ins_code, /* First insertion code (may be skipped if it is empty ) */
        lastAtom->auth_comp_id,
        /* Second altId placeholder */
        lastAtom->auth_seq_id,        /* Second author residue number */
        lastAtom->pdbx_PDB_ins_code   /* Second insertion code placeholder (may be skipped if it is empty) */
    );

    /* Write out the step name in pieces */
    str = malloc(strLen + strLenModels + 1);
    strPtr = str;
    /* PDB ID */
    strPtr += sprintf(
        strPtr,
        "%s",
        lwrEntryId
    );
    /* Model tag, if there are multiple models */
    if (hasMultipleModels)
        strPtr += sprintf(strPtr, "-m%d", firstAtom->pdbx_PDB_model_num);
    /* PDB ID, Chain, first base name */
    strPtr += sprintf(
        strPtr,
        "_%s_%s",
        firstAtom->auth_asym_id, firstAtom->auth_comp_id
    );
    /* altId, if any */
    if (altIdFirst != '\0') {
        *strPtr++ = '.';
        *strPtr++ = altIdFirst;
    }
    /* First author residue number */
    strPtr += sprintf(strPtr, "_%d", firstAtom->auth_seq_id);
    /* First insertion code */
    if (strcmp(firstAtom->pdbx_PDB_ins_code, LLKA_NO_INSCODE)) {
        strPtr += sprintf(
            strPtr,
            ".%s",
            firstAtom->pdbx_PDB_ins_code
        );
    }

    /* Second base */
    strPtr += sprintf(strPtr, "_%s", lastAtom->auth_comp_id);
    /* Second altId, if any */
    if (altIdSecond != '\0') {
        *strPtr++ = '.';
        *strPtr++ = altIdSecond;
    }
    /* Second author residue number */
    strPtr += sprintf(strPtr, "_%d", lastAtom->auth_seq_id);
    /* Second insertion code */
    if (strcmp(lastAtom->pdbx_PDB_ins_code, LLKA_NO_INSCODE)) {
        strPtr += sprintf(
            strPtr,
            ".%s",
            lastAtom->pdbx_PDB_ins_code
        );
    }
    *strPtr = '\0';

    return str;
}

int structureHasMultipleModels(const LLKA_Structures *steps)
{
    int32_t modelNo;
    size_t idx;

    if (steps->nStrus == 0)
        return 0;
    if (steps->strus[0].nAtoms == 0)
        return 0;

    modelNo = steps->strus[0].atoms[0].pdbx_PDB_model_num;
    for (idx = 1; idx < steps->nStrus; idx++) {
        const LLKA_Structure *stru = &steps->strus[idx];
        if (stru->nAtoms == 0)
            continue;
        if (stru->atoms[0].pdbx_PDB_model_num != modelNo)
            return 1;
    }

    return 0;
}
