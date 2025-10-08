#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <json.hpp>
#include <CLI11.hpp>

#include <ntc_constants.h>
#include <llka_structure.h>
#include <llka_connectivity_similarity.h>
#include <llka_ntc.h>
#include <ntc.hpp>
#include <llka_classification.h>
#include <llka_minicif.h>
#include <llka_nucleotide.h>
#include <llka_resource_loaders.h>
#include <llka_util.h>

/* We expect _USE_MATH_DEFINES to be defined */
#include <math.h>

#define CLUSTERS_FILE LLKA_PathLiteral("clusters.csv")
#define CONFALS_FILE LLKA_PathLiteral("confals.csv")
#define GOLDEN_STEPS_FILE LLKA_PathLiteral("golden_steps.csv")
#define NU_ANGLES_FILE LLKA_PathLiteral("nu_angles.csv")
#define CONFAL_PERCENTILES_FILE LLKA_PathLiteral("confal_percentiles.csv")

#define MAX_CLOSE_ENOUGH_RMSD 0.5

#define SZ_DETAILS 255

#define SIZEOF_ARRAY(array) (sizeof(array) / sizeof(array[0]))

#define FULL_PATH_SIZE 1024

const LLKA_NtC NTCS[] = { LLKA_AA00, LLKA_AA02, LLKA_AA03, LLKA_AA04, LLKA_AA08, LLKA_AA09, LLKA_AA01, LLKA_AA05, LLKA_AA06, LLKA_AA10, LLKA_AA11, LLKA_AA07, LLKA_AA12, LLKA_AA13,
    LLKA_AB01, LLKA_AB02, LLKA_AB03, LLKA_AB04, LLKA_AB05,
    LLKA_BA01, LLKA_BA05, LLKA_BA09, LLKA_BA08, LLKA_BA10, LLKA_BA13, LLKA_BA16, LLKA_BA17,
    LLKA_BB00, LLKA_BB01, LLKA_BB17, LLKA_BB02, LLKA_BB03, LLKA_BB11, LLKA_BB16, LLKA_BB04, LLKA_BB05, LLKA_BB07, LLKA_BB08, LLKA_BB10, LLKA_BB12, LLKA_BB13, LLKA_BB14, LLKA_BB15, LLKA_BB20,
    LLKA_IC01, LLKA_IC02, LLKA_IC03, LLKA_IC04, LLKA_IC05, LLKA_IC06, LLKA_IC07,
    LLKA_OP01, LLKA_OP02, LLKA_OP03, LLKA_OP04, LLKA_OP05, LLKA_OP06, LLKA_OP07, LLKA_OP08, LLKA_OP09, LLKA_OP10, LLKA_OP11, LLKA_OP12, LLKA_OP13, LLKA_OP14, LLKA_OP15, LLKA_OP16, LLKA_OP17, LLKA_OP18, LLKA_OP19, LLKA_OP20, LLKA_OP21, LLKA_OP22, LLKA_OP23, LLKA_OP24, LLKA_OP25, LLKA_OP26, LLKA_OP27, LLKA_OP28, LLKA_OP29, LLKA_OP30, LLKA_OP31, LLKA_OPS1, LLKA_OP1S,
    LLKA_AAS1, LLKA_AB1S, LLKA_AB2S,
    LLKA_BB1S, LLKA_BB2S, LLKA_BBS1,
    LLKA_ZZ01, LLKA_ZZ02, LLKA_ZZ1S, LLKA_ZZ2S, LLKA_ZZS1, LLKA_ZZS2, LLKA_INVALID_NTC};

const int DISTANCE_PRECISION = 0.001;
const int EULER_PRECISION = 0.1;
const int RMSD_PRECISION = 0.001;

typedef struct LLKA_SortedStep{
    bool havePrevStep = false;
    bool haveNextStep = false;
    LLKA_Structure currStep;
    LLKA_Structure prevStep;
    LLKA_Structure nextStep;
} LLKA_SortedStep;

typedef struct LLKA_SortedSteps{
    LLKA_SortedStep *steps;
    size_t nSteps;
} LLKA_SortedSteps;

typedef struct LLKA_AnalyzedStep{
    const char *auth_asym_id_1;
    const char *label_comp_id_1;
    int32_t PDB_ins_code_1;
    int32_t auth_seq_id_1;
    char label_alt_id_1;
    const char *auth_asym_id_2;
    const char *label_comp_id_2;
    int32_t PDB_ins_code_2;
    int32_t auth_seq_id_2;
    char label_alt_id_2;
} LLKA_AnalyzedStep;

typedef struct LLKA_AnalyzedSteps{
    LLKA_AnalyzedStep* steps;
    size_t nSteps;
} LLKA_AnalyzedSteps;

char* create_full_path(const char *filename) {
    const char *prefix_path = getenv("DNATCO_ASSETS_PATH");

    /* Allocate memory for the full path */
    static char full_path[FULL_PATH_SIZE];

    if (prefix_path != NULL) {
        /* Calculate the length of the prefix path and filename */
        size_t prefix_length = strlen(prefix_path);
        size_t filename_length = strlen(filename);

        /* Check if the combined length exceeds the buffer size */
        if (prefix_length + filename_length + 2 > FULL_PATH_SIZE) {
            fprintf(stderr, "The combined length of the prefix and filename exceeds the buffer size.\n");
            exit(1);  /* Exit if the buffer size is not sufficient */
        }

        /* Check if the last character of the prefix path is a slash */
        if (prefix_path[prefix_length - 1] != '/' && prefix_path[prefix_length - 1] != '\\') {
            snprintf(full_path, sizeof(full_path), "%s/%s", prefix_path, filename);
        } else {
            snprintf(full_path, sizeof(full_path), "%s%s", prefix_path, filename);
        }

    } else {
        strncpy(full_path, filename, sizeof(full_path) - 1);
        full_path[sizeof(full_path) - 1] = '\0'; /* Ensure null-termination */
    }

    return full_path;
}

/* DNATCO step name functions signatures */
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
    str = (char*)malloc(strLen + strLenModels + 1);
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


static
void prepareExtenedCifCategories(
    LLKA_CifDataCategory *overall,
    LLKA_CifDataCategory *ntc,
    LLKA_CifDataCategory *summary,
    LLKA_CifDataCategory *parameters,
    LLKA_CifDataCategory *sugar
)
{
    const char *overall_keywords[] = {
        "entry_id",
        "confal_score",
        "confal_percentile",
        "ntc_version",
        "cana_version",
        "num_steps",
        "num_classified",
        "num_unclassified",
        "num_unclassified_rmsd_close"
    };
    const char *ntc_keywords[] = {
        "id",
        "name",
        "PDB_model_number",
        "label_entity_id_1",
        "label_asym_id_1",
        "label_seq_id_1",
        "label_comp_id_1",
        "label_alt_id_1",
        "label_entity_id_2",
        "label_asym_id_2",
        "label_seq_id_2",
        "label_comp_id_2",
        "label_alt_id_2",
        "auth_asym_id_1",
        "auth_seq_id_1",
        "auth_asym_id_2",
        "auth_seq_id_2",
        "PDB_ins_code_1",
        "PDB_ins_code_2"
    };
    const char *summary_keywords[] = {
        "step_id",
        "assigned_CANA",
        "assigned_NtC",
        "confal_score",
        "euclidean_distance_NtC_ideal",
        "cartesian_rmsd_closest_NtC_representative",
        "closest_CANA",
        "closest_NtC",
        "closest_step_golden"
    };
    const char *parameters_keywords[] = {
        "step_id", "tor_delta_1", "tor_epsilon_1",
        "tor_zeta_1", "tor_alpha_2", "tor_beta_2",
        "tor_gamma_2", "tor_delta_2", "tor_chi_1",
        "tor_chi_2", "dist_NN", "dist_CC",
        "tor_NCCN", "diff_tor_delta_1", "diff_tor_epsilon_1",
        "diff_tor_zeta_1", "diff_tor_alpha_2", "diff_tor_beta_2",
        "diff_tor_gamma_2", "diff_tor_delta_2", "diff_tor_chi_1",
        "diff_tor_chi_2", "diff_dist_NN", "diff_dist_CC",
        "diff_tor_NCCN", "confal_tor_delta_1", "confal_tor_epsilon_1",
        "confal_tor_zeta_1", "confal_tor_alpha_2", "confal_tor_beta_2",
        "confal_tor_gamma_2", "confal_tor_delta_2", "confal_tor_chi_1",
        "confal_tor_chi_2", "confal_dist_NN", "confal_dist_CC",
        "confal_tor_NCCN", "details"
    };
    const char *sugar_keywords[] = {
        "step_id", "P_1", "tau_1",
        "Pn_1", "P_2", "tau_2",
        "Pn_2", "nu_1_1", "nu_1_2",
        "nu_1_3", "nu_1_4", "nu_1_5",
        "nu_2_1", "nu_2_2", "nu_2_3",
        "nu_2_4", "nu_2_5", "diff_nu_1_1",
        "diff_nu_1_2", "diff_nu_1_3", "diff_nu_1_4",
        "diff_nu_1_5", "diff_nu_2_1", "diff_nu_2_2",
        "diff_nu_2_3", "diff_nu_2_4", "diff_nu_2_5"
    };
    size_t idx;

    /* Note that the functions below may return NULL if a new entry cannot be added.
     * Properly written code should check the return values. */

    for (idx = 0; idx < SIZEOF_ARRAY(overall_keywords); idx++)
        LLKA_cifDataCategory_addItem(overall, overall_keywords[idx]);

    for (idx = 0; idx < SIZEOF_ARRAY(ntc_keywords); idx++)
        LLKA_cifDataCategory_addItem(ntc, ntc_keywords[idx]);

    for (idx = 0; idx < SIZEOF_ARRAY(summary_keywords); idx++)
        LLKA_cifDataCategory_addItem(summary, summary_keywords[idx]);

    for (idx = 0; idx < SIZEOF_ARRAY(parameters_keywords); idx++)
        LLKA_cifDataCategory_addItem(parameters, parameters_keywords[idx]);

    for (idx = 0; idx < SIZEOF_ARRAY(sugar_keywords); idx++)
        LLKA_cifDataCategory_addItem(sugar, sugar_keywords[idx]);
}

static
void extendCifData(LLKA_CifData *cifData, const LLKA_ClassifiedSteps *classifiedSteps, const LLKA_AverageConfal *avgConfal, const char *entryId, const LLKA_Structures *steps, int hasMultipleModels)
{
    /* Get the first Cif block. We expect the structure to be in the first block
     * and therefore we want to extend the first block */
    LLKA_CifDataBlock *block = &cifData->blocks[0];

    /* Add new categories to the block. Note that the functions below may return NULL
     * if a new category cannot be added. Properly written code should check the return value. */
    LLKA_CifDataCategory *NdbStructNtCOverall = LLKA_cifDataBlock_addCategory(block, "ndb_struct_ntc_overall");
    LLKA_CifDataCategory *NdbStructNtC = LLKA_cifDataBlock_addCategory(block, "ndb_struct_ntc_step");
    LLKA_CifDataCategory *NdbStructStepSummary = LLKA_cifDataBlock_addCategory(block, "ndb_struct_ntc_step_summary");
    LLKA_CifDataCategory *NdbStructStepParameters = LLKA_cifDataBlock_addCategory(block, "ndb_struct_ntc_step_parameters");
    LLKA_CifDataCategory *NdbStructStepSugar = LLKA_cifDataBlock_addCategory(block, "ndb_struct_sugar_step_parameters");
    /* Helper variables */
    int32_t assignedStepsCount = 0;
    int32_t unassignedButCloseEnoughStepsCount = 0;
    int32_t totalStepsCount = 0; /* Different than the number of steps in the LLKA_ClassifiedSteps object because some steps there might be bogus */
    LLKA_CifDataItem *item;
    LLKA_RetCode tRet;
    size_t idx;
    LLKA_CifDataValue cifValue;
    char *strbuf;
    char *detailsBuf;

    strbuf = (char*) malloc(32); /* 32 characters is enough for everything except details */
    /* MiniCif makes internal copies of data values so it is safe to create
     * one data value object and reuse it */
    cifValue.text = strbuf; /* Work around the fact that "text" is const char * but we need a writable buffer */
    cifValue.state = LLKA_MINICIF_VALUE_SET;

    detailsBuf = (char*) malloc(SZ_DETAILS+1);

    /* Add entries to categories */
    prepareExtenedCifCategories(NdbStructNtCOverall, NdbStructNtC, NdbStructStepSummary, NdbStructStepParameters, NdbStructStepSugar);

    /* Walk the list of assigned steps and extend the Cif data */
    for (idx = 0; idx < classifiedSteps->nAttemptedSteps; idx++) {
        const LLKA_AttemptedClassifiedStep *attempted = &classifiedSteps->attemptedSteps[idx];
        const LLKA_ClassifiedStep *cs = &attempted->step;
        const LLKA_Structure *struStep = &steps->strus[idx];
        const LLKA_Atom *firstAtom = &struStep->atoms[0];
        const LLKA_Atom *lastAtom = &struStep->atoms[struStep->nAtoms - 1];
        char *DNATCOStepName;

        item = NdbStructStepSummary->firstItem;

        /* Attempted step could not have been classified, ignore it */
        if (attempted->status != LLKA_OK)
            continue;

        totalStepsCount++;
        if (cs->violations == LLKA_CLASSIFICATION_OK)
            assignedStepsCount++;
        else if (cs->violations & LLKA_CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH)
            unassignedButCloseEnoughStepsCount++;

        /* Write the values.
         * Here we are taking advantage of the fact that the order of entries
         * in a category matches the order in which they were added. As such,
         * we can simply call _nextItem() to jump to the next entry.
         * In cases where the addition order is not known or if some entries
         * were deleted, use _findEntry() function to get the corresponding
         * entry.
         */

        /* Write step summary row */

        snprintf(strbuf, 32, "%d", totalStepsCount);
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = LLKA_CANAToName(cs->assignedCANA);
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = LLKA_NtCToName(cs->assignedNtC);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%d", (int)(cs->confalScore.total + 0.5));
        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = strbuf;
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->euclideanDistanceNtCIdeal);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.3f", cs->rmsdToClosestNtC);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = LLKA_CANAToName(cs->closestCANA);
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = LLKA_NtCToName(cs->closestNtC);
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        if (strlen(cs->closestGoldenStep) == 0) {
            cifValue.state = LLKA_MINICIF_VALUE_NONE;
            LLKA_cifDataItem_addValue(item, &cifValue);
            cifValue.state = LLKA_MINICIF_VALUE_SET;
        } else {
            cifValue.text = cs->closestGoldenStep;
            LLKA_cifDataItem_addValue(item, &cifValue);
        }

        /* Write NtC row */

        item = NdbStructNtC->firstItem;
        snprintf(strbuf, 32, "%d", totalStepsCount);
        cifValue.text = strbuf;
        LLKA_cifDataItem_addValue(item, &cifValue);

        DNATCOStepName = deriveDNATCOStepName(entryId,  hasMultipleModels, struStep);
        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = DNATCOStepName;
        LLKA_cifDataItem_addValue(item, &cifValue);
        free(DNATCOStepName);

        snprintf(strbuf, 32, "%d", firstAtom->pdbx_PDB_model_num);
        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = strbuf;
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = firstAtom->label_entity_id;
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = firstAtom->label_asym_id;
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%d", firstAtom->label_seq_id);
        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = strbuf;
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = firstAtom->auth_comp_id;
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        if (firstAtom->label_alt_id == LLKA_NO_ALTID) {
            cifValue.state = LLKA_MINICIF_VALUE_NONE;
            LLKA_cifDataItem_addValue(item, &cifValue);
            cifValue.state = LLKA_MINICIF_VALUE_SET;
        } else {
            snprintf(strbuf, 32, "%c", firstAtom->label_alt_id);
            cifValue.text = strbuf;
            LLKA_cifDataItem_addValue(item, &cifValue);
        }

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = lastAtom->label_entity_id;
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = lastAtom->label_asym_id;
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%d", lastAtom->label_seq_id);
        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = strbuf;
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = lastAtom->auth_comp_id;
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        if (lastAtom->label_alt_id == LLKA_NO_ALTID) {
            cifValue.state = LLKA_MINICIF_VALUE_NONE;
            LLKA_cifDataItem_addValue(item, &cifValue);
            cifValue.state = LLKA_MINICIF_VALUE_SET;
        } else {
            snprintf(strbuf, 32, "%c", lastAtom->label_alt_id);
            cifValue.text = strbuf;
            LLKA_cifDataItem_addValue(item, &cifValue);
        }

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = firstAtom->auth_asym_id;
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%d", firstAtom->auth_seq_id);
        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = strbuf;
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = lastAtom->auth_asym_id;
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%d", lastAtom->auth_seq_id);
        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = strbuf;
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        if (strcmp(firstAtom->pdbx_PDB_ins_code, LLKA_NO_INSCODE) == 0) {
            cifValue.state = LLKA_MINICIF_VALUE_NONE;
            LLKA_cifDataItem_addValue(item, &cifValue);
            cifValue.state = LLKA_MINICIF_VALUE_SET;
        } else {
            size_t len = strlen(firstAtom->pdbx_PDB_ins_code);
            char *xbuf = (char*) malloc(len + 1);
            snprintf(xbuf, len + 1, "%s", firstAtom->pdbx_PDB_ins_code);
            cifValue.text = xbuf;
            LLKA_cifDataItem_addValue(item, &cifValue);
            free(xbuf);
        }

        item = LLKA_cifDataCategory_nextItem(item);
        if (strcmp(lastAtom->pdbx_PDB_ins_code, LLKA_NO_INSCODE) == 0) {
            cifValue.state = LLKA_MINICIF_VALUE_NONE;
            LLKA_cifDataItem_addValue(item, &cifValue);
            cifValue.state = LLKA_MINICIF_VALUE_SET;
        } else {
            size_t len = strlen(lastAtom->pdbx_PDB_ins_code);
            char *xbuf = (char*) malloc(len + 1);
            snprintf(xbuf, len + 1, "%s", lastAtom->pdbx_PDB_ins_code);
            cifValue.text = xbuf;
            LLKA_cifDataItem_addValue(item, &cifValue);
            free(xbuf);
        }

        /* Write step parameters row */

        item = NdbStructStepParameters->firstItem;

        cifValue.text = strbuf;

        snprintf(strbuf, 32, "%d", totalStepsCount);
        LLKA_cifDataItem_addValue(item, &cifValue);

        /** Actual step metrics */
        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.delta_1)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.epsilon_1)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.zeta_1)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.alpha_2)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.beta_2)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.gamma_2)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.delta_2)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.chi_1)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.chi_2)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.2f", cs->metrics.NN);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.2f", cs->metrics.CC);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->metrics.mu)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        /** Metrics differences from averages of the closest NtC class */
        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.delta_1));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.epsilon_1));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.zeta_1));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.alpha_2));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.beta_2));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.gamma_2));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.delta_2));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.chi_1));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.chi_2));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.2f", cs->differencesFromNtCAverages.NN);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.2f", cs->differencesFromNtCAverages.CC);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->differencesFromNtCAverages.mu));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        /** Confals */
        snprintf(strbuf, 32, "%.1f", cs->confalScore.delta_1);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.epsilon_1);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.zeta_1);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.alpha_2);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.beta_2);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.gamma_2);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.delta_2);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.chi_1);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.chi_2);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.NN);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.CC);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", cs->confalScore.mu);
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        /** Details */
        item = LLKA_cifDataCategory_nextItem(item);
        {
            const char *TORSION_NAMES[] = { "d1", "e1", "z1", "a2", "b2", "g2", "d2", "ch1", "ch2" };
            int haveDetails = 0;
            size_t idx;

            detailsBuf[0] = '\0';

            if (cs->violations & LLKA_CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT) {
                for (idx = 0; idx < 9; idx++) {
                    if (cs->violatingTorsionsAverage & (1 << idx)) {
                        if (haveDetails) strncat(detailsBuf, ";", SZ_DETAILS);
                        strncat(detailsBuf, "cAn", SZ_DETAILS);
                        strncat(detailsBuf, TORSION_NAMES[idx], SZ_DETAILS);
                        haveDetails = 1;
                    }
                }
            }

            if (cs->violations & LLKA_CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT) {
                for (idx = 0; idx < 9; idx++) {
                    if (cs->violatingTorsionsAverage & (1 << idx)) {
                        if (haveDetails) strncat(detailsBuf, ";", SZ_DETAILS);
                        strncat(detailsBuf, "cNn", SZ_DETAILS);
                        strncat(detailsBuf, TORSION_NAMES[idx], SZ_DETAILS);
                        haveDetails = 1;
                    }
                }
            }

            if (cs->violations & LLKA_CLASSIFICATION_E_NN_TOO_LOW || cs->violations & LLKA_CLASSIFICATION_E_NN_TOO_HIGH) {
                if (haveDetails) strncat(detailsBuf, ";", SZ_DETAILS);
                strncat(detailsBuf, "cNN", SZ_DETAILS);
                haveDetails = 1;
            }

            if (cs->violations & LLKA_CLASSIFICATION_E_CC_TOO_LOW || cs->violations & LLKA_CLASSIFICATION_E_CC_TOO_HIGH) {
                if (haveDetails) strncat(detailsBuf, ";", SZ_DETAILS);
                strncat(detailsBuf, "cCC", SZ_DETAILS);
                haveDetails = 1;
            }

            if (cs->violations & LLKA_CLASSIFICATION_E_MU_TOO_LOW || cs->violations & LLKA_CLASSIFICATION_E_MU_TOO_HIGH) {
                if (haveDetails) strncat(detailsBuf, ";", SZ_DETAILS);
                strncat(detailsBuf, "cmu", SZ_DETAILS);
                haveDetails = 1;
            }

            if (cs->violations & LLKA_CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH) {
                if (haveDetails) strncat(detailsBuf, ";", SZ_DETAILS);
                strncat(detailsBuf, "cMB", SZ_DETAILS);
                haveDetails = 1;
            }

            if (cs->violations & LLKA_CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT) {
                if (haveDetails) strncat(detailsBuf, ";", SZ_DETAILS);
                strncat(detailsBuf, "cP", SZ_DETAILS);
                haveDetails = 1;
            }

            if (cs->violations & LLKA_CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT) {
                if (haveDetails) strncat(detailsBuf, ";", SZ_DETAILS);
                strncat(detailsBuf, "cP1", SZ_DETAILS);
                haveDetails = 1;
            }

            if (cs->violations & LLKA_CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED) {
                if (haveDetails) strncat(detailsBuf, ";", SZ_DETAILS);
                strncat(detailsBuf, "DELTA", SZ_DETAILS);
                haveDetails = 1;
            }

            if (haveDetails) {
                cifValue.text = detailsBuf;
                LLKA_cifDataItem_addValue(item, &cifValue);
            } else {
                cifValue.state = LLKA_MINICIF_VALUE_NONE;
                LLKA_cifDataItem_addValue(item, &cifValue);
                cifValue.state = LLKA_MINICIF_VALUE_SET;
            }
        }

        /* Write sugar parameters row */

        item = NdbStructStepSugar->firstItem;

        cifValue.text = strbuf;

        snprintf(strbuf, 32, "%d", totalStepsCount);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->ribosePseudorotation_1)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->tau_1)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = LLKA_sugarPuckerToName(cs->sugarPucker_1, LLKA_SPN_VERY_TERSE);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->ribosePseudorotation_2)));
        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = strbuf;
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->tau_2)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = LLKA_sugarPuckerToName(cs->sugarPucker_2, LLKA_SPN_VERY_TERSE);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_1.nu_0)));
        item = LLKA_cifDataCategory_nextItem(item);
        cifValue.text = strbuf;
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_1.nu_1)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_1.nu_2)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_1.nu_3)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_1.nu_4)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_2.nu_0)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_2.nu_1)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_2.nu_2)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_2.nu_3)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs->nuAngles_2.nu_4)));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_1.nu_0));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_1.nu_1));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_1.nu_2));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_1.nu_3));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_1.nu_4));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_2.nu_0));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_2.nu_1));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_2.nu_2));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_2.nu_3));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);

        snprintf(strbuf, 32, "%.1f", LLKA_rad2deg(cs->nuAngleDifferences_2.nu_4));
        item = LLKA_cifDataCategory_nextItem(item);
        LLKA_cifDataItem_addValue(item, &cifValue);
    }

    item = NdbStructNtCOverall->firstItem;
    cifValue.text = entryId;
    LLKA_cifDataItem_addValue(item, &cifValue);

    snprintf(strbuf, 32, "%.1f", avgConfal->score);
    item = LLKA_cifDataCategory_nextItem(item);
    cifValue.text = strbuf;
    LLKA_cifDataItem_addValue(item, &cifValue);

    snprintf(strbuf, 32, "%d", (int)avgConfal->percentile);
    item = LLKA_cifDataCategory_nextItem(item);
    LLKA_cifDataItem_addValue(item, &cifValue);

    strncpy(strbuf, LLKA_INTERNAL_NTC_VERSION, 32);
    item = LLKA_cifDataCategory_nextItem(item);
    LLKA_cifDataItem_addValue(item, &cifValue);

    strncpy(strbuf, LLKA_INTERNAL_CANA_VERSION, 32);
    item = LLKA_cifDataCategory_nextItem(item);
    LLKA_cifDataItem_addValue(item, &cifValue);

    snprintf(strbuf, 32, "%d", totalStepsCount);
    item = LLKA_cifDataCategory_nextItem(item);
    LLKA_cifDataItem_addValue(item, &cifValue);

    snprintf(strbuf, 32, "%d", assignedStepsCount);
    item = LLKA_cifDataCategory_nextItem(item);
    LLKA_cifDataItem_addValue(item, &cifValue);

    snprintf(strbuf, 32, "%d", totalStepsCount - assignedStepsCount);
    item = LLKA_cifDataCategory_nextItem(item);
    LLKA_cifDataItem_addValue(item, &cifValue);

    snprintf(strbuf, 32, "%d", unassignedButCloseEnoughStepsCount);
    item = LLKA_cifDataCategory_nextItem(item);
    LLKA_cifDataItem_addValue(item, &cifValue);

    /* Addition of entries or values to entries set a "tainted" flag on Cif data.
     * Before the data can be written to a string, it must be "detainted".
     * Only valid Cif data will pass a validation check that will clear the flag.
     * All functions that require valid Cif data check the "tainted" flag and
     * return an error if they get tainted data.
     */
    tRet = LLKA_cifData_detaint(cifData);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Extended Cif data is invalid: %s\n", LLKA_errorToString(tRet));
        abort();
    }

    free(strbuf);
    free(detailsBuf);
}

static
LLKA_ClassificationContext * initializeClassificationContext(void)
{
    LLKA_RetCode tRet;
    LLKA_ClassificationContext *ctx;
    LLKA_ClassificationLimits limits;
    LLKA_Resource goldenSteps;
    LLKA_Resource clusters;
    LLKA_Resource confals;
    LLKA_Resource nuAngles;
    LLKA_Resource confalPercentiles;

    goldenSteps.type = LLKA_RES_GOLDEN_STEPS;
    clusters.type = LLKA_RES_CLUSTERS;
    confals.type = LLKA_RES_CONFALS;
    nuAngles.type = LLKA_RES_AVERAGE_NU_ANGLES;
    confalPercentiles.type = LLKA_RES_CONFAL_PERCENTILES;

    limits.averageNeighborsTorsionCutoff = LLKA_deg2rad(28.0);
    limits.nearestNeighborTorsionsCutoff = LLKA_deg2rad(28.0);
    limits.totalDistanceCutoff = LLKA_deg2rad(60.0);
    limits.pseudorotationCutoff = LLKA_deg2rad(72.0);
    limits.minimumClusterVotes = 0.001111;
    limits.minimumNearestNeighbors = 7;
    limits.numberOfUsedNearestNeighbors = 11;

    /* Load classification parameters from files.
     * There are currently four parameter files */

    tRet = LLKA_loadResourceFile(create_full_path(GOLDEN_STEPS_FILE), &goldenSteps);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load golden steps: %s\nThe file is expected in cwd or you can set DNATCO_ASSETS_PATH shell variable.\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    tRet = LLKA_loadResourceFile(create_full_path(CLUSTERS_FILE), &clusters);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load clusters: %s\nThe file is expected in cwd or you can set DNATCO_ASSETS_PATH shell variable.\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    tRet = LLKA_loadResourceFile(create_full_path(CONFALS_FILE), &confals);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load confals: %s\nThe file is expected in cwd or you can set DNATCO_ASSETS_PATH shell variable.\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    tRet = LLKA_loadResourceFile(create_full_path(NU_ANGLES_FILE), &nuAngles);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load average nu angles: %s\nThe file is expected in cwd or you can set DNATCO_ASSETS_PATH shell variable.\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    tRet = LLKA_loadResourceFile(create_full_path(CONFAL_PERCENTILES_FILE), &confalPercentiles);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load confal percentiles: %s\nThe file is expected in cwd or you can set DNATCO_ASSETS_PATH shell variable.\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    /* Try to initialize the context */
    tRet = LLKA_initializeClassificationContext(
        clusters.data.clusters, clusters.count,
        goldenSteps.data.goldenSteps, goldenSteps.count,
        confals.data.confals, confals.count,
        nuAngles.data.clusterNuAngles, nuAngles.count,
        confalPercentiles.data.confalPercentiles, confalPercentiles.count,
        &limits,
        MAX_CLOSE_ENOUGH_RMSD,
        &ctx
    );
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to initialize classification context: %s\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    /* Classification context keeps internal copies of all of the parameters.
     * Delete the resources we used for initialization. */
    LLKA_destroyResource(&goldenSteps);
    LLKA_destroyResource(&clusters);
    LLKA_destroyResource(&confals);
    LLKA_destroyResource(&nuAngles);
    LLKA_destroyResource(&confalPercentiles);

    return ctx;
}

static
void classifyStructure(const LLKA_Structures *steps, const LLKA_ClassificationContext *ctx, LLKA_ClassifiedSteps *classifiedSteps, LLKA_AverageConfal *avgConfal)
{
    LLKA_RetCode tRet;

    /* Classify all steps */
    tRet = LLKA_classifyStepsMultiple(steps, ctx, classifiedSteps);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to classify steps: %s\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    *avgConfal = LLKA_averageConfalAttempted(classifiedSteps, ctx);
}

static
int makeStructureFromCif(const LLKA_PathChar *path, LLKA_ImportedStructure *importedStru)
{
    LLKA_RetCode tRet;
    char *error;
    tRet = LLKA_cifFileToStructure(path, importedStru, &error, LLKA_MINICIF_GET_CIFDATA);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Cannot read CIF file: %s\n", LLKA_errorToString(tRet));

        /* If the CIF file looks malformed, display more information about why the parsing failed */
        if (error) {
            fprintf(stderr, "CIF parsing error: %s\n", error);
            LLKA_destroyString(error);
        }

        return 0;
    }

    return 1;
}

bool isDigitFunction(std::string &input){
    if(input.empty()){
        return false;
    }
    for (size_t i = 0; i < input.size(); i++){
        if(!std::isdigit(input[i])){
            return false;
        }
    }
    return true;
}

double roundTo(double number, double precision){
    return std::round(number / precision) * precision;
}

bool fileExists(std::string name){
    std::ifstream jsonFile(name);
    if (jsonFile.is_open()){
        jsonFile.close();
        return true;
    } else {
        return false;
    }
}

nlohmann::json createJson(){
    return nlohmann::json::object();
}

LLKA_RetCode saveJson(nlohmann::ordered_json jsonData, std::string fileName){
//    std::ofstream outputFile(fileName+".json");
    std::ofstream outputFile(fileName);
    if(outputFile.is_open()){
        outputFile << jsonData.dump(0);
        outputFile.close();
        std::cout << "Data was succesfully written to: "+fileName << std::endl;
    }else{
        std::cerr << "No file named: "<<fileName <<", found! " << std::endl;
        return LLKA_E_NO_FILE;
    }
    return LLKA_OK;
}

LLKA_RetCode writeConnectivityToJson(LLKA_Connectivity *conn, size_t step, LLKA_NtC currNtC, LLKA_NtC secondNtC, bool forward, nlohmann::json &jsonData){
    size_t currId = step+1;
    if(!forward){
        jsonData["Connectivity"][std::to_string(currId)]["PrevConnectivity"][LLKA_NtCToName(secondNtC)][LLKA_NtCToName(currNtC)] = {{"C5pDis", roundTo(conn->C5PrimeDistance, 0.001)}, {"O3PDis", roundTo(conn->O3PrimeDistance, 0.001)}};
    }else{
        jsonData["Connectivity"][std::to_string(currId)]["NextConnectivity"][LLKA_NtCToName(currNtC)][LLKA_NtCToName(secondNtC)] = {{"C5pDis", roundTo(conn->C5PrimeDistance, 0.001)}, {"O3PDis", roundTo(conn->O3PrimeDistance, 0.001)}};
    }
    return LLKA_OK;
}

LLKA_RetCode writeSimilaritiesToJson(LLKA_Similarities &sims, nlohmann::json &jsonData, int stepName, double cutRMSD){
    bool isSetCutRMSD = cutRMSD != -1;
    size_t numElements = sizeof(NTCS) / sizeof(NTCS[0]);
    int amountOfWrittenSimilarity = 0;
    

    //Set first similarity to the smallest
    LLKA_Similarity smallest = sims.similars[0];
    int ntcsPosition = 0;

    for(size_t i = 0; i < numElements; i++){
        LLKA_Similarity sim = sims.similars[i];
        if(LLKA_INVALID_NTC == NTCS[i]){
            continue;
        }
        if(isSetCutRMSD){
            if(sim.rmsd > cutRMSD){
                if(smallest.rmsd > sim.rmsd){
                    smallest = sim;
                    ntcsPosition = i;
                }
                continue;
            }
        }
        amountOfWrittenSimilarity++;
        jsonData["Similarity"][std::to_string(stepName)][LLKA_NtCToName(NTCS[i])] = {{"rmsd", roundTo(sim.rmsd, 0.001)}, {"ED", roundTo(sim.euclideanDistance, 0.1)}};
    }
    if(amountOfWrittenSimilarity == 0){
        jsonData["Similarity"][std::to_string(stepName)][LLKA_NtCToName(NTCS[ntcsPosition])] = {{"rmsd", roundTo(smallest.rmsd, 0.001)}, {"ED", roundTo(smallest.euclideanDistance, 0.1)}};
    }

    return LLKA_OK;
}

LLKA_RetCode writeSimilaritiesToJson(LLKA_Similarity &sim, nlohmann::json &jsonData, LLKA_NtC ntcs ,int stepName){
    jsonData["Similarity"][std::to_string(stepName)][LLKA_NtCToName(ntcs)] = {{"rmsd", roundTo(sim.rmsd, 0.001)}, {"ED", roundTo(sim.euclideanDistance, 0.1)}};
    return LLKA_OK;
}

/*
    LLKA_Structure &firstStep       First step
    LLKA_NtC firstNtC               NtC for first step
    LLKA_Structure &secondStep      Second step
    LLKA_NtC secondNtC              Ntc for second step
    int currentStepId               Id of current step
    bool isNextConnectivity         Is a previsou or next connectivity?
    nlohmann::json &jsonData        Json data where result is saved
*/
LLKA_RetCode calculateConnectivitiesStructureSingle(
    LLKA_Structure &firstStep,
    LLKA_NtC firstNtC,
    LLKA_Structure &secondStep,
    LLKA_NtC secondNtC,
    LLKA_Connectivity &connsData
){
    LLKA_RetCode tRet;        

    tRet = LLKA_measureStepConnectivityNtCs(&firstStep, firstNtC, &secondStep, secondNtC, &connsData);
    if (tRet != LLKA_OK){
        std::cerr << "Error: " << LLKA_errorToString(tRet) << std::endl;
        return tRet;
    }
    return LLKA_OK;
}

bool cutOffDistanceCheck(double cutOffDistance, LLKA_Connectivity conns){
    double average = (conns.C5PrimeDistance + conns.O3PrimeDistance) / 2;
    if(average > cutOffDistance){
        return false;
    }else{
        return true;
    }
}

// Main function to calculate connectivity
LLKA_RetCode calculateConnectivity(
    LLKA_SortedSteps &steps,
    LLKA_NtC ntcs,
    size_t stepId,
    double cutDistance,
    std::string fileName,
    std::string prevNtC,
    std::string nextNtC
){
    int def = -1;
    bool isSetNtCs = ntcs != LLKA_INVALID_NTC;
    bool isSetStepId = (stepId != def);
    bool isSetCutDistance = (cutDistance != def);
    bool isSetNextNtCs;
    bool isSetPrevNtCs;
    LLKA_NtC prevNtCCon;
    LLKA_NtC nextNtCCon;
    LLKA_RetCode tRet;

    if(prevNtC == "default"){
        isSetPrevNtCs = false;
        prevNtCCon = LLKA_INVALID_NTC;
    }else if(prevNtC == "true"){
        isSetPrevNtCs = true;
        prevNtCCon = LLKA_INVALID_NTC;
    }else{
        isSetPrevNtCs = true;
        prevNtCCon = LLKA_nameToNtC(prevNtC.c_str());
    }

    if(nextNtC == "default"){
        isSetNextNtCs = false;
        nextNtCCon = LLKA_INVALID_NTC;
    }else if(nextNtC == "true"){
        isSetNextNtCs = true;
        nextNtCCon = LLKA_INVALID_NTC;
    }else{
        isSetNextNtCs = true;
        nextNtCCon = LLKA_nameToNtC(nextNtC.c_str());
    }    

    nlohmann::json jsonData = createJson();

    // Main process of calculating connectivity 
    
    for(size_t i = 0; i <= steps.nSteps-1; i++){
        //Check if stepIdIsSet, if is true then skip all others calculation
        if(isSetStepId){
            if(i != stepId){
                continue;
            } 
        }
        // For each NtC
        for(int o = 0; o < 96; o++){
            if(isSetNtCs){
                if(o != ntcs){
                    continue;
                }
            }
            // For each next connectivity
            int smallestIndxData = 0;
            LLKA_Connectivity smallestData;
            int totalWritten = 0;
            // Next connectivity working
            for(int a = 0; a < 96; a++){
                if(isSetPrevNtCs && (!isSetNextNtCs)){
                    break;
                }
                if(nextNtCCon != LLKA_INVALID_NTC && nextNtCCon != NTCS[a]){
                    continue;
                }
                if(!(steps.steps[i].haveNextStep)){
                    break;
                }
                LLKA_Connectivity connData;

                tRet = calculateConnectivitiesStructureSingle(steps.steps[i].currStep, NTCS[o], steps.steps[i].nextStep, NTCS[a], connData);
                // If calculation failed, return error code
                if(tRet != LLKA_OK){
                  return tRet;
                }
                if(isSetCutDistance){
                    if(a == 0){
                        smallestData = connData;
                        smallestIndxData = a;
                    }else{
                        if(!(totalWritten > 0)){
                            double average = (connData.C5PrimeDistance + connData.O3PrimeDistance) / 2;
                            double smallAverage = (smallestData.C5PrimeDistance + smallestData.O3PrimeDistance)/2;
                            if(average < smallAverage){
                                smallestData = connData;
                                smallestIndxData = a;
                            }
                        }
                    }
                    if(cutOffDistanceCheck(cutDistance, connData)){
                        tRet = writeConnectivityToJson(&connData, i, NTCS[o], NTCS[a], true, jsonData);
                        if(tRet != LLKA_OK){
                          return tRet;
                        }
                        totalWritten++;
                    }
                }else{
                    tRet = writeConnectivityToJson(&connData, i, NTCS[o], NTCS[a], true, jsonData);
                    if(tRet != LLKA_OK){
                      return tRet;
                    }
                    totalWritten++;
                }
            }
            if(totalWritten == 0 && !(isSetPrevNtCs && (!isSetNextNtCs))){
                writeConnectivityToJson(&smallestData, i, NTCS[o], NTCS[smallestIndxData], true, jsonData);    
            }
            // For each prev connectivity
            smallestIndxData = 0;
            totalWritten = 0;
            for(int b = 0; b < 96; b++){
                if(isSetNextNtCs && (!isSetPrevNtCs)){
                    break;
                }
                if(prevNtCCon != LLKA_INVALID_NTC && prevNtCCon != NTCS[b]){
                    continue;
                }
                if(!(steps.steps[i].havePrevStep)){
                    break;
                }
                LLKA_Connectivity connData;
                tRet = calculateConnectivitiesStructureSingle(steps.steps[i].prevStep, NTCS[b], steps.steps[i].currStep, NTCS[o], connData);
                if(tRet != LLKA_OK){
                  return tRet;
                }
                if(isSetCutDistance){
                    if(b == 0){
                        smallestData = connData;
                        smallestIndxData = b;
                    }else{
                        if(!(totalWritten > 0)){
                            double average = (connData.C5PrimeDistance + connData.O3PrimeDistance) / 2;
                            double smallAverage = (smallestData.C5PrimeDistance + smallestData.O3PrimeDistance)/2;
                            if(average < smallAverage){
                                smallestData = connData;
                                smallestIndxData = b;
                            }
                        }
                    }
                    if(cutOffDistanceCheck(cutDistance, connData)){
                        tRet = writeConnectivityToJson(&connData, i, NTCS[b], NTCS[o], false, jsonData);
                        if(tRet != LLKA_OK){
                          return tRet;
                        }
                        totalWritten++;
                    }
                }else{
                    tRet = writeConnectivityToJson(&connData, i, NTCS[b], NTCS[o], false, jsonData);
                    if(tRet != LLKA_OK){
                      return tRet;
                    }
                    totalWritten++;
                }
            }
            if(totalWritten == 0){
                tRet = writeConnectivityToJson(&smallestData, i, NTCS[smallestIndxData], NTCS[o], false, jsonData);
                if(tRet != LLKA_OK){
                  return tRet;
                }
            }

        }
    }
    saveJson(jsonData, fileName);
    return LLKA_OK;

}

LLKA_RetCode calculateSimilarities(
    LLKA_Structure currentStruc,
    LLKA_NtC ntc,
    LLKA_Similarity &simData
){

    LLKA_RetCode tRet;
    tRet = LLKA_measureStepSimilarityNtC(&currentStruc, ntc, &simData);

    if(tRet != LLKA_OK){
        fprintf(stderr, "Failed to calculate similarity: %s\n", LLKA_errorToString(tRet));
        return tRet;
    }
    return LLKA_OK;
    }

// Main function to calculate Similarity
LLKA_RetCode calculateSimilarity(
    LLKA_Structures &steps,
    LLKA_NtC ntcs,
    unsigned int stepId,
    double cutRMSD,
    std::string fileName
){
    int def = -1;
    bool isSetNtCs = ntcs != LLKA_INVALID_NTC;
    bool isSetStepId = (stepId != def);
    bool isSetCutDistance = (cutRMSD != def);
    LLKA_RetCode tRet;

    nlohmann::json jsonData = createJson();

    //Calculating similarity
    for(size_t i = 0; i<steps.nStrus; i++){
        //Controling if stepId is selected
        if(isSetStepId){
            if(i != stepId){
                continue;
            }
        }
        LLKA_Similarity smallestData;
        int smallestIndxData = 0;
        int totalWritten = 0;
        for(int a = 0; a < 96; a++){
            if(isSetNtCs){
                if(ntcs != NTCS[a]){
                    continue;
                }
            }
            LLKA_Similarity simData;
            tRet = calculateSimilarities(steps.strus[i], NTCS[a],simData);
            if(tRet != LLKA_OK){
              return tRet;
            }
            if(isSetCutDistance){
                if(a == 0){
                    smallestData = simData;
                    smallestIndxData = a;
                }else{
                    if(!(totalWritten > 0)){
                        if(simData.rmsd < smallestData.rmsd){
                            smallestData = simData;
                            smallestIndxData = a;
                        }
                    }
                }
                if(simData.rmsd < cutRMSD){
                    tRet = writeSimilaritiesToJson(simData, jsonData, NTCS[a], i+1);
                    if(tRet != LLKA_OK){
                      return tRet;
                    }
                    totalWritten++;
                }
            }else{
                tRet = writeSimilaritiesToJson(simData, jsonData, NTCS[a], i+1);
                if(tRet != LLKA_OK){
                  return tRet;
                }
                totalWritten++;
            }

        }
        if(totalWritten == 0){
            tRet = writeSimilaritiesToJson(smallestData, jsonData, NTCS[smallestIndxData], i+1);
            if(tRet != LLKA_OK){
              return tRet;
            }
        }        
    }
    //Saving
    tRet = saveJson(jsonData, fileName);
    if(tRet != LLKA_OK){
      return tRet;
    }
    return LLKA_OK;
}

LLKA_RetCode testSimilarity(LLKA_Structure steps, LLKA_NtC ntc){
    LLKA_Similarity testSimilarity;
    calculateSimilarities(steps, ntc, testSimilarity);
    std::cout << "RMSD: " << testSimilarity.rmsd << " " << "ED: " << testSimilarity.euclideanDistance << std::endl;
    return LLKA_OK;
}

LLKA_RetCode testConnectivityNext(LLKA_Structure stepsFirst, LLKA_Structure stepsSecond, LLKA_NtC ntcFirst, LLKA_NtC ntcSecond){
    LLKA_Connectivity connectivity;
    calculateConnectivitiesStructureSingle(stepsFirst, ntcFirst, stepsSecond, ntcSecond, connectivity);
    std::cout << "C5Distance: " << connectivity.C5PrimeDistance << " " << "O3Distance: " << connectivity.O3PrimeDistance << std::endl;
    return LLKA_OK;
}

LLKA_RetCode testConnectivityPrev(LLKA_Structure stepsFirst, LLKA_Structure stepsSecond, LLKA_NtC ntcFirst, LLKA_NtC ntcSecond){
    LLKA_Connectivity connectivity;
    calculateConnectivitiesStructureSingle(stepsSecond, ntcSecond, stepsFirst, ntcFirst, connectivity);
    std::cout << "C5Distance: " << connectivity.C5PrimeDistance << " " << "O3Distance: " << connectivity.O3PrimeDistance << std::endl;
    return LLKA_OK;
}

LLKA_RetCode analyzeStep(LLKA_Structure &step, LLKA_AnalyzedStep &result){
    for (int j = 0; j < step.nAtoms-1; j++){
        if(j == 0){

            result.auth_asym_id_1 = step.atoms[j].auth_asym_id;
            result.label_comp_id_1 = step.atoms[j].label_comp_id;
            result.PDB_ins_code_1 = step.atoms[j].pdbx_PDB_model_num;
            result.auth_seq_id_1 = step.atoms[j].auth_seq_id;
            result.label_alt_id_1 = step.atoms[j].label_alt_id;
            continue;
        }
        if( strcmp(result.auth_asym_id_1, step.atoms[j].auth_asym_id) == 0 && strcmp(result.label_comp_id_1, step.atoms[j].label_comp_id) == 0 &&
        result.PDB_ins_code_1 == step.atoms[j].pdbx_PDB_model_num && result.auth_seq_id_1 == step.atoms[j].auth_seq_id &&
        result.label_alt_id_1 == step.atoms[j].label_alt_id)
        {
            continue;
        }
        else
        {
            result.auth_asym_id_2 = step.atoms[j].auth_asym_id;
            result.label_comp_id_2 = step.atoms[j].label_comp_id;
            result.PDB_ins_code_2 = step.atoms[j].pdbx_PDB_model_num;
            result.auth_seq_id_2 = step.atoms[j].auth_seq_id;
            result.label_alt_id_2 = step.atoms[j].label_alt_id;
            break;
        }
    }
    return LLKA_OK;
}

LLKA_RetCode sortSteps(LLKA_Structures &steps, LLKA_SortedSteps &sortedSteps) {
    LLKA_RetCode tRet;
    sortedSteps.nSteps = steps.nStrus;
    sortedSteps.steps = new LLKA_SortedStep[steps.nStrus];

    LLKA_AnalyzedSteps allSteps;
    allSteps.nSteps = steps.nStrus;
    allSteps.steps = new LLKA_AnalyzedStep[steps.nStrus];

    for (size_t s = 0; s < steps.nStrus; s++) {
        LLKA_AnalyzedStep analyzedStep;
        tRet = analyzeStep(steps.strus[s], analyzedStep);
        if(tRet != LLKA_OK){
          return tRet;
        }
        allSteps.steps[s] = analyzedStep;
    }

    for (size_t i = 0; i < steps.nStrus; i++) {
        LLKA_SortedStep currSortStep{};

        currSortStep.currStep = steps.strus[i];
        currSortStep.havePrevStep = false;
        currSortStep.haveNextStep = false;

        LLKA_AnalyzedStep currStepInformation = allSteps.steps[i];

        for (size_t j = 0; j < steps.nStrus; j++) {
            if (i == j) continue;
            LLKA_AnalyzedStep prevStepInformation = allSteps.steps[j];
            if (strcmp(currStepInformation.auth_asym_id_1, prevStepInformation.auth_asym_id_2) == 0 &&
                strcmp(currStepInformation.label_comp_id_1, prevStepInformation.label_comp_id_2) == 0 &&
                currStepInformation.PDB_ins_code_1 == prevStepInformation.PDB_ins_code_2 &&
                currStepInformation.auth_seq_id_1 == prevStepInformation.auth_seq_id_2 &&
                currStepInformation.label_alt_id_1 == prevStepInformation.label_alt_id_2)
            {
                currSortStep.havePrevStep = true;
                currSortStep.prevStep = steps.strus[j];
                break;
            }
        }

        for (size_t j = 0; j < steps.nStrus; j++) {
            if (i == j) continue;

            LLKA_AnalyzedStep nextStepInformation = allSteps.steps[j];
            if (strcmp(currStepInformation.auth_asym_id_2, nextStepInformation.auth_asym_id_1) == 0 &&
                strcmp(currStepInformation.label_comp_id_2, nextStepInformation.label_comp_id_1) == 0 &&
                currStepInformation.PDB_ins_code_2 == nextStepInformation.PDB_ins_code_1 &&
                currStepInformation.auth_seq_id_2 == nextStepInformation.auth_seq_id_1 &&
                currStepInformation.label_alt_id_2 == nextStepInformation.label_alt_id_1)
            {
                currSortStep.haveNextStep = true;
                currSortStep.nextStep = steps.strus[j];
                break;
            }
        }
        sortedSteps.steps[i] = currSortStep;
    }
    delete[] allSteps.steps;
    return LLKA_OK;
}

std::string getStepName(LLKA_ImportedStructure importedStru, int idOfStep, LLKA_Structures &steps){
    return deriveDNATCOStepName(importedStru.entry.id, structureHasMultipleModels(&steps),&steps.strus[idOfStep]);
}

void creatingStepIdFile(LLKA_ImportedStructure importedStru, LLKA_Structures &steps, std::string prefix){
    nlohmann::ordered_json information = createJson();
    for(size_t i = 0; i < steps.nStrus; i++){
        information["steps"][std::to_string(i+1)] = getStepName(importedStru, i,steps);
    }
    saveJson(information, prefix + "_step_information.json");
}

LLKA_RetCode createAndWriteCifFile(std::string cifString, std::string fileName){
    //std::string fileNameEdited = fileName + ".cif";
    std::ofstream outputFile(fileName);
    if(!outputFile){
        return LLKA_E_CANNOT_READ_FILE;
    }
    outputFile << cifString;
    outputFile.close();
    std::cout << "Data succesfully written to "+fileName << std::endl;
    return LLKA_OK;
}


#ifdef LLKA_PLATFORM_WIN32
int wmain(int argc, wchar_t *argv[])
#else
int main(int argc, char *argv[])
#endif /* LLKA_PLATFORM_WIN32 */
{
    LLKA_RetCode tRet;
    LLKA_Structures steps;
    LLKA_ClassificationContext *ctx;
    char *cifString;
    LLKA_ImportedStructure importedStru = {};
    LLKA_ClassifiedSteps classifiedSteps = {};
    LLKA_AverageConfal avgConfal = {};
    LLKA_SortedSteps sortedSteps = {};
    nlohmann::json jsonFile;

    // Section to setup flags or input options
    // TODO: Write more description
    CLI::App app{"Calculation of Similarity & Connectivity, Expand mmcif file"};

    argv = app.ensure_utf8(argv);

    // ==============================================
    // Initialization default parms
    // DO NOT CHANGE THIS
    // int
    size_t def = -1;
    size_t idOfStep = -1;
    // double
    double cutOffRmsd = -1.0;
    double cutOffDis = -1.0;
    // bool
    bool version = false;
    bool test = false;
    bool overwrite = false;
    // string
    std::string currNtC = "default";
    std::string prevNtCCon = "default";
    std::string nextNtCCon = "default";
    std::string similarity = "default";
    std::string connectivity = "default";
    std::string filename = "default";
    std::string outputCif = "default";
    std::string outputPrefix = "default";

    // Read inputs
    app.add_option("-x,--step", idOfStep, "Enter id of the step you want to analyze [1..N] (default: process all)");
    app.add_option("-r,--rmsd", cutOffRmsd, "Limit Similarity output to RMSD in Å (default: all)");
    app.add_option("-d,--distance_cutoff", cutOffDis, "Limit Connectivity output to (prev/next)steps within given distance");
    app.add_option("-n,--ntc", currNtC, "Return data only for given NtC type");
    app.add_option("-i,--mmcif-in", filename, "Input mmcif file name")
                                            ->check(CLI::ExistingFile)
                                            ->required();
    auto* ocif = app.add_flag("-o,--mmcif-out", outputCif, "Extendend mmCIF output file name")
                                            ->take_last()
                                            ->expected(0,1);
    app.add_option("-p,--prefix", outputPrefix, "Prefix for output files (default: entry.id content)");
    app.add_flag("-v,--version", version, "Show application version");
    app.add_flag("-t,--test", test, "Some nice description.. :)")->group("");
    app.add_flag("-f,--force_overwrite", overwrite, "Force overwrite of existing file");
    auto* sims = app.add_flag("-s,--similarity", similarity, "Calculate similarity, output filename is optional")
                                            ->take_last()
                                            ->expected(0,1);
    auto* cons = app.add_flag("-c,--connectivity", connectivity, "Calculate connectivity, output filename is optional")
                                            ->take_last()
                                            ->expected(0,1);
    app.add_flag("-b,--prev-ntc", prevNtCCon, "Limit Connectivity output for the previous step to given NtC type")
                                            ->take_last()
                                            ->expected(0,1);
    app.add_flag("-a,--next-ntc", nextNtCCon, "Limit Connectivity output for the next step to given NtC type")
                                            ->take_last()
                                            ->expected(0,1);                                            

    CLI11_PARSE(app, argc, argv);
    // End of input flags
    // ==============================================

    // Check if user choose to show version
    if(version){
        printf("Similarity & Connectivity version 0.5.0\n");
        printf("DNATCO v5.0, NtC v" LLKA_INTERNAL_NTC_VERSION ", CANA v" LLKA_INTERNAL_CANA_VERSION "\n");
        return EXIT_SUCCESS;
    }

    // Read structure from Cif File
    if (!makeStructureFromCif(filename.c_str(), &importedStru))
        return EXIT_FAILURE;

    tRet = LLKA_splitStructureToDinucleotideSteps(&importedStru.structure, &steps);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to split the structure into steps: %s\n", LLKA_errorToString(tRet));
        return EXIT_FAILURE;
    }

    ctx = initializeClassificationContext();

    classifyStructure(&steps, ctx, &classifiedSteps, &avgConfal);

    // Controling inputs
    // ==============================================
    LLKA_NtC currNtCIn = LLKA_INVALID_NTC;
    if(currNtC != "default"){
        currNtCIn = LLKA_nameToNtC(currNtC.c_str());
        if (currNtCIn == LLKA_INVALID_NTC){
            fprintf(stderr, "You entered wrong ntc Global!");
            return EXIT_FAILURE;
        }
    }

    LLKA_NtC nextNtCIn = LLKA_INVALID_NTC;
    if(nextNtCCon != "true" && nextNtCCon != "default"){
        nextNtCIn = LLKA_nameToNtC(nextNtCCon.c_str());
        if (nextNtCIn == LLKA_INVALID_NTC){
            fprintf(stderr, "You entered wrong ntc Next!");
            return EXIT_FAILURE;
        }
    }
    
    LLKA_NtC prevNtCIn = LLKA_INVALID_NTC;
    if(prevNtCCon != "true" && prevNtCCon != "default"){
        prevNtCIn = LLKA_nameToNtC(prevNtCCon.c_str());
        if (prevNtCIn == LLKA_INVALID_NTC){
            fprintf(stderr, "You entered wrong ntc Prev!");
            return EXIT_FAILURE;
        }
    }

    // Control StepId
    if(idOfStep != def){
        if(idOfStep > steps.nStrus){
            fprintf(stderr, "The step id you entered is larger than number of steps in the structure.");
            return EXIT_FAILURE;
        }else if(idOfStep<=0){
            fprintf(stderr, "The step id should be positive number.");
            return EXIT_FAILURE;
        }
        else{
            idOfStep = idOfStep-1;
        }
    }
    // End of Controling statement
    // ==============================================

    // Start of Test statement
    // ==============================================

    if(test){
        if (idOfStep == def || currNtC == "default"){
            return EXIT_FAILURE;
        }
        if(sims->count()>0){
            testSimilarity(steps.strus[idOfStep], currNtCIn);
            return EXIT_SUCCESS;
        }
        if(prevNtCCon == "default" || nextNtCCon == "default"){
            return EXIT_FAILURE;
        }
        if(cons->count()>0){
            testConnectivityPrev(steps.strus[idOfStep], steps.strus[idOfStep-1],currNtCIn, prevNtCIn);
            testConnectivityNext(steps.strus[idOfStep], steps.strus[idOfStep+1],currNtCIn, nextNtCIn);
            return EXIT_SUCCESS;
        }
        return EXIT_SUCCESS;
    }

    // End of Test statement
    // ==============================================

    std::string prefix;
    if (outputPrefix == "default"){
        prefix = std::string(importedStru.entry.id);
    } else {
        prefix = outputPrefix;
    }

    // Generating information file
    creatingStepIdFile(importedStru, steps, prefix);

    // Setting name
    if(similarity=="true"){
        similarity = prefix + "_similarity.json";
    }
    if(connectivity == "true"){
        connectivity = prefix + "_connectivity.json";
    }
    if(outputCif == "true"){
        outputCif = prefix + "_extended.cif";
    }

    // Exists output file?
    // ==============================================
    if(sims->count()>0){
        if(fileExists(similarity)){
            if(overwrite){
                std::ofstream file(similarity, std::ios::trunc);
                if(!file) {
                  std::cerr << "Error opening file " << similarity << std::endl;
                  file.close();
                  return EXIT_FAILURE;
                }
                std::cout << "File \"" << similarity << "\" was successfully cleared." << std::endl;
                file.close();
            }else{
                std::cerr << "File \"" << similarity << "\" already exists! Data not written, use -f to overwrite the file or choose different output file name." << std::endl;
                return EXIT_SUCCESS;
            }
        }
    }

    if(cons->count()>0){
        if(fileExists(connectivity)){
            if(overwrite){
                std::ofstream file(connectivity, std::ios::trunc);
                if(!file) {
                    std::cerr << "Error opening file " << connectivity << std::endl;
                    file.close();
                    return EXIT_FAILURE;
                }
                std::cout << "File \"" << connectivity << "\" was successfully cleared." << std::endl;
                file.close();
            }else{
                std::cerr << "File \"" << connectivity << "\" already exists! Data not written, use -f to overwrite the file or choose different output file name." << std::endl;
                return EXIT_SUCCESS;
            }
        }
    }

    if(ocif->count()>0){
        if(fileExists(outputCif)){
            if(overwrite){
                std::ofstream file(outputCif, std::ios::trunc);
                if(!file) {
                    std::cerr << "Error opening file " << outputCif << std::endl;
                    file.close();
                    return EXIT_FAILURE;
                }
                std::cout << "File \"" << outputCif << "\" was successfully cleared." << std::endl;
                file.close();
            }else{
                std::cerr << "File \"" << outputCif << "\" already exists! Data not written, use -f to overwrite the file or choose different output file name." << std::endl;
                return EXIT_SUCCESS;
            }
        }
    }
    // ==============================================

    sortSteps(steps, sortedSteps);

    // Main logic
    if(sims->count() > 0){
        tRet = calculateSimilarity(steps, currNtCIn, idOfStep, cutOffRmsd, similarity);
        if (tRet != LLKA_OK) {
            fprintf(stderr, "Cannot calculate similarity: %s\n", LLKA_errorToString(tRet));
            goto out;
        }
    }
    if(cons->count() > 0){
        tRet= calculateConnectivity(sortedSteps, currNtCIn, idOfStep, cutOffDis, connectivity, prevNtCCon, nextNtCCon);
        if (tRet != LLKA_OK) {
            fprintf(stderr, "Cannot calculate connectivity: %s\n", LLKA_errorToString(tRet));
            goto out;
        }
    }
    if(ocif->count() > 0){
        extendCifData(importedStru.cifData, &classifiedSteps, &avgConfal, importedStru.entry.id, &steps, structureHasMultipleModels(&steps));
        tRet = LLKA_cifDataToString(importedStru.cifData, LLKA_TRUE, &cifString);
        if (tRet != LLKA_OK) {
            fprintf(stderr, "Cannot write out CifData as Cif string: %s\n", LLKA_errorToString(tRet));
            goto out;
        } else {
            createAndWriteCifFile(std::string(cifString), outputCif);
            LLKA_destroyString(cifString);
        }
    }

out:
    LLKA_destroyClassifiedSteps(&classifiedSteps);
    LLKA_destroyClassificationContext(ctx);

    LLKA_destroyStructures(&steps);
    LLKA_destroyImportedStructure(&importedStru);

    return EXIT_SUCCESS;
}

