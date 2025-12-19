/* vim: set sw=4 ts=4 sts=4 expandtab : */

/****************************************************************

  This example demonstrates how to use the MiniCif module
  to extend a Cif file with NtC assignment data. The example
  is based on  the "simple_NtC_assignment" code. For details
  how to use this library to do the NtC assigment, see
  the "simple_NtC_assignment" example code first.

  Output of this example code creates a Cif string whose
  contents mostly match the output of DNATCO 4.1.

****************************************************************/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
char * deriveDNATCOStepName(const char *entryId, int hasMultipleModels, const LLKA_Structure *step);
int structureHasMultipleModels(const LLKA_Structures *steps);

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

    strbuf = malloc(32); /* 32 characters is enough for everything except details */
    /* MiniCif makes internal copies of data values so it is safe to create
     * one data value object and reuse it */
    cifValue.text = strbuf; /* Work around the fact that "text" is const char * but we need a writable buffer */
    cifValue.state = LLKA_MINICIF_VALUE_SET;

    detailsBuf = malloc(SZ_DETAILS+1);

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
            char *xbuf = malloc(len + 1);
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
            char *xbuf = malloc(len + 1);
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
    LLKA_ImportedStructure importedStru = {0};
    LLKA_ClassifiedSteps classifiedSteps = {0};
    LLKA_AverageConfal avgConfal = {0};

    /* Read structure from Cif and assign NtCs */

    if (argc > 1 && strcmp(argv[1], "--version") == 0) {
        printf("classify_and_write_cif version 0.0.2\n");
        printf("DNATCO v5.0, NtC v" LLKA_INTERNAL_NTC_VERSION ", CANA v" LLKA_INTERNAL_CANA_VERSION "\n");
        return 0;
    }

    if (argc < 2) {
        fprintf(stderr, "No input CIF file\n");
        return EXIT_FAILURE;
    }

    if (!makeStructureFromCif(argv[1], &importedStru))
        return EXIT_FAILURE;

    tRet = LLKA_splitStructureToDinucleotideSteps(&importedStru.structure, &steps);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to split the structure into steps: %s\n", LLKA_errorToString(tRet));
        return EXIT_FAILURE;
    }

    ctx = initializeClassificationContext();

    classifyStructure(&steps, ctx, &classifiedSteps, &avgConfal);
    /* NtC assignment is done. Extend Cif data with NtC assignment results */

    extendCifData(importedStru.cifData, &classifiedSteps, &avgConfal, importedStru.entry.id, &steps, structureHasMultipleModels(&steps));

    /* Write out Cif data as a string */
    tRet = LLKA_cifDataToString(importedStru.cifData, LLKA_TRUE, &cifString);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Cannot write out CifData as Cif string: %s\n", LLKA_errorToString(tRet));
        goto out;
    } else {
        fprintf(stdout, "%s\n", cifString);
        LLKA_destroyString(cifString);
    }

out:
    LLKA_destroyClassifiedSteps(&classifiedSteps);
    LLKA_destroyClassificationContext(ctx);

    LLKA_destroyStructures(&steps);
    LLKA_destroyImportedStructure(&importedStru);

    return EXIT_SUCCESS;
}
