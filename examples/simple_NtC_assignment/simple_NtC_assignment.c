/* vim: set sw=4 ts=4 sts=4 expandtab : */

/****************************************************************

  This example demostrates how to use the DNATCO nucleic acid
  processing library to assign NtC classes to a nucleic acid
  structure.

  The example accepts a mmCif file with the structure as
  a command line parameter. The code shown here:

  1) Reads the mmCif file and creates an internal representation
     of the structure.
  2) Splits the structure into a list of dinucletide steps.
  3) Classifies each step, meaning that it attempts to assign
     NtC class to each step and calculates the accompanying
     information.
  4) Prints out the result.

  Please refer to the documentation in the header files for
  more information about how the individual parts of the
  DNATCO nucleic acid processing library are supposed to
  be used.

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

/* DNATCO step name function signatures */
char * deriveDNATCOStepName(const char *entryId, int hasMultipleModels, const LLKA_Structure *step);
int structureHasMultipleModels(const LLKA_Structures *steps);

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

    tRet = LLKA_loadResourceFile(GOLDEN_STEPS_FILE, &goldenSteps);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load golden steps: %s\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    tRet = LLKA_loadResourceFile(CLUSTERS_FILE, &clusters);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load clusters: %s\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    tRet = LLKA_loadResourceFile(CONFALS_FILE, &confals);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load confals: %s\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    tRet = LLKA_loadResourceFile(NU_ANGLES_FILE, &nuAngles);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load average nu angles: %s\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    tRet = LLKA_loadResourceFile(CONFAL_PERCENTILES_FILE, &confalPercentiles);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to load confal percentiles: %s\n", LLKA_errorToString(tRet));
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
void annotateStructure(const char *entryId, const LLKA_Structures *steps)
{
    LLKA_RetCode tRet;
    LLKA_ClassificationContext *ctx;
    LLKA_ClassifiedSteps classifiedSteps;
    size_t idx;
    int hasMultipleModels = structureHasMultipleModels(steps); /* Needed to derive DNATCO step name */

    /* First we need to initialize the classification context.
     * The classification process is highly parametrized and requires external data to operate */
    ctx = initializeClassificationContext();

    /* Classify all steps */
    tRet = LLKA_classifyStepsMultiple(steps, ctx, &classifiedSteps);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to classify steps: %s\n", LLKA_errorToString(tRet));
        exit(EXIT_FAILURE);
    }

    /* We can use indices of the steps in the "Structures" object to access the
     * classified steps. The library guarantees that the index of a step in the input
     * will correspond to the index of the classified step in the output. */
    for (idx = 0; idx < steps->nStrus; idx++) {
        /* Create shorthand aliases for brevity */
        const LLKA_AttemptedClassifiedStep *attemptedStep = &classifiedSteps.attemptedSteps[idx];
        const LLKA_ClassifiedStep *cStep = &attemptedStep->step;
        char *DNATCOName;

        /* The step was rejected by the classification routine.
         * The structure may not be a valid NtC step or it has unacceptable geometry */
        if (attemptedStep->status != LLKA_OK) {
            printf("WARN: Could not classify step %zu: %s\n", idx, LLKA_errorToString(attemptedStep->status));
            continue;
        }

        DNATCOName = deriveDNATCOStepName(entryId, hasMultipleModels, &steps->strus[idx]);

        /* Print basic step classification info */
        printf(
            "* Step %zu\n"
            "  DNATCO name:                       %s\n"
            "  Assigned NtC:  %-4s, Closest NtC:  %-4s\n"
            "  Assigned CANA: %-4s, Closest CANA: %-4s\n"
            "  Golden step DNATCO name:           %s\n"
            "  RMSD to closest NtC:               %g\n",
            idx + 1,
            DNATCOName,
            LLKA_NtCToName(cStep->assignedNtC), LLKA_NtCToName(cStep->closestNtC),
            LLKA_CANAToName(cStep->assignedCANA), LLKA_CANAToName(cStep->closestCANA),
            cStep->closestGoldenStep,
            cStep->rmsdToClosestNtC
        );
        free(DNATCOName);

        if (cStep->violations != LLKA_CLASSIFICATION_OK) {
            /* If a step could not have been fully classified, show why. */

            int bit;

            printf("  Classification violations:         ");
            for (bit = 0; bit < 32; bit++) {
                int viol = (1 << bit);
                if (cStep->violations & viol)
                    printf("%s ", LLKA_classificationViolationToName(viol));
            }
            putchar('\n');
        }
        /* Step has been fully classified, show confal scores. */
        printf(
            "  Confal scores:\n"
            "  %-6s %d\n"
            "  %-6s %d\n"
            "  %-6s %d\n"
            "  %-6s %d\n"
            "  %-6s %d\n"
            "  %-6s %d\n"
            "  %-6s %d\n"
            "  %-6s %d\n"
            "  %-6s %d\n"
            "  %-5s %d\n"
            "  %-5s %d\n"
            "  %-6s %d\n"
            "  %s %d\n",
            LLKA_dinucleotideTorsionName(LLKA_TOR_DELTA_1, LLKA_TRUE), (int)(cStep->confalScore.delta_1 + 0.5),
            LLKA_dinucleotideTorsionName(LLKA_TOR_EPSILON_1, LLKA_TRUE), (int)(cStep->confalScore.epsilon_1 + 0.5),
            LLKA_dinucleotideTorsionName(LLKA_TOR_ZETA_1, LLKA_TRUE), (int)(cStep->confalScore.zeta_1 + 0.5),
            LLKA_dinucleotideTorsionName(LLKA_TOR_ALPHA_2, LLKA_TRUE), (int)(cStep->confalScore.alpha_2 + 0.5),
            LLKA_dinucleotideTorsionName(LLKA_TOR_BETA_2, LLKA_TRUE), (int)(cStep->confalScore.beta_2 + 0.5),
            LLKA_dinucleotideTorsionName(LLKA_TOR_GAMMA_2, LLKA_TRUE), (int)(cStep->confalScore.gamma_2 + 0.5),
            LLKA_dinucleotideTorsionName(LLKA_TOR_DELTA_2, LLKA_TRUE), (int)(cStep->confalScore.delta_2 + 0.5),
            LLKA_dinucleotideTorsionName(LLKA_TOR_CHI_1, LLKA_TRUE), (int)(cStep->confalScore.chi_1 + 0.5),
            LLKA_dinucleotideTorsionName(LLKA_TOR_CHI_2, LLKA_TRUE), (int)(cStep->confalScore.chi_2 + 0.5),
            LLKA_crossResidueMetricName(LLKA_XR_DIST_CC, LLKA_TRUE), (int)(cStep->confalScore.CC + 0.5),
            LLKA_crossResidueMetricName(LLKA_XR_DIST_NN, LLKA_TRUE), (int)(cStep->confalScore.NN + 0.5),
            LLKA_crossResidueMetricName(LLKA_XR_TOR_MU, LLKA_TRUE), (int)(cStep->confalScore.mu + 0.5),
            "total", (int)(cStep->confalScore.total + 0.5)
        );
        putchar('\n');
    }

    /* Delete classification data */
    LLKA_destroyClassifiedSteps(&classifiedSteps);
    /* Delete classification context */
    LLKA_destroyClassificationContext(ctx);
}

static
int makeStructureFromCif(const LLKA_PathChar *path, LLKA_ImportedStructure *importedStru)
{
    LLKA_RetCode tRet;
    char *error;

    tRet = LLKA_cifFileToStructure(path, importedStru, &error, LLKA_FALSE);
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
    LLKA_ImportedStructure importedStru = {0};

    if (argc < 2) {
        fprintf(stderr, "No input CIF file\n");
        return EXIT_FAILURE;
    }

    /* Read CIF file and create LLKA_Structure object based on its contents */
    if (!makeStructureFromCif(argv[1], &importedStru))
        return EXIT_FAILURE;

    /* "stru" object contains the entire structure. In order to do the NtC assignment,
     * the structure has to be split up to dinucleotide steps. */
    tRet = LLKA_splitStructureToDinucleotideSteps(&importedStru.structure, &steps);
    if (tRet != LLKA_OK) {
        fprintf(stderr, "Failed to split the structure into steps: %s\n", LLKA_errorToString(tRet));
        return EXIT_FAILURE;
    }

    /* Pass the steps on for annotation */
    annotateStructure(importedStru.entry.id, &steps);

    /* Destroy the steps and the "total" structure */
    LLKA_destroyStructures(&steps);
    LLKA_destroyImportedStructure(&importedStru);

    return EXIT_SUCCESS;
}
