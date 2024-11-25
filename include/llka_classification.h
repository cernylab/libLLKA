/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_CLASSIFICATION_H
#define _LLKA_CLASSIFICATION_H

#include "llka_ntc.h"
#include "llka_nucleotide.h"

/*!
 * List of possible violations that prevented assigment of NtC/CANA classes to step
 */
enum {
    LLKA_CLASSIFICATION_OK = 0,
    LLKA_CLASSIFICATION_E_SCORE_TOO_LOW = (1 << 0),
    LLKA_CLASSIFICATION_E_NOT_ENOUGH_NEAREST_NEIGHBORS = (1 << 1),
    LLKA_CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT = (1 << 2),
    LLKA_CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT = (1 << 3),
    LLKA_CLASSIFICATION_E_CC_TOO_LOW = (1 << 4),
    LLKA_CLASSIFICATION_E_CC_TOO_HIGH = (1 << 5),
    LLKA_CLASSIFICATION_E_NN_TOO_LOW = (1 << 6),
    LLKA_CLASSIFICATION_E_NN_TOO_HIGH = (1 << 7),
    LLKA_CLASSIFICATION_E_MU_TOO_LOW = (1 << 8),
    LLKA_CLASSIFICATION_E_MU_TOO_HIGH = (1 << 9),
    LLKA_CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH = (1 << 10),
    LLKA_CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT = (1 << 11),
    LLKA_CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT = (1 << 12),
    LLKA_CLASSIFICATION_E_BEST_CLUSTER_DOES_NOT_HAVE_ENOUGH_VOTES = (1 << 13),
    LLKA_CLASSIFICATION_E_DELTA_TORSION_ANGLE_REJECTED = (1 << 14),   /*!< delta_1 or delta_2 backbone torsion exceeded tolerances. */
    LLKA_CLASSIFICATION_E_WRONG_METRICS = (1 << 15),   /*!< Calculated metrics contains invalid data. If this flag is set, it will be the only flag set. */
    LLKA_CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH = (1 << 16) /*!< Step is unassigned but its RMSD against the best matching reference step is sufficiently low.
                                                                        The "close enough" RMSD value is set in \p LLKA_ClassificationContext. */
};

/*!
 * List of torsions that exceeded tolerance
 */
enum {
    LLKA_CLASSIFICATION_E_TORSION_DELTA_1   = (1 << 0),
    LLKA_CLASSIFICATION_E_TORSION_EPSILON_1 = (1 << 1),
    LLKA_CLASSIFICATION_E_TORSION_ZETA_1    = (1 << 2),
    LLKA_CLASSIFICATION_E_TORSION_ALPHA_2   = (1 << 3),
    LLKA_CLASSIFICATION_E_TORSION_BETA_2    = (1 << 4),
    LLKA_CLASSIFICATION_E_TORSION_GAMMA_2   = (1 << 5),
    LLKA_CLASSIFICATION_E_TORSION_DELTA_2   = (1 << 6),
    LLKA_CLASSIFICATION_E_TORSION_CHI_1     = (1 << 7),
    LLKA_CLASSIFICATION_E_TORSION_CHI_2     = (1 << 8)
};

/*!
 * Average confal score of multiple steps and the statistical percentile of the average score
 */
typedef struct LLKA_AverageConfal {
    double score;         /*!< Average score */
    double percentile;    /*!< Statistical percentile of the average score */
} LLKA_AverageConfal;
LLKA_IS_POD(LLKA_AverageConfal)

/*!
 * Tolerances and options of the classification logic
 */
typedef struct LLKA_ClassificationLimits {
    size_t minimumNearestNeighbors;          /*!< Mimimum number of golden steps required to attempt to classify a step. Must be positive. */
    size_t numberOfUsedNearestNeighbors;     /*!< Number of the best matching golden steps used to determine NtC/CANA classes. Must be positive and greater than \p minimumNearestNeighbors */
    double averageNeighborsTorsionCutoff;    /*!< Maximum allowed difference between averaged nearest neighbors torsion angles and the classifiee. Must be positive. In radians.  */
    double nearestNeighborTorsionsCutoff;    /*!< Maximum allowed difference between the nearest neighbor golden step and the classifiee. Must be positive. In radians. */
    double totalDistanceCutoff;              /*!< Maximum allowed total distance between differences of the first seven backbone torsions of the assigned cluster and the classifiee. Must be positive. In radians. */
    double pseudorotationCutoff;             /*!< Maximum allowed difference between the ribose pseudorotaion value of the assigned cluster and the classifiee. In radians. */
    double minimumClusterVotes;              /*!< Minimum number of "votes" the best matching clusted must have to be considered as classified. */
} LLKA_ClassificationLimits;
LLKA_IS_POD(LLKA_ClassificationLimits)

/*!
 * Tolerances of classification metrics
 */
typedef struct LLKA_ClassificationMetric {
    double deviation;    /*!< Standard devitation. Used to derive tolerances. */
    double minValue;     /*!< Minimum value. If the value is an angle, it shall be in radians in [0; 2PI] range. Must be smaller than \p meanValue and \p maxValue. */
    double meanValue;    /*!< Mean value. If the value is an angle, it shall be in radians in [0; 2PI] range. Must be between \p minValue and \p maxValue. */
    double maxValue;     /*!< Maximum value. If the value is an angle, it shall be in radians in [0; 2PI] range. Must be greater than \p minValue and \p meanValue. */
} LLKA_ClassificationMetric;
LLKA_IS_POD(LLKA_ClassificationMetric)

/*!
 * Torsion angles of the (deoxy)ribose bonds as classification metrics.
 */
typedef struct LLKA_NuAnglesMetrics {
    LLKA_ClassificationMetric nu_0;    /*!< C4'-O4'-C1'-C2' */
    LLKA_ClassificationMetric nu_1;    /*!< O4'-C1'-C2'-C3' */
    LLKA_ClassificationMetric nu_2;    /*!< C1'-C2'-C3'-C4' */
    LLKA_ClassificationMetric nu_3;    /*!< C2'-C3'-C4'-O4' */
    LLKA_ClassificationMetric nu_4;    /*!< C3'-C4'-O4'-C1' */
} LLKA_NuAnglesMetrics;
LLKA_IS_POD(LLKA_NuAnglesMetrics)

/*!
 * Used to pass nu angles for the whole cluster
 */
typedef struct LLKA_ClusterNuAngles {
    LLKA_NuAnglesMetrics firstNucleotide;     /*!< Nu angles of first nucleotide in step */
    LLKA_NuAnglesMetrics secondNucleotide;    /*!< Nu angles of second nucleotide in step */
    int32_t clusterNumber;                    /*!< Number of the associated classification cluster */
} LLKA_ClusterNuAngles;
LLKA_IS_POD(LLKA_ClusterNuAngles)

/*!
 * Classification cluster that corresponds to one specific NtC and CANA class.
 */
typedef struct LLKA_ClassificationCluster {
    LLKA_ClassificationMetric delta_1;                   /*!< delta_1 backbone torsion angle classificaion metric. In radians. */
    LLKA_ClassificationMetric epsilon_1;                 /*!< epsilon_1 backbone torsion angle classificaion metric. In radians. */
    LLKA_ClassificationMetric zeta_1;                    /*!< zeta_1 backbone torsion angle classificaion metric. In radians. */
    LLKA_ClassificationMetric alpha_2;                   /*!< alpha_2 backbone torsion angle classificaion metric. In radians. */
    LLKA_ClassificationMetric beta_2;                    /*!< beta_2 backbone torsion angle classificaion metric. In radians. */
    LLKA_ClassificationMetric gamma_2;                   /*!< gamma_2 backbone torsion angle classificaion metric. In radians. */
    LLKA_ClassificationMetric delta_2;                   /*!< delta_2 backbone torsion angle classificaion metric. In radians. */
    LLKA_ClassificationMetric chi_1;                     /*!< chi_1 base torsion angle classificaion metric. In radians. */
    LLKA_ClassificationMetric chi_2;                     /*!< chi_2 base torsion angle classificaion metric. In radians. */
    LLKA_ClassificationMetric CC;                        /*!< Distance between the N[1|9] atoms of the first and second nucleotide classificaion metric. In Ångströms. */
    LLKA_ClassificationMetric NN;                        /*!< Distance between the C1' atoms of the first and second nucleotide classification metric. In Ångströms. */
    LLKA_ClassificationMetric mu;                        /*!< epsilon_1 backbone torsion angle classificaion metric. In radians. */
    LLKA_NuAnglesMetrics nusFirst;                       /*!< nusFirst nu angles on the first (deoxy)ribose classification metrics. */
    LLKA_NuAnglesMetrics nusSecond;                      /*!< nusSecond nu angles on the first (deoxy)ribose classification metrics. */
    double ribosePseudorotation_1;                       /*!< Ribose pseudorotation for the first nucleotide */
    double ribosePseudorotation_2;                       /*!< Ribose pseudorotation for the second nucleotide */
    int32_t number;                                      /*!< Cluster number, must be unique */
    LLKA_NtC NtC;                                        /*!< Associated NtC class */
    LLKA_CANA CANA;                                      /*!< Associated CANA class */
} LLKA_ClassificationCluster;
LLKA_IS_POD(LLKA_ClassificationCluster)

/*!
 * Y value of the confal cumulative probability function.
 * Wrapped in a struct so that in can be used in the resource loader module.
 */
typedef struct LLKA_ConfalPercentile {
    double value;
} LLKA_ConfalPercentile;
LLKA_IS_POD(LLKA_ConfalPercentile)

/*!
 * Confal scores assigned to a step.
 * Confal scores are assigned only to steps that have been fully classified.
 */
typedef struct LLKA_ConfalScore {
    double delta_1;      /*!< delta_1 backbone torsion angle score */
    double epsilon_1;    /*!< epsilon_1 backbone torsion angle score */
    double zeta_1;       /*!< zeta_1 backbone torsion angle score */
    double alpha_2;      /*!< alpha_2 backbone torsion angle score */
    double beta_2;       /*!< beta_2 backbone torsion angle score */
    double gamma_2;      /*!< gamma_2 backbone torsion angle score */
    double delta_2;      /*!< delta_2 backbone torsion angle score */
    double chi_1;        /*!< chi_1 base torsion angle score */
    double chi_2;        /*!< chi_2 base torsion angle score */
    double CC;           /*!< CC dsitance angle score */
    double NN;           /*!< NN distance angle score */
    double mu;           /*!< mu backbone torsion angle score */
    double total;        /*!< total score */
} LLKA_ConfalScore;
LLKA_IS_POD(LLKA_ConfalScore)

/*!
 * Results of an attempt to classify a step.
 */
typedef struct LLKA_ClassifiedStep {
    LLKA_NtC assignedNtC;                           /*!< NtC class of the assigned cluster. Set only if there are no classification violations, otherwise \p LLKA_INVALID_NTC */
    LLKA_CANA assignedCANA;                         /*!< CANA class of the assigned cluster. Set only if there are no classification violations, otherwise \p LLKA_INVALID_CANA */
    LLKA_NtC closestNtC;                            /*!< NtC class of the closest matching golden step. Set even if there are classification violations. */
    LLKA_CANA closestCANA;                          /*!< CANA class of the closest matching golden step. Set even if there are classification violations. */
    LLKA_ConfalScore confalScore;                   /*!< Confal score. Calculated only if there are no violations. */
    double euclideanDistanceNtCIdeal;               /*!< Euclidean distance between the classifiee and mean values of the assigned cluster. */
    LLKA_StepMetrics metrics;                       /*!< Actual metrics of the step */
    LLKA_StepMetrics differencesFromNtCAverages;    /*!< Differences between the classifiee and mean values of the assigned cluster. */
    LLKA_NuAngles nuAngles_1;                       /*!< Nu angles of the first nucleotide of the classifiee. */
    LLKA_NuAngles nuAngles_2;                       /*!< Nu angles of the second nucleotide of the classifiee. */
    double ribosePseudorotation_1;                  /*!< Ribose pseudorotation of the first nucleotide of the classifiee. */
    double ribosePseudorotation_2;                  /*!< Ribose pseudorotation of the second nucleotide of the classifiee. */
    double tau_1;
    double tau_2;
    LLKA_SugarPucker sugarPucker_1;                 /*!< Type of sugar pucker of the first nucleotide */
    LLKA_SugarPucker sugarPucker_2;                 /*!< Type of sugar pucker of the second nucleotide */
    LLKA_NuAngles nuAngleDifferences_1;             /*!< First nucleotide nu angles differences between the assigned cluster mean values and the classifiee. In radians. */
    LLKA_NuAngles nuAngleDifferences_2;             /*!< Second nucleotide nu angles differences between the assigned cluster mean values and the classifiee. In radians. */
    double rmsdToClosestNtC;                        /*!< Root Mean Square Distance between the assigned NtC representative and the classifiee.
                                                         \p closestNtC representative of there are classification violations */
    const char * closestGoldenStep;                 /*!< Name of the golden step with the shortest overall euclidean distance to the classifiee.
                                                         Note that this is a pointer to the \p ClassificationContext object and will become invalid if that object is destroyed. */
    int32_t violations;                             /*!< Bitfield of violations that prevented the step from being fully classified. Set to \p LLKA_CLASSIFICATION_OK if there were no violations. */
    int16_t violatingTorsionsAverage;               /*!< Bitfield indicating which torsions exceeded tolerance measured against average values of all nearest neighbors.
                                                         This field is valid only of the \p LLKA_CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT is set in \p violations. */
    int16_t violatingTorsionsNearest;               /*!< Bitfield indicating which torsions exceeded tolerance measured against the nearest neighbor.
                                                         This field is valid only of the \p LLKA_CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT is set in \p violations. */
} LLKA_ClassifiedStep;
LLKA_IS_POD(LLKA_ClassifiedStep)

/*!
 * Result of an attempt to classify a step combined with an error code.
 * Used by \p LLKA_classifyStepsMultiple()
 */
typedef struct LLKA_AttemptedClassifiedStep {
    LLKA_ClassifiedStep step;    /*!< Result of the attempt to classify a step */
    LLKA_RetCode status;         /*!< Return code of the attempt. Contains the corresponding error code if no attempt to classify a step could have been made */
} LLKA_AttemptedClassifiedStep;
LLKA_IS_POD(LLKA_AttemptedClassifiedStep)

/*!
 * Mutiple attempted classified steps
 */
typedef struct LLKA_ClassifiedSteps {
    LLKA_AttemptedClassifiedStep *attemptedSteps;    /*!< Array of attemped classified steps */
    size_t nAttemptedSteps;                         /*!< Number of items in \p attemptedSteps array */
} LLKA_ClassifiedSteps;
LLKA_IS_POD(LLKA_ClassifiedSteps)

/*!
 * Confal definition.
 */
typedef struct LLKA_Confal {
    double delta_1;           /*!< delta_1 backbone torsion angle */
    double epsilon_1;         /*!< epsilon_1 backbone torsion angle */
    double zeta_1;            /*!< zeta_1 backbone torsion angle */
    double alpha_2;           /*!< alpha_2 backbone torsion angle */
    double beta_2;            /*!< beta_2 backbone torsion angle */
    double gamma_2;           /*!< gamma_2 backbone torsion angle */
    double delta_2;           /*!< delta_2 backbone torsion angle */
    double chi_1;             /*!< chi_1 base torsion angle */
    double chi_2;             /*!< chi_2 base torsion angle */
    double CC;                /*!< CC distance */
    double NN;                /*!< NN distance torsion angle */
    double mu;                /*!< mu cross-residue "torsion" angle */
    LLKA_NuAngles nusFirst;   /*!< nusFirst nu angles of the first (deoxy)ribose */
    LLKA_NuAngles nusSecond;  /*!< nusSecond nu angles of the first (deoxy)ribose */
    int32_t clusterNumber;    /*!< Number of the associated cluster */
} LLKA_Confal;
LLKA_IS_POD(LLKA_Confal)

/*!
 * Golden step definition
 */
typedef struct LLKA_GoldenStep {
    /* Sugar pucker geometries */
    LLKA_SugarPucker pucker_1;    /*!< Sugar pucker type of the first nucleotide */
    LLKA_SugarPucker pucker_2;    /*!< Sugar pucker type of the second nucleotide */
    LLKA_NuAngles nuAngles_1;     /*!< Nu angles of the first nucleotide */
    LLKA_NuAngles nuAngles_2;     /*!< Nu angles of the second nucleotide */
    /* NtC step metrics */
    LLKA_StepMetrics metrics;     /*!< NtC metrics of the golden step */
    const char *name;             /*!< DNATCO name of the step. */
    size_t clusterIdx;            /*!< Used internally by \p ClassificationContext */
    int32_t clusterNumber;        /*!< Name of the cluster associated with the golden step */
} LLKA_GoldenStep;
LLKA_IS_POD(LLKA_GoldenStep)

/* Opaque type - not to be accessed from the outside */
typedef struct LLKA_ClassificationContext LLKA_ClassificationContext;

LLKA_BEGIN_API_FUNCTIONS

/*!
 * Calculates average confal score for the given steps.
 *
 * @param[in] classifiedSteps Classified steps to calculate the average for
 * @param[in] nClassifiedSteps Number of steps in the array
 * @param[in] ctx Classification context
 *
 * @returns LLKA_AverageConfal object with the results.
 */
LLKA_API LLKA_AverageConfal LLKA_CC LLKA_averageConfal(const LLKA_ClassifiedStep *classifiedSteps, size_t nClassifiedSteps, const LLKA_ClassificationContext *ctx);

/*!
 * Calculates average confal score for the given steps.
 *
 * @param[in] classifiedSteps Classified steps to calculate the average for
 * @param[in] nClassifiedSteps Number of steps in the array
 * @param[in] ctx Classification context
 *
 * @returns LLKA_AverageConfal object with the results.
 */
LLKA_API LLKA_AverageConfal LLKA_CC LLKA_averageConfalAttempted(const LLKA_ClassifiedSteps *attemptedClassifiedSteps, const LLKA_ClassificationContext *ctx);

/*!
 * Retrieves classification cluster for the given NtC.
 *
 * @param[in] ntc NtC for which to get the corresponding cluster
 * @param[in] ctx Classification context object
 * @param[out] cluster Retrieved cluster
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_INVALID_ARGUMENT Cluster for the given NtC was not found in the passed classification context
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_classificationClusterForNtC(LLKA_NtC ntc, const LLKA_ClassificationContext *ctx, LLKA_ClassificationCluster *cluster);

/*!
 * Translates classification violation flag to string.
 *
 * @param[in] violation Violation flat to translate. If the value contains multiple violation flags the result is undefined.
 *
 * @returns String representation of violation flag
 */
LLKA_API const char * LLKA_CC LLKA_classificationViolationToName(int32_t violation);

/*!
 * Attempts to classify a step.
 *
 * @param[in] stru Step to classify. Structure must be a valid step.
 * @param[in] ctx Classification context.
 * @param[out] classifiedStep classification Result of the classification.
 *
 * @returns LLKA_OK Success, corresponding error code if an attempt to classify the structure
 *          could not have been made, eg. if the passed structure was not a step.
 *          Note that even if the function returns LLKA_OK, it is necessary to check
 *          the \p violations field of \p LLKA_ClassifiedStep to check whether the step was fully classified.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_classifyStep(const LLKA_Structure *stru, const LLKA_ClassificationContext *ctx, LLKA_ClassifiedStep *classifiedStep);

/*!
 * Attempts to classify multiple steps,
 *
 * @param[in] strus Array of structures to attempt to classify.
 * @param[in] ctx Classification context.
 * @param[out] classifiedSteps Results of classification.
 *                             Order of the results in the \p attemptedSteps array will be the same
 *                             as the order of steps in the input \p strus array.
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_INVALID_ARGUMENT Array of structures to classify was empty.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_classifyStepsMultiple(const LLKA_Structures *strus, const LLKA_ClassificationContext *ctx, LLKA_ClassifiedSteps *classifiedSteps);

/*!
 * Retrieves confal for the given NtC.
 *
 * @param[in] ntc NtC for which to get the corresponding confal
 * @param[in] ctx Classification context object
 * @param[out] confal Retrieved cluster
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_INVALID_ARGUMENT Confal for the given NtC was not found in the passed classification context
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_confalForNtC(LLKA_NtC ntc, const LLKA_ClassificationContext *ctx, LLKA_Confal *confal);

/*!
 * Gets the corresponding percentile for the given confal score.
 *
 * @param[in] confalScore Confal score to get the percentile for.
 * @param[in] ctx Classification context.
 *
 * @returns Confal percentile or -1 if the confal score is not within [0; 100]
 */
LLKA_API double LLKA_CC LLKA_confalPercentile(double confalScore, const LLKA_ClassificationContext *ctx);

/*!
 * Destroys classification context.
 * Note that some data returned in \p LLKA_ClassifiedStep may point to data in \p LLKA_ClassificationContext. The context
 * can therefore be safely destroyed only of there are no classification results hanging around.
 *
 * @param[in] ctx Context to destroy.
 */
LLKA_API void LLKA_CC LLKA_destroyClassificationContext(LLKA_ClassificationContext *ctx);

/*!
 * Destroys results set by \p LLKA_classifyStepsMultiple()
 *
 * @param[in] classifiedSteps Results to destroy.
 */
LLKA_API void LLKA_CC LLKA_destroyClassifiedSteps(LLKA_ClassifiedSteps *classifiedSteps);

/*!
 * Initializes classification context.
 *
 * Classification context creates internal copies of all data that are used to initialize it.
 *
 * @param[in] clusters Array of of classification clusters.
 * @param[in] nClusters Number of golden steps.
 * @param[in] goldenSteps Array of golden steps.
 * @param[in] nClusterElements Number of golden steps. Cluster numbers of golden steps must match to the clusters.
 * @param[in] confals Array of confal parameters.
 * @param[in] nConfals Number of confals. Number of confals must be the same as the number of steps and the cluster numbers must match.
 * @param[in] clusterNuAngles Array of cluster nu angles. Number of cluster nu angles must be the same as the number of clusters and the cluster numbers must match.
 * @param[in] limits Classifier tolerances and options.
 * @param[in] maxCloseEnoughRmsd Maximum RMSD value of an unassigned step at which the step is considered as "unassigned but close"
 * @param[out] ctx Classification context to be initialized.
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_INVALID_ARGUMENT Some argument was invalid.
 * @retval LLKA_E_MISMATCHING_SIZES Sizes of clusters, confals or cluster nu angles did not match.
 * @retval LLKA_E_BAD_CLASSIFICATION_CLUSTERS Some classification clusters have invalid definitions or duplicit cluster numbers.
 * @retval LLKA_E_BAD_GOLDEN_STEPS Some golden steps have invalid definitions or unknown cluster numbers.
 * @retval LLKA_E_BAD_CONFALS Some confals have invalid definitions or unknown cluster numbers.
 * @retval LLKA_E_BAD_AVERAGE_NU_ANGLES Some cluster nu angles have invalid definitions or unknown cluster numbers.
 * @retval LLKA_E_BAD_CLASSIFICATION_LIMITS Classification limits have invalid values.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_initializeClassificationContext(
    const LLKA_ClassificationCluster *clusters, size_t nClusters,
    const LLKA_GoldenStep *goldenSteps, size_t nGoldenSteps,
    const LLKA_Confal *confals, size_t nConfals,
    const LLKA_ClusterNuAngles *clusterNuAngles, size_t nClusterNuAngles,
    const LLKA_ConfalPercentile *confalPercentiles, size_t nConfalPercentiles,
    const LLKA_ClassificationLimits *limits,
    double maxCloseEnoughRmsd,
    LLKA_ClassificationContext **ctx
);

LLKA_END_API_FUNCTIONS

#endif /* _LLKA_CLASSIFICATION_H */
