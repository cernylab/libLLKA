// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_CONNECTIVITY_SIMILARITY_H
#define _LLKA_CONNECTIVITY_SIMILARITY_H

#include "llka_ntc.h"

/*!
 * Connectivity measure of two consecutive NtC steps.
 *
 * Connectivity is measured as the distance between the C5' atom of the \a second residue of the \a first step and the C5' atom of the \a first residue of the \a second step
 * and the O3' atom of the \a second residue of the \a first step and the O3' atom of the \a first residue of the \a second step.
 */
typedef struct LLKA_Connectivity {
    double C5PrimeDistance;    /*!< Distance between C5' on the second residue of the first step and C5' on the first residue of the second step. */
    double O3PrimeDistance;    /*!< Distance between O3' on the second residue of the first step and O3' on the first residue of the second step. */
} LLKA_Connectivity;
LLKA_IS_POD(LLKA_Connectivity)

/*!
 * List of connectivity measures.
 */
typedef struct LLKA_Connectivities {
    LLKA_Connectivity *conns;    /*!< Array of \p LLKA_Connectivity s */
    size_t nConns;               /*!< Number of elements in the conns array */
} LLKA_Connectivities;
LLKA_IS_POD(LLKA_Connectivities)

/*!
 * Similarity measure of two NtC steps.
 *
 * Similarity is measured in two dimensions. First dimension is the RMSD between the measured and reference dinucleotide.
 * Second dimension is a (weighted) difference between the 12 parameters of the NtC step metrics. Reference metrics
 * is a series of precalculated averages for each NtC.
 */
typedef struct LLKA_Similarity {
    double rmsd;                 /*!< RMSD of the steps after Kabsch superposition */
    double euclideanDistance;    /*!< Sum of the euclidean distances between the torsions or distances */
} LLKA_Similarity;
LLKA_IS_POD(LLKA_Similarity)

/*!
 * List of Similarities.
 */
typedef struct LLKA_Similarities {
    LLKA_Similarity *similars;    /*!< Array of \p LLKA_Similarity s */
    size_t nSimilars;             /*!< Number of elements in the simils array */
} LLKA_Similarities;
LLKA_IS_POD(LLKA_Similarities)

LLKA_BEGIN_API_FUNCTIONS

/*
 * Measures connectivity between two "synthetic" steps whose geometry corresponds exactly to concrete NtCs.
 *
 * This function takes two NtC steps from a real structure and two LLKA_NtC items. Steps from a real structure
 * is used to determine positions of the "synthetic" steps in space. "Synthetic" steps are first superposed onto the
 * real steps and then the connectivity is calculated.
 *
 * First step is the step with lower \a label_seq_id values.
 *
 * @param[in] positionFirst Structure used to determine position of the first "synthetic" step. Usually taken from a real structure. Structure must contain all atoms of the extended backbone.
 * @param[in] ntcFirst NtC identifier of the first "synthetic" step.
 * @param[in] positionSecond Structure used to determine position ofthe second "synthetic" step. Usuall taken from a real structure. Structure must contain all atoms of the extended backbone.
 * @param[in] ntcSecond NtC identifier of the second "synthetic" step.
 * @param[out] result Calculated connectivity
 *
 * @return LLKA_OK on success, LLKA_E_ error code if the connectivity cannot be measured or the input parameters are invalid.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityNtCs(const LLKA_Structure *positionFirst, LLKA_NtC ntcFirst, const LLKA_Structure *positionSecond, LLKA_NtC ntcSecond, LLKA_Connectivity *result);

/*
 * Measures connectivity between two "synthetic" steps whose geometry corresponds exactly to concrete NtCs.
 *
 * This function does the same thing as \p LLKA_measureStepConnectivityNtCs but allows to pass multiple NtCs for the first step.
 * This has the benefit of reduced overhead compared to calling \p LLKA_measureStepConnectivityNtCs for all NtCs individually.
 *
 * @param[in] positionFirst Structure used to determine position of the first "synthetic" step. Usually taken from a real structure. Structure must contain all atoms of the extended backbone.
 * @param[in] ntcsFirst Array of NtC identifiers of the first "synthetic" steps. The last element of the array must be LLKA_INVALID_NTC
 * @param[in] positionSecond Structure used to determine position ofthe second "synthetic" step. Usuall taken from a real structure. Structure must contain all atoms of the extended backbone.
 * @param[in] ntcSecond NtC identifier of the second "synthetic" step.
 * @param[out] results List of calculated connectivity. The order is the same as the order of \p ntcsFirst
 *
 * @return LLKA_OK on success, LLKA_E_ error code if the connectivity cannot be measured or the input parameters are invalid.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityNtCsMultipleFirst(const LLKA_Structure *positionFirst, const LLKA_NtC ntcsFirst[], const LLKA_Structure *positionSecond, LLKA_NtC ntcSecond, LLKA_Connectivities *results);

/*
 * Measures connectivity between two "synthetic" steps whose geometry corresponds exactly to concrete NtCs.
 *
 * This function does the same thing as \p LLKA_measureStepConnectivityNtCs but allows to pass multiple NtCs for the second step.
 * This has the benefit of reduced overhead compared to calling \p LLKA_measureStepConnectivityNtCs for all NtCs individually.
 *
 * @param[in] positionFirst Structure used to determine position of the first "synthetic" step. Usually taken from a real structure. Structure must contain all atoms of the extended backbone.
 * @param[in] ntcFirst NtC identifier of the first "synthetic" step.
 * @param[in] positionSecond Structure used to determine position ofthe second "synthetic" step. Usually taken from a real structure. Structure must contain all atoms of the extended backbone.
 * @param[in] ntcsSecond Array of NtC identifiers of the second "synthetic" steps. The last element of the array must be LLKA_INVALID_NTC
 * @param[out] results List of calculated connectivity. The order is the same as the order of \p ntcsFirst
 *
 * @return LLKA_OK on success, LLKA_E_ error code if the connectivity cannot be measured or the input parameters are invalid.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityNtCsMultipleSecond(const LLKA_Structure *positionFirst, LLKA_NtC ntcFirst, const LLKA_Structure *positionSecond, const LLKA_NtC ntcSecond[], LLKA_Connectivities *results);

/*!
 * Measures connectivity between two dinucleotides as if they were two consecutive steps in a molecule.
 *
 * First step is the step with lower \a label_seq_id values.
 *
 * @param[in] positionFirst Structure that the position of the first step in the molecule will be calculated from. This structure must
 *                          contain atoms of the extended backbone, otherwise the function will fail. In most cases,
 *                          this parameter will be an actual step in a molecule.
 * @param[in] dinuFirst Structure to be used as the first step. This structure must contain atoms of the extended
 *                      backbone, otherwise the function will fail. In most cases, this parameter will be a reference
 *                      NtC dinucleotide.
 * @param[in] positionSecond Structure that the position of the second step in the molecule will be calculated from. This structure
 *                           must contain atoms of the extended backbone, otherwise the function will fail. In most
 *                           cases, this parameter will be an actual step in a molecule.
 * @param[in] dinuSecond Structure to be used as the second step. This structure must contain atoms of the extended
 *                       backbone, otherwise the function will fail. In most cases, this parameter will be a reference
 *                       NtC dinucleotide.
 * @param[out] result The computed connectivity.
 *
 * @return LLKA_OK on success, LLKA_E_ if the connectivity cannot be measured.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityStructures(const LLKA_Structure *positionFirst, const LLKA_Structure *dinuFirst, const LLKA_Structure *positionSecond, const LLKA_Structure *dinuSecond, LLKA_Connectivity *result);

/*!
 * Measures connectivity of one structure against multiple reference structures.
 *
 * @param[in] positionFirst Structure that the position of the first step in the molecule will be calculated from. This structure must
 *                          contain atoms of the extended backbone, otherwise the function will fail. In most cases,
 *                          this parameter will be an actual step in a molecule.
 * @param[in] dinuFirst Structure to be used as the first step. This structure must contain atoms of the extended
 *                      backbone, otherwise the function will fail. In most cases, this parameter will be a reference
 *                      NtC dinucleotide.
 * @param[in] positionSecond Structure that the position of the second step in the molecule will be calculated from. This structure
 *                           must contain atoms of the extended backbone, otherwise the function will fail. In most
 *                           cases, this parameter will be an actual step in a molecule.
 * @param[in] dinusSecond Multiple structures to be used as the second steps. These structures must contain atoms of the extended
 *                        backbone, otherwise the function will fail. In most cases, this parameter will be referennces to
 *                        NtC dinucleotides.
 * @param[out] results Computed connectivities. Order of the results is the same as the order of structures in \p dinusSecond.
 *
 * @return LLKA_OK on success, LLKA_E_ if the connectivity cannot be measured.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_measureStepConnectivityStructuresMultiple(const LLKA_Structure *positionFirst, const LLKA_Structure *dinuFirst, const LLKA_Structure *positionSecond, const LLKA_Structures *dinusSecond, LLKA_Connectivities *results);

/*!
 * Measures similarity between given a dinucleotide and a reference NtC. Input structure must be a valid NtC step.
 *
 * @param[in] stepStru Structure to be measured
 * @param[in] ntc ID of the NtC to use as reference
 * @param[out] result The computed similarity.
 *
 * @return LLKA_OK on success, LLKA_E_ error code if the similarity cannot be measured.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_measureStepSimilarityNtC(const LLKA_Structure *stepStru, LLKA_NtC ntc, LLKA_Similarity *result);

LLKA_API LLKA_RetCode LLKA_CC LLKA_measureStepSimilarityNtCMultiple(const LLKA_Structure *stepStru, const LLKA_NtC ntcs[], LLKA_Similarities *results);

/*!
 * Measures similarity between two dinucleotides. Both dinuclecotides must be valid NtC steps.
 *
 * For description how the similarity is measured, see LLKA_measureStepSimilarityNtC().
 * The difference between this function and LLKA_measureStepSimilarityNtC() is that the reference NtC metrics
 * is calculated directly from the reference structure instead of being taken from precalculated data.
 *
 * @param[in] stepStru Structure to be measured
 * @param[in] refStru Structure to be used as reference
 * @param[out] result Calculated similarity
 *
 * @return LLKA_OK on success, LLKA_E_ error code if the similarity cannot be measured.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_measureStepSimilarityStructure(const LLKA_Structure *stepStru, const LLKA_Structure *refStru, LLKA_Similarity *result);

LLKA_END_API_FUNCTIONS

#endif // _LLKA_CONNECTIVITY_SIMILARITY_H
