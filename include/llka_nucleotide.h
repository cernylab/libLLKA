/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_NUCLEOTIDE_H
#define _LLKA_NUCLEOTIDE_H

#include "llka_structure.h"

/*!
 * Torsion angles of the (deoxy)ribose bonds as scalars.
 */
typedef struct LLKA_NuAngles {
    double nu_0;    /*!< C4'-O4'-C1'-C2' */
    double nu_1;    /*!< O4'-C1'-C2'-C3' */
    double nu_2;    /*!< C1'-C2'-C3'-C4' */
    double nu_3;    /*!< C2'-C3'-C4'-O4' */
    double nu_4;    /*!< C3'-C4'-O4'-C1' */
} LLKA_NuAngles;
LLKA_IS_POD(LLKA_NuAngles)

/*!
 * All recognized sugar pucker geometries.
 */
typedef enum LLKA_SugarPucker {
    LLKA_C3_ENDO = 0,
    LLKA_C4_EXO,
    LLKA_O4_ENDO,
    LLKA_C1_EXO,
    LLKA_C2_ENDO,
    LLKA_C3_EXO,
    LLKA_C4_ENDO,
    LLKA_O4_EXO,
    LLKA_C1_ENDO,
    LLKA_C2_EXO,
    LLKA_INVALID_SUGAR_PUCKER = 0x7FFFFFFF
} LLKA_SugarPucker;

typedef struct LLKA_RiboseMetrics {
    LLKA_NuAngles nus;          /*!< Torsion angles of the (deoxy)ribose bonds */
    double P;                   /*!< Pseudorotation */
    double tMax;                /*!< tauMax */
    LLKA_SugarPucker pucker;    /*!< Sugar pucker conformation */
} LLKA_RiboseMetrics;
LLKA_IS_POD(LLKA_RiboseMetrics)

/*!
 * Levels of brevity of sugar pucker name
 */
typedef enum LLKA_SugarPuckerNameBrevity {
    LLKA_SPN_VERY_TERSE,    /*!< C1exo, C1end */
    LLKA_SPN_TERSE,         /*!< C1exo, C1endo */
    LLKA_SPN_FANCY          /*!< "C1' exo", "C1' endo" */
} LLKA_SugarPuckerNameBrevity;

LLKA_BEGIN_API_FUNCTIONS

/*!
 * Extracts a single nucleotide from a structure as LLKA_Structure.
 *
 * Note that this function performs exhaustive search on the entire structure.
 * Calling this function in a loop may lead to suboptimal performance.
 *
 * @param[in] stru Structure to extract from.
 * @param[in] pdbx_PDB_model_num Number of the model to extract from.
 * @param[in] label_asym_id Chain to extract from.
 * @param[in] label_seq_id Residue number to extract from.
 *
 * @retval LLKA_OK Success.
 */
LLKA_API LLKA_Structure LLKA_CC LLKA_extractNucleotide(const LLKA_Structure *stru, int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id);

/*!
 * Extracts a single nucleotide from a structure as LLKA_StructureView.
 *
 * Note that this function performs exhaustive search on the entire structure.
 * Calling this function in a loop may lead to suboptimal performance.
 *
 * @param[in] stru Structure to extract from.
 * @param[in] pdbx_PDB_model_num Number of the model to extract from.
 * @param[in] label_asym_id Chain to extract from.
 * @param[in] label_seq_id Residue number to extract from.
 *
 * @retval LLKA_OK Success.
 */
LLKA_API LLKA_StructureView LLKA_CC LLKA_extractNucleotideView(const LLKA_Structure *stru, int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id);

/*!
 * Extracts ribose ring from a structure as LLKA_Structure.
 * The structure must contain exactly one ribose ring, otherwise the result is undefined.
 *
 * @param[in] stru Structure to extract from.
 * @param[out] riboseStru Extracted ribose structure.
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_MISSING_ATOMS Structure does not contain all atoms of ribose ring.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_extractRibose(const LLKA_Structure *stru, LLKA_Structure *riboseStru);

/*!
 * Extracts ribose ring from a structure as LLKA_StructureView.
 * The structure must contain exactly one ribose ring, otherwise the result is undefined.
 *
 * @param[in] stru Structure to extract from.
 * @param[out] riboseStru Extracted ribose structure view.
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_MISSING_ATOMS Structure does not contain all atoms of ribose ring.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_extractRiboseView(const LLKA_Structure *stru, LLKA_StructureView *riboseView);

/*!
 * Checks if a compound name is recognized as nucleotide.
 *
 * @param[in] compId Name of the compound
 *
 * @returns LLKA_TRUE if the compound name is recognized as nucleotide, LLKA_FALSE otherwise
 */
LLKA_API LLKA_Bool LLKA_CC LLKA_isNucleotideCompound(const char *compId);

/*!
 * Returns sugar pucker type from the corresponding name.
 *
 * Recognized sample inputs are "C1end", "C1endo", "C1'end", "C1'endo", "C1' endo"...
 *
 * @param[in] name Sugar pucker name
 * @returns Sugar pucker type of LLKA_INVALID_SUGAR_PUCKER if the name is not recognized
 */
LLKA_API LLKA_SugarPucker LLKA_CC LLKA_nameToSugarPucker(const char *name);

/*!
 * Calculates metrics for the ribose core of a nucleotide.
 * The input structure must contain exactly one ribose core, otherwise the result is undefined.
 *
 * @param[in] stru Structure with the ribose core
 * @param[out] metrics Calculated ribose metrics
 *
 * @retval LLKA_OK Success
 * @rerval LLKA_E_INVALID_ARGUMENT Input structure does not seem to contain a ribose core
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_riboseMetrics(const LLKA_Structure *stru, LLKA_RiboseMetrics *metrics);

/*!
 * Calculates metrics for the ribose core of a nucleotide view.
 * The input structure must contain exactly one ribose core, otherwise the result is undefined.
 *
 * @param[in] stru Structure view with the ribose core
 * @param[out] metrics Calculated ribose metrics
 *
 * @retval LLKA_OK Success
 * @rerval LLKA_E_INVALID_ARGUMENT Input structure does not seem to contain a ribose core
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_riboseMetricsView(const LLKA_Structure *view, LLKA_RiboseMetrics *metrics);

/*!
 * Determines the sugar pucker of the ribose ring contained in the input structure.
 * The input structure must contain exactly one ribose ring, otherwise the results are undefined.
 *
 * @param[in] stru Structure to process
 * @param[out] pucker Determined sugar pucker
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_INVALID_ARGUMENT Input structure does not contain a ribose ring
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_sugarPucker(const LLKA_Structure *stru, LLKA_SugarPucker *pucker);

/*!
 * Determines the sugar pucker of the ribose ring contained in the input structure view.
 * The input structure must contain exactly one ribose ring, otherwise the results are undefined.
 *
 * @param[in] stru Structure to process
 * @param[out] pucker Determined sugar pucker
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_INVALID_ARGUMENT Input structure does not contain a ribose ring
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_sugarPuckerView(const LLKA_StructureView *view, LLKA_SugarPucker *pucker);

/*!
 * Returns name of the sugar pucker type.
 *
 * @param[in] pucker Sugar pucker type
 * @param[in] brevity Level of brevity of the returned name
 *
 * @returns Sugar pucker name or an empty string if the type is LLKA_INVALID_SUGAR_PUCKER
 */
LLKA_API const char * LLKA_CC LLKA_sugarPuckerToName(LLKA_SugarPucker pucker, LLKA_SugarPuckerNameBrevity brevity);

LLKA_END_API_FUNCTIONS

#endif /* _LLKA_NUCLEOTIDE_H */
