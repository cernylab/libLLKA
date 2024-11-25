/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_RESOURCE_LOADERS_H
#define _LLKA_RESOURCE_LOADERS_H

#include "llka_classification.h"

/*!
 * List of kinds of resources that can be loaded by libLLKA
 */
typedef enum LLKA_ResourceType {
    LLKA_RES_AVERAGE_NU_ANGLES,    /*!< Array of average Nu angles, used by classification context */
    LLKA_RES_CLUSTERS,             /*!< Array of classification clusters, used by classification context */
    LLKA_RES_CONFALS,              /*!< Array of confals, used by classification context */
    LLKA_RES_GOLDEN_STEPS,         /*!< Array of golden steps, used by classification context */
    LLKA_RES_CONFAL_PERCENTILES    /*!< Array of Y values of a cumulative probability function of confal scores. Values are wrapped in a struct for technical reasons */
} LLKA_ResourceType;

typedef union LLKA_UResource {
    LLKA_ClusterNuAngles *clusterNuAngles;
    LLKA_ClassificationCluster *clusters;
    LLKA_Confal *confals;
    LLKA_GoldenStep *goldenSteps;
    LLKA_ConfalPercentile *confalPercentiles;
} LLKA_UResource;

typedef struct LLKA_Resource {
    LLKA_ResourceType type;    /*!< Type of the resource */
    LLKA_UResource data;       /*!< Union of all resource data that can be loaded by libLLKA */
    size_t count;              /*!< If the resource data is an array, count indicates the number of items in the array. Otherwise it is set to 1 */
} LLKA_Resource;
LLKA_IS_POD(LLKA_Resource)

LLKA_BEGIN_API_FUNCTIONS

/*!
 * Frees LLKA_Resource along with the resource data
 *
 * @param[in] resource Resource to destroy
 */
LLKA_API void LLKA_CC LLKA_destroyResource(LLKA_Resource *resource);

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
/*!
 * Loads resource from file of a text string.
 *
 * @param[in] path Path to the file to read data from.
 * @param[in, out] resource Loaded resource data. The \p type member must be set to the type of resource to load.
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_NO_FILE File not found.
 * @retval LLKA_E_CANNOT_READ_FILE File exists but it cannot be read
 * @retval LLKA_E_BAD_DATA Input data is malformed.
 * @retval LLKA_E_MISMATCHING_DATA Input data is malformed.
 * @retval LLKA_E_MISMATCHING_SIZES Input data is malformed.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_loadResourceFile(const LLKA_PathChar *path, LLKA_Resource *resource);
#endif /* LLKA_FILESYSTEM_ACCESS_DISABLED */

/*!
 * Loads resource from file of a text string.
 *
 * @param[in] text Text string to process.
 * @param[in, out] resource Loaded resource data. The \p type member must be set to the type of resource to load.
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_BAD_DATA Input data is malformed.
 * @retval LLKA_E_MISMATCHING_DATA Input data is malformed.
 * @retval LLKA_E_MISMATCHING_SIZES Input data is malformed.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_loadResourceText(const char *text, LLKA_Resource *resource);

LLKA_END_API_FUNCTIONS

#endif /* _LLKA_RESOURCE_LOADERS_H */
