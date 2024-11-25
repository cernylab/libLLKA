// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_SUPERPOSITION_H
#define _LLKA_SUPERPOSITION_H

#include "llka_structure.h"

LLKA_BEGIN_API_FUNCTIONS

/*!
 * Applies transformation matrix on an array of points.
 *
 * @param[in,out] what Points to be transformed
 * @param[in] matrix Transformation matrix
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_INVALID_ARGUMENT Wrong transformation matrix
 */
LLKA_API LLKA_RetCode LLKA_applyTransformationPoints(LLKA_Points *what, const LLKA_Matrix *matrix);

/*!
 * Applies transformation matrix on atoms of a structure.
 *
 * @param[in,out] what Structure to be transformed
 * @param[in] matrix Transformation matrix
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_INVALID_ARGUMENT Wrong transformation matrix
 */
LLKA_API LLKA_RetCode LLKA_applyTransformationStructure(LLKA_Structure *what, const LLKA_Matrix *matrix);

/*!
 * Calculates the centroid of the given set of points.
 *
 * @param[in] points Set of points to calculate the centroid for
 *
 * @return Calculated centroid
 */
LLKA_API LLKA_Point LLKA_CC LLKA_centroidPoints(const LLKA_Points *points);

/*!
 * Calculates the centroid of a structure.
 *
 * @param[in] stru The structure to calculate the centroid for.
 *
 * @return Calculated centroid
 */
LLKA_API LLKA_Point LLKA_CC LLKA_centroidStructure(const LLKA_Structure *stru);

/*!
 * Calculates the Root Mean Square Distance (RMSD) between two sets of points.
 * The sets of points must have the same size. If the two sets of points are not ordered
 * in the same way the calculated RMSD value will be incorrect.
 *
 * @param[in] a First set of points
 * @param[in] b Second set of points
 * @param[out] rmsd Calculated RMSD
 *
 * @retval LLKA_OK RMSD was successfully calculated
 * @retval LLKA_E_MISMATCHING_SIZES Passed sets of points do not have the same size
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_rmsdPoints(const LLKA_Points *a, const LLKA_Points *b, double *rmsd);

/*!
 * Calculates the Root Mean Square Distance (RMSD) between two structures.
 * The structures must have the same number of atoms and the atoms must be ordered in
 * the same way. Otherwise the calculated RMSD value will be incorrect.
 *
 * @param[in] a First structure
 * @param[in] b Second structure
 * @param[out] rmsd Calculated RMSD
 *
 * @retval LLKA_OK RMSD was successfully calculated
 * @retval LLKA_E_MISMATCHING_SIZES Passed sets of points do not have the same size
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_rmsdStructures(const LLKA_Structure *a, const LLKA_Structure *b, double *rmsd);

/*!
 * Applies the Kabsch algorithm to superpose two sets of points onto each other.
 * Input sets of points are expected to be ordered and of the same size as the Kabsch algorithm requires.
 *
 * @param[in,out] what Set of points to be superposed onto the target
 * @param[in] onto Target of the superposition
 * @param[out] rmsd RMSD of the two sets of points after superposition.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_superposePoints(LLKA_Points *what, const LLKA_Points *onto, double *rmsd);

/*!
 * Applies the Kabsch algoritm to superpose two structures onto each other.
 * The structures must contain the same number of atoms ordered in the same way.
 *
 * @param[in,out] what Structure to be superposed on the target
 * @param[in] onto Target of the superposition
 * @param[out] rmsd RMSD of the two structures after superposition
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_MISMATCHING_SIZES Superposition objects do not have the same number of elements
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_superposeStructures(LLKA_Structure *what, const LLKA_Structure *onto, double *rmsd);

/*!
 * Applies the Kabsch algoritm to superpose two structures onto each other.
 * The structures must contain the same number of atoms ordered in the same way.
 *
 * @param[in,out] what Structure to be superposed on the target
 * @param[in] onto Target of the superposition
 * @param[out] rmsd RMSD of the two structures after superposition
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_MISMATCHING_SIZES Superposition objects do not have the same number of elements
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_superposeStructuresView(LLKA_Structure *what, const LLKA_StructureView *onto, double *rmsd);


/*!
 * Applies the Kabsch algorithm to calculate transformation matrix for superposition of two sets of points onto each other.
 * Input sets of points are expected to be ordered and of the same size as the Kabsch algorithm requires.
 *
 * @param[in] what Set of points to be superposed onto the target
 * @param[in] onto Target of the superposition
 * @param[out] matrix Transformation matrix of the superposition
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_MISMATCHING_SIZES Passed sets of points do not have the same size
 * @retval LLKA_E_MISMATCHING_SIZES Passed sets of points are empty
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_superpositionMatrixPoints(const LLKA_Points *what, const LLKA_Points *onto, LLKA_Matrix *matrix);


/*!
 * Applies the Kabsch algoritm to calculate transformation matrix for superposition of two structures onto each other.
 * The structures must contain the same number of atoms ordered in the same way.
 *
 * @param[in,out] what Structure to be superposed on the target
 * @param[in] onto Target of the superposition
 * @param[out] matrix Transformation matrix of the superposition
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_MISMATCHING_SIZES Superposition objects do not have the same number of elements
 * @retval LLKA_E_MISMATCHING_SIZES Superposition objects are empty
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_superpositionMatrixStructures(const LLKA_Structure *what, const LLKA_Structure *onto, LLKA_Matrix *matrix);

/*!
 * Applies the Kabsch algoritm to calculate transformation matrix for superposition of two structures onto each other.
 * The structures must contain the same number of atoms ordered in the same way.
 *
 * @param[in,out] what Structure to be superposed on the target
 * @param[in] onto Target of the superposition
 * @param[out] matrix Transformation matrix of the superposition
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_MISMATCHING_SIZES Superposition objects do not have the same number of elements
 * @retval LLKA_E_MISMATCHING_SIZES Superposition objects are empty
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_superpositionMatrixStructureViews(const LLKA_StructureView *what, const LLKA_StructureView *onto, LLKA_Matrix *matrix);

LLKA_END_API_FUNCTIONS

#endif // _LLKA_SUPERPOSITION_H
