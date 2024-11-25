/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_MAIN_H
#define _LLKA_MAIN_H

#include <stddef.h>
#include <stdint.h>

#include "llka_module.h"

#ifdef __cplusplus
    #include <cassert>
#else
    #include <assert.h>
#endif /* __cplusplus */

#define ENUM_FORCE_INT32_SIZE(ename) , _ECHMET_ST_ENUM_PLACEHOLDER ## ename = (0x7FFFFFFF)
#define ENUM_FORCE_INT32_SIZE_ITEM(ename) _ECHMET_ST_ENUM_PLACEHOLDER ## ename

#ifdef __cplusplus
    #include <type_traits>
#endif /* __cplusplus */

#ifdef __cplusplus
    #define LLKA_BEGIN_API_FUNCTIONS extern "C" {
    #define LLKA_END_API_FUNCTIONS }

    #define LLKA_IS_POD(theType) static_assert(std::is_trivial<theType>::value && std::is_standard_layout<theType>::value, "Type " #theType " is not a POD type");
#else
    #define LLKA_BEGIN_API_FUNCTIONS
    #define LLKA_END_API_FUNCTIONS

    #define LLKA_IS_POD(theType)
#endif /* __cplusplus */

typedef enum LLKA_RetCode {
    LLKA_OK                            = 0x0,     /*!< Success */
    LLKA_E_INVALID_ARGUMENT            = 0x1,     /*!< Invalid argument passed to a function */
    LLKA_E_MISMATCHING_SIZES           = 0x2,     /*!< Function expected arguments of the same size or length but the actual arguments' sizes did not match */
    LLKA_E_NOT_IMPLEMENTED             = 0x3,     /*!< Function is currently not implemented */
    LLKA_E_MULTIPLE_ALT_IDS            = 0x4,     /*!< Structure was expected to contain only atoms with the same alternate positions */
    LLKA_E_MISSING_ATOMS               = 0x5,     /*!< Structure was expected to contain some specific atoms */
    LLKA_E_MISMATCHING_DATA            = 0x6,     /*!< Data that were supposed to match did not match */
    LLKA_E_BAD_CLASSIFICATION_CLUSTERS = 0x7,     /*!< Classification cluster definitions passed to initialize classification context were invalid */
    LLKA_E_BAD_GOLDEN_STEPS            = 0x8,     /*!< Golden steps definitions passed to initialize classification context were invalid */
    LLKA_E_BAD_CONFALS                 = 0x9,     /*!< Confals definitions passed to initialize classification context were invalid */
    LLKA_E_BAD_AVERAGE_NU_ANGLES       = 0x10,    /*!< Nu angles definitions passed to initialize classification context were invalid */
    LLKA_E_BAD_CLASSIFICATION_LIMITS   = 0x11,    /*!< Classification tolerances and control parameters passed to initialize classification context were invalid */
    LLKA_E_NO_FILE                     = 0x12,    /*!< File does not exist */
    LLKA_E_CANNOT_READ_FILE            = 0x13,    /*!< File exists but cannot be read */
    LLKA_E_BAD_DATA                    = 0x14,    /*!< Data contains wrong or unexpected values */
    LLKA_E_NO_DATA                     = 0x15,    /*!< Data block is empty when it should not be */
    LLKA_E_NOTHING_TO_CLASSIFY         = 0x16     /*!< Empty data was passed to classification module
                                                       The usual cause of this error is when at attempt is made to classify a structure that does not contain any nucleic acid residues */
    ENUM_FORCE_INT32_SIZE(LLKA_RetCode)
} LLKA_RetCode;

typedef unsigned char LLKA_Bool;
#define LLKA_TRUE (unsigned char)1
#define LLKA_FALSE (unsigned char)0

/*!
 * Column-major ordered matrix of doubles
 */
typedef struct LLKA_Matrix {
    double *data;   /*!< Pointer to matrix data in column-major ordering */
    size_t nCols;   /*!< Number of columns */
    size_t nRows;   /*!< Number of rows */
} LLKA_Matrix;
LLKA_IS_POD(LLKA_Matrix)

/*!
 * Position in 3D space in cartesian coordinates
 */
typedef struct LLKA_Point {
    double x;   /*!< x coordinate */
    double y;   /*!< y coordinate */
    double z;   /*!< z coordinate */
} LLKA_Point;
LLKA_IS_POD(LLKA_Point)

/*!
 * Set of cartesian coorinates in 3D space
 */
typedef struct LLKA_Points {
    union {
        LLKA_Point *points; /*!< Points in 3D space as an array of <tt>LLKA_Point</tt> */
        double *raw;        /*!< Points in 3D space as a raw array of <tt>double</tt>s */
    }
#ifdef __cplusplus
    ;
#else
    data;   /* ANSI C does not allow anonymous unions. */
#endif
    size_t nPoints;         /*!< Number of points in the set */
} LLKA_Points;
LLKA_IS_POD(LLKA_Points)

LLKA_BEGIN_API_FUNCTIONS

/*!
 * Gets element of matrix
 *
 * @param[in] row Row coordinate
 * @param[in] col Column coordinate
 * @param[in] matrix The matrix
 *
 * @returns Matrix element at (row, column)
 */
LLKA_STATIC_INLINE
double LLKA_matrixGet(size_t row, size_t col, const LLKA_Matrix *matrix)
{
    assert(row < matrix->nRows && col < matrix->nCols);
    return matrix->data[matrix->nRows * col + row];
}

/*!
 * Sets element of matrix
 *
 * @param[in] value Value to set
 * @param[in] row Row coordinate
 * @param[in] col Column coordinate
 * @param[in,out] matrix The matrix
 */
LLKA_STATIC_INLINE
void LLKA_matrixSet(double v, size_t row, size_t col, const LLKA_Matrix *matrix)
{
    assert(row < matrix->nRows && col < matrix->nCols);
    matrix->data[matrix->nRows * col + row] = v;
}

/*!
 * Destroys matrix
 *
 * @param matrix The matrix
 */
LLKA_API void LLKA_CC LLKA_destroyMatrix(const LLKA_Matrix *matrix);

/*!
 * Destroys \p LLKA_Points structure
 *
 * @param[in] points LLKA_Points to destroy
 */
LLKA_API void LLKA_CC LLKA_destroyPoints(const LLKA_Points *points);

/*!
 * Destroys a string created by libLLKA
 *
 * @param[in] str String to destroy
 */
LLKA_API void LLKA_CC LLKA_destroyString(const char *str);

/*!
 * Returns a string representation of \p LLKA_RetCode.
 *
 * @param[in] tRet LLKA_RetCode
 *
 * @return String representation of the error code or <tt>"Unknown error"</tt> if the \p tRet value is not recognized
 */
LLKA_API const char * LLKA_CC LLKA_errorToString(LLKA_RetCode tRet);

/*!
 * Initializes matrix
 *
 * @param[in] rows Number of rows
 * @param[in] cols Number of columns
 * @param[out] matrix The matrix to initialize
 */
LLKA_API void LLKA_CC LLKA_initMatrix(size_t rows, size_t columns, LLKA_Matrix *matrix);

/*!
 * Creates \p LLKA_Points structure from array of \p LLKA_Point s. The returned \p LLKA_Points contains copies of the input points.
 *
 * @param[in] points Array of \p LLKA_Point s
 * @param[ing nPoints Number of elements in the \p points array
 *
 * @return LLKA_Points structure.
 */
LLKA_API LLKA_Points LLKA_CC LLKA_makePoints(const LLKA_Point *points, size_t nPoints);

/*!
 * Creates \p LLKA_Points structure with a specified number of points. Coordinates of all points are initially set to zero.
 *
 * @param[ing nPoints Number of points to create
 *
 * @return LLKA_Points structure.
 */
LLKA_API LLKA_Points LLKA_CC LLKA_zeroPoints(size_t nPoints);



LLKA_END_API_FUNCTIONS

#endif /* _LLKA_MAIN_H */
