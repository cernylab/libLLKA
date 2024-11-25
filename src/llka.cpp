/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_main.h>

#include <cstring>

#define LLKA_RETCODE_CASE(code) case code: return #code

void LLKA_CC LLKA_destroyMatrix(const LLKA_Matrix *matrix)
{
    delete [] matrix->data;
}

void LLKA_CC LLKA_initMatrix(size_t rows, size_t cols, LLKA_Matrix *matrix)
{
    matrix->data = new double[rows * cols];
    matrix->nRows = rows;
    matrix->nCols = cols;
}

void LLKA_CC LLKA_destroyPoints(const LLKA_Points *points)
{
    delete [] points->points;
}

void LLKA_CC LLKA_destroyString(const char *str)
{
    delete [] str;
}

const char * LLKA_CC LLKA_errorToString(LLKA_RetCode tRet)
{
    switch (tRet) {
    LLKA_RETCODE_CASE(LLKA_OK);
    LLKA_RETCODE_CASE(LLKA_E_INVALID_ARGUMENT);
    LLKA_RETCODE_CASE(LLKA_E_MISMATCHING_SIZES);
    LLKA_RETCODE_CASE(LLKA_E_NOT_IMPLEMENTED);
    LLKA_RETCODE_CASE(LLKA_E_MULTIPLE_ALT_IDS);
    LLKA_RETCODE_CASE(LLKA_E_MISSING_ATOMS);
    LLKA_RETCODE_CASE(LLKA_E_MISMATCHING_DATA);
    LLKA_RETCODE_CASE(LLKA_E_BAD_CLASSIFICATION_CLUSTERS);
    LLKA_RETCODE_CASE(LLKA_E_BAD_GOLDEN_STEPS);
    LLKA_RETCODE_CASE(LLKA_E_BAD_CONFALS);
    LLKA_RETCODE_CASE(LLKA_E_BAD_AVERAGE_NU_ANGLES);
    LLKA_RETCODE_CASE(LLKA_E_BAD_CLASSIFICATION_LIMITS);
    LLKA_RETCODE_CASE(LLKA_E_NO_FILE);
    LLKA_RETCODE_CASE(LLKA_E_CANNOT_READ_FILE);
    LLKA_RETCODE_CASE(LLKA_E_BAD_DATA);
    LLKA_RETCODE_CASE(LLKA_E_NO_DATA);
    LLKA_RETCODE_CASE(LLKA_E_NOTHING_TO_CLASSIFY);
    case ENUM_FORCE_INT32_SIZE_ITEM(LLKA_RetCode): return "Unknown return code";
    }

    return "Unknown return code";
}

LLKA_Points LLKA_CC LLKA_makePoints(const LLKA_Point *points, size_t nPoints)
{
    LLKA_Points pts{
        .points = nullptr,
        .nPoints = nPoints
    };

    if (nPoints < 1)
        return pts;

    // TODO: We need to allocate aligned and padded memory blocks here
    pts.points = new LLKA_Point[nPoints];
    for (size_t idx = 0; idx < nPoints; idx++)
        pts.points[idx] = points[idx];

    return pts;
}


LLKA_Points LLKA_CC LLKA_zeroPoints(size_t nPoints)
{
    if (nPoints == 0) {
        return {
            .points = nullptr,
            .nPoints = 0
        };
    }

    // TODO: We need to allocate aligned and padded memory blocks here
    LLKA_Points pts{
        .points = new LLKA_Point[nPoints],
        .nPoints = nPoints
    };
    std::memset(pts.raw, 0, sizeof(LLKA_Point) * nPoints);

    return pts;
}
