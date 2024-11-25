/* vim : sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_TRACING_H
#define _LLKA_TRACING_H

#include "llka_main.h"

typedef struct LLKA_TracepointInfo {
    int32_t TPID;               /*!< ID of the tracepoint */
    const char *description;    /*!< Description of what kind of information or event is logged by the tracepoint */
} LLKA_TracepointInfo;
LLKA_IS_POD(LLKA_TracepointInfo)

typedef struct LLKA_TracepointInfos {
    LLKA_TracepointInfo *infos;    /*!< Array of \p LLKA_TracepointInfo s*/
    size_t nInfos;                 /*!< Length of the array */
} LLKA_TracepointInfos;
LLKA_IS_POD(LLKA_TracepointInfos)

LLKA_BEGIN_API_FUNCTIONS

/*!
 * Destroys trace string
 *
 * @param[in] trace Trace string to destroy
 */
LLKA_API void LLKA_CC LLKA_destroyTrace(const char *trace);

/*!
 * Destroys LLKA_TracepointInfos
 *
 * @param[in] infos LLKA_TracepointInfos object to destroy
 */
LLKA_API void LLKA_CC LLKA_destroyTracepointInfo(LLKA_TracepointInfos *infos);

/*!
 * Sets all tracepoints to the given state.
 *
 * @param[in] state If \p true all tracepoints will be enabled and vice versa.
 */
LLKA_API void LLKA_CC LLKA_toggleAllTracepoints(LLKA_Bool state);

/*!
 * Set state of one tracepoint.
 *
 * @param[in] TPID Internal ID of the tracepoint to set.
 * @param[in] state If \p true the tracepoint will be enabled and vice versa.
 */
LLKA_API void LLKA_CC LLKA_toggleTracepoint(int32_t TPID, LLKA_Bool state);

/*!
 * Returns the complete trace.
 *
 * @param[in] dontClear If \p true the trace log will not be cleared.
 *
 * @return String containing the whole trace.
 */
LLKA_API const char * LLKA_CC LLKA_trace(LLKA_Bool dontClear);

/*!
 * Returns information about available tracepoints.
 *
 * @retval Pointer to a vector of all available tracepoints. May be \p NULL if
 *         the operation fails or if no tracepoins are available.
 */
LLKA_API LLKA_TracepointInfos LLKA_CC LLKA_tracepointInfo();

/*!
 * Returns the state of a given tracepoint.
 *
 * @param[in] TPID Internal ID of the tracepoint whose state is requested.
 *
 * @retval true if the tracepoint is enabled and vice versa.
 */
LLKA_API LLKA_Bool LLKA_CC LLKA_tracepointState(int32_t TPID);

LLKA_END_API_FUNCTIONS

#endif /* _LLKA_TRACING_H */
