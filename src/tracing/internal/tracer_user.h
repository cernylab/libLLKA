#ifndef _ECHMET_TRACER_USER_H
#define _ECHMET_TRACER_USER_H

#include "tracer_base.h"

#ifndef ECHMET_TRACER_DISABLE_TRACING

/*!
 * \def ECHMET_S(str)
 * Converts \str to C string
 */
#define ECHMET_S(str) \
	std::string{str}.c_str()

/*!
 * \def ECHMET_TRACE(TracerClass, TPID, Args...)
 * Logs data for the given \TracerClass and tracepoint
 *
 * @param TracerClass Tracer class
 * @param ID of the tracepoint to log
 * @param Args... Arguments for the logging function
 */
#define ECHMET_TRACE(TracerClass, TPID, ...) \
	::ECHMET::_ECHMET_TRACE<TracerClass, TracerClass::TPID>(__VA_ARGS__)

#define ECHMET_TRACE_NOARGS(TracerClass, TPID) \
	::ECHMET::_ECHMET_TRACE<TracerClass, TracerClass::TPID>()

#define ECHMET_TRACE_T1(TracerClass, TPID, T1, ...) \
	::ECHMET::_ECHMET_TRACE_T1<TracerClass, TracerClass::TPID, T1>(__VA_ARGS__)

#define ECHMET_TRACE_T2(TracerClass, TPID, T1, T2, ...) \
	::ECHMET::_ECHMET_TRACE_T2<TracerClass, TracerClass::TPID, T1, T2>(__VA_ARGS__)

#define ECHMET_TRACE_T3(TracerClass, TPID, T1, T2, T3, ...) \
	::ECHMET::_ECHMET_TRACE_T3<TracerClass, TracerClass::TPID, T1, T2, T3>(__VA_ARGS__)

#define ECHMET_TRACE_T4(TracerClass, TPID, T1, T2, T3, T4, ...) \
	::ECHMET::_ECHMET_TRACE_T4<TracerClass, TracerClass::TPID, T1, T2, T3, T4>(__VA_ARGS__)

#define ECHMET_TRACE_T5(TracerClass, TPID, T1, T2, T3, T4, T5, ...) \
	::ECHMET::_ECHMET_TRACE_T5<TracerClass, TracerClass::TPID, T1, T2, T3, T4, T5>(__VA_ARGS__)

/*!
 * \def ECHMET_TRACER_LOG(TracerClass)
 * Returns the complete log from a given \TracerClass
 *
 * @param TracerClass Tracer class
 */
#define ECHMET_TRACER_LOG(TracerClass) \
	TRACER_INSTANCE<TracerClass>().logged()

#else

inline constinit std::string __empty_trace_string{};

#define ECHMET_S(str) __empty_trace_string.c_str()
#define ECHMET_TRACE(TraceClass, TPID, ...)
#define ECHMET_TRACE_NOARGS(TraceClass, TPID)
#define ECHMET_TRACE_T1(TraceClass, TPID, ...)
#define ECHMET_TRACE_T2(TraceClass, TPID, ...)
#define ECHMET_TRACE_T3(TraceClass, TPID, ...)
#define ECHMET_TRACE_T4(TraceClass, TPID, ...)
#define ECHMET_TRACE_T5(TraceClass, TPID, ...)
#define ECHMET_TRACER_LOG(TracerClass) std::string{}

#endif // ECHMET_TRACER_DISABLE_TRACING

#endif // _ECHMET_TRACER_USER_H
