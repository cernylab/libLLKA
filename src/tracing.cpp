/* vim: sw=4 ts=4 sts=4 expandtab */

#include <llka_tracing.h>

#include "tracing/internal/tracer_util.h"
#include "tracing/llka_tracer_impl.h"

#include <cstring>

#ifndef LLKA_DISABLE_TRACING

ECHMET_MAKE_TRACER_NOINLINE(LLKATracing)

#else

ECHMET_MAKE_TRACER_NOINLINE(ECHMET::__DUMMY_TRACER_CLASS)

#endif // LLKA_DISABLE_TRACING

LLKA_API void LLKA_CC LLKA_destroyTrace(const char *trace)
{
    delete [] trace;
}

void LLKA_CC LLKA_destroyTracepointInfo(LLKA_TracepointInfos *infos)
{
    for (size_t idx = 0; idx < infos->nInfos; idx++) {
        auto &info = infos->infos[idx];
        delete [] info.description;
    }

    delete [] infos->infos;
}

void LLKA_CC LLKA_toggleAllTracepoints(LLKA_Bool state)
{
	if (state)
		ECHMET::TRACER_INSTANCE<LLKATracing>().enableAllTracepoints();
	else
		ECHMET::TRACER_INSTANCE<LLKATracing>().disableAllTracepoints();
}

void LLKA_CC LLKA_toggleTracepoint(int32_t TPID, LLKA_Bool state)
{
	if (state)
		ECHMET::TRACER_INSTANCE<LLKATracing>().enableTracepoint(TPID);
	else
		ECHMET::TRACER_INSTANCE<LLKATracing>().disableTracepoint(TPID);
}

const char * LLKA_CC LLKA_trace(LLKA_Bool dontClear)
{
#ifdef ECHMET_TRACER_DISABLE_TRACING
	(void)dontClear;
	return nullptr;
#else
	const auto &trace = ECHMET::TRACER_INSTANCE<LLKATracing>().logged(dontClear);

    char *rawTrace = new char[trace.length() + 1];
    std::strcpy(rawTrace, trace.c_str());

    return rawTrace;
#endif // ECHMET_TRACER_DISABLE_TRACING
}

LLKA_TracepointInfos LLKA_CC LLKA_tracepointInfo()
{
#ifdef ECHMET_TRACER_DISABLE_TRACING
    static const LLKA_TracepointInfos empty = {};
	return empty;
#else // ECHMET_TRACER_DISABLE_TRACING

	auto &tracerInstance = ECHMET::TRACER_INSTANCE<LLKATracing>();
	auto tracepoints = tracerInstance.tracepoints();

    LLKA_TracepointInfos infos{};
    infos.infos = new LLKA_TracepointInfo[tracepoints.size()];
    infos.nInfos = tracepoints.size();

    for (size_t idx = 0; idx < infos.nInfos; idx++) {
        const auto &[TPID, desc] = tracepoints[idx];
        auto &info = infos.infos[idx];

        auto rawDesc = new char[desc.length() + 1];
        std::strcpy(rawDesc, desc.c_str());

        info.TPID = TPID;
        info.description = rawDesc;
    }

	return infos;
#endif // ECHMET_TRACER_DISABLE_TRACING
}

LLKA_Bool LLKA_CC LLKA_tracepointState(int32_t TPID)
{
	return ECHMET::TRACER_INSTANCE<LLKATracing>().isTracepointEnabled(TPID);
}
