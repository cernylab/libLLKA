#ifndef _ECHMET_TRACER_TYPES_H
#define _ECHMET_TRACER_TYPES_H

#include <string>

namespace ECHMET {

#ifdef ECHMET_TRACER_DISABLE_TRACING

enum class __DUMMY_TRACER_CLASS : bool {
	NONE
};

#endif // ECHMET_TRACER_DISABLE_TRACING

typedef int32_t TPIDInt;

template <typename TracepointIDs>
class Tracepoint
{
	static_assert(std::is_enum<TracepointIDs>::value, "TracepointIDs is not an enum");

public:
	const TracepointIDs ID;
	const std::string description;
};

} // namespace ECHMET

#endif // _ECHMET_TRACER_TYPES_H
