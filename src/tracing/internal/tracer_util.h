#ifndef _ECHMET_TRACER_UTIL_H
#define _ECHMET_TRACER_UTIL_H

#include "tracer_types.h"
#include <map>
#include <tuple>
#include <vector>

namespace ECHMET {

/*!
 * Returns human-readable information about a tracepoint
 *
 * @tparam TracerClass Tracer class
 * @tparam TPID ID of the requested tracepoint
 * @return Tracepoint object
 */
template <typename TracerClass, TracerClass TPID>
Tracepoint<TracerClass> TRACEPOINT_INFO();

/*!
 * Returns ID of the first tracepoint for a given tracer
 *
 * @tparam TracerClass Tracer class
 * @return ID of the first tracepoint
 */
template <typename TracerClass>
inline
constexpr TracerClass FIRST_TRACEPOINT_ID();

/*!
 * Returns ID of the next tracepoint
 *
 * @tparam TracerClass Tracer class
 * @tparam TPID ID of the tracepoint to get the next tracepoint for
 * @return ID of the next tracepoint
 */
template <typename TracerClass, TracerClass TPID>
inline
constexpr TracerClass NEXT_TRACEPOINT_ID()
{
	typedef typename std::underlying_type<TracerClass>::type UType;
	return static_cast<TracerClass>(static_cast<UType>(TPID) + 1);
}

template <typename TracerClass>
inline
constexpr TracerClass LAST_TRACEPOINT_ID();

template <typename TracerClass>
struct TUTYPE {
	typedef typename std::underlying_type<TracerClass>::type type;
};

template <typename TracerClass, typename UType = typename std::underlying_type<TracerClass>::type>
inline
constexpr UType TUTYPE_CAST(const TracerClass &TPID)
{
	return static_cast<UType>(TPID);
}

template <typename RTPID, typename TracerClass>
inline
bool IS_TPID_VALID(const RTPID &tpid)
{
	static_assert(std::is_same<RTPID, typename TUTYPE<TracerClass>::type>::value, "Incompatible raw tracepoint IDs");

	const typename TUTYPE<TracerClass>::type first = TUTYPE_CAST(FIRST_TRACEPOINT_ID<TracerClass>());
	const typename TUTYPE<TracerClass>::type last = TUTYPE_CAST(LAST_TRACEPOINT_ID<TracerClass>());

	return (tpid >= first && tpid < last);
}

/*!
 * Builds a vector with information about all tracepoints available for
 * a given tracer. This function is recursive, each call appends a new
 * entry to the vector.
 *
 * @tparam TracerClass Tracer class
 * @tparam TPID ID of the currently being added tracepoint
 * @return Vector of Tracepoint<TracerClass> objects
 */
template <typename TracerClass, TracerClass TPID>
inline
void TRACEPOINT_INFO_BUILD(std::vector<std::tuple<TPIDInt, std::string>> &tracepointInfoVec)
{
#ifndef ECHMET_TRACER_DISABLE_TRACING
	const auto tpinfo = TRACEPOINT_INFO<TracerClass, TPID>();
	tracepointInfoVec.emplace_back(static_cast<TPIDInt>(tpinfo.ID), tpinfo.description);
	TRACEPOINT_INFO_BUILD<TracerClass, NEXT_TRACEPOINT_ID<TracerClass, TPID>()>(tracepointInfoVec);
#else
	(void)tracepointInfoVec;
	return;
#endif // TRACER_DISABLE_TRACING
}

template <typename TracerClass, TracerClass TPID>
inline
void TOGGLE_ALL_TRACEPOINTS(const bool state, std::map<TracerClass, bool> &enabledTracepoints)
{
#ifndef ECHMET_TRACER_DISABLE_TRACING
	enabledTracepoints[TPID] = state;
	TOGGLE_ALL_TRACEPOINTS<TracerClass, NEXT_TRACEPOINT_ID<TracerClass, TPID>()>(state, enabledTracepoints);
#else
	(void)enabledTracepoints;
	return;
#endif // ECHMET_TRACER_DISABLE_TRACING
}

/*!
 * Logging functor.
 *
 * @tparam Tracer class
 * @tparam TPID ID of the tracepoint
 * @tparam TFArgs Optional template arguments template of the logging function
 */
template <typename TracerClass, TracerClass TPID, typename... TFArgs>
class TracepointLogger;

} // namespace ECHMET

/*!
 * \def ECHMET_MAKE_TRACEPOINT(TracerClass, TPID, description)
 * Defines functions necessary to query information about tracepoints for the given \TracerClass
 *
 * @param TracerClass Tracer class
 * @param TPID ID of the tracepoint
 * @param description Human-readable description of the tracepoint
 */
#define ECHMET_MAKE_TRACEPOINT(TracerClass, TPID, description) \
	template <> \
	inline \
	Tracepoint<TracerClass> TRACEPOINT_INFO<TracerClass, TracerClass::TPID>() { return Tracepoint<TracerClass>{TracerClass::TPID, description}; }

/*!
 * \def ECHMET_MAKE_TRACEPOINT_NOINLINE(TracerClass, TPID, description)
 * Defines functions necessary to query information about tracepoints for the given \TracerClass
 * Use this one to define tracepoints in auxiliary compilation units to prevent linking issues.
 *
 * @param TracerClass Tracer class
 * @param TPID ID of the tracepoint
 * @param description Human-readable description of the tracepoint
 */
#define ECHMET_MAKE_TRACEPOINT_NOINLINE(TracerClass, TPID, description) \
	template <> \
	Tracepoint<TracerClass> TRACEPOINT_INFO<TracerClass, TracerClass::TPID>() { return Tracepoint<TracerClass>{TracerClass::TPID, description}; }

/*!
 * \def ECHMET_MAKE_LOGGER(TracerClass, TPID, Args...)
 * Defines logging functor for a given tracer and its tracepoint
 *
 * @param TracerClass Tracer class
 * @param TPID ID of the tracepoint whose logging function is being declared
 * @Args... Argument template of the logging function
 */
#define ECHMET_BEGIN_MAKE_LOGGER(TracerClass, TPID, ...) \
	template <> \
	class TracepointLogger<::TracerClass, ::TracerClass::TPID> { \
	public: \
		static std::string call(__VA_ARGS__)

#define ECHMET_BEGIN_MAKE_LOGGER_NOARGS(TracerClass, TPID) \
	template <> \
	class TracepointLogger<::TracerClass, ::TracerClass::TPID> { \
	public: \
		static std::string call()

#define ECHMET_BEGIN_MAKE_LOGGER_T1(TracerClass, TPID, ...) \
	template <typename T1> \
	class TracepointLogger<::TracerClass, ::TracerClass::TPID, T1> { \
	public: \
		static std::string call(__VA_ARGS__)

#define ECHMET_BEGIN_MAKE_LOGGER_T2(TracerClass, TPID, ...) \
	template <typename T1, typename T2> \
	class TracepointLogger<::TracerClass, ::TracerClass::TPID, T1, T2> { \
	public: \
		static std::string call(__VA_ARGS__)

#define ECHMET_BEGIN_MAKE_LOGGER_T3(TracerClass, TPID, ...) \
	template <typename T1, typename T2, typename T3> \
	class TracepointLogger<::TracerClass, ::TracerClass::TPID, T1, T2, T3> { \
	public: \
		static std::string call(__VA_ARGS__)

#define ECHMET_BEGIN_MAKE_LOGGER_T4(TracerClass, TPID, ...) \
	template <typename T1, typename T2, typename T3, typename T4> \
	class TracepointLogger<::TracerClass, ::TracerClass::TPID, T1, T2, T3, T4> { \
	public: \
		static std::string call(__VA_ARGS__)

#define ECHMET_BEGIN_MAKE_LOGGER_T5(TracerClass, TPID, ...) \
	template <typename T1, typename T2, typename T3, typename T4, typename T5> \
	class TracepointLogger<::TracerClass, ::TracerClass::TPID, T1, T2, T3, T4, T5> { \
	public: \
		static std::string call(__VA_ARGS__)

/* Five template parameters ought to be enough for everyone! */

#define ECHMET_END_MAKE_LOGGER };

#define ECHMET_LOGGER_ADD_OVERLOAD(...) \
	static std::string call(__VA_ARGS__)

#ifndef _ECHMET_TRACER_IMPL_SECTION
	#ifndef ECHMET_TRACER_DISABLE_TRACING
	/*!
	 * \def ECHMET_MAKE_TRAECPOINT_IDS(TracerClass, first, last)
	 * Defines functions necessary to build a list of tracepoints for the given \TracerClass
	 *
	 * @param TracerClass Tracer class
	 * @param first ID of the first tracepoint
	 * @param last ID of the last tracepoint
	 */
	#define ECHMET_MAKE_TRACEPOINT_IDS(TracerClass, first, last) \
		namespace ECHMET { \
			template <> \
			inline \
			constexpr ::TracerClass FIRST_TRACEPOINT_ID<::TracerClass>() { return ::TracerClass::first; } \
			template <> \
			inline \
			constexpr ::TracerClass LAST_TRACEPOINT_ID<::TracerClass>() { return ::TracerClass::last; } \
			template <> \
			inline \
			void TRACEPOINT_INFO_BUILD<::TracerClass, ::TracerClass::last>(std::vector<std::tuple<TPIDInt, std::string>> &tracepointInfoVec) \
			{ \
				(void)tracepointInfoVec; \
				return; /* No-op for the last dummy tracepoint */ \
			} \
			template <> \
			inline \
			void TOGGLE_ALL_TRACEPOINTS<::TracerClass, ::TracerClass::last>(const bool state, std::map<TracerClass, bool> &enabledTracepoints) \
			{ \
				enabledTracepoints[TracerClass::last] = state; \
			} \
			template <> \
			inline \
			bool IS_TPID_VALID<::TracerClass, ::TracerClass>(const ::TracerClass &) { return true; } \
		} // namespace ECHMET
	#else
	#define ECHMET_MAKE_TRACEPOINT_IDS(TracerClass, first, last)
	namespace ECHMET {
		template <>
		constexpr __DUMMY_TRACER_CLASS FIRST_TRACEPOINT_ID<__DUMMY_TRACER_CLASS>() { return __DUMMY_TRACER_CLASS::NONE; }
		template <>
		inline
		void TRACEPOINT_INFO_BUILD<__DUMMY_TRACER_CLASS, __DUMMY_TRACER_CLASS::NONE>(std::vector<std::tuple<TPIDInt, std::string>> &)
		{
			return;
		}
		template <>
		inline
		void TOGGLE_ALL_TRACEPOINTS<__DUMMY_TRACER_CLASS, __DUMMY_TRACER_CLASS::NONE>(const bool, std::map<__DUMMY_TRACER_CLASS, bool> &)
		{
			return;
		}
	} // namespace ECHMET
	#endif // TRACER_DISABLE_TRACING
#else

#define ECHMET_MAKE_TRACEPOINT_IDS(TracerClass, first, last)

#endif // _ECHMET_TRACER_IMPL_SECTION

#endif // _ECHMET_TRACER_UTIL_H
