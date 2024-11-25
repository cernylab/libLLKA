/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "llka_cpp.h"

#include <stdexcept>

#ifdef LLKA_PLATFORM_EMSCRIPTEN
#include <emscripten/bind.h>
#endif // LLKA_PLATFORM_EMSCRIPTEN

namespace DNATCOExtras {

class LLKA_CPP_API Error : public std::exception {
public:
    Error(LLKA_RetCode tRet) noexcept;
    auto retcode() const noexcept -> LLKA_RetCode;
    auto what() const noexcept -> const char *;

private:
    LLKA_RetCode m_tRet;
};

auto addDNATCOCategoriesToCif(
    LLKA::CifData cifData,
    const LLKA::AttemptedClassifiedSteps &attemptedSteps,
    const LLKA_AverageConfal &avgConfal,
    const LLKA::Structures &steps,
    const std::string &entryId
) -> LLKA::CifData;

} // namespace DNATCOExtras

#ifdef LLKA_PLATFORM_EMSCRIPTEN

EMSCRIPTEN_BINDINGS(DNATCOExtras) {
    emscripten::function("addDNATCOCategoriesToCif", &DNATCOExtras::addDNATCOCategoriesToCif);
}

#endif // LLKA_PLATFORM_EMSCRIPTEN
