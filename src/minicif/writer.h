/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_minicif.h>

#include <string>

namespace LLKAInternal::MiniCif {

auto dataToString(const LLKA_CifData &cifData, const bool pretty) -> std::string;

} // namespace LLKAInternal::MiniCif
