/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_MINICIF_CATEGORIES_ENTRY_H
#define _LLKA_MINICIF_CATEGORIES_ENTRY_H

#include <llka_minicif.h>

#include "../schema.hpp"

namespace LLKAInternal::MiniCif::Categories {

template <typename _Type, _Type LLKA_StructureEntry::* Field>
struct EntrySetter {
	using Type = _Type;
	static constexpr auto set(LLKA_StructureEntry &e, const Type &v) { e.*Field = v; }
};

using Entry = Schema<
    "entry",
    LLKA_StructureEntry,
    TagValue<"id", EntrySetter<const char *, &LLKA_StructureEntry::id>>
>;

} // namespace LLKAInternal::MiniCif::Categories

#endif // _LLKA_MINICIF_CATEGORIES_ENTRY_H
