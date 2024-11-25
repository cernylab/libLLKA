/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_MINICIF_CATEGORIES_ATOM_SITE_H
#define _LLKA_MINICIF_CATEGORIES_ATOM_SITE_H

#include <llka_structure.h>

#include "../schema.hpp"

namespace LLKAInternal::MiniCif::Categories {

template <typename _Type, _Type LLKA_Atom::* Field>
struct LLKA_AtomSetter {
	using Type = _Type;
	static constexpr auto set(LLKA_Atom &atom, const Type &v) { atom.*Field = v; }
};

template <double LLKA_Point::* Field>
struct LLKA_AtomCoordsSetter {
	using Type = double;
	static constexpr auto set(LLKA_Atom &atom, const double v) { atom.coords.*Field = v; }
};

using AtomSite = Schema<
	"atom_site",
	LLKA_Atom,
	//TagValue<"group_PDB", LLKA_StructureSetter<typename _Type, _Type LLKA_Structure::*Field>
	TagValue<"id", LLKA_AtomSetter<uint32_t, &LLKA_Atom::id>>,
	TagValue<"type_symbol", LLKA_AtomSetter<const char *, &LLKA_Atom::type_symbol>>,
	TagValue<"label_atom_id", LLKA_AtomSetter<const char *, &LLKA_Atom::label_atom_id>>,
	TagValue<"label_entity_id", LLKA_AtomSetter<const char *, &LLKA_Atom::label_entity_id>>,
	TagValue<"label_comp_id", LLKA_AtomSetter<const char *, &LLKA_Atom::label_comp_id>>,
	TagValue<"label_asym_id", LLKA_AtomSetter<const char *, &LLKA_Atom::label_asym_id>>,
	TagValue<"auth_atom_id", LLKA_AtomSetter<const char *, &LLKA_Atom::auth_atom_id>, DefaultValue<const char *, nullptr>>, // Fix up later
	TagValue<"auth_comp_id", LLKA_AtomSetter<const char *, &LLKA_Atom::auth_comp_id>, DefaultValue<const char *, nullptr>>, // Fix up later
	TagValue<"auth_asym_id", LLKA_AtomSetter<const char *, &LLKA_Atom::auth_asym_id>, DefaultValue<const char *, nullptr>>, // Fix up later
	TagValue<"cartn_x", LLKA_AtomCoordsSetter<&LLKA_Point::x>>,
	TagValue<"cartn_y", LLKA_AtomCoordsSetter<&LLKA_Point::y>>,
	TagValue<"cartn_z", LLKA_AtomCoordsSetter<&LLKA_Point::z>>,
	TagValue<"label_seq_id", LLKA_AtomSetter<int32_t, &LLKA_Atom::label_seq_id>, DefaultValue<int32_t, 0>>,
	TagValue<"auth_seq_id", LLKA_AtomSetter<int32_t, &LLKA_Atom::auth_seq_id>>,
	TagValue<"pdbx_pdb_model_num", LLKA_AtomSetter<int32_t, &LLKA_Atom::pdbx_PDB_model_num>, DefaultValue<int32_t, 1>>,
	TagValue<"pdbx_pdb_ins_code", LLKA_AtomSetter<const char *, &LLKA_Atom::pdbx_PDB_ins_code>, DefaultValue<const char *, nullptr>>, // Fix up later
	TagValue<"label_alt_id", LLKA_AtomSetter<char, &LLKA_Atom::label_alt_id>, DefaultValue<char, LLKA_NO_ALTID>>
>;

} // namespace LLKAInternal::MiniCif::Categories

#endif // _LLKA_MINICIF_CATEGORIES_ATOM_SITE_H
