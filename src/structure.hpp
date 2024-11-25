/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _STRUCTURE_HPP
#define _STRUCTURE_HPP

#include <llka_structure.h>

#include <cassert>

namespace LLKAInternal {

inline
auto initStructureView(LLKA_StructureView &view, size_t capacity)
{
    view.atoms = new const LLKA_Atom *[capacity];
    view.nAtoms = 0;
    view.capacity = capacity;
}

inline
auto makeStructureView(size_t capacity) -> LLKA_StructureView
{
    LLKA_StructureView view;

    view.atoms = new const LLKA_Atom *[capacity];
    view.nAtoms = 0;
    view.capacity = capacity;

    return view;
}

inline
auto structureViewToStructure(const LLKA_StructureView &view, LLKA_Structure &stru)
{
    assert(view.nAtoms > 0);

    stru.atoms = new LLKA_Atom[view.nAtoms];
    stru.nAtoms = view.nAtoms;

    for (size_t idx = 0; idx < view.nAtoms; idx++)
        LLKA_duplicateAtom(view.atoms[idx], &stru.atoms[idx]);
}

} // namespace LLKAInternal

#endif // _STRUCTURE_HPP
