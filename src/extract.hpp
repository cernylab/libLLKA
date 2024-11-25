/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _EXTRACT_HPP
#define _EXTRACT_HPP

#include <extract.h>
#include <structure.hpp>

#include <util/templates.hpp>

#include <cstring>

static
auto _isAtomCandidate(const LLKA_Atom *atom, const LLKAInternal::AtomToExtract &ate)
{
    bool isCandidate = (
        atom->pdbx_PDB_model_num == ate.modelNum &&
        atom->label_seq_id == ate.seqId &&
        (std::strcmp(atom->label_asym_id, ate.asymId.c_str()) == 0) &&
        (std::strcmp(atom->label_comp_id, ate.compId.c_str()) == 0) &&
        (std::strcmp(atom->label_atom_id, ate.name.c_str()) == 0)
    );

    return isCandidate;
}

namespace LLKAInternal {

inline
auto getMatchingAtom(const LLKA_Structure *stru, const AtomToExtract &ate) -> LLKA_Atom *
{
    for (size_t idx = 0; idx < stru->nAtoms; idx++) {
        auto atom = getAtomPtr(*stru, idx);
        bool isCandidate = _isAtomCandidate(atom, ate);

        if (!isCandidate)
            continue;

        if (ate.altId == LLKA_NO_ALTID)
            return atom;
        else if (atom->label_alt_id == LLKA_NO_ALTID || atom->label_alt_id == ate.altId)
            return atom;
    }

    return nullptr;
}

inline
auto getMatchingAtom(const LLKA_StructureView *view, const AtomToExtract &ate) -> const LLKA_Atom *
{
    for (size_t idx = 0; idx < view->nAtoms; idx++) {
        const auto atom = getAtomPtr(*view, idx);
        bool isCandidate = _isAtomCandidate(atom, ate);

        if (!isCandidate)
            continue;

        if (ate.altId == LLKA_NO_ALTID)
            return atom;
        else if (atom->label_alt_id == LLKA_NO_ALTID || atom->label_alt_id == ate.altId)
            return atom;
    }

    return nullptr;
}

template <typename StructureType> requires LLKAStructureType<StructureType>
inline
auto extractAtoms(const StructureType *stru, const std::vector<AtomToExtract> &toExtract, LLKA_StructureView *extracted) -> LLKA_RetCode
{
    assert(toExtract.size() <= extracted->capacity);

    size_t atomIdx = 0;
    // This is crazy slow but we have to account for the fact that the atoms
    // may be in arbitrary order in the input structure.
    // We also want to return the atoms in the order specified by the "toExtract" vector
    for (const auto &ate : toExtract) {
        const auto atom = getMatchingAtom(stru, ate);
        if (!atom)
            return LLKA_E_MISSING_ATOMS;

        extracted->atoms[atomIdx++] = atom;
    }

    extracted->nAtoms = toExtract.size();

    return LLKA_OK;
}

} // namespace LLKAInternal

#endif // _EXTRACT_HPP
