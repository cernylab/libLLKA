/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _RESIDUES_H
#define _RESIDUES_H

#include <cassert>
#include <set>

#include "standard_residues.h"
#include "non_standard_residues.h"

namespace LLKAInternal {

static
auto initializeAllKnownResidues()
{
    std::set<LLKABones::ANString> all;

    // Standard residues are just a vector of strings
    for (const auto &name : KNOWN_STANDARD_RESIDUES) {
        assert(!all.contains(name));
        all.insert(name);
    }

    // Non-standard residues are a map
    for (const auto &[key, _] : KNOWN_NON_STANDARD_RESIDUES) {
        assert(!all.contains(key));
        all.insert(key);
    }

    return all;
}

inline const std::set<LLKABones::ANString> ALL_KNOWN_RESIDUES = initializeAllKnownResidues();

inline
bool isKnownResidue(const LLKABones::ANString &name)
{
    // Check priority residues first. ALL_KNOWN_RESIDUES contains them too
    // but the lookup would take more time because we'd have to do a search in a set
    if (std::find(PRIORITY_RESIDUES.cbegin(), PRIORITY_RESIDUES.cend(), name) != PRIORITY_RESIDUES.cend())
        return true;

    return ALL_KNOWN_RESIDUES.contains(name);
}

inline
bool isStandardResidue(const LLKABones::ANString &name)
{
    // Check priority residues first. ALL_KNOWN_RESIDUES contains them too
    // but the lookup would take more time because we'd have to do a search in a set
    if (std::find(PRIORITY_RESIDUES.cbegin(), PRIORITY_RESIDUES.cend(), name) != PRIORITY_RESIDUES.cend())
        return true;

    return std::find(KNOWN_STANDARD_RESIDUES.cbegin(), KNOWN_STANDARD_RESIDUES.cend(), name) != KNOWN_STANDARD_RESIDUES.cend();
}

inline
bool findBaseKind(const LLKABones::ANString &name, LLKA_BaseKind &kind)
{
    auto itStd = STANDARD_RESIDUES_BASES.find(name);
    if (itStd != STANDARD_RESIDUES_BASES.cend()) {
        kind = itStd->second;
        return true;
    }

    auto itNonStd = KNOWN_NON_STANDARD_RESIDUES.find(name);
    if (itNonStd != KNOWN_NON_STANDARD_RESIDUES.cend()) {
        kind = LLKA_NON_STANDARD_BASE;
        return true;
    }

    return false;
}

inline
const LLKABones::Bone & findBone(const LLKABones::ANString &name)
{
    if (isStandardResidue(name)) {
        const auto base = STANDARD_RESIDUES_BASES.at(name);
        switch (base) {
        case LLKA_PURINE:
            return LLKABones::StandardPurineBone;
        case LLKA_PYRIMIDINE:
            return LLKABones::StandardPyrimidineBone;
        default:
            assert(false);
        }
    }

    auto nonStdId = KNOWN_NON_STANDARD_RESIDUES.find(name);
    if (nonStdId == KNOWN_NON_STANDARD_RESIDUES.cend()) {
        return LLKABones::StandardPurineBone; // Fall back to this because we don't know what else to do.
    }

    return nonStdId->second;
}

} // namespace LLKAInternal

#endif // _RESIDUES_H
