/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_STRUCTURE_UTIL_H
#define _LLKA_STRUCTURE_UTIL_H

#include <llka_structure.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <numeric>
#include <tuple>
#include <vector>

namespace LLKAInternal {

auto dinucleotideToSteps(const LLKA_Structure &firstNucl, const LLKA_Structure &secondNucl) -> std::vector<LLKA_Structure>;

inline
auto isSameChain(const LLKA_Atom &atom, int32_t pdbx_PDB_model_num, const char *label_asym_id)
{
    return
        atom.pdbx_PDB_model_num == pdbx_PDB_model_num
        &&
        (std::strcmp(atom.label_asym_id, label_asym_id) == 0);
}

inline
auto isSameChain(const LLKA_Atom &lhs, const LLKA_Atom &rhs)
{
    return
        lhs.pdbx_PDB_model_num == rhs.pdbx_PDB_model_num
        &&
        (std::strcmp(lhs.label_asym_id, rhs.label_asym_id) == 0);
}

inline
auto isSameResidue(const LLKA_Atom &atom, int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id, char label_alt_id)
{
    bool maybeSame =
        atom.pdbx_PDB_model_num == pdbx_PDB_model_num
        &&
        atom.label_seq_id == label_seq_id
        &&
        (std::strcmp(atom.label_asym_id, label_asym_id) == 0);

    maybeSame &= (label_alt_id == LLKA_NO_ALTID) ? true : (atom.label_alt_id == label_alt_id);

    return maybeSame;
}

inline
auto isSameResidue(const LLKA_Atom &lhs, const LLKA_Atom &rhs)
{
    bool maybeSame =
        lhs.pdbx_PDB_model_num == rhs.pdbx_PDB_model_num
        &&
        lhs.label_seq_id == rhs.label_seq_id
        &&
        (std::strcmp(lhs.label_asym_id, rhs.label_asym_id) == 0);

    maybeSame &= (lhs.label_alt_id == LLKA_NO_ALTID || rhs.label_alt_id == LLKA_NO_ALTID) || (lhs.label_alt_id == rhs.label_alt_id);

    return maybeSame;
}

enum class StrucutreMergePolicy {
    SHALLOW,
    MOVE
};

// PERF: Consider making this a function with variable number of arguments to avoid creating
// any intermediate data structures altogether
template <size_t NStructures, StrucutreMergePolicy Policy>
inline
auto mergeStructures(const std::array<LLKA_Structure *, NStructures> &toMerge)
{
    LLKA_Structure merged{};

    const auto N = std::accumulate(
        toMerge.cbegin(),
        toMerge.cend(),
        0ULL,
        [](size_t acc, const LLKA_Structure *s) { return acc + s->nAtoms; }
    );

    merged.atoms = new LLKA_Atom[N];
    merged.nAtoms = N;

    size_t idx = 0;
    for (auto &s : toMerge) {
        const size_t NA = s->nAtoms;
        std::copy_n(s->atoms, NA, merged.atoms + idx);
        idx += NA;

        if constexpr (Policy == StrucutreMergePolicy::MOVE) {
            delete [] s->atoms;
            s->atoms = nullptr;
        }
    }

    assert(idx == N);

    return merged;
}

template <size_t NStructures>
inline
auto mergeStructuresMove(const std::array<LLKA_Structure *, NStructures> &toMerge)
{
    return mergeStructures<NStructures, StrucutreMergePolicy::MOVE>(toMerge);
}

template <size_t NStructures>
inline
auto mergeStructuresShallow(const std::array<const LLKA_Structure *, NStructures> &toMerge)
{
    return mergeStructures<NStructures, StrucutreMergePolicy::SHALLOW>(toMerge);
}

auto splitByAltIds(const LLKA_Structure &stru) -> std::tuple<std::vector<LLKA_Structure>, std::vector<char>>;

} // namespace LLKAInternal

#endif // _LLKA_STRUCTURE_UTIL_H

