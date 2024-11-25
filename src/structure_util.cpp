/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "structure_util.hpp"
#include "util/elementaries.h"
#include "util/geometry.h"

#include <cassert>

namespace LLKAInternal {

static
auto getAtomByName(const LLKA_Structure &stru, const char *name) -> LLKA_Atom *
{
    for (size_t idx = 0; idx < stru.nAtoms; idx++) {
        auto atom = &stru.atoms[idx];
        if (std::strcmp(atom->auth_atom_id, name) == 0)
            return atom;
    }

    return nullptr;
}

auto dinucleotideToSteps(const LLKA_Structure &firstNucl, const LLKA_Structure &secondNucl) -> std::vector<LLKA_Structure>
{
    auto [ splittedFirst, altIdsFirst ] = splitByAltIds(firstNucl);
    auto [ splittedSecond, altIdsSecond ] = splitByAltIds(secondNucl);

    assert(splittedFirst.size() > 0);
    assert(splittedSecond.size() > 0);

    std::vector<LLKA_Structure> steps{};

    const double MAX_O3_P_DISTANCE_ANGSTROMS = 1.9;
    if (altIdsFirst.empty() && altIdsSecond.empty()) [[ likely ]] {
        auto atomO3 = getAtomByName(splittedFirst.front(), "O3'");
        auto atomP = getAtomByName(splittedSecond.front(), "P");
        if (atomO3 == nullptr || atomP == nullptr) [[ unlikely ]] {
            LLKA_destroyStructure(&splittedFirst[0]);
            LLKA_destroyStructure(&splittedSecond[0]);

            return {}; // Weird stucture, ignore it
        }

        // Check that the O3' -> P distance is reasonable
        auto dist = spatialDistance<double>(atomO3->coords, atomP->coords);
        if (dist <= MAX_O3_P_DISTANCE_ANGSTROMS) [[ likely ]]
            steps.emplace_back(mergeStructuresMove<2>({ &splittedFirst[0], &splittedSecond[0] }));
        else {
            LLKA_destroyStructure(&splittedFirst[0]);
            LLKA_destroyStructure(&splittedSecond[0]);
        }
    } else {
        steps.reserve(splittedFirst.size() * splittedSecond.size());

        for (auto &sf : splittedFirst) {
            auto atomO3 = getAtomByName(sf, "O3'");
            if (atomO3 == nullptr) [[ unlikely ]]
                continue; // Weird nucleotide, just skip it

            for (auto &ss : splittedSecond) {
                auto atomP = getAtomByName(ss, "P");

                if (atomP == nullptr) [[ unlikely ]]
                    continue; // Another weird nucleotide, skip it

                // Check that the O3' -> P distance is reasonable
                auto dist = spatialDistance<double>(atomO3->coords, atomP->coords);
                if (dist <= MAX_O3_P_DISTANCE_ANGSTROMS) [[ likely ]] {
                    // We need to duplicate the second structure because we may need it more than once
                    LLKA_Structure sfCopy = LLKA_makeStructure(sf.atoms, sf.nAtoms);
                    LLKA_Structure ssCopy = LLKA_makeStructure(ss.atoms, ss.nAtoms);
                    steps.emplace_back(mergeStructuresMove<2>({ &sfCopy, &ssCopy }));
                }
            }
        }

        for (auto &s : splittedFirst)
            LLKA_destroyStructure(&s);
        for (auto &s : splittedSecond)
            LLKA_destroyStructure(&s);
    }

    return steps;
}

auto splitByAltIds(const LLKA_Structure &stru) -> std::tuple<std::vector<LLKA_Structure>, std::vector<char>>
{
    std::vector<char> altIds;
    std::vector<LLKA_Structure> strus;

    for (size_t idx = 0; idx < stru.nAtoms; idx++) {
        const auto altId = stru.atoms[idx].label_alt_id;

        if (altId != LLKA_NO_ALTID && !contains(altIds, altId))
            altIds.push_back(altId);
    }

    if (altIds.empty()) {
        strus.resize(1);
        LLKA_duplicateStructure(&stru, &strus[0]);
    } else {
        const auto N = altIds.size();
        strus.resize(N);

        for (size_t sdx = 0; sdx < N; sdx++) {
            const auto altId = altIds[sdx];
            auto &splitStru = strus[sdx];
            splitStru.nAtoms = 0;

            for (size_t idx = 0; idx < stru.nAtoms; idx++) {
                const auto &atom = stru.atoms[idx];

                if (atom.label_alt_id == LLKA_NO_ALTID || atom.label_alt_id == altId)
                    LLKA_appendAtom(&atom, &splitStru);
            }
        }
    }

    assert(!strus.empty());

    return { std::move(strus), std::move(altIds) };
}

} // namespace LLKAInternal
