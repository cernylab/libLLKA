/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "assignment.h"

#include <algorithm>
#include <cctype>
#include <tuple>

static
auto nucleotideName(const std::string &base, int32_t seqId, const std::string &insCode, char altId)
{
    std::string altIdStr = altId != LLKA_NO_ALTID ? "." + std::string{altId} : "";
    std::string insCodeStr = insCode != LLKA_NO_INSCODE ? "." + insCode : "";
    return base + altIdStr + "_" + std::to_string(seqId) + insCodeStr;
}

static
auto stepAltIds(const LLKA::Structure &lStep)
{
    auto it = lStep.cbegin();

    auto firstSeqId = it->label_seq_id;
    char altIdFirst = it->label_alt_id;
    char altIdSecond = '\0';

    // Look for alt id
    while (it != lStep.cend() && altIdFirst == '\0' && it->label_seq_id == firstSeqId) {
        altIdFirst = it->label_alt_id;
        it++;
    }

    // Scroll to the second nucleotide
    while (it != lStep.cend() && it->label_seq_id == firstSeqId) it++;
    if (it == lStep.cend())
        return std::make_tuple(altIdFirst, altIdSecond); // No second residue? Odd, but pretend that we didn't see that.

    while (it != lStep.cend() && altIdSecond == '\0') {
        altIdSecond = it->label_alt_id;
        it++;
    }

    return std::make_tuple(altIdFirst, altIdSecond);
}

LoadedStructure::LoadedStructure(const bool hasMultipleModels, std::string id, std::vector<Step> steps) noexcept :
    hasMultipleModels{hasMultipleModels},
    id{std::move(id)},
    steps{std::move(steps)}
{
}

Step::Step(const bool multipleModels, const std::string &id, const LLKA::Structure &structure) noexcept :
    state{NotYetClassified},
    structure{structure}
{
    name = DNATCOName(multipleModels, id, structure);
}

Step::Step(const bool multipleModels, const std::string &id, const LLKA::Structure &structure, const LLKA::ClassifiedStep &classification) noexcept :
    state{Classified},
    structure{structure},
    classification{classification}
{
    name = DNATCOName(multipleModels, id, structure);
}

auto Step::DNATCOName(const bool multipleModels, const std::string &id, const LLKA::Structure &structure) -> std::string
{
    auto [ altIdFirst, altIdSecond ] = stepAltIds(structure);

    auto chain = structure.front().auth_asym_id;
    auto baseFirst = structure.front().auth_comp_id;
    auto seqIdFirst = structure.front().auth_seq_id;
    auto baseSecond = structure.back().auth_comp_id;
    auto seqIdSecond = structure.back().auth_seq_id;
    auto insCodeFirst = structure.front().pdbx_PDB_ins_code;
    auto insCodeSecond = structure.back().pdbx_PDB_ins_code;
    auto modelPrefix = multipleModels ? "-m" + std::to_string(structure.front().pdbx_PDB_model_num) : "";

    return id + modelPrefix + "_" + chain + "_" + nucleotideName(baseFirst, seqIdFirst, insCodeFirst, altIdFirst) + "_" + nucleotideName(baseSecond, seqIdSecond, insCodeSecond, altIdSecond);
}

auto loadStructure(const std::filesystem::path &path) -> LoadedStructure
{
    auto resStru = LLKA::cifToStructure(path);

    if (!resStru.isSuccess()) {
        const auto &fail = resStru.failure();

        throw AssignmentError{fail.tRet, fail.error};
    }

    const auto &importedStru = resStru.success();

    auto resSteps = LLKA::splitStructureToDinucleotideSteps(importedStru.structure);
    if (!resSteps.isSuccess()) {
        const auto &fail = resSteps.failure();

        throw AssignmentError{fail, ""};
    }

    const auto &lSteps = resSteps.success();
    bool hasMultipleModels = false;
    if (!lSteps.empty()) {
        auto modelNo = lSteps.front().front().pdbx_PDB_model_num;

        for (size_t idx = 1; idx < lSteps.size(); idx++) {
            if (lSteps[idx].front().pdbx_PDB_model_num != modelNo) {
                hasMultipleModels = true;
                break;
            }
        }
    }

    std::string lwrId{importedStru.id};
    std::transform(lwrId.begin(), lwrId.end(), lwrId.begin(), ::tolower);

    std::vector<Step> steps;
    for (const auto &ls : lSteps)
        steps.emplace_back(hasMultipleModels, lwrId, ls);

    return { hasMultipleModels, lwrId, steps };
}
