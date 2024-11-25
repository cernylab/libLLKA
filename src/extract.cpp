/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <extract.hpp>

#include <cassert>
#include <memory>

namespace LLKAInternal {

auto extractAtoms(const LLKA_Structure *stru, const std::vector<AtomToExtract> &toExtract, LLKA_Structure *extracted) -> LLKA_RetCode
{
    const size_t nAtoms = toExtract.size();
    auto view = makeStructureView(nAtoms);

    auto tRet = extractAtoms(stru, toExtract, &view);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructureView(&view);
        return tRet;
    }

    structureViewToStructure(view, *extracted);
    LLKA_destroyStructureView(&view);

    return LLKA_OK;
}

auto getAllMatchingAtoms(const LLKA_Structure *stru, const AtomToExtract &ate) -> std::vector<LLKA_Atom *>
{
    std::vector<LLKA_Atom *> matching{};
    matching.reserve(stru->nAtoms / 2 + 1);

    for (size_t idx = 0; idx < stru->nAtoms; idx++) {
        const auto atom = &stru->atoms[idx];
        bool isCandidate = (
            atom->pdbx_PDB_model_num == ate.modelNum &&
            atom->label_seq_id == ate.seqId &&
            (std::strcmp(atom->label_asym_id, ate.asymId.c_str()) == 0)
        );
        if (!ate.compId.empty())
            isCandidate &= std::strcmp(atom->label_comp_id, ate.compId.c_str()) == 0;
        if (!ate.name.empty())
            isCandidate &= std::strcmp(atom->label_atom_id, ate.name.c_str()) == 0;

        if (!isCandidate)
            continue;

        if (ate.altId == LLKA_NO_ALTID)
            matching.push_back(atom);
        else if (atom->label_alt_id == LLKA_NO_ALTID || atom->label_alt_id == ate.altId)
            matching.push_back(atom);
    }

    return matching;
}

} // namespace LLKAInternal
