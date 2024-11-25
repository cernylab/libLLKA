// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _LLKA_EXTRACT_H
#define _LLKA_EXTRACT_H

#include <llka_structure.h>

#include <cstdint>
#include <string>
#include <vector>

namespace LLKAInternal {

class AtomToExtract {
public:
    AtomToExtract() :
        modelNum{-1},
        seqId{-1},
        asymId{""},
        name{""},
        altId{LLKA_NO_ALTID}
    {}

    AtomToExtract(int32_t modelNum, std::string asymId, int32_t seqId, std::string compId, std::string name = {}, char altId = LLKA_NO_ALTID) :
        modelNum{modelNum},
        seqId{seqId},
        asymId{std::move(asymId)},
        compId{std::move(compId)},
        name{std::move(name)},
        altId{altId}
    {}

    int32_t modelNum;
    int32_t seqId;
    std::string asymId;
    std::string compId;
    std::string name;
    char altId;
};

auto extractAtoms(const LLKA_Structure *stru, const std::vector<AtomToExtract> &toExtract, LLKA_Structure *extracted) -> LLKA_RetCode;
auto getAllMatchingAtoms(const LLKA_Structure *stru, const AtomToExtract &ate) -> std::vector<LLKA_Atom *>;

} // namespace LLKAInternal

#endif // _LLKA_EXTRACT_H
