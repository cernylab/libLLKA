/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "ntc_bones.hpp"

#include <map>
#include <mutex>

namespace LLKABones {

static std::map<std::string, const BoneBackboneQuads *> quadsCache;
static std::mutex quadsCacheLock;

auto nonStandardBackboneQuads(const Bone &first, const Bone &second) -> const BoneBackboneQuads &
{
    std::string cacheKey = std::string{first.name} + std::string{second.name};

    quadsCacheLock.lock();
    auto it = quadsCache.find(cacheKey);
    if (it != quadsCache.cend()) {
        const BoneBackboneQuads *quads = it->second;
        quadsCacheLock.unlock();
        return *quads;
    }
    quadsCacheLock.unlock();

    // We have non-standard backbone. Create a custom list of backbone quads
    auto quads = new BoneBackboneQuads{};
    for (size_t ctr = 0; ctr < quads->size(); ctr++) {
        size_t a = (ctr + 2) % 7;
        size_t b = (ctr + 3) % 7;
        size_t c = (ctr + 4) % 7;
        size_t d = (ctr + 5) % 7;

        // Read atom names from the "secondResidue backbone" becase we have to walk the elementaries
        // step backbone. We need to start from what would be the C5' atom on a standard backbone
        // amd wrap around what would be O3',
        const auto &nameA = (ctr > 3 ? first : second).secondResidue[a];
        const auto &nameB = (ctr > 2 ? first : second).secondResidue[b];
        const auto &nameC = (ctr > 1 ? first : second).secondResidue[c];
        const auto &nameD = (ctr > 0 ? first : second).secondResidue[d];

        (*quads)[ctr][0] = { nameA, 1 + (ctr > 3) };
        (*quads)[ctr][1] = { nameB, 1 + (ctr > 2) };
        (*quads)[ctr][2] = { nameC, 1 + (ctr > 1) };
        (*quads)[ctr][3] = { nameD, 1 + (ctr > 0) };
    }

    quadsCacheLock.lock();
    it = quadsCache.find(cacheKey);
    if (it != quadsCache.cend()) {
        delete quads;
        auto _quads = it->second;
        quadsCacheLock.unlock();
        return *_quads;
    }

    quadsCache[cacheKey] = quads;
    quadsCacheLock.unlock();

    return *quads;
}

} // namespace LLKABones
