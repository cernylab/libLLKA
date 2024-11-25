/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _NTC_BONES_HPP
#define _NTC_BONES_HPP

#include <array>
#include <string>
#include <span>
#include <tuple>

namespace LLKABones {

// TODO: std::string wastes a lot of memory when we use it just for atom names.
// Use some more memory-efficient custom structure
using ANString = std::string;

using Quad = std::array<ANString, 4>;
using TaggedQuad = std::array<std::tuple<ANString, int>, 4>;
using BoneFirstResidue = std::array<ANString, 7>;
using BoneSecondResidue = std::array<ANString, 9>;
using BoneBase = std::array<ANString, 2>;
using BoneBackboneQuads = std::array<TaggedQuad, 7>;

/*
 * Collection of identifiers of atoms that are needed to calculate the necessary NtC metrics.
 * This information needs to be parametrized because some odd residues may have different bases
 * or, in a less common case, even different atoms on the backbone.
 */
class Bone {
public:
    // List of atom names when a given residue is used as first or second residue in a step

    // List of atom names when the residue is first in step
    const BoneFirstResidue &firstResidue;
    // List of atom names when the residue is second in step
    const BoneSecondResidue &secondResidue;

    // Base atoms to use for calculations of distances and pseudotorsions.
    const BoneBase &base;
    // Four atoms used to calculate chi torsions
    const Quad &baseQuad;

    const bool hasStandardBackbone;
    const bool hasStandardBase;
    const char name[4]; // Used only for non-standard bones
};

inline const BoneFirstResidue STANDARD_FIRST_RESIDUE = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue STANDARD_SECOND_RESIDUE = { "P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase STANDARD_PURINE = { "N9", "C4" };
inline const BoneBase STANDARD_PYRIMIDINE = { "N1", "C2" };
inline const BoneBackboneQuads STANDARD_BACKBONE_QUADS = {{
    {{  // Delta 1
        { "C5'", 1 },
        { "C4'", 1 },
        { "C3'", 1 },
        { "O3'", 1 }
    }},
    {{ // Epsilon 1
        { "C4'", 1 },
        { "C3'", 1 },
        { "O3'", 1 },
        { "P",   2 }
    }},
    {{  // Zeta 1
        { "C3'", 1 },
        { "O3'", 1 },
        { "P",   2 },
        { "O5'", 2 }
    }},
    {{  // Alpha 2
        { "O3'", 1 },
        { "P",   2 },
        { "O5'", 2 },
        { "C5'", 2 }
    }},
    {{  // Beta 2
        { "P",   2 },
        { "O5'", 2 },
        { "C5'", 2 },
        { "C4'", 2 }
    }},
    {{  // Gamma 2
        { "O5'", 2 },
        { "C5'", 2 },
        { "C4'", 2 },
        { "C3'", 2 }
    }},
    {{  // Delta 2
        { "C5'", 2 },
        { "C4'", 2 },
        { "C3'", 2 },
        { "O3'", 2 }
    }}
}};
inline const Quad STANDARD_PURINE_QUAD = { "O4'", "C1'", "N9", "C4" };
inline const Quad STANDARD_PYRIMIDINE_QUAD = { "O4'", "C1'", "N1", "C2" };

inline constexpr Bone StandardPurineBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = STANDARD_PURINE,
    .baseQuad = STANDARD_PURINE_QUAD,

    .hasStandardBackbone = true,
    .hasStandardBase = true,
    .name = { '\0' }
};

inline constexpr Bone StandardPyrimidineBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = STANDARD_PYRIMIDINE,
    .baseQuad = STANDARD_PYRIMIDINE_QUAD,

    .hasStandardBackbone = true,
    .hasStandardBase = true,
    .name = { '\0' }
};

inline constexpr size_t NUM_EXTENDED_ATOMS = STANDARD_FIRST_RESIDUE.size() + STANDARD_SECOND_RESIDUE.size() + 2; // Backbone atoms with two pairs of base atoms but without O4'
inline constexpr size_t NUM_PLAIN_ATOMS = STANDARD_FIRST_RESIDUE.size() + STANDARD_SECOND_RESIDUE.size() - 4;    // Backbone atoms without C1' and O4'
inline constexpr size_t NUM_TORSIONS_ATOMS = STANDARD_FIRST_RESIDUE.size() + STANDARD_SECOND_RESIDUE.size() + 4; // Everything and two pairs of base atoms
inline constexpr size_t NUM_FIRST_EXTENDED_BACKONE_ATOMS = STANDARD_FIRST_RESIDUE.size() - 1;
inline constexpr size_t NUM_SECOND_EXTENDED_BACKONE_ATOMS = STANDARD_SECOND_RESIDUE.size() - 1;

auto nonStandardBackboneQuads(const Bone &first, const Bone &second) -> const BoneBackboneQuads &;

inline
auto BACKBONE_QUADS(const Bone &first, const Bone &second) -> const BoneBackboneQuads &
{
    if (first.hasStandardBackbone && second.hasStandardBackbone) [[ likely ]]
        return STANDARD_BACKBONE_QUADS;

    return nonStandardBackboneQuads(first, second);
}

template <size_t N>
inline
constexpr auto EXTENDED_BACKBONE(const std::array<ANString, N> &atoms)
{
    static_assert(N > 1);
    return std::span(atoms.cbegin(), N - 1);
}

template <size_t N>
inline
constexpr auto PLAIN_BACKBONE(const std::array<ANString, N> &atoms)
{
    static_assert(N > 2);
    return std::span(atoms.cbegin(), N - 2);
}

template <size_t N>
inline
constexpr auto TORSIONS_BACKBONE(const std::array<ANString, N> &atoms)
{
    return atoms;
}

} // LLKABones

#endif // _NTC_BONES_HPP
