/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _NON_STANDARD_RESIDUES_H
#define _NON_STANDARD_RESIDUES_H

#include "ntc_bones.hpp"

#include <map>

/*
 * IMPORTANT NOTE:
 * Contents of this file are mostly auto-generated. The logic used to find what atoms to use instead of
 * the standard atoms is rather crude and it is expected that some atoms will be wrong.
 */

namespace LLKABones {

inline const BoneFirstResidue _05ABoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue _05ABoneSecondResidue = { "C3", "N2", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase _05ABase = { "N1", "C6" };
inline const Quad _05ABaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _05ABone = {
    .firstResidue = _05ABoneFirstResidue,
    .secondResidue = _05ABoneSecondResidue,
    .base = _05ABase,

    .baseQuad = _05ABaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { '0', '5', 'A', '\0' }
};

inline const BoneFirstResidue _05HBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue _05HBoneSecondResidue = { "C71", "N5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase _05HBase = { "N1", "C2" };
inline const Quad _05HBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone _05HBone = {
    .firstResidue = _05HBoneFirstResidue,
    .secondResidue = _05HBoneSecondResidue,
    .base = _05HBase,

    .baseQuad = _05HBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { '0', '5', 'H', '\0' }
};

inline const BoneFirstResidue _05KBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue _05KBoneSecondResidue = { "C71", "N5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase _05KBase = { "N1", "C2" };
inline const Quad _05KBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone _05KBone = {
    .firstResidue = _05KBoneFirstResidue,
    .secondResidue = _05KBoneSecondResidue,
    .base = _05KBase,

    .baseQuad = _05KBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { '0', '5', 'K', '\0' }
};

inline const BoneBase _0AUBase = { "N1", "C6" };
inline const Quad _0AUBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _0AUBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _0AUBase,
    .baseQuad = _0AUBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '0', 'A', 'U', '\0' }
};

inline const BoneBase _0UBase = { "N1", "C6" };
inline const Quad _0UBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _0UBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _0UBase,
    .baseQuad = _0UBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '0', 'U', '\0' }
};

inline const BoneBase _0U1Base = { "N1", "C6" };
inline const Quad _0U1BaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _0U1Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _0U1Base,
    .baseQuad = _0U1BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '0', 'U', '1', '\0' }
};

inline const BoneFirstResidue _128BoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue _128BoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase _128Base = { "N9", "C4" };
inline const Quad _128BaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone _128Bone = {
    .firstResidue = _128BoneFirstResidue,
    .secondResidue = _128BoneSecondResidue,
    .base = _128Base,

    .baseQuad = _128BaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { '1', '2', '8', '\0' }
};

inline const BoneBase _1RNBase = { "N1", "C6" };
inline const Quad _1RNBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _1RNBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _1RNBase,
    .baseQuad = _1RNBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '1', 'R', 'N', '\0' }
};

inline const BoneBase _1TLBase = { "N1", "C6" };
inline const Quad _1TLBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _1TLBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _1TLBase,
    .baseQuad = _1TLBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '1', 'T', 'L', '\0' }
};

inline const BoneBase _1W5Base = { "C1", "C2" };
inline const Quad _1W5BaseQuad = { "O4'", "C1'", "C1", "C2" };
inline constexpr Bone _1W5Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _1W5Base,
    .baseQuad = _1W5BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '1', 'W', '5', '\0' }
};

inline const BoneBase _2DFBase = { "N1", "O2" };
inline const Quad _2DFBaseQuad = { "O4'", "C1'", "N1", "O2" };
inline constexpr Bone _2DFBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _2DFBase,
    .baseQuad = _2DFBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '2', 'D', 'F', '\0' }
};

inline const BoneBase _2LABase = { "N9", "C8" };
inline const Quad _2LABaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone _2LABone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _2LABase,
    .baseQuad = _2LABaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '2', 'L', 'A', '\0' }
};

inline const BoneBase _2OMBase = { "N1", "C6" };
inline const Quad _2OMBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _2OMBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _2OMBase,
    .baseQuad = _2OMBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '2', 'O', 'M', '\0' }
};

inline const BoneBase _3MUBase = { "N1", "C6" };
inline const Quad _3MUBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _3MUBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _3MUBase,
    .baseQuad = _3MUBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '3', 'M', 'U', '\0' }
};

inline const BoneBase _3TDBase = { "C5", "C4" };
inline const Quad _3TDBaseQuad = { "O4'", "C1'", "C5", "C4" };
inline constexpr Bone _3TDBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _3TDBase,
    .baseQuad = _3TDBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '3', 'T', 'D', '\0' }
};

inline const BoneBase _4ENBase = { "N8", "N9" };
inline const Quad _4ENBaseQuad = { "O4'", "C1'", "N8", "N9" };
inline constexpr Bone _4ENBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _4ENBase,
    .baseQuad = _4ENBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '4', 'E', 'N', '\0' }
};

inline const BoneBase _4MFBase = { "N1", "C7A" };
inline const Quad _4MFBaseQuad = { "O4'", "C1'", "N1", "C7A" };
inline constexpr Bone _4MFBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _4MFBase,
    .baseQuad = _4MFBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '4', 'M', 'F', '\0' }
};

inline const BoneBase _56BBase = { "N9", "C8" };
inline const Quad _56BBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone _56BBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _56BBase,
    .baseQuad = _56BBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '5', '6', 'B', '\0' }
};

inline const BoneFirstResidue _5FABoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue _5FABoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase _5FABase = { "N9", "C4" };
inline const Quad _5FABaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone _5FABone = {
    .firstResidue = _5FABoneFirstResidue,
    .secondResidue = _5FABoneSecondResidue,
    .base = _5FABase,

    .baseQuad = _5FABaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { '5', 'F', 'A', '\0' }
};

inline const BoneBase _5NCBase = { "N1", "C6" };
inline const Quad _5NCBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _5NCBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _5NCBase,
    .baseQuad = _5NCBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '5', 'N', 'C', '\0' }
};

inline const BoneFirstResidue _5UABoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue _5UABoneSecondResidue = { "C6'", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase _5UABase = { "N9", "C4" };
inline const Quad _5UABaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone _5UABone = {
    .firstResidue = _5UABoneFirstResidue,
    .secondResidue = _5UABoneSecondResidue,
    .base = _5UABase,

    .baseQuad = _5UABaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { '5', 'U', 'A', '\0' }
};

inline const BoneBase _64PBase = { "N1", "C6" };
inline const Quad _64PBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _64PBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _64PBase,
    .baseQuad = _64PBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '6', '4', 'P', '\0' }
};

inline const BoneBase _64TBase = { "N1", "O2" };
inline const Quad _64TBaseQuad = { "O4'", "C1'", "N1", "O2" };
inline constexpr Bone _64TBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _64TBase,
    .baseQuad = _64TBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '6', '4', 'T', '\0' }
};

inline const BoneBase _6FKBase = { "N9", "C8" };
inline const Quad _6FKBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone _6FKBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _6FKBase,
    .baseQuad = _6FKBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '6', 'F', 'K', '\0' }
};

inline const BoneFirstResidue _6FMBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue _6FMBoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase _6FMBase = { "N9", "C4" };
inline const Quad _6FMBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone _6FMBone = {
    .firstResidue = _6FMBoneFirstResidue,
    .secondResidue = _6FMBoneSecondResidue,
    .base = _6FMBase,

    .baseQuad = _6FMBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { '6', 'F', 'M', '\0' }
};

inline const BoneBase _6FUBase = { "N1", "C6" };
inline const Quad _6FUBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _6FUBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _6FUBase,
    .baseQuad = _6FUBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '6', 'F', 'U', '\0' }
};

inline const BoneBase _6HABase = { "N9", "C4" };
inline const Quad _6HABaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone _6HABone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _6HABase,
    .baseQuad = _6HABaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '6', 'H', 'A', '\0' }
};

inline const BoneBase _6HCBase = { "N1", "C6" };
inline const Quad _6HCBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _6HCBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _6HCBase,
    .baseQuad = _6HCBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '6', 'H', 'C', '\0' }
};

inline const BoneBase _6HGBase = { "N9", "C4" };
inline const Quad _6HGBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone _6HGBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _6HGBase,
    .baseQuad = _6HGBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '6', 'H', 'G', '\0' }
};

inline const BoneBase _6HTBase = { "N1", "C6" };
inline const Quad _6HTBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _6HTBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _6HTBase,
    .baseQuad = _6HTBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '6', 'H', 'T', '\0' }
};

inline const BoneBase _6MIBase = { "N1M", "C8A" };
inline const Quad _6MIBaseQuad = { "O4'", "C1'", "N1M", "C8A" };
inline constexpr Bone _6MIBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _6MIBase,
    .baseQuad = _6MIBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '6', 'M', 'I', '\0' }
};

inline const BoneBase _7ATBase = { "N9", "N8" };
inline const Quad _7ATBaseQuad = { "O4'", "C1'", "N9", "N8" };
inline constexpr Bone _7ATBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _7ATBase,
    .baseQuad = _7ATBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '7', 'A', 'T', '\0' }
};

inline const BoneBase _7MGBase = { "N9", "C8" };
inline const Quad _7MGBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone _7MGBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _7MGBase,
    .baseQuad = _7MGBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '7', 'M', 'G', '\0' }
};

inline const BoneBase _7SNBase = { "N9", "C8" };
inline const Quad _7SNBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone _7SNBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _7SNBase,
    .baseQuad = _7SNBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '7', 'S', 'N', '\0' }
};

inline const BoneBase _8AABase = { "N9", "C8" };
inline const Quad _8AABaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone _8AABone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _8AABase,
    .baseQuad = _8AABaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '8', 'A', 'A', '\0' }
};

inline const BoneBase _8AGBase = { "N9", "C8" };
inline const Quad _8AGBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone _8AGBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _8AGBase,
    .baseQuad = _8AGBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '8', 'A', 'G', '\0' }
};

inline const BoneBase _8AZBase = { "N9", "N8" };
inline const Quad _8AZBaseQuad = { "O4'", "C1'", "N9", "N8" };
inline constexpr Bone _8AZBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _8AZBase,
    .baseQuad = _8AZBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '8', 'A', 'Z', '\0' }
};

inline const BoneBase _8MGBase = { "N9", "C8" };
inline const Quad _8MGBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone _8MGBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _8MGBase,
    .baseQuad = _8MGBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '8', 'M', 'G', '\0' }
};

inline const BoneBase _8PYBase = { "N9", "C8" };
inline const Quad _8PYBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone _8PYBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _8PYBase,
    .baseQuad = _8PYBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '8', 'P', 'Y', '\0' }
};

inline const BoneBase _8ROBase = { "N1", "C6" };
inline const Quad _8ROBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _8ROBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _8ROBase,
    .baseQuad = _8ROBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '8', 'R', 'O', '\0' }
};

inline const BoneBase _8YNBase = { "C1", "C2" };
inline const Quad _8YNBaseQuad = { "O4'", "C1'", "C1", "C2" };
inline constexpr Bone _8YNBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _8YNBase,
    .baseQuad = _8YNBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '8', 'Y', 'N', '\0' }
};

inline const BoneBase _93DBase = { "C1", "C6" };
inline const Quad _93DBaseQuad = { "O4'", "C1'", "C1", "C6" };
inline constexpr Bone _93DBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _93DBase,
    .baseQuad = _93DBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '9', '3', 'D', '\0' }
};

inline const BoneBase _9V9Base = { "N1", "C6" };
inline const Quad _9V9BaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone _9V9Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = _9V9Base,
    .baseQuad = _9V9BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { '9', 'V', '9', '\0' }
};

inline const BoneFirstResidue A1PBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue A1PBoneSecondResidue = { "O2P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase A1PBase = { "N9", "C4" };
inline const Quad A1PBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone A1PBone = {
    .firstResidue = A1PBoneFirstResidue,
    .secondResidue = A1PBoneSecondResidue,
    .base = A1PBase,

    .baseQuad = A1PBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'A', '1', 'P', '\0' }
};

inline const BoneFirstResidue A3PBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue A3PBoneSecondResidue = { "P2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase A3PBase = { "N9", "C4" };
inline const Quad A3PBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone A3PBone = {
    .firstResidue = A3PBoneFirstResidue,
    .secondResidue = A3PBoneSecondResidue,
    .base = A3PBase,

    .baseQuad = A3PBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'A', '3', 'P', '\0' }
};

inline const BoneBase A6CBase = { "N1", "C6" };
inline const Quad A6CBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone A6CBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = A6CBase,
    .baseQuad = A6CBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'A', '6', 'C', '\0' }
};

inline const BoneBase A6UBase = { "N1", "C6" };
inline const Quad A6UBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone A6UBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = A6UBase,
    .baseQuad = A6UBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'A', '6', 'U', '\0' }
};

inline const BoneBase A7CBase = { "N9", "N8" };
inline const Quad A7CBaseQuad = { "O4'", "C1'", "N9", "N8" };
inline constexpr Bone A7CBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = A7CBase,
    .baseQuad = A7CBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'A', '7', 'C', '\0' }
};

inline const BoneBase A7EBase = { "N9", "N8" };
inline const Quad A7EBaseQuad = { "O4'", "C1'", "N9", "N8" };
inline constexpr Bone A7EBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = A7EBase,
    .baseQuad = A7EBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'A', '7', 'E', '\0' }
};

inline const BoneFirstResidue AD2BoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue AD2BoneSecondResidue = { "P1", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase AD2Base = { "N9", "C4" };
inline const Quad AD2BaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone AD2Bone = {
    .firstResidue = AD2BoneFirstResidue,
    .secondResidue = AD2BoneSecondResidue,
    .base = AD2Base,

    .baseQuad = AD2BaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'A', 'D', '2', '\0' }
};

inline const BoneFirstResidue ADXBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue ADXBoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase ADXBase = { "N9", "C4" };
inline const Quad ADXBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone ADXBone = {
    .firstResidue = ADXBoneFirstResidue,
    .secondResidue = ADXBoneSecondResidue,
    .base = ADXBase,

    .baseQuad = ADXBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'A', 'D', 'X', '\0' }
};

inline const BoneBase B8HBase = { "C5", "C4" };
inline const Quad B8HBaseQuad = { "O4'", "C1'", "C5", "C4" };
inline constexpr Bone B8HBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = B8HBase,
    .baseQuad = B8HBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'B', '8', 'H', '\0' }
};

inline const BoneBase B8KBase = { "N9", "C8" };
inline const Quad B8KBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone B8KBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = B8KBase,
    .baseQuad = B8KBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'B', '8', 'K', '\0' }
};

inline const BoneBase B8NBase = { "C5", "C4" };
inline const Quad B8NBaseQuad = { "O4'", "C1'", "C5", "C4" };
inline constexpr Bone B8NBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = B8NBase,
    .baseQuad = B8NBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'B', '8', 'N', '\0' }
};

inline const BoneBase B8QBase = { "N1", "C6" };
inline const Quad B8QBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone B8QBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = B8QBase,
    .baseQuad = B8QBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'B', '8', 'Q', '\0' }
};

inline const BoneBase B9HBase = { "N1", "C6" };
inline const Quad B9HBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone B9HBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = B9HBase,
    .baseQuad = B9HBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'B', '9', 'H', '\0' }
};

inline const BoneBase BGHBase = { "N9", "C8" };
inline const Quad BGHBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone BGHBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = BGHBase,
    .baseQuad = BGHBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'B', 'G', 'H', '\0' }
};

inline const BoneBase BGMBase = { "N9", "C8" };
inline const Quad BGMBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone BGMBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = BGMBase,
    .baseQuad = BGMBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'B', 'G', 'M', '\0' }
};

inline const BoneBase BMNBase = { "C1", "C6" };
inline const Quad BMNBaseQuad = { "O4'", "C1'", "C1", "C6" };
inline constexpr Bone BMNBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = BMNBase,
    .baseQuad = BMNBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'B', 'M', 'N', '\0' }
};

inline const BoneBase BMQBase = { "N1", "C6" };
inline const Quad BMQBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone BMQBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = BMQBase,
    .baseQuad = BMQBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'B', 'M', 'Q', '\0' }
};

inline const BoneBase C36Base = { "N1", "C6" };
inline const Quad C36BaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone C36Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = C36Base,
    .baseQuad = C36BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'C', '3', '6', '\0' }
};

inline const BoneBase C4JBase = { "C5", "C4" };
inline const Quad C4JBaseQuad = { "O4'", "C1'", "C5", "C4" };
inline constexpr Bone C4JBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = C4JBase,
    .baseQuad = C4JBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'C', '4', 'J', '\0' }
};

inline const BoneBase CARBase = { "N1", "C6" };
inline const Quad CARBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone CARBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = CARBase,
    .baseQuad = CARBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'C', 'A', 'R', '\0' }
};

inline const BoneBase CDWBase = { "N1", "C6" };
inline const Quad CDWBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone CDWBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = CDWBase,
    .baseQuad = CDWBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'C', 'D', 'W', '\0' }
};

inline const BoneBase CGYBase = { "C1", "C2" };
inline const Quad CGYBaseQuad = { "O4'", "C1'", "C1", "C2" };
inline constexpr Bone CGYBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = CGYBase,
    .baseQuad = CGYBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'C', 'G', 'Y', '\0' }
};

inline const BoneBase CJ1Base = { "N9", "C8" };
inline const Quad CJ1BaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone CJ1Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = CJ1Base,
    .baseQuad = CJ1BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'C', 'J', '1', '\0' }
};

inline const BoneFirstResidue CSFBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue CSFBoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase CSFBase = { "N1", "C2" };
inline const Quad CSFBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone CSFBone = {
    .firstResidue = CSFBoneFirstResidue,
    .secondResidue = CSFBoneSecondResidue,
    .base = CSFBase,

    .baseQuad = CSFBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'C', 'S', 'F', '\0' }
};

inline const BoneBase CSMBase = { "N1", "C6" };
inline const Quad CSMBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone CSMBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = CSMBase,
    .baseQuad = CSMBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'C', 'S', 'M', '\0' }
};

inline const BoneBase CTGBase = { "N1", "C2" };
inline const Quad CTGBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone CTGBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = CTGBase,
    .baseQuad = CTGBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'C', 'T', 'G', '\0' }
};

inline const BoneBase CVCBase = { "N9", "C4" };
inline const Quad CVCBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone CVCBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = CVCBase,
    .baseQuad = CVCBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'C', 'V', 'C', '\0' }
};

inline const BoneBase D3Base = { "N1A", "C5A" };
inline const Quad D3BaseQuad = { "O4'", "C1'", "N1A", "C5A" };
inline constexpr Bone D3Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = D3Base,
    .baseQuad = D3BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'D', '3', '\0' }
};

inline const BoneBase D33Base = { "N1", "C4" };
inline const Quad D33BaseQuad = { "O4'", "C1'", "N1", "C4" };
inline constexpr Bone D33Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = D33Base,
    .baseQuad = D33BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'D', '3', '3', '\0' }
};

inline const BoneBase D3NBase = { "N1", "C6" };
inline const Quad D3NBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone D3NBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = D3NBase,
    .baseQuad = D3NBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'D', '3', 'N', '\0' }
};

inline const BoneBase DFTBase = { "C1", "C6" };
inline const Quad DFTBaseQuad = { "O4'", "C1'", "C1", "C6" };
inline constexpr Bone DFTBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = DFTBase,
    .baseQuad = DFTBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'D', 'F', 'T', '\0' }
};

inline const BoneFirstResidue DGIBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue DGIBoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase DGIBase = { "N9", "C4" };
inline const Quad DGIBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone DGIBone = {
    .firstResidue = DGIBoneFirstResidue,
    .secondResidue = DGIBoneSecondResidue,
    .base = DGIBase,

    .baseQuad = DGIBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'D', 'G', 'I', '\0' }
};

inline const BoneBase DPYBase = { "C1", "C2" };
inline const Quad DPYBaseQuad = { "O4'", "C1'", "C1", "C2" };
inline constexpr Bone DPYBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = DPYBase,
    .baseQuad = DPYBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'D', 'P', 'Y', '\0' }
};

inline const BoneBase DRPBase = { "C1", "C2" };
inline const Quad DRPBaseQuad = { "O4'", "C1'", "C1", "C2" };
inline constexpr Bone DRPBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = DRPBase,
    .baseQuad = DRPBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'D', 'R', 'P', '\0' }
};

inline const BoneBase DZBase = { "C1", "C6" };
inline const Quad DZBaseQuad = { "O4'", "C1'", "C1", "C6" };
inline constexpr Bone DZBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = DZBase,
    .baseQuad = DZBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'D', 'Z', '\0' }
};

inline const BoneBase E3CBase = { "N1", "C6" };
inline const Quad E3CBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone E3CBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = E3CBase,
    .baseQuad = E3CBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'E', '3', 'C', '\0' }
};

inline const BoneBase E7GBase = { "N9", "C8" };
inline const Quad E7GBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone E7GBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = E7GBase,
    .baseQuad = E7GBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'E', '7', 'G', '\0' }
};

inline const BoneBase EDCBase = { "N1", "C6" };
inline const Quad EDCBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone EDCBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = EDCBase,
    .baseQuad = EDCBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'E', 'D', 'C', '\0' }
};

inline const BoneFirstResidue ENPBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue ENPBoneSecondResidue = { "P1", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase ENPBase = { "N9", "C4" };
inline const Quad ENPBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone ENPBone = {
    .firstResidue = ENPBoneFirstResidue,
    .secondResidue = ENPBoneSecondResidue,
    .base = ENPBase,

    .baseQuad = ENPBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'E', 'N', 'P', '\0' }
};

inline const BoneBase EW3Base = { "N1", "C6" };
inline const Quad EW3BaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone EW3Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = EW3Base,
    .baseQuad = EW3BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'E', 'W', '3', '\0' }
};

inline const BoneBase F2TBase = { "N1", "C6" };
inline const Quad F2TBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone F2TBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = F2TBase,
    .baseQuad = F2TBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'F', '2', 'T', '\0' }
};

inline const BoneBase F3HBase = { "N1", "C6" };
inline const Quad F3HBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone F3HBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = F3HBase,
    .baseQuad = F3HBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'F', '3', 'H', '\0' }
};

inline const BoneFirstResidue FAGBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue FAGBoneSecondResidue = { "O3P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase FAGBase = { "N9", "C4" };
inline const Quad FAGBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone FAGBone = {
    .firstResidue = FAGBoneFirstResidue,
    .secondResidue = FAGBoneSecondResidue,
    .base = FAGBase,

    .baseQuad = FAGBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'F', 'A', 'G', '\0' }
};

inline const BoneBase FFDBase = { "C6", "C5" };
inline const Quad FFDBaseQuad = { "O4'", "C1'", "C6", "C5" };
inline constexpr Bone FFDBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = FFDBase,
    .baseQuad = FFDBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'F', 'F', 'D', '\0' }
};

inline const BoneBase FHUBase = { "C5", "F5" };
inline const Quad FHUBaseQuad = { "O4'", "C1'", "C5", "F5" };
inline constexpr Bone FHUBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = FHUBase,
    .baseQuad = FHUBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'F', 'H', 'U', '\0' }
};

inline const BoneFirstResidue G4PBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue G4PBoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase G4PBase = { "N9", "C4" };
inline const Quad G4PBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone G4PBone = {
    .firstResidue = G4PBoneFirstResidue,
    .secondResidue = G4PBoneSecondResidue,
    .base = G4PBase,

    .baseQuad = G4PBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'G', '4', 'P', '\0' }
};

inline const BoneFirstResidue GMXBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue GMXBoneSecondResidue = { "OP3", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase GMXBase = { "N9", "C4" };
inline const Quad GMXBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone GMXBone = {
    .firstResidue = GMXBoneFirstResidue,
    .secondResidue = GMXBoneSecondResidue,
    .base = GMXBase,

    .baseQuad = GMXBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'G', 'M', 'X', '\0' }
};

inline const BoneBase GN7Base = { "N7", "C4" };
inline const Quad GN7BaseQuad = { "O4'", "C1'", "N7", "C4" };
inline constexpr Bone GN7Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = GN7Base,
    .baseQuad = GN7BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'G', 'N', '7', '\0' }
};

inline const BoneBase I2TBase = { "C5", "C4" };
inline const Quad I2TBaseQuad = { "O4'", "C1'", "C5", "C4" };
inline constexpr Bone I2TBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = I2TBase,
    .baseQuad = I2TBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'I', '2', 'T', '\0' }
};

inline const BoneBase ICBase = { "N1", "C6" };
inline const Quad ICBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone ICBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = ICBase,
    .baseQuad = ICBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'I', 'C', '\0' }
};

inline const BoneBase IMCBase = { "N1", "C6" };
inline const Quad IMCBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone IMCBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = IMCBase,
    .baseQuad = IMCBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'I', 'M', 'C', '\0' }
};

inline const BoneFirstResidue IOOBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C2", "O4'" };
inline const BoneSecondResidue IOOBoneSecondResidue = { "P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2", "O4'" };
inline const BoneBase IOOBase = { "N9", "C1'" };
inline const Quad IOOBaseQuad = { "O4'", "C1'", "N9", "C1'" };
inline constexpr Bone IOOBone = {
    .firstResidue = IOOBoneFirstResidue,
    .secondResidue = IOOBoneSecondResidue,
    .base = IOOBase,

    .baseQuad = IOOBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'I', 'O', 'O', '\0' }
};

inline const BoneBase IRNBase = { "N1", "C4" };
inline const Quad IRNBaseQuad = { "O4'", "C1'", "N1", "C4" };
inline constexpr Bone IRNBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = IRNBase,
    .baseQuad = IRNBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'I', 'R', 'N', '\0' }
};

inline const BoneBase J4TBase = { "N1", "C6" };
inline const Quad J4TBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone J4TBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = J4TBase,
    .baseQuad = J4TBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'J', '4', 'T', '\0' }
};

inline const BoneBase JLNBase = { "N1", "C4" };
inline const Quad JLNBaseQuad = { "O4'", "C1'", "N1", "C4" };
inline constexpr Bone JLNBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = JLNBase,
    .baseQuad = JLNBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'J', 'L', 'N', '\0' }
};

inline const BoneBase JMHBase = { "N1", "C6" };
inline const Quad JMHBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone JMHBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = JMHBase,
    .baseQuad = JMHBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'J', 'M', 'H', '\0' }
};

inline const BoneBase JSPBase = { "C1", "C2" };
inline const Quad JSPBaseQuad = { "O4'", "C1'", "C1", "C2" };
inline constexpr Bone JSPBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = JSPBase,
    .baseQuad = JSPBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'J', 'S', 'P', '\0' }
};

inline const BoneBase LCCBase = { "N1", "C6" };
inline const Quad LCCBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone LCCBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = LCCBase,
    .baseQuad = LCCBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'L', 'C', 'C', '\0' }
};

inline const BoneBase LHOBase = { "N1", "C4" };
inline const Quad LHOBaseQuad = { "O4'", "C1'", "N1", "C4" };
inline constexpr Bone LHOBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = LHOBase,
    .baseQuad = LHOBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'L', 'H', 'O', '\0' }
};

inline const BoneFirstResidue LMSBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue LMSBoneSecondResidue = { "S", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase LMSBase = { "N9", "C4" };
inline const Quad LMSBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone LMSBone = {
    .firstResidue = LMSBoneFirstResidue,
    .secondResidue = LMSBoneSecondResidue,
    .base = LMSBase,

    .baseQuad = LMSBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'L', 'M', 'S', '\0' }
};

inline const BoneBase MBZBase = { "N1", "C9" };
inline const Quad MBZBaseQuad = { "O4'", "C1'", "N1", "C9" };
inline constexpr Bone MBZBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MBZBase,
    .baseQuad = MBZBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'B', 'Z', '\0' }
};

inline const BoneBase MDJBase = { "N1", "C2" };
inline const Quad MDJBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone MDJBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MDJBase,
    .baseQuad = MDJBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'D', 'J', '\0' }
};

inline const BoneBase MDKBase = { "N1", "C2" };
inline const Quad MDKBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone MDKBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MDKBase,
    .baseQuad = MDKBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'D', 'K', '\0' }
};

inline const BoneBase MDQBase = { "N1", "C6" };
inline const Quad MDQBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone MDQBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MDQBase,
    .baseQuad = MDQBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'D', 'Q', '\0' }
};

inline const BoneBase MDUBase = { "N1", "C6" };
inline const Quad MDUBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone MDUBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MDUBase,
    .baseQuad = MDUBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'D', 'U', '\0' }
};

inline const BoneBase ME6Base = { "N1", "C6" };
inline const Quad ME6BaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone ME6Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = ME6Base,
    .baseQuad = ME6BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'E', '6', '\0' }
};

inline const BoneBase MFOBase = { "N9", "C8" };
inline const Quad MFOBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone MFOBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MFOBase,
    .baseQuad = MFOBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'F', 'O', '\0' }
};

inline const BoneBase MFTBase = { "N1", "C6" };
inline const Quad MFTBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone MFTBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MFTBase,
    .baseQuad = MFTBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'F', 'T', '\0' }
};

inline const BoneFirstResidue MGQBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue MGQBoneSecondResidue = { "PBE", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase MGQBase = { "N9", "C4" };
inline const Quad MGQBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone MGQBone = {
    .firstResidue = MGQBoneFirstResidue,
    .secondResidue = MGQBoneSecondResidue,
    .base = MGQBase,

    .baseQuad = MGQBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'M', 'G', 'Q', '\0' }
};

inline const BoneBase MHGBase = { "N9", "C8" };
inline const Quad MHGBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone MHGBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MHGBase,
    .baseQuad = MHGBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'H', 'G', '\0' }
};

inline const BoneBase MM7Base = { "C1", "C2" };
inline const Quad MM7BaseQuad = { "O4'", "C1'", "C1", "C2" };
inline constexpr Bone MM7Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MM7Base,
    .baseQuad = MM7BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'M', '7', '\0' }
};

inline const BoneFirstResidue MMTBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue MMTBoneSecondResidue = { "NP", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase MMTBase = { "N1", "C2" };
inline const Quad MMTBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone MMTBone = {
    .firstResidue = MMTBoneFirstResidue,
    .secondResidue = MMTBoneSecondResidue,
    .base = MMTBase,

    .baseQuad = MMTBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'M', 'M', 'T', '\0' }
};

inline const BoneBase MTRBase = { "C1", "C6" };
inline const Quad MTRBaseQuad = { "O4'", "C1'", "C1", "C6" };
inline constexpr Bone MTRBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MTRBase,
    .baseQuad = MTRBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'T', 'R', '\0' }
};

inline const BoneBase MTUBase = { "N9", "C8" };
inline const Quad MTUBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone MTUBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = MTUBase,
    .baseQuad = MTUBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'M', 'T', 'U', '\0' }
};

inline const BoneBase N5IBase = { "NE1", "CE2" };
inline const Quad N5IBaseQuad = { "O4'", "C1'", "NE1", "CE2" };
inline constexpr Bone N5IBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = N5IBase,
    .baseQuad = N5IBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'N', '5', 'I', '\0' }
};

inline const BoneBase N6GBase = { "N9", "C8" };
inline const Quad N6GBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone N6GBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = N6GBase,
    .baseQuad = N6GBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'N', '6', 'G', '\0' }
};

inline const BoneBase NCUBase = { "N1", "C6" };
inline const Quad NCUBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone NCUBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = NCUBase,
    .baseQuad = NCUBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'N', 'C', 'U', '\0' }
};

inline const BoneBase NF2Base = { "C1", "C2" };
inline const Quad NF2BaseQuad = { "O4'", "C1'", "C1", "C2" };
inline constexpr Bone NF2Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = NF2Base,
    .baseQuad = NF2BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'N', 'F', '2', '\0' }
};

inline const BoneBase NP3Base = { "N1", "C2" };
inline const Quad NP3BaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone NP3Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = NP3Base,
    .baseQuad = NP3BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'N', 'P', '3', '\0' }
};

inline const BoneBase NTTBase = { "N1", "C6" };
inline const Quad NTTBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone NTTBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = NTTBase,
    .baseQuad = NTTBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'N', 'T', 'T', '\0' }
};

inline const BoneFirstResidue OADBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue OADBoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase OADBase = { "N9", "C4" };
inline const Quad OADBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone OADBone = {
    .firstResidue = OADBoneFirstResidue,
    .secondResidue = OADBoneSecondResidue,
    .base = OADBase,

    .baseQuad = OADBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'O', 'A', 'D', '\0' }
};

inline const BoneFirstResidue OKQBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue OKQBoneSecondResidue = { "P1", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase OKQBase = { "N1", "C2" };
inline const Quad OKQBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone OKQBone = {
    .firstResidue = OKQBoneFirstResidue,
    .secondResidue = OKQBoneSecondResidue,
    .base = OKQBase,

    .baseQuad = OKQBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'O', 'K', 'Q', '\0' }
};

inline const BoneBase ONEBase = { "N1", "C6" };
inline const Quad ONEBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone ONEBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = ONEBase,
    .baseQuad = ONEBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'O', 'N', 'E', '\0' }
};

inline const BoneBase P2UBase = { "C5", "C4" };
inline const Quad P2UBaseQuad = { "O4'", "C1'", "C5", "C4" };
inline constexpr Bone P2UBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = P2UBase,
    .baseQuad = P2UBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'P', '2', 'U', '\0' }
};

inline const BoneBase P7GBase = { "N9", "C8" };
inline const Quad P7GBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone P7GBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = P7GBase,
    .baseQuad = P7GBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'P', '7', 'G', '\0' }
};

inline const BoneBase PBTBase = { "N1", "O2" };
inline const Quad PBTBaseQuad = { "O4'", "C1'", "N1", "O2" };
inline constexpr Bone PBTBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = PBTBase,
    .baseQuad = PBTBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'P', 'B', 'T', '\0' }
};

inline const BoneBase PSUBase = { "C5", "C4" };
inline const Quad PSUBaseQuad = { "O4'", "C1'", "C5", "C4" };
inline constexpr Bone PSUBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = PSUBase,
    .baseQuad = PSUBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'P', 'S', 'U', '\0' }
};

inline const BoneBase PYYBase = { "C1", "C2" };
inline const Quad PYYBaseQuad = { "O4'", "C1'", "C1", "C2" };
inline constexpr Bone PYYBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = PYYBase,
    .baseQuad = PYYBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'P', 'Y', 'Y', '\0' }
};

inline const BoneBase QBTBase = { "N1", "C6" };
inline const Quad QBTBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone QBTBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = QBTBase,
    .baseQuad = QBTBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'Q', 'B', 'T', '\0' }
};

inline const BoneBase QCKBase = { "N1", "C6" };
inline const Quad QCKBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone QCKBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = QCKBase,
    .baseQuad = QCKBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'Q', 'C', 'K', '\0' }
};

inline const BoneBase RCEBase = { "N1", "C6" };
inline const Quad RCEBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone RCEBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = RCEBase,
    .baseQuad = RCEBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'R', 'C', 'E', '\0' }
};

inline const BoneFirstResidue RIABoneFirstResidue = { "C5'", "C4'", "O1'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue RIABoneSecondResidue = { "P'", "O5'", "C5'", "C4'", "O1'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase RIABase = { "O2A", "O3A" };
inline const Quad RIABaseQuad = { "O4'", "C1'", "O2A", "O3A" };
inline constexpr Bone RIABone = {
    .firstResidue = RIABoneFirstResidue,
    .secondResidue = RIABoneSecondResidue,
    .base = RIABase,

    .baseQuad = RIABaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'R', 'I', 'A', '\0' }
};

inline const BoneFirstResidue RTPBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue RTPBoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase RTPBase = { "N1", "C5" };
inline const Quad RTPBaseQuad = { "O4'", "C1'", "N1", "C5" };
inline constexpr Bone RTPBone = {
    .firstResidue = RTPBoneFirstResidue,
    .secondResidue = RTPBoneSecondResidue,
    .base = RTPBase,

    .baseQuad = RTPBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'R', 'T', 'P', '\0' }
};

inline const BoneBase SAYBase = { "CAA", "CAF" };
inline const Quad SAYBaseQuad = { "O4'", "C1'", "CAA", "CAF" };
inline constexpr Bone SAYBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = SAYBase,
    .baseQuad = SAYBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'S', 'A', 'Y', '\0' }
};

inline const BoneBase T0TBase = { "C1", "C6" };
inline const Quad T0TBaseQuad = { "O4'", "C1'", "C1", "C6" };
inline constexpr Bone T0TBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = T0TBase,
    .baseQuad = T0TBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'T', '0', 'T', '\0' }
};

inline const BoneBase TDYBase = { "N1", "C6" };
inline const Quad TDYBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone TDYBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = TDYBase,
    .baseQuad = TDYBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'T', 'D', 'Y', '\0' }
};

inline const BoneFirstResidue TFFBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue TFFBoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase TFFBase = { "N1", "C2" };
inline const Quad TFFBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone TFFBone = {
    .firstResidue = TFFBoneFirstResidue,
    .secondResidue = TFFBoneSecondResidue,
    .base = TFFBase,

    .baseQuad = TFFBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'T', 'F', 'F', '\0' }
};

inline const BoneFirstResidue THPBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue THPBoneSecondResidue = { "P2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase THPBase = { "N1", "C2" };
inline const Quad THPBaseQuad = { "O4'", "C1'", "N1", "C2" };
inline constexpr Bone THPBone = {
    .firstResidue = THPBoneFirstResidue,
    .secondResidue = THPBoneSecondResidue,
    .base = THPBase,

    .baseQuad = THPBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'T', 'H', 'P', '\0' }
};

inline const BoneBase TLBBase = { "N1", "C6" };
inline const Quad TLBBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone TLBBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = TLBBase,
    .baseQuad = TLBBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'T', 'L', 'B', '\0' }
};

inline const BoneBase TLCBase = { "N1", "C6" };
inline const Quad TLCBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone TLCBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = TLCBase,
    .baseQuad = TLCBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'T', 'L', 'C', '\0' }
};

inline const BoneBase TLNBase = { "N1", "C6" };
inline const Quad TLNBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone TLNBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = TLNBase,
    .baseQuad = TLNBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'T', 'L', 'N', '\0' }
};

inline const BoneFirstResidue TPGBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue TPGBoneSecondResidue = { "PAT", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase TPGBase = { "N9", "C4" };
inline const Quad TPGBaseQuad = { "O4'", "C1'", "N9", "C4" };
inline constexpr Bone TPGBone = {
    .firstResidue = TPGBoneFirstResidue,
    .secondResidue = TPGBoneSecondResidue,
    .base = TPGBase,

    .baseQuad = TPGBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'T', 'P', 'G', '\0' }
};

inline const BoneBase U23Base = { "N1", "C6" };
inline const Quad U23BaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone U23Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = U23Base,
    .baseQuad = U23BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', '2', '3', '\0' }
};

inline const BoneFirstResidue UBDBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue UBDBoneSecondResidue = { "P1", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase UBDBase = { "N1", "C6" };
inline const Quad UBDBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone UBDBone = {
    .firstResidue = UBDBoneFirstResidue,
    .secondResidue = UBDBoneSecondResidue,
    .base = UBDBase,

    .baseQuad = UBDBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'U', 'B', 'D', '\0' }
};

inline const BoneFirstResidue UDPBoneFirstResidue = { "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneSecondResidue UDPBoneSecondResidue = { "PA", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C1'", "O4'" };
inline const BoneBase UDPBase = { "N1", "C6" };
inline const Quad UDPBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone UDPBone = {
    .firstResidue = UDPBoneFirstResidue,
    .secondResidue = UDPBoneSecondResidue,
    .base = UDPBase,

    .baseQuad = UDPBaseQuad,

    .hasStandardBackbone = false,
    .hasStandardBase = false,
    .name = { 'U', 'D', 'P', '\0' }
};

inline const BoneBase UF2Base = { "N1", "C6" };
inline const Quad UF2BaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone UF2Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = UF2Base,
    .baseQuad = UF2BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'F', '2', '\0' }
};

inline const BoneBase UMSBase = { "N1", "C6" };
inline const Quad UMSBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone UMSBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = UMSBase,
    .baseQuad = UMSBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'M', 'S', '\0' }
};

inline const BoneBase UMXBase = { "N1", "C6" };
inline const Quad UMXBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone UMXBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = UMXBase,
    .baseQuad = UMXBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'M', 'X', '\0' }
};

inline const BoneBase UOBBase = { "N1", "C6" };
inline const Quad UOBBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone UOBBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = UOBBase,
    .baseQuad = UOBBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'O', 'B', '\0' }
};

inline const BoneBase UR3Base = { "N1", "C6" };
inline const Quad UR3BaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone UR3Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = UR3Base,
    .baseQuad = UR3BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'R', '3', '\0' }
};

inline const BoneBase URXBase = { "N1", "C6" };
inline const Quad URXBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone URXBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = URXBase,
    .baseQuad = URXBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'R', 'X', '\0' }
};

inline const BoneBase US4Base = { "N1", "C6" };
inline const Quad US4BaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone US4Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = US4Base,
    .baseQuad = US4BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'S', '4', '\0' }
};

inline const BoneBase USMBase = { "N1", "C6" };
inline const Quad USMBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone USMBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = USMBase,
    .baseQuad = USMBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'S', 'M', '\0' }
};

inline const BoneBase UVXBase = { "N1", "C6" };
inline const Quad UVXBaseQuad = { "O4'", "C1'", "N1", "C6" };
inline constexpr Bone UVXBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = UVXBase,
    .baseQuad = UVXBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'V', 'X', '\0' }
};

inline const BoneBase UY1Base = { "C5", "C4" };
inline const Quad UY1BaseQuad = { "O4'", "C1'", "C5", "C4" };
inline constexpr Bone UY1Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = UY1Base,
    .baseQuad = UY1BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'U', 'Y', '1', '\0' }
};

inline const BoneBase WC7Base = { "N1", "C4" };
inline const Quad WC7BaseQuad = { "O4'", "C1'", "N1", "C4" };
inline constexpr Bone WC7Bone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = WC7Base,
    .baseQuad = WC7BaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'W', 'C', '7', '\0' }
};

inline const BoneBase XAEBase = { "N9", "C8" };
inline const Quad XAEBaseQuad = { "O4'", "C1'", "N9", "C8" };
inline constexpr Bone XAEBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = XAEBase,
    .baseQuad = XAEBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'X', 'A', 'E', '\0' }
};

inline const BoneBase XCSBase = { "C8", "C6" };
inline const Quad XCSBaseQuad = { "O4'", "C1'", "C8", "C6" };
inline constexpr Bone XCSBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = XCSBase,
    .baseQuad = XCSBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'X', 'C', 'S', '\0' }
};

inline const BoneBase XTYBase = { "C8", "C6" };
inline const Quad XTYBaseQuad = { "O4'", "C1'", "C8", "C6" };
inline constexpr Bone XTYBone = {
    .firstResidue = STANDARD_FIRST_RESIDUE,
    .secondResidue = STANDARD_SECOND_RESIDUE,
    .base = XTYBase,
    .baseQuad = XTYBaseQuad,

    .hasStandardBackbone = true,
    .hasStandardBase = false,
    .name = { 'X', 'T', 'Y', '\0' }
};


} // namespace LLKABones

namespace LLKAInternal {

inline const std::map<LLKABones::ANString, const LLKABones::Bone &> KNOWN_NON_STANDARD_RESIDUES = {
    { "05A", LLKABones::_05ABone },
    { "05H", LLKABones::_05HBone },
    { "05K", LLKABones::_05KBone },
    { "0AU", LLKABones::_0AUBone },
    {  "0U",  LLKABones::_0UBone },
    { "0U1", LLKABones::_0U1Bone },
    { "128", LLKABones::_128Bone },
    { "1RN", LLKABones::_1RNBone },
    { "1TL", LLKABones::_1TLBone },
    { "1W5", LLKABones::_1W5Bone },
    { "2DF", LLKABones::_2DFBone },
    { "2LA", LLKABones::_2LABone },
    { "2OM", LLKABones::_2OMBone },
    { "3MU", LLKABones::_3MUBone },
    { "3TD", LLKABones::_3TDBone },
    { "4EN", LLKABones::_4ENBone },
    { "4MF", LLKABones::_4MFBone },
    { "56B", LLKABones::_56BBone },
    { "5FA", LLKABones::_5FABone },
    { "5NC", LLKABones::_5NCBone },
    { "5UA", LLKABones::_5UABone },
    { "64P", LLKABones::_64PBone },
    { "64T", LLKABones::_64TBone },
    { "6FK", LLKABones::_6FKBone },
    { "6FM", LLKABones::_6FMBone },
    { "6FU", LLKABones::_6FUBone },
    { "6HA", LLKABones::_6HABone },
    { "6HC", LLKABones::_6HCBone },
    { "6HG", LLKABones::_6HGBone },
    { "6HT", LLKABones::_6HTBone },
    { "6MI", LLKABones::_6MIBone },
    { "7AT", LLKABones::_7ATBone },
    { "7MG", LLKABones::_7MGBone },
    { "7SN", LLKABones::_7SNBone },
    { "8AA", LLKABones::_8AABone },
    { "8AG", LLKABones::_8AGBone },
    { "8AZ", LLKABones::_8AZBone },
    { "8MG", LLKABones::_8MGBone },
    { "8PY", LLKABones::_8PYBone },
    { "8RO", LLKABones::_8ROBone },
    { "8YN", LLKABones::_8YNBone },
    { "93D", LLKABones::_93DBone },
    { "9V9", LLKABones::_9V9Bone },
    { "A1P", LLKABones::A1PBone },
    { "A3P", LLKABones::A3PBone },
    { "A6C", LLKABones::A6CBone },
    { "A6U", LLKABones::A6UBone },
    { "A7C", LLKABones::A7CBone },
    { "A7E", LLKABones::A7EBone },
    { "AD2", LLKABones::AD2Bone },
    { "ADX", LLKABones::ADXBone },
    { "B8H", LLKABones::B8HBone },
    { "B8K", LLKABones::B8KBone },
    { "B8N", LLKABones::B8NBone },
    { "B8Q", LLKABones::B8QBone },
    { "B9H", LLKABones::B9HBone },
    { "BGH", LLKABones::BGHBone },
    { "BGM", LLKABones::BGMBone },
    { "BMN", LLKABones::BMNBone },
    { "BMQ", LLKABones::BMQBone },
    { "C36", LLKABones::C36Bone },
    { "C4J", LLKABones::C4JBone },
    { "CAR", LLKABones::CARBone },
    { "CDW", LLKABones::CDWBone },
    { "CGY", LLKABones::CGYBone },
    { "CJ1", LLKABones::CJ1Bone },
    { "CSF", LLKABones::CSFBone },
    { "CSM", LLKABones::CSMBone },
    { "CTG", LLKABones::CTGBone },
    { "CVC", LLKABones::CVCBone },
    {  "D3",  LLKABones::D3Bone },
    { "D33", LLKABones::D33Bone },
    { "D3N", LLKABones::D3NBone },
    { "DFT", LLKABones::DFTBone },
    { "DGI", LLKABones::DGIBone },
    { "DPY", LLKABones::DPYBone },
    { "DRP", LLKABones::DRPBone },
    {  "DZ",  LLKABones::DZBone },
    { "E3C", LLKABones::E3CBone },
    { "E7G", LLKABones::E7GBone },
    { "EDC", LLKABones::EDCBone },
    { "ENP", LLKABones::ENPBone },
    { "EW3", LLKABones::EW3Bone },
    { "F2T", LLKABones::F2TBone },
    { "F3H", LLKABones::F3HBone },
    { "FAG", LLKABones::FAGBone },
    { "FFD", LLKABones::FFDBone },
    { "FHU", LLKABones::FHUBone },
    { "G4P", LLKABones::G4PBone },
    { "GMX", LLKABones::GMXBone },
    { "GN7", LLKABones::GN7Bone },
    { "I2T", LLKABones::I2TBone },
    {  "IC",  LLKABones::ICBone },
    { "IMC", LLKABones::IMCBone },
    { "IOO", LLKABones::IOOBone },
    { "IRN", LLKABones::IRNBone },
    { "J4T", LLKABones::J4TBone },
    { "JLN", LLKABones::JLNBone },
    { "JMH", LLKABones::JMHBone },
    { "JSP", LLKABones::JSPBone },
    { "LCC", LLKABones::LCCBone },
    { "LHO", LLKABones::LHOBone },
    { "LMS", LLKABones::LMSBone },
    { "MBZ", LLKABones::MBZBone },
    { "MDJ", LLKABones::MDJBone },
    { "MDK", LLKABones::MDKBone },
    { "MDQ", LLKABones::MDQBone },
    { "MDU", LLKABones::MDUBone },
    { "ME6", LLKABones::ME6Bone },
    { "MFO", LLKABones::MFOBone },
    { "MFT", LLKABones::MFTBone },
    { "MGQ", LLKABones::MGQBone },
    { "MHG", LLKABones::MHGBone },
    { "MM7", LLKABones::MM7Bone },
    { "MMT", LLKABones::MMTBone },
    { "MTR", LLKABones::MTRBone },
    { "MTU", LLKABones::MTUBone },
    { "N5I", LLKABones::N5IBone },
    { "N6G", LLKABones::N6GBone },
    { "NCU", LLKABones::NCUBone },
    { "NF2", LLKABones::NF2Bone },
    { "NP3", LLKABones::NP3Bone },
    { "NTT", LLKABones::NTTBone },
    { "OAD", LLKABones::OADBone },
    { "OKQ", LLKABones::OKQBone },
    { "ONE", LLKABones::ONEBone },
    { "P2U", LLKABones::P2UBone },
    { "P7G", LLKABones::P7GBone },
    { "PBT", LLKABones::PBTBone },
    { "PSU", LLKABones::PSUBone },
    { "PYY", LLKABones::PYYBone },
    { "QBT", LLKABones::QBTBone },
    { "QCK", LLKABones::QCKBone },
    { "RCE", LLKABones::RCEBone },
    { "RIA", LLKABones::RIABone },
    { "RTP", LLKABones::RTPBone },
    { "SAY", LLKABones::SAYBone },
    { "T0T", LLKABones::T0TBone },
    { "TDY", LLKABones::TDYBone },
    { "TFF", LLKABones::TFFBone },
    { "THP", LLKABones::THPBone },
    { "TLB", LLKABones::TLBBone },
    { "TLC", LLKABones::TLCBone },
    { "TLN", LLKABones::TLNBone },
    { "TPG", LLKABones::TPGBone },
    { "U23", LLKABones::U23Bone },
    { "UBD", LLKABones::UBDBone },
    { "UDP", LLKABones::UDPBone },
    { "UF2", LLKABones::UF2Bone },
    { "UMS", LLKABones::UMSBone },
    { "UMX", LLKABones::UMXBone },
    { "UOB", LLKABones::UOBBone },
    { "UR3", LLKABones::UR3Bone },
    { "URX", LLKABones::URXBone },
    { "US4", LLKABones::US4Bone },
    { "USM", LLKABones::USMBone },
    { "UVX", LLKABones::UVXBone },
    { "UY1", LLKABones::UY1Bone },
    { "WC7", LLKABones::WC7Bone },
    { "XAE", LLKABones::XAEBone },
    { "XCS", LLKABones::XCSBone },
    { "XTY", LLKABones::XTYBone },
};

} // namespace LLKAInternal

#endif // _NON_STANDARD_RESIDUES_H
