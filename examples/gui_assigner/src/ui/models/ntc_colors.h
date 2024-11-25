/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _NTC_COLORS_H
#define _NTC_COLORS_H

#include <llka_cpp.h>

#include <QColor>

#include <map>

inline constexpr auto COLOR_A = QColor{0xFF, 0xC1, 0xC1};
inline constexpr auto COLOR_B = QColor{0xC8, 0xCF, 0xFF};
inline constexpr auto COLOR_BII = QColor{0x00, 0x59, 0xDA};
inline constexpr auto COLOR_miB = QColor{0x3B, 0xE8, 0xFB};
inline constexpr auto COLOR_Z = QColor{0x01, 0xF6, 0x0E};
inline constexpr auto COLOR_IC = QColor{0xFA, 0x5C, 0xFB};
inline constexpr auto COLOR_OPN = QColor{0xE9, 0x00, 0x00};
inline constexpr auto COLOR_SYN = QColor{0xFF, 0xFF, 0x01};
inline constexpr auto COLOR_N = QColor{0xF2, 0xF2, 0xF2};

inline const std::map<LLKA_NtC, std::pair<QColor, QColor>> NTC_COLORS{
    { LLKA_INVALID_NTC, { COLOR_N, COLOR_N } },
    { LLKA_AA00, { COLOR_A, COLOR_A } },
    { LLKA_AA02, { COLOR_A, COLOR_A } },
    { LLKA_AA03, { COLOR_A, COLOR_A } },
    { LLKA_AA04, { COLOR_A, COLOR_A } },
    { LLKA_AA08, { COLOR_A, COLOR_A } },
    { LLKA_AA09, { COLOR_A, COLOR_A } },
    { LLKA_AA01, { COLOR_A, COLOR_A } },
    { LLKA_AA05, { COLOR_A, COLOR_A } },
    { LLKA_AA06, { COLOR_A, COLOR_A } },
    { LLKA_AA10, { COLOR_A, COLOR_A } },
    { LLKA_AA11, { COLOR_A, COLOR_A } },
    { LLKA_AA07, { COLOR_A, COLOR_A } },
    { LLKA_AA12, { COLOR_A, COLOR_A } },
    { LLKA_AA13, { COLOR_A, COLOR_A } },
    { LLKA_AB01, { COLOR_A, COLOR_B } },
    { LLKA_AB02, { COLOR_A, COLOR_B } },
    { LLKA_AB03, { COLOR_A, COLOR_B } },
    { LLKA_AB04, { COLOR_A, COLOR_B } },
    { LLKA_AB05, { COLOR_A, COLOR_B } },
    { LLKA_BA01, { COLOR_B, COLOR_A } },
    { LLKA_BA05, { COLOR_B, COLOR_A } },
    { LLKA_BA09, { COLOR_B, COLOR_A } },
    { LLKA_BA08, { COLOR_BII, COLOR_A } },
    { LLKA_BA10, { COLOR_B, COLOR_A } },
    { LLKA_BA13, { COLOR_BII, COLOR_A } },
    { LLKA_BA16, { COLOR_BII, COLOR_A } },
    { LLKA_BA17, { COLOR_BII, COLOR_A } },
    { LLKA_BB00, { COLOR_B, COLOR_B } },
    { LLKA_BB01, { COLOR_B, COLOR_B } },
    { LLKA_BB17, { COLOR_B, COLOR_B } },
    { LLKA_BB02, { COLOR_B, COLOR_B } },
    { LLKA_BB03, { COLOR_B, COLOR_B } },
    { LLKA_BB11, { COLOR_B, COLOR_B } },
    { LLKA_BB16, { COLOR_B, COLOR_B } },
    { LLKA_BB04, { COLOR_B, COLOR_BII } },
    { LLKA_BB05, { COLOR_B, COLOR_BII } },
    { LLKA_BB07, { COLOR_BII, COLOR_BII } },
    { LLKA_BB08, { COLOR_BII, COLOR_BII } },
    { LLKA_BB10, { COLOR_miB, COLOR_miB } },
    { LLKA_BB12, { COLOR_miB, COLOR_miB } },
    { LLKA_BB13, { COLOR_miB, COLOR_miB } },
    { LLKA_BB14, { COLOR_miB, COLOR_miB } },
    { LLKA_BB15, { COLOR_miB, COLOR_miB } },
    { LLKA_BB20, { COLOR_miB, COLOR_miB } },
    { LLKA_IC01, { COLOR_IC, COLOR_IC } },
    { LLKA_IC02, { COLOR_IC, COLOR_IC } },
    { LLKA_IC03, { COLOR_IC, COLOR_IC } },
    { LLKA_IC04, { COLOR_IC, COLOR_IC } },
    { LLKA_IC05, { COLOR_IC, COLOR_IC } },
    { LLKA_IC06, { COLOR_IC, COLOR_IC } },
    { LLKA_IC07, { COLOR_IC, COLOR_IC } },
    { LLKA_OP01, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP02, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP03, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP04, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP05, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP06, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP07, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP08, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP09, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP10, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP11, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP12, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP13, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP14, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP15, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP16, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP17, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP18, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP19, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP20, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP21, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP22, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP23, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP24, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP25, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP26, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP27, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP28, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP29, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP30, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP31, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OPS1, { COLOR_OPN, COLOR_OPN } },
    { LLKA_OP1S, { COLOR_OPN, COLOR_SYN } },
    { LLKA_AAS1, { COLOR_SYN, COLOR_A } },
    { LLKA_AB1S, { COLOR_A, COLOR_SYN } },
    { LLKA_AB2S, { COLOR_A, COLOR_SYN } },
    { LLKA_BB1S, { COLOR_B, COLOR_SYN } },
    { LLKA_BB2S, { COLOR_B, COLOR_SYN } },
    { LLKA_BBS1, { COLOR_SYN, COLOR_B } },
    { LLKA_ZZ01, { COLOR_Z, COLOR_Z } },
    { LLKA_ZZ02, { COLOR_Z, COLOR_Z } },
    { LLKA_ZZ1S, { COLOR_Z, COLOR_SYN } },
    { LLKA_ZZ2S, { COLOR_Z, COLOR_SYN } },
    { LLKA_ZZS1, { COLOR_SYN, COLOR_Z } },
    { LLKA_ZZS2, { COLOR_SYN, COLOR_Z } }
};

#endif // _NTC_COLORS_H
