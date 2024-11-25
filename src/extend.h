/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_EXTEND_H
#define _LLKA_EXTEND_H

#include <llka_structure.h>

namespace LLKAInternal {

enum ExtendDirection {
    EXT_DIR_FWD  = (1 << 0),
    EXT_DIR_BACK = (1 << 1)
};

class ExtendSource {
public:
    const LLKA_Structure *stru;
    size_t fromIndex;
};

LLKA_Structure extendToResidue(const ExtendSource &src, const int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id, char label_alt_id, int direction = EXT_DIR_FWD | EXT_DIR_BACK);

} // namespace LLKAInternal

#endif // _LLKA_EXTEND_H
