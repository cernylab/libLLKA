/* vim: set sw=4 ts=4 sts=4 expandtab : */


#include "extend.h"


#include "structure_util.hpp"


#include <cassert>
#include <cstring>
#include <vector>



namespace LLKAInternal {

/*  
    TODO: @sromdani Check this if is correct
    src                     ->      Source atoms
    pbdx_PBD_model_number   ->      Number of model in Protein Data Bank
    label_asym_id           ->      Identification asymetric value
    label_seq_id            ->      Identification sequence residue
    label_alt_id            ->      Identification alternative identifier
    direction               ->      Direction of extension
*/
LLKA_Structure extendToResidue(const ExtendSource &src, const int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id, char label_alt_id, int direction)
{

    assert(direction > 0);
    std::vector<LLKA_Atom *> atoms{ &src.stru->atoms[src.fromIndex] };

    if (direction & EXT_DIR_BACK) {


        for (size_t idx = src.fromIndex - 1; idx < src.stru->nAtoms; idx--) {
            const auto atom = &src.stru->atoms[idx];
            if (!isSameResidue(*atom, pdbx_PDB_model_num, label_asym_id, label_seq_id, label_alt_id))
                break;
            atoms.push_back(atom);
        }
    }


    if (direction & EXT_DIR_FWD) {
        for (size_t idx = src.fromIndex + 1; idx < src.stru->nAtoms; idx++) {
            const auto atom = &src.stru->atoms[idx];
            if (!isSameResidue(*atom, pdbx_PDB_model_num, label_asym_id, label_seq_id, label_alt_id))
                break;
            atoms.push_back(atom);
        }
    }
    return LLKA_makeStructureFromPtrs(atoms.data(), atoms.size());
}

} // namespace LLKAInternal
