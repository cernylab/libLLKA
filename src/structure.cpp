// vim: set sw=4 ts=4 sts=4 expandtab :

#include <llka_structure.h>
#include <llka_nucleotide.h>

#include "extend.h"
#include "nucleotide.hpp"
#include "structure_util.hpp"
#include "util/elementaries.h"

#include <cassert>
#include <cstring>

static
auto atomMatchesCriteria(
    const LLKA_Atom *atom,
    const char *label_atom_id,
    const char *label_comp_id,
    const char *label_asym_id,
    int32_t label_seq_id,
    char label_alt_id,
    const char *pdbx_PDB_ins_code,
    int32_t pdbx_PDB_model_num
) {
    assert(pdbx_PDB_ins_code != nullptr);

    if (label_atom_id && std::strcmp(label_atom_id, atom->label_atom_id))
        return false;

    if (label_comp_id && std::strcmp(label_comp_id, atom->label_comp_id))
        return false;

    if (label_asym_id && std::strcmp(label_asym_id, atom->label_asym_id))
        return false;

    if (label_seq_id >= 0 && atom->label_seq_id != label_seq_id)
        return false;

    if (label_alt_id != LLKA_NO_ALTID && atom->label_alt_id != label_alt_id)
        return false;

    if (std::strcmp(pdbx_PDB_ins_code, LLKA_NO_INSCODE) && std::strcmp(pdbx_PDB_ins_code, atom->pdbx_PDB_ins_code))
        return false;

    return atom->pdbx_PDB_model_num == pdbx_PDB_model_num;
}

void LLKA_CC LLKA_appendAtom(const LLKA_Atom *atom, LLKA_Structure *stru)
{
    auto newAtoms = new LLKA_Atom[stru->nAtoms + 1];
    if (stru->nAtoms > 0)
        std::memcpy(newAtoms, stru->atoms, stru->nAtoms * sizeof(LLKA_Atom));

    LLKA_duplicateAtom(atom, &newAtoms[stru->nAtoms]);
    delete[] stru->atoms;
    stru->atoms = newAtoms;
    stru->nAtoms++;
}

void LLKA_CC LLKA_appendAtomFromParams(
    uint32_t id,
    const char *type_symbol,
    const char *label_atom_id,
    const char *label_entity_id,
    const char *label_comp_id,
    const char *label_asym_id,
    const char *auth_atom_id,
    const char *auth_comp_id,
    const char* auth_asym_id,
    int32_t label_seq_id,
    char label_alt_id,
    int32_t auth_seq_id,
    const char *pdbx_PDB_ins_code,
    int32_t pdbx_PDB_model_num,
    const LLKA_Point *coords,
    LLKA_Structure *stru
)
{
    auto newAtoms = new LLKA_Atom[stru->nAtoms + 1];
    if (stru->nAtoms > 0)
        std::memcpy(newAtoms, stru->atoms, stru->nAtoms * sizeof(LLKA_Atom));

    LLKA_initAtom(
        id,
        type_symbol,
        label_atom_id, label_entity_id, label_comp_id, label_asym_id,
        auth_atom_id, auth_comp_id, auth_asym_id,
        label_seq_id, label_alt_id, auth_seq_id,
        pdbx_PDB_ins_code,
        pdbx_PDB_model_num,
        coords,
        &newAtoms[stru->nAtoms]
    );

    delete[] stru->atoms;
    stru->atoms = newAtoms;
    stru->nAtoms++;
}

LLKA_API void LLKA_CC LLKA_appendStructure(const LLKA_Structure *appendend, LLKA_Structure *appendee)
{
    const size_t nNewAtoms = appendee->nAtoms + appendend->nAtoms;

    auto newAtoms = new LLKA_Atom[nNewAtoms];
    if (appendee->nAtoms > 0)
        std::memcpy(newAtoms, appendee->atoms, appendee->nAtoms * sizeof(LLKA_Atom));

    for (size_t idx = appendee->nAtoms; idx < nNewAtoms; idx++)
        LLKA_duplicateAtom(&appendend->atoms[idx - appendee->nAtoms], &newAtoms[idx]);

    delete [] appendee->atoms;
    appendee->atoms = newAtoms;
    appendee->nAtoms = nNewAtoms;
}

LLKA_RetCode LLKA_CC LLKA_baseKind(const char *compId, LLKA_BaseKind *kind)
{
    if (LLKAInternal::findBaseKind(compId, *kind))
        return LLKA_OK;

    return LLKA_E_INVALID_ARGUMENT;
}

LLKA_Bool LLKA_CC LLKA_compareAtoms(const LLKA_Atom *a, const LLKA_Atom *b, LLKA_Bool ignoreId)
{
    bool matches =
        (a->id == b->id || ignoreId) &&
        a->label_seq_id == b->label_seq_id &&
        a->auth_seq_id == b->auth_seq_id &&
        a->label_alt_id == b->label_alt_id &&
        a->pdbx_PDB_model_num == b->pdbx_PDB_model_num;
    if (!matches)
        return LLKA_FALSE;

    // SOMETHING TO CONSIDER: These exact FP comparisons are perhaps not ideal.
    matches =
        a->coords.x == b->coords.x &&
        a->coords.y == b->coords.y &&
        a->coords.z == b->coords.z;
    if (!matches)
        return LLKA_FALSE;

    if (std::strcmp(a->pdbx_PDB_ins_code, b->pdbx_PDB_ins_code))
        return LLKA_FALSE;
    if (std::strcmp(a->type_symbol, b->type_symbol))
        return LLKA_FALSE;
    if (std::strcmp(a->label_atom_id, b->label_atom_id))
        return LLKA_FALSE;
    if (std::strcmp(a->label_entity_id, b->label_entity_id))
        return LLKA_FALSE;
    if (std::strcmp(a->label_comp_id, b->label_comp_id))
        return LLKA_FALSE;
    if (std::strcmp(a->label_asym_id, b->label_asym_id))
        return LLKA_FALSE;

    // We are not comparing "auth" identifiers because
    // two atoms with the same "label" but different "auth"
    // identifiers would be really weird and probably indicative
    // of malformed data.

    return LLKA_TRUE;
}

void LLKA_CC LLKA_destroyAlternatePositions(const LLKA_AlternatePositions *alts)
{
    delete [] alts->positions;
}

void LLKA_CC LLKA_destroyAtom(const LLKA_Atom *atom)
{
    LLKAInternal::destroyString(atom->type_symbol);
    LLKAInternal::destroyString(atom->label_atom_id);
    LLKAInternal::destroyString(atom->label_entity_id);
    LLKAInternal::destroyString(atom->label_comp_id);
    LLKAInternal::destroyString(atom->label_asym_id),
    LLKAInternal::destroyString(atom->auth_atom_id);
    LLKAInternal::destroyString(atom->auth_comp_id);
    LLKAInternal::destroyString(atom->auth_asym_id);
    LLKAInternal::destroyString(atom->pdbx_PDB_ins_code);
}

void LLKA_CC LLKA_destroyStructure(const LLKA_Structure *stru)
{
    for (size_t idx = 0; idx < stru->nAtoms; idx++)
        LLKA_destroyAtom(&stru->atoms[idx]);

    delete [] stru->atoms;
}

void LLKA_CC LLKA_destroyStructureView(const LLKA_StructureView *view)
{
    delete [] view->atoms;
}

void LLKA_CC LLKA_destroyStructures(const LLKA_Structures *strus)
{
    for (size_t idx = 0; idx < strus->nStrus; idx++)
        LLKA_destroyStructure(&strus->strus[idx]);

    delete [] strus->strus;
}

void LLKA_CC LLKA_duplicateAtom(const LLKA_Atom *source, LLKA_Atom *target)
{
    target->id = source->id;
    target->type_symbol = LLKAInternal::duplicateString(source->type_symbol);
    target->label_entity_id = LLKAInternal::duplicateString(source->label_entity_id);
    target->label_atom_id = LLKAInternal::duplicateString(source->label_atom_id);
    target->label_comp_id = LLKAInternal::duplicateString(source->label_comp_id);
    target->label_asym_id = LLKAInternal::duplicateString(source->label_asym_id);
    target->auth_atom_id = LLKAInternal::duplicateString(source->auth_atom_id);
    target->auth_comp_id = LLKAInternal::duplicateString(source->auth_comp_id);
    target->auth_asym_id = LLKAInternal::duplicateString(source->auth_asym_id);
    target->coords = source->coords;
    target->label_seq_id = source->label_seq_id;
    target->auth_seq_id = source->auth_seq_id;
    target->pdbx_PDB_ins_code = LLKAInternal::duplicateString(source->pdbx_PDB_ins_code);
    target->pdbx_PDB_model_num = source->pdbx_PDB_model_num;
    target->label_alt_id = source->label_alt_id;
}

void LLKA_CC LLKA_duplicateStructure(const LLKA_Structure *source, LLKA_Structure *target)
{
    target->atoms = nullptr;
    target->nAtoms = source->nAtoms;

    if (target->nAtoms  < 1)
        return;

    target->atoms = new LLKA_Atom[target->nAtoms];
    for (size_t idx = 0; idx < target->nAtoms; idx++)
        LLKA_duplicateAtom(&source->atoms[idx], &target->atoms[idx]);
}

LLKA_Atom LLKA_CC LLKA_emptyAtom(void)
{
    LLKA_Atom empty{
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        { 0, 0, 0 },
        0,
        0,
        0,
        0,
        LLKA_NO_INSCODE,
        LLKA_NO_ALTID
    };

    return empty;
}

LLKA_Atom * LLKA_findAtom(
    LLKA_Structure *stru,
    const char *label_atom_id,
    const char *label_comp_id,
    const char *label_asym_id,
    int32_t label_seq_id,
    char label_alt_id,
    const char *pdbx_PDB_ins_code,
    int32_t pdbx_PDB_model_num
) {
    for (size_t idx = 0; idx < stru->nAtoms; idx++) {
        LLKA_Atom *atom = &stru->atoms[idx];

        if (atomMatchesCriteria(
            atom,
            label_atom_id, label_comp_id, label_asym_id, label_seq_id, label_alt_id, pdbx_PDB_ins_code, pdbx_PDB_model_num
        ))
            return atom;
    }

    return nullptr;
}

const LLKA_Atom * LLKA_findAtomView(
    LLKA_Structure *view,
    const char *label_atom_id,
    const char *label_comp_id,
    const char *label_asym_id,
    int32_t label_seq_id,
    char label_alt_id,
    const char *pdbx_PDB_ins_code,
    int32_t pdbx_PDB_model_num
) {
    for (size_t idx = 0; idx < view->nAtoms; idx++) {
        const LLKA_Atom *atom = &view->atoms[idx];

        if (atomMatchesCriteria(
            atom,
            label_atom_id, label_comp_id, label_asym_id, label_seq_id, label_alt_id, pdbx_PDB_ins_code, pdbx_PDB_model_num
        ))
            return atom;
    }

    return nullptr;
}


void LLKA_CC LLKA_initAtom(
    uint32_t id,
    const char *type_symbol,
    const char *label_atom_id,
    const char *label_entity_id,
    const char *label_comp_id,
    const char *label_asym_id,
    const char *auth_atom_id,
    const char *auth_comp_id,
    const char *auth_asym_id,
    int32_t label_seq_id,
    char label_alt_id,
    int32_t auth_seq_id,
    const char *pdbx_PDB_ins_code,
    int32_t pdbx_PDB_model_num,
    const LLKA_Point *coords,
    LLKA_Atom *atom
)
{
    atom->id = id;
    atom->type_symbol = LLKAInternal::duplicateString(type_symbol);
    atom->label_atom_id = LLKAInternal::duplicateString(label_atom_id);
    atom->label_entity_id = LLKAInternal::duplicateString(label_entity_id);
    atom->label_comp_id = LLKAInternal::duplicateString(label_comp_id);
    atom->label_asym_id = LLKAInternal::duplicateString(label_asym_id);
    atom->auth_atom_id = auth_atom_id ? LLKAInternal::duplicateString(auth_atom_id) : LLKAInternal::duplicateString(label_atom_id);
    atom->auth_comp_id = auth_comp_id ? LLKAInternal::duplicateString(auth_comp_id) : LLKAInternal::duplicateString(label_comp_id);
    atom->auth_asym_id = auth_asym_id ? LLKAInternal::duplicateString(auth_asym_id) : LLKAInternal::duplicateString(label_asym_id);
    atom->coords = *coords;
    atom->label_seq_id = label_seq_id;
    atom->auth_seq_id = auth_seq_id;
    atom->pdbx_PDB_ins_code = LLKAInternal::duplicateString(pdbx_PDB_ins_code),
    atom->pdbx_PDB_model_num = pdbx_PDB_model_num,
    atom->label_alt_id = label_alt_id;
}

void LLKA_CC LLKA_initStructure(const LLKA_Atom *atoms, size_t nAtoms, LLKA_Structure *stru)
{
    if (nAtoms < 1) {
        stru->atoms = nullptr;
        stru->nAtoms = 0;
        return;
    }

    stru->atoms  = new LLKA_Atom[nAtoms];
    for (size_t idx = 0; idx < nAtoms; idx++)
        LLKA_duplicateAtom(&atoms[idx], &stru->atoms[idx]);

    stru->nAtoms = nAtoms;
}

LLKA_Atom LLKA_CC LLKA_makeAtom(
    uint32_t id,
    const char *type_symbol,
    const char *label_atom_id,
    const char *label_entity_id,
    const char *label_comp_id,
    const char *label_asym_id,
    const char *auth_atom_id,
    const char *auth_comp_id,
    const char *auth_asym_id,
    int32_t label_seq_id,
    char label_alt_id,
    int32_t auth_seq_id,
    const char *pdbx_PDB_ins_code,
    int32_t pdbx_PDB_model_num,
    const LLKA_Point *coords
)
{
    LLKA_Atom atom{};

    LLKA_initAtom(
        id,
        type_symbol,
        label_atom_id,
        label_entity_id,
        label_comp_id,
        label_asym_id,
        auth_atom_id,
        auth_comp_id,
        auth_asym_id,
        label_seq_id,
        label_alt_id,
        auth_seq_id,
        pdbx_PDB_ins_code,
        pdbx_PDB_model_num,
        coords,
        &atom
    );

    return atom;
}

LLKA_Structure LLKA_CC LLKA_makeStructure(const LLKA_Atom *atoms, size_t nAtoms)
{
    LLKA_Structure stru{
        nullptr,
        nAtoms
    };

    if (nAtoms < 1)
        return stru;

    stru.atoms  = new LLKA_Atom[stru.nAtoms];
    for (size_t idx = 0; idx < nAtoms; idx++)
        LLKA_duplicateAtom(&atoms[idx], &stru.atoms[idx]);

    return stru;
}

LLKA_Structure LLKA_CC LLKA_makeStructureFromPtrs(const LLKA_Atom *const *atoms, size_t nAtoms)
{
    LLKA_Structure stru{
        nullptr,
        nAtoms
    };

    if (nAtoms < 1)
        return stru;

    stru.atoms  = new LLKA_Atom[stru.nAtoms];
    for (size_t idx = 0; idx < nAtoms; idx++)
        LLKA_duplicateAtom(atoms[idx], &stru.atoms[idx]);

    return stru;
}

void LLKA_CC LLKA_removeAtomById(uint32_t id, LLKA_Structure *stru)
{
    size_t idx = 0;
    for (; idx < stru->nAtoms; idx++) {
        if (stru->atoms[idx].id == id)
            break;
    }
    if (idx == stru->nAtoms)
        return;

    auto newAtoms = new LLKA_Atom[stru->nAtoms - 1];
    if (idx > 0)
        std::memcpy(newAtoms, stru->atoms, idx * sizeof(LLKA_Atom));
    if (idx < stru->nAtoms - 1)
        std::memcpy(newAtoms + idx, stru->atoms + idx + 1, (stru->nAtoms - idx - 1) * sizeof(LLKA_Atom));

    LLKA_destroyAtom(&stru->atoms[idx]);
    delete [] stru->atoms;
    stru->atoms = newAtoms;
    stru->nAtoms--;
}

LLKA_API LLKA_Structures LLKA_CC LLKA_splitByAltIds(const LLKA_Structure *stru, LLKA_AlternatePositions *alts)
{
    auto [ strus, altIds ] = LLKAInternal::splitByAltIds(*stru);
    const auto N = strus.size();

    LLKA_Structures retStrus {
        new LLKA_Structure[N],
        N
    };
    std::copy_n(strus.cbegin(), N, retStrus.strus);

    if (!altIds.empty()) {
        char *pos = new char[N];
        std::copy_n(altIds.cbegin(), N, pos);

        alts->positions = pos;
        alts->nPositions = N;
    } else {
        alts->positions = nullptr;
        alts->nPositions = 0;
    }

    return retStrus;
}

LLKA_RetCode LLKA_CC LLKA_splitStructureToDinucleotideSteps(const LLKA_Structure *stru, LLKA_Structures *steps)
{
    std::vector<LLKA_Structure> allSteps;

    size_t idx = 0;
    while (idx < stru->nAtoms) {
        const auto &atom = stru->atoms[idx];

        if (!LLKAInternal::isNucleotideCompound(atom.label_comp_id)) {
            idx++;
            continue;
        }

        // Notice how we instruct the extendToResidue() function to only look forward.
        // This is becasue the loop above looks for the first viable atom. Going backwards from
        // what has to be a first atom in a residue makes no sense.
        // Besides being a small performance tweak at also allows us to deal with structures
        // with microheterogetnity. Microheterogenity means that one residue in a structure
        // can correspond to two different compounds. Look at 6r93 as an example.
        // The loop above would correctly skip an MHET compoud if it has an unknown base.
        // However, a bidirectional residue extension would put those skipped portion back
        // into the returned structure and we do not want that.
        // Note that this logic can still fail if an MHET structure intertwines known and unknown
        // bases instead of putting the in contiguous blocks. To deal with that we would have to
        // filter out any unknown bases after extending to residue.
        // That would be rather slow and it probably is not necessary for any currently deposited structures.

        LLKA_Structure firstNucl = LLKAInternal::extendToResidue({ stru, idx }, atom.pdbx_PDB_model_num, atom.label_asym_id, atom.label_seq_id, LLKA_NO_ALTID, LLKAInternal::EXT_DIR_FWD);
        idx += firstNucl.nAtoms;

        if (idx >= stru->nAtoms) {
            LLKA_destroyStructure(&firstNucl);
            break;
        }

        const auto &atom2 = stru->atoms[idx];
        if (!LLKAInternal::isSameChain(atom, atom2)) {
            LLKA_destroyStructure(&firstNucl);
            continue;
        }

        LLKA_Structure secondNucl = LLKAInternal::extendToResidue({ stru, idx }, atom2.pdbx_PDB_model_num, atom2.label_asym_id, atom2.label_seq_id, LLKA_NO_ALTID, LLKAInternal::EXT_DIR_FWD);

        auto _steps = LLKAInternal::dinucleotideToSteps(firstNucl, secondNucl);
        allSteps.insert(allSteps.end(), _steps.cbegin(), _steps.cend());

        LLKA_destroyStructure(&firstNucl);
        LLKA_destroyStructure(&secondNucl);
    }

    const size_t N = allSteps.size();
    steps->strus = new LLKA_Structure[N];
    steps->nStrus = N;

    std::copy_n(allSteps.cbegin(), N, steps->strus);

    return LLKA_OK;
}
