/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_STRUCTURE_H
#define _LLKA_STRUCTURE_H

#include "llka_main.h"

#define LLKA_NO_ALTID '\0'
#define LLKA_NO_INSCODE ""

typedef enum LLKA_BaseKind {
    LLKA_PURINE,
    LLKA_PYRIMIDINE,
    LLKA_NON_STANDARD_BASE
    ENUM_FORCE_INT32_SIZE(LLKA_BaseKind)
} LLKA_BaseKind;

typedef struct LLKA_AlternatePositions {
    const char *positions;
    size_t nPositions;
} LLKA_AlternatePositions;
LLKA_IS_POD(LLKA_AlternatePositions)

/*!
 * Defines an atom.
 *
 * With the exception of <tt>cootds</tt>, fields in the structure
 * correspond to the definition of <tt>atom_site</tt> category.
 * See "https://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Categories/atom_site.html" for details.
 */
typedef struct LLKA_Atom {
    const char *type_symbol;
    const char *label_atom_id;
    const char *label_entity_id;
    const char *label_comp_id;
    const char *label_asym_id;
    const char *auth_atom_id;
    const char *auth_comp_id;
    const char *auth_asym_id;
    LLKA_Point coords;
    uint32_t id;
    int32_t label_seq_id;
    int32_t auth_seq_id;
    int32_t pdbx_PDB_model_num;
    const char *pdbx_PDB_ins_code;
    char label_alt_id;
} LLKA_Atom;
LLKA_IS_POD(LLKA_Atom)

/*!
 * A structure
 *
 * Structure is a set of atoms.
 */
typedef struct LLKA_Structure {
    LLKA_Atom *atoms;    /*!< Array of \p LLKA_Atom s */
    size_t nAtoms;       /*!< Number of atoms in the set */
} LLKA_Structure;
LLKA_IS_POD(LLKA_Structure)

/*!
 * An immutable view on a section of a \p LLKA_Structure object.
 *
 * The source \p LLKA_Structure object must remain valid througout the lifetime of the view.
 */
typedef struct LLKA_StructureView {
    const LLKA_Atom **atoms;    /*!< Array of pointers to <tt>LLKA_Atom</tt>s contained in the source <tt>LLKA_Structure</tt>.
                                     If the source structure becomes invalid, these become dangling pointers. */
    size_t nAtoms;              /*!< Number of atoms in the view */
    size_t capacity;            /*!< Actual capacity of the atoms array.
                                     May be bigger that the number of viewed atoms. Used internally by the library . */
} LLKA_StructureView;
LLKA_IS_POD(LLKA_StructureView)

/*!
 * A set of structures
 */
typedef struct LLKA_Structures {
    LLKA_Structure *strus;    /*!< Array of \p LLKA_Structure s */
    size_t nStrus;            /*!< Number of structures in the set */
} LLKA_Structures;
LLKA_IS_POD(LLKA_Structures)

LLKA_BEGIN_API_FUNCTIONS

/*!
 * Appends atom to a structure.
 *
 * Contents of the input <tt>LLKA_Atom</tt> are deep-copied to the <tt>LLKA_Structure</tt> object.
 *
 * @param[in] atom LLKA_Atom to append.
 * @param[in,out] stru LLKA_Structure to append the atom to.
 */
LLKA_API void LLKA_CC LLKA_appendAtom(const LLKA_Atom *atom, LLKA_Structure *stru);

/*!
 * Appends atom to a structure.
 *
 * Instead of taking and copying and exiting atom, this function takes parameters
 * of the appended atom and constructs the atom in-place.
 *
 * @param[in] id id
 * @param[in] type_symbol type_symbol
 * @param[in] label_atom_id label_atom_id
 * @param[in] label_entity_id label_entity_id
 * @param[in] label_comp_id label_comp_id
 * @param[in] label_asym_id label_asym_id
 * @param[in] auth_atom_id auth_atom_id This may be <tt>NULL</tt> if no author atom id is provided
 * @param[in] auth_asym_id auth_comp_id This may be <tt>NULL</tt> if no author compound id is provided
 * @param[in] auth_asym_id auth_asym_id This may be <tt>NULL</tt> if no author assembly id is provided
 * @param[in] label_seq_id label_seq_id
 * @param[in] label_alt_id label_alt_id
 * @param[in] auth_seq_id auth_seq_id
 * @param[in] pdbx_PDB_model_num pdbx_PDB_model_num
 * @param[in] coords coords
 * @param[in,out] stru LLKA_Structure to append the atom to.
 */
LLKA_API void LLKA_CC LLKA_appendAtomFromParams(
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
    LLKA_Structure *stru
);

LLKA_API void LLKA_CC LLKA_appendStructure(const LLKA_Structure *appendend, LLKA_Structure *appendee);

/*!
 * Returns whether a residue is either purine or pyrimidine base
 *
 * @param[in] compId Name of the residue
 * @param[out] kind Determined kind of the residue.
 *
 * @retval LLKA_OK Passed residue name was recognized and determined to be either purine or pyrimidine base
 * @retval LLKA_E_INVALID_ARGUMENT Passed residue name was not recognized. Value of \p kind is undefined.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_baseKind(const char *compId, LLKA_BaseKind *kind);

/*!
 * Compares two \p LLKA_Atom s
 *
 * @param[in] a First atom.
 * @param[in] b Second atom.
 * @param[in] ignoreId If \p LLKA_TRUE, value of \p LLKA_Atom\::id is ignored in comparison
 *
 * @retval LLKA_TRUE The atoms' properties match
 * @retval LLKA_FALSE otherwise.
 */
LLKA_API LLKA_Bool LLKA_CC LLKA_compareAtoms(const LLKA_Atom *a, const LLKA_Atom *b, LLKA_Bool ignoreId);

/*!
 * Releases all resources claimed by \p LLKA_AlternatePositions
 *
 * @param[in] alts The object to destroy
 */
LLKA_API void LLKA_CC LLKA_destroyAlternatePositions(const LLKA_AlternatePositions *alts);

/*!
 * Releases all resources claimed by \p LLKA_Atom
 *
 * @param[in] atom The atom to destroy
 */
LLKA_API void LLKA_CC LLKA_destroyAtom(const LLKA_Atom *atom);

/*!
 * Releases all resources claimed by \p LLKA_Structure , including the atoms it is made of.
 *
 * @param[in] stru \p LLKA_Structure to release.
 */
LLKA_API void LLKA_CC LLKA_destroyStructure(const LLKA_Structure *stru);

/*!
 * Releases all resources claimed by \p LLKA_StructureView.
 * Note that viewer atoms are *not* owned by \p LLKA_StructureView
 *
 * @param[in] stru \p LLKA_Structure to release.
 */
LLKA_API void LLKA_CC LLKA_destroyStructureView(const LLKA_StructureView *stru);

/*!
 * Releases all resources claimed by a set of \p LLKA_Structure s, including the structures and atoms it is made of.
 * Note that this function can be safely used only on LLKA_Structures objects that were initialized by libLLKA.
 *
 * @param[in] strus The set of structures to destroy
 */
LLKA_API void LLKA_CC LLKA_destroyStructures(const LLKA_Structures *strus);

/*!
 * Makes a deep copy of \p LLKA_Atom
 *
 * @param[in] source Atom to copy from.
 * @param[out] target Atom to copy to.
 */
LLKA_API void LLKA_CC LLKA_duplicateAtom(const LLKA_Atom *source, LLKA_Atom *target);

/*! Makes a deep copy of \p LLKA_Structure
 *
 * @param[in] source Structure to copy from.
 * @param[out] target Structore to copy to.
 */
LLKA_API void LLKA_CC LLKA_duplicateStructure(const LLKA_Structure *source, LLKA_Structure *target);

/*!
 * Returns a fully intialized \p LLKA_Atom with default attributes.
 *
 * @return An empty \p LLKA_Atom
 */
LLKA_API LLKA_Atom LLKA_CC LLKA_emptyAtom(void);

/*!
 * Searches for an atom with the given attributes in \p stru.
 *
 * @param[in] stru The structure to search for the atom in.
 * @param[in] label_atom_id label_atom_id. If the value is <tt>NULL</tt>, the parameter is disregarded.
 * @param[in] label_comp_id label_comp_id. If the value is <tt>NULL</tt>, the parameter is disregarded.
 * @param[in] label_asym_id label_asym_id.If the value is <tt>NULL</tt>, the parameter is disregarded.
 * @param[in] label_seq_id label_seq_id. If the value is -1, the parameter is disregarded.
 * @param[in] label_alt_id label_alt_id. If the value is <tt>LLKA_NO_ALTID</tt>, the parameter is disregarded.
 * @param[in] pdbx_PDB_ins_code pdbx_PDB_ins_code. If the values is <tt>LLKA_NO_INSCODE</tt>, the parameter is disregarded.
 * @param[in] pdbx_PDB_model_num pdbx_PDB_model_num. Set to 1 for structures that contain just one model.
 *
 * @return A pointer to the first atom that matches the criteria, or <tt>NULL</tt> if no matching atom was found.
 */
LLKA_API LLKA_Atom * LLKA_findAtom(
    LLKA_Structure *stru,
    const char *label_atom_id,
    const char *label_comp_id,
    const char *label_asym_id,
    int32_t label_seq_id,
    char label_alt_id,
    const char *pdbx_PDB_ins_code,
    int32_t pdbx_PDB_model_num
);

/*!
 * Searches for an atom with the given attributes in \p view.
 *
 * @param[in] stru The structure to search for the atom in.
 * @param[in] label_atom_id label_atom_id. If the value is <tt>NULL</tt>, the parameter is disregarded.
 * @param[in] label_comp_id label_comp_id. If the value is <tt>NULL</tt>, the parameter is disregarded.
 * @param[in] label_asym_id label_asym_id.If the value is <tt>NULL</tt>, the parameter is disregarded.
 * @param[in] label_seq_id label_seq_id. If the value is -1, the parameter is disregarded.
 * @param[in] label_alt_id label_alt_id. If the value is <tt>LLKA_NO_ALTID</tt>, the parameter is disregarded.
 * @param[in] pdbx_PDB_ins_code pdbx_PDB_ins_code. If the values is <tt>LLKA_NO_INSCODE</tt>, the parameter is disregarded.
 * @param[in] pdbx_PDB_model_num pdbx_PDB_model_num. Set to 1 for structures that contain just one model.
 *
 * @return A pointer to the first atom that matches the criteria, or <tt>NULL</tt> if no matching atom was found.
 */
LLKA_API const LLKA_Atom * LLKA_findAtomView(
    LLKA_StructureView *view,
    const char *label_atom_id,
    const char *label_comp_id,
    const char *label_asym_id,
    int32_t label_seq_id,
    char label_alt_id,
    const char *pdbx_PDB_ins_code,
    int32_t pdbx_PDB_model_num
);

/*!
 * Initializes \p LLKA_Atom with the given attributes.
 *
 * The initialized \p LLKA_Atom contains copies of the passed attributes.
 *
 * @param[in] id id
 * @param[in] type_symbol type_symbol
 * @param[in] label_atom_id label_atom_id
 * @param[in] label_entity_id label_entity_id
 * @param[in] label_comp_id label_comp_id
 * @param[in] label_asym_id label_asym_id
 * @param[in] auth_atom_id auth_atom_id This may be <tt>NULL</tt> if no author atom id is provided
 * @param[in] auth_asym_id auth_comp_id This may be <tt>NULL</tt> if no author compound id is provided
 * @param[in] auth_asym_id auth_asym_id This may be <tt>NULL</tt> if no author assembly id is provided
 * @param[in] label_seq_id label_seq_id
 * @param[in] label_alt_id label_alt_id
 * @param[in] auth_seq_id auth_seq_id
 * @param[in] pdbx_PDB_model_num pdbx_PDB_model_num
 * @param[in] coords coords
 * @param[out] atom The atom to initialize
 */
LLKA_API void LLKA_CC LLKA_initAtom(
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
    LLKA_Atom *atom
);

LLKA_API void LLKA_CC LLKA_initStructure(const LLKA_Atom *atoms, size_t nAtoms, LLKA_Structure *stru);

/*!
 * Creates \p LLKA_Atom with the given a attributes.
 *
 * The created \p LLKA_Atom contains copies of the passed attributes.
 *
 * @param[in] id id
 * @param[in] type_symbol type_symbol
 * @param[in] label_atom_id label_atom_id
 * @param[in] label_entity_id label_entity_id
 * @param[in] label_comp_id label_comp_id
 * @param[in] label_asym_id label_asym_id
 * @param[in] auth_atom_id auth_atom_id This may be <tt>NULL</tt> if no author atom id is provided
 * @param[in] auth_comp_id auth_comp_id This may be <tt>NULL</tt> if no author compound id is provided
 * @param[in] auth_asym_id auth_asym_id This may be <tt>NULL</tt> if no author assembly id is provided
 * @param[in] label_seq_id label_seq_id
 * @param[in] label_alt_id label_alt_id
 * @param[in] auth_seq_id auth_seq_id
 * @param[in] pdbx_PDB_model_num pdbx_PDB_model_num
 * @param[in] coords coords
 *
 * @return Created LLKA_Atom
 */
LLKA_API LLKA_Atom LLKA_CC LLKA_makeAtom(
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
    const LLKA_Point *coords
);

/*!
 * Creates new \p LLKA_Structure from an array of atoms.
 *
 * The returned structure will contain deep copies of the input atoms.
 *
 * @param[in] atoms Atoms of the substructure
 * @param[in] nAtoms Number of input atoms
 *
 * @return Created LLKA_Structure
 */
LLKA_API LLKA_Structure LLKA_CC LLKA_makeStructure(const LLKA_Atom *atoms, size_t nAtoms);

/*!
 * Creates new \p LLKA_Structure from an array of pointers to atoms.
 *
 * The returned structure will contain deep copies of the input atoms.
 *
 * @param[in] atoms Pointers to atoms of the substructure
 * @param[in] nAtoms Number of input atoms
 *
 * @return Created LLKA_Structure
 */
LLKA_API LLKA_Structure LLKA_CC LLKA_makeStructureFromPtrs(const LLKA_Atom *const *atoms, size_t nAtoms);

/*!
 * Removes atom with a given id from the structure.
 *
 * @param[in] id Id of the atom to be removed.
 * @param[in,out] stru Structure to remove the atom from.
 */
LLKA_API void LLKA_CC LLKA_removeAtomById(uint32_t id, LLKA_Structure *stru);

/*!
 * Splits the structure into multiple substructures.
 * Each substructure will contain only atom with one alternate position ID.
 *
 * @param[in] stru structure to split.
 * @param[out] alts List of alternate position IDs found in the structure as an array of chars (not zero-terminated).
 *                  The list is ordered in the same way as the returned substructures, e. g. if the second element
 *                  of this list has value "B", the second substructure will contain only atoms with alternate position "B".
 * @returns \p LLKA_Structures with the splitted substructures
 */
LLKA_API LLKA_Structures LLKA_CC LLKA_splitByAltIds(const LLKA_Structure *stru, LLKA_AlternatePositions *alts);

/*!
 * Splits structure into a list of dinucleotide steps. A dinucleotide step contains two nucleotides that spatially
 * follow one another in the structure. If there are alternate positions, only atoms with matching alternate position IDs
 * are allowed to form a step.
 *
 * Note that the structure must be normalized (atoms must be grouped by model, chain and sequence id) in order for this
 * function to operate correctly.
 *
 * @param[in] stru Structure to be split into steps.
 * @param[out] steps LLKA_Structures object with the results. If the function does not return LLKA_OK,
 *                   content of this parameter is not modified.
 * @retval LLKA_OK Success.
 * @retval LLKA_E_BAD_DATA Structure cannot be split into steps because it containts data that this function cannot process.
 *                         The structure may not be normalized or it specifies alternate positions in an unusual way.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_splitStructureToDinucleotideSteps(const LLKA_Structure *stru, LLKA_Structures *steps);

LLKA_END_API_FUNCTIONS

#endif /* _LLKA_STRUCTURE_H */
