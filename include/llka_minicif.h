/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_MINICIF_H
#define _LLKA_MINICIF_H

#include "llka_structure.h"

/*! Flags that control import behavior. */
enum {
    LLKA_MINICIF_NORMALIZE   = (1 << 0),    /*!< Attempt to sort imported structure by model, chain and sequence id */
    LLKA_MINICIF_GET_CIFDATA = (1 << 1)     /*!< Export complete processed Cif data */
};

typedef struct LLKA_CifData LLKA_CifData;
typedef struct LLKA_CifDataPrivate LLKA_CifDataPrivate;
typedef struct LLKA_CifDataBlockPrivate LLKA_CifDataBlockPrivate;
typedef struct LLKA_CifDataCategoryPrivate LLKA_CifDataCategoryPrivate;
typedef struct LLKA_CifDataItemPrivate LLKA_CifDataItemPrivate;

typedef enum LLKA_CifDataValueState {
    LLKA_MINICIF_VALUE_SET  = 0,
    LLKA_MINICIF_VALUE_NONE = 1,
    LLKA_MINICIF_VALUE_UNKW = 2
    ENUM_FORCE_INT32_SIZE(LLKA_CifDataValueState)
} LLKA_CifDataValueState;

typedef struct LLKA_CifDataValue {
    const char *text;
    LLKA_CifDataValueState state;
} LLKA_CifDataValue;
LLKA_IS_POD(LLKA_CifDataValue)

typedef struct LLKA_CifDataItem {
    const char *keyword;
    LLKA_CifDataValue *values;
    size_t nValues;
    LLKA_CifDataItemPrivate *p;
} LLKA_CifDataItem;
LLKA_IS_POD(LLKA_CifDataItem)

typedef struct LLKA_CifDataCategory {
    const char *name;
    LLKA_CifDataItem *firstItem;
    LLKA_CifDataCategoryPrivate *p;
} LLKA_CifDataCategory;
LLKA_IS_POD(LLKA_CifDataCategory)

typedef struct LLKA_CifDataBlock {
    const char *name;
    LLKA_CifDataCategory *firstCategory;
    LLKA_CifDataBlockPrivate *p;
} LLKA_CifDataBlock;
LLKA_IS_POD(LLKA_CifDataBlock)

struct LLKA_CifData {
    LLKA_CifDataBlock *blocks;
    size_t nBlocks;
    LLKA_CifDataPrivate *p;
};
LLKA_IS_POD(LLKA_CifData)

typedef struct LLKA_StructureEntry {
    const char *id;
} LLKA_StructureEntry;
LLKA_IS_POD(LLKA_StructureEntry)

typedef struct LLKA_ImportedStructure {
    LLKA_StructureEntry entry;
    LLKA_Structure structure;
    LLKA_CifData *cifData;
} LLKA_ImportedStructure;
LLKA_IS_POD(LLKA_ImportedStructure)

LLKA_BEGIN_API_FUNCTIONS

/*!
 * Adds a new block to CifData.
 *
 * If you keep references to blocks in CifData, mind that these references might be
 * invalidated when this function returns.
 *
 * @param[in] data CifData to add the block to.
 * @param[in] name Name of the block, may be an empty string.
 */
LLKA_API void LLKA_CC LLKA_cifData_addBlock(LLKA_CifData *data, const char *name);

/*!
 * Validates contents of CifData and clears the taint flag if the data is valid.
 *
 * @param[in] cifData CifData to detaint.
 *
 * @retval LLKA_OK Data is valid. Taint flag cleared.
 * @retval LLKA_E_MISMATCHING_SIZES Number of values in entries of one category does not match.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifData_detaint(LLKA_CifData *cifData);

/*!
 * Removes a block from CifData.
 *
 * If you keep references to blocks in CifData, mind that these references might be
 * invalidated if this function succeeds.
 *
 * @param[in] data CifData to remove the block from.
 * @param[in] idx Index of the block in range from [0;nBlocks]
 *
 * @retval LLKA_OK Block removed.
 * @retval LLKA_E_INVALID_ARGUMENT Invalid block index.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifData_deleteBlock(LLKA_CifData *data, size_t idx);

/*!
 * Creates empty CifData object.
 *
 * @returns Empty CifData object
 */
LLKA_API LLKA_CifData * LLKA_CC LLKA_cifData_empty(void);

/*!
 * Checks if the CifData object is tainted.
 *
 * @param[in] data CifData to check
 *
 * @retval LLKA_TRUE Data is tainted.
 * @retval LLKA_TRUE Data is not tainted.
 */
LLKA_API LLKA_Bool LLKA_CC LLKA_cifData_isTainted(const LLKA_CifData *data);

/*!
 * Adds a new category to a block.
 *
 * Pointer to the category object will remain valid until the category or its containing block is removed.
 * Adding or deleting other categories will not invalidate the pointer.
 *
 * @param[in] block Block to add the category to.
 * @param[in] name Name of the category. The name must not be empty or the same an already existing category.
 *
 * @returns Pointer to the new category object or NULL if the function fails.
 */
LLKA_API LLKA_CifDataCategory * LLKA_CC LLKA_cifDataBlock_addCategory(LLKA_CifDataBlock *block, const char *name);

/*!
 * Deletes a category from a block.
 *
 * @param[in] block Block to remove the category from.
 * @param[in] name Name of the category.
 *
 * @retval LLKA_OK Category deleted.
 * @retval LLKA_E_INVALID_ARGUMENT Category of the given name was not found.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifDataBlock_deleteCategory(LLKA_CifDataBlock *block, const char *name);

/*!
 * Finds a category in a block by name.
 *
 * @param[in] block Block to search ing
 * @param[in] name Name of the category to find.
 *
 * @returns Pointer to the catgory object or NULL if a category with a matching name was not found.
 */
LLKA_API LLKA_CifDataCategory * LLKA_CC LLKA_cifDataBlock_findCategory(const LLKA_CifDataBlock *block, const char *name);

/*!
 * Initializes an empty \p CifDataBlock object.
 *
 * This function is useful for code that preallocates an array for \p CifDataBlock s that it needs to initialize.
 * Pointers to \p name and \p p must be set to zero, otherwise the function will fail and not modify the object.
 *
 * @param[in,out] block The block to initialize
 * @param[in] name Name of the block. May be an empty string.
 * @param[in] data The CifData object that contains the block.
 *
 * @retval LLKA_OK Success
 * @retval LLKA_E_INVALID_ARGUMENT Passed block was not zeroed.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifDataBlock_initialize(LLKA_CifDataBlock *block, const char *name, LLKA_CifData *data);

/*!
 * Jumps to the next category in a block.
 *
 * If the order in which all categories were added to a block is known,
 * this function may be used to quicky iterate over all categories.
 * This function shall be used with care if categories are added conditionally
 * of if some categories are deleted.
 *
 * @param[in] cat Current category.
 *
 * @retval Pointer to the next category or NULL if there is no next category.
 */
LLKA_API LLKA_CifDataCategory * LLKA_CC LLKA_cifDataBlock_nextCategory(const LLKA_CifDataCategory *cat);


/*!
 * Jumps to the previous category in a block.
 *
 * If the order in which all categories were added to a block is known,
 * this function may be used to quicky iterate over all categories.
 * This function shall be used with care if categories are added conditionally
 * of if some categories are deleted.
 *
 * @param[in] cat Current category.
 *
 * @retval Pointer to the previous category or NULL if there is no previous category.
 */
LLKA_API LLKA_CifDataCategory * LLKA_CC LLKA_cifDataBlock_previousCategory(const LLKA_CifDataCategory *cat);

/*!
 * Add an item to category.
 *
 * Pointer to the item object will remain valid until the item or its containing category is removed.
 * Adding or deleting other entries will not invalidate the pointer.
 *
 * Adding an item taints CifData.
 *
 * @param[in] cat Category to add the item to.
 * @param[in] keyword Keyword of the item, The keyword must not be an empty string or same as a keyword of already exiting item in the category.
 *
 * @returns Pointer to the item object or NULL if the function fails.
 */
LLKA_API LLKA_CifDataItem * LLKA_CC LLKA_cifDataCategory_addItem(LLKA_CifDataCategory *cat, const char *keyword);

/*!
 * Deletes an item from a category.
 *
 * @param[in] cat Category to delete the item from.
 * @param[in] keyword Keyword of the item to delete.
 *
 * @retval LLKA_OK Item deleted.
 * @retval LLKA_E_INVALID_ARGUMENT Item with the given keyword was not found.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifDataCategory_deleteItem(LLKA_CifDataCategory *cat, const char *keyword);

/*!
 * Finds an item in a category by name.
 *
 * @param[in] cat Category to search in.
 * @param[in] keyword Keyword to look for.
 *
 * @returns Pointer to the item object or NULL if no item with a matching keyword was found.
 */
LLKA_API LLKA_CifDataItem * LLKA_CC LLKA_cifDataCategory_findItem(const LLKA_CifDataCategory *cat, const char *keyword);

/*!
 * Checks if a category is a loop or a key-value kind.
 *
 * @param[in] cat Category to check.
 *
 * @retval LLKA_TRUE Category is a loop
 * @retval LLKA_FALSE Category is a key-value
 */
LLKA_API LLKA_Bool LLKA_CC LLKA_cifDataCategory_isLoop(LLKA_CifDataCategory *cat);

/*!
 * Jumps to the next item in a category,
 *
 * If the order in which all entries were added to a category is known,
 * this function may be used to quicky iterate over all entries.
 * This function shall be used with care if entries are added conditionally
 * of if some entries are deleted.
 *
 * @param[in] item Current item.
 *
 * @returns Pointer to the next item or NULL if there is no next item.
 */
LLKA_API LLKA_CifDataItem * LLKA_CC LLKA_cifDataCategory_nextItem(const LLKA_CifDataItem *item);

/*!
 * Jumps to the previous item in a category,
 *
 * If the order in which all entries were added to a category is known,
 * this function may be used to quicky iterate over all entries.
 * This function shall be used with care if entries are added conditionally
 * of if some entries are deleted.
 *
 * @param[in] item Current item.
 *
 * @returns Pointer to the previous item or NULL if there is no previous item.
 */
LLKA_API LLKA_CifDataItem * LLKA_CC LLKA_cifDataCategory_previousItem(const LLKA_CifDataItem *item);

/*!
 * Appends new value to the list of values.
 *
 * MiniCif will cretae an internal copy of the value.
 *
 * Adding a value taints CifData.
 *
 * @param[in] item Item to append the value to.
 * @param[in] value Value to append.
 */
LLKA_API void LLKA_CC LLKA_cifDataItem_addValue(LLKA_CifDataItem *item, const LLKA_CifDataValue *value);

/*!
 * Appends new values to the list of values.
 *
 * MiniCif will create an internal copy of the values.
 *
 * Adding values taints CifData.
 *
 * @param[in] item Item to append the value to.
 * @param[in] values Array of values to append.
 * @param[in] nValues Length of the array of values to append.
 */
LLKA_API void LLKA_CC LLKA_cifDataItem_addValues(LLKA_CifDataItem *item, const LLKA_CifDataValue *values, size_t nValues);

/*!
 * Sets new values for an item.
 *
 * MiniCif will create an internal copy of the values.
 *
 * Setting the values taints CifData.
 *
 * @param[in] imte Item to set data for.
 * @param[in] values An array of strings to use as values.
 * @param[in] nValues Length of the values array.
 */
LLKA_API void LLKA_CC LLKA_cifDataItem_setValues(LLKA_CifDataItem *item, const LLKA_CifDataValue *values, size_t nValues);

/*!
 * Converts processed data into a (mm)Cif string.
 *
 * @param[in] cifData Cif data to be written out as (mm)Cif string.
 * @param[in] pretty Produce a neatly padded output. Pretty output is larger and slower to generate
 *                   but easier for humans to read.
 * @param[out] cifString Pointer to the output string;
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_BAD_DATA Cif data contains invalid content.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifDataToString(const LLKA_CifData *cifData, LLKA_Bool pretty, char **cifString);

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
/*!
 * Creates LLKA_CifData from CIF file.
 *
 * @param[in] path Path to the CIF file.
 * @param[out] data LLKA_CifData to be created from the CIF file.
 * @param[out] error Contains brief error description in case the CIF file cannot be parsed. Set to NULL of the file is parsed successfully.
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_NO_FILE CIF file does not exist.
 * @retval LLKA_E_CANNOT_READ_FILE CIF file exists but cannot be read.
 * @retval LLKA_E_BAD_DATA CIF file is malformed. Check the \p error string for more information.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifFileToData(const LLKA_PathChar *path, LLKA_CifData **data, char **error);
#endif /* LLKA_FILESYSTEM_ACCESS_DISABLED */

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
/*!
 * Creates LLKA_Structure from CIF file.
 *
 * @param[in] path Path to the CIF file.
 * @param[out] stru LLKA_Structure to be created from the CIF file.
 * @param[out] error Contains brief error description in case the CIF file cannot be parsed. Set to NULL of the file is parsed successfully.
 * @param[in] options Additional import flags.
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_NO_FILE CIF file does not exist.
 * @retval LLKA_E_NO_DATA CIF file is empty.
 * @retval LLKA_E_CANNOT_READ_FILE CIF file exists but cannot be read.
 * @retval LLKA_E_BAD_DATA CIF file is malformed. Check the \p error string for more information.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifFileToStructure(const LLKA_PathChar *path, LLKA_ImportedStructure *importedStru, char **error, int32_t options);
#endif /* LLKA_FILESYSTEM_ACCESS_DISABLED */

/*!
 * Creates LLKA_Structure from a CIF-formatted raw text.
 *
 * @param[in] text Text to parse as CIF data.
 * @param[out] data LLKA_CifData to be created from the CIF file.
 * @param[out] error Contains brief error description in case the CIF file cannot be parsed. Set to NULL of the file is parsed successfully.
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_NO_FILE CIF file does not exist.
 * @retval LLKA_E_NO_DATA CIF file is empty.
 * @retval LLKA_E_CANNOT_READ_FILE CIF file exists but cannot be read.
 * @retval LLKA_E_BAD_DATA CIF file is malformed. Check the \p error string for more information.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifTextToData(const char *text, LLKA_CifData **data, char **error);

/*!
 * Creates LLKA_Structure from a CIF-formatted raw text.
 *
 * @param[in] text Text to parse as CIF data.
 * @param[out] stru LLKA_Structure to be created from the CIF file.
 * @param[out] error Contains brief error description in case the CIF file cannot be parsed. Set to NULL of the file is parsed successfully.
 * @param[in] options Additional import flags.
 *
 * @retval LLKA_OK Success.
 * @retval LLKA_E_NO_FILE CIF file does not exist.
 * @retval LLKA_E_CANNOT_READ_FILE CIF file exists but cannot be read.
 * @retval LLKA_E_BAD_DATA CIF file is malformed. Check the \p error string for more information.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_cifTextToStructure(const char *text, LLKA_ImportedStructure *importedStru, char **error, int32_t options);

/*!
 * Destroys Cif data.
 *
 * Note that \p LLKA_destroyImportedStructure() calls this function internally so there is
 * no need to explicitly destroy Cif data returned in as part of LLKA_ImportedStructure.
 */
LLKA_API void LLKA_CC LLKA_destroyCifData(LLKA_CifData *cifData);

/*!
 * Destroys structure imported from Cif
 *
 * @param[in] importedStru Structure to destroy
 */
LLKA_API void LLKA_CC LLKA_destroyImportedStructure(LLKA_ImportedStructure *importedStru);

LLKA_END_API_FUNCTIONS

#endif /* _LLKA_MINICIF_H */
