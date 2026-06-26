/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_minicif.h>

#include "effedup.hpp"

static
auto test_ok()
{
    LLKA_ImportedStructure importedStru{};
    char *error;

    auto tRet = LLKA_cifFileToStructure(LLKA_PathLiteral("./1BNA.cif"), &importedStru, &error, 0);
    EFF_expect(tRet, LLKA_OK, "unexpected return value");
    EFF_expect(error, nullptr, "error should have been set to nullptr");

    EFF_expect(importedStru.structure.nAtoms, 566UL, "wrong number of atoms in structure");
    EFF_expect(importedStru.structure.atoms[165].label_seq_id, 9, "unexpected label_seq_id");

    LLKA_destroyImportedStructure(&importedStru);
}

static
auto test_broken()
{
    LLKA_ImportedStructure importedStru{};
    char *error;

    auto tRet = LLKA_cifFileToStructure(LLKA_PathLiteral("./1BNA_broken.cif"), &importedStru, &error, 0);
    EFF_expect(tRet, LLKA_E_BAD_DATA, "unexpected return value");
    EFF_expect((const char *)error, "Malformed loop on line 1089, unexpected token kind 8", "unexpected error message");

    LLKA_destroyString(error);
}

static
auto test_manipulate_categories()
{
    auto data = LLKA_cifData_empty();

    LLKA_cifData_addBlock(data, "");
    auto block = &data->blocks[0];

    // Add three categories
    auto c = LLKA_cifDataBlock_addCategory(block, "first");
    if (c == nullptr)
        EFF_fail("Failed to add first category");

    c = LLKA_cifDataBlock_addCategory(block, "second");
    if (c == nullptr)
        EFF_fail("Failed to add second category");

    c = LLKA_cifDataBlock_addCategory(block, "third");
    if (c == nullptr)
        EFF_fail("Failed to add third category");

    EFF_expect(
        LLKA_cifDataBlock_addCategory(block, "third"),
        nullptr,
        "Failed to prevent creation of duplicate category"
    );

    // Check the order
    auto cat = block->firstCategory;
    EFF_expect(cat->name, "first", "wrong first category name");

    cat = LLKA_cifDataBlock_nextCategory(cat);
    EFF_expect(cat->name, "second", "wrong second category name");

    cat = LLKA_cifDataBlock_nextCategory(cat);
    EFF_expect(cat->name, "third", "wrong third category name");

    cat = LLKA_cifDataBlock_nextCategory(cat);
    EFF_expect(cat, nullptr, "Last category does not look as last");

    // Delete second category and check the order again
    auto tRet = LLKA_cifDataBlock_deleteCategory(block, "second");
    EFF_expect(tRet, LLKA_OK, std::string{"Failed to delete Cif category: "} + LLKA_errorToString(tRet));

    cat = block->firstCategory;
    EFF_expect(cat->name, "first", "wrong first category name");

    cat = LLKA_cifDataBlock_nextCategory(cat);
    EFF_expect(cat->name, "third", "wrong third category name");

    // Delete first category and check that the third category is set as the first and only category
    tRet = LLKA_cifDataBlock_deleteCategory(block, "first");
    EFF_expect(tRet, LLKA_OK, std::string{"Failed to delete Cif category: "} + LLKA_errorToString(tRet));

    cat = block->firstCategory;
    EFF_expect(cat->name, "third", "wrong third category name");

    EFF_expect(LLKA_cifDataBlock_nextCategory(cat), nullptr, "Category is not the only left category (has next check)");
    EFF_expect(LLKA_cifDataBlock_previousCategory(cat), nullptr, "Category is not the only left category (has previous check)");

    LLKA_destroyCifData(data);
}

static
auto test_manipulate_entries()
{
    auto data = LLKA_cifData_empty();

    LLKA_cifData_addBlock(data, "");
    auto block = &data->blocks[0];

    auto cat = LLKA_cifDataBlock_addCategory(block, "xyz");

    // Add three entries
    auto itemA = LLKA_cifDataCategory_addItem(cat, "a");
    if (itemA == nullptr)
        EFF_fail("Cannot add itemA");

    auto itemB = LLKA_cifDataCategory_addItem(cat, "b");
    if (itemB == nullptr)
        EFF_fail("Cannot add itemB");

    auto itemC = LLKA_cifDataCategory_addItem(cat, "c");
    if (itemC == nullptr)
        EFF_fail("Cannot add itemC");

    EFF_expect(
        LLKA_cifDataCategory_addItem(cat, "b"),
        nullptr,
        "Failed to prevent creation of duplicate item"
    );

    // Check order
    auto item = cat->firstItem;
    EFF_expect(item->keyword, "a", "Wrong keyword");

    item = LLKA_cifDataCategory_nextItem(item);
    EFF_expect(item->keyword, "b", "Wrong keyword");

    item = LLKA_cifDataCategory_nextItem(item);
    EFF_expect(item->keyword, "c", "Wrong keyword");

    item = LLKA_cifDataCategory_nextItem(item);
    EFF_expect(item, nullptr, "Last item does not look as last");

    // Delete item "c" and check again
    auto tRet = LLKA_cifDataCategory_deleteItem(cat, "c");
    EFF_expect(tRet,  LLKA_OK, "Failed to delete item");

    item = cat->firstItem;
    EFF_expect(item->keyword, "a", "Wrong keyword");

    item = LLKA_cifDataCategory_nextItem(item);
    EFF_expect(item->keyword, "b", "Wrong keyword");

    EFF_expect(LLKA_cifDataCategory_nextItem(item), nullptr, "Last item does not look like last (has next check)");

    item = LLKA_cifDataCategory_previousItem(item);
    EFF_expect(item, itemA, "Entry should have a previous item but it does not");

    // Add some values and check tainting
    auto valuesA = new LLKA_CifDataValue[2];
    valuesA[0].text = "Run to the hills";
    valuesA[0].state = LLKA_MINICIF_VALUE_SET;
    valuesA[1].text = "Run for your life";
    valuesA[1].state = LLKA_MINICIF_VALUE_SET;
    LLKA_cifDataItem_setValues(itemA, valuesA, 2);
    delete [] valuesA;

    LLKA_CifDataValue val{ .text = "I am a man who walks alone", .state = LLKA_MINICIF_VALUE_SET };
    LLKA_cifDataItem_addValue(itemB, &val);

    LLKA_CifDataValue val2{ .text = "when I'm walking a dark road", .state = LLKA_MINICIF_VALUE_SET };
    LLKA_cifDataItem_addValue(itemB, &val2);

    EFF_expect(LLKA_cifData_isTainted(data), LLKA_TRUE, "Data is not marked as tainted when it should be");
    tRet = LLKA_cifData_detaint(data);
    EFF_expect(tRet, LLKA_OK, "Failed to detaint valid Cif data");

    LLKA_CifDataValue val3{ .text = "at night or strolling through a park", .state = LLKA_MINICIF_VALUE_SET };
    LLKA_cifDataItem_addValue(itemB, &val3);
    tRet = LLKA_cifData_detaint(data);
    EFF_expect(tRet, LLKA_E_MISMATCHING_SIZES, "Incorrectly detainted invalid data");

    LLKA_destroyCifData(data);
}

static
auto test_loop_not_loop()
{
    auto data = LLKA_cifData_empty();

    LLKA_cifData_addBlock(data, "");
    auto block = &data->blocks[0];

    auto loopCat = LLKA_cifDataBlock_addCategory(block, "xyz");
    auto notLoopCat = LLKA_cifDataBlock_addCategory(block, "abc");

    auto loopEntry = LLKA_cifDataCategory_addItem(loopCat, "loop_item");
    auto notLoopEntry = LLKA_cifDataCategory_addItem(notLoopCat, "not_loop_item");

    LLKA_CifDataValue val{ .text = "Osculum", .state = LLKA_MINICIF_VALUE_SET };
    LLKA_cifDataItem_addValue(loopEntry, &val);
    LLKA_CifDataValue val2{ .text = "obscenum", .state = LLKA_MINICIF_VALUE_SET };
    LLKA_cifDataItem_addValue(loopEntry, &val2);

    LLKA_CifDataValue val3{ .text = "Spillways", .state = LLKA_MINICIF_VALUE_SET };
    LLKA_cifDataItem_addValue(notLoopEntry, &val3);

    LLKA_cifData_detaint(data);

    EFF_expect(LLKA_cifDataCategory_isLoop(loopCat), LLKA_TRUE, "Category is not a loop but it should be");
    EFF_expect(LLKA_cifDataCategory_isLoop(notLoopCat), LLKA_FALSE, "Category is a loop but it should not be");

    LLKA_destroyCifData(data);
}

static
auto GET_AND_CHECK_TEXT_CIF_VALUE(const LLKA_CifDataCategory *itemC, const char *name, const char *text)
{
    auto item = LLKA_cifDataCategory_findItem(itemC, name);
    if (item == nullptr)
        EFF_fail("cannot find item \"" + std::string{name} + "\"");
    EFF_expect(item->nValues, 1ULL, "unexpected number of values");
    EFF_expect(item->values[0].state, LLKA_MINICIF_VALUE_SET, "unexpected state of CIF value");
    EFF_expect(item->values[0].text, text, "unexpected text of value");
}
static
auto test_quote_in_quotes()
{
    LLKA_CifData *data;
    char *error;

    auto tRet = LLKA_cifFileToData(LLKA_PathLiteral("./quote_in_quotes.cif"), &data, &error);
    EFF_expect(tRet, LLKA_OK, "unexpected return value");
    EFF_expect(data->nBlocks, 1ULL, "unexpected number of blocks");

    auto block = &data->blocks[0];
    auto itemC = LLKA_cifDataBlock_findCategory(block, "some");
    if (itemC == nullptr)
        EFF_fail("cannot find the expected category");

    GET_AND_CHECK_TEXT_CIF_VALUE(itemC, "double_quote_string", "\"ABC\"");
    GET_AND_CHECK_TEXT_CIF_VALUE(itemC, "single_quote_string", "'ABC'");
    GET_AND_CHECK_TEXT_CIF_VALUE(itemC, "single_quote_in_double_quote_string", "\"AB'C\"");
    GET_AND_CHECK_TEXT_CIF_VALUE(itemC, "double_quote_in_single_quote_string", "'AB\"C'");
    GET_AND_CHECK_TEXT_CIF_VALUE(itemC, "single_quote_in_single_quote_string", "'AB'C'");
    GET_AND_CHECK_TEXT_CIF_VALUE(itemC, "single_quote_in_single_quote_string_terminating", "'ABC''");
    GET_AND_CHECK_TEXT_CIF_VALUE(itemC, "double_quote_in_double_quote_string_terminating", "\"ABC\"\"");


    LLKA_destroyCifData(data);
}

static
auto test_unterminated_quote()
{
    LLKA_CifData *data;
    char *error;

    auto tRet = LLKA_cifFileToData(LLKA_PathLiteral("./unterminated_quote.cif"), &data, &error);
    EFF_expect(tRet, LLKA_E_BAD_DATA, "unexpected return value");
    if (error == nullptr)
        EFF_fail("expected to have an error string");
    EFF_expect(error, "Unterminated quoted token that begins on line 3", "unexpected error condition");

    LLKA_destroyString(error);
}

static
auto test_no_entry_category()
{
    static const char *cifText =
        "data_\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.segment_id\n"
        "_atom_site.cartn_x\n"
        "_atom_site.cartn_y\n"
        "_atom_site.cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.cartn_x_esd\n"
        "_atom_site.cartn_y_esd\n"
        "_atom_site.cartn_z_esd\n"
        "_atom_site.occupancy_esd\n"
        "_atom_site.B_iso_or_equiv_esd\n"
        "_atom_site.pdbx_formal_charge\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1  P P     . DA A 1 99 . . 8.258  14.006  12.104 1 30 ? ? ? ? ? ? 1 DT A P     1\n"
        "#\n"
        "";

    LLKA_ImportedStructure importedStru{};
    char *error;

    auto tRet = LLKA_cifTextToStructure(cifText, &importedStru, &error, LLKA_MINICIF_ALLOW_NO_ENTRY_CATEGORY);
    EFF_expect(tRet, LLKA_OK, "unexpected return value");
    EFF_expect(error, nullptr, "error should have been set to nullptr");

    EFF_expect(importedStru.structure.nAtoms, 1ULL, "wrong number of atoms in structure");
    EFF_expect(importedStru.structure.atoms[0].label_seq_id, 99, "unexpected label_seq_id");
    EFF_expect(importedStru.entry.id, "", "unexpected entry id");

    LLKA_destroyImportedStructure(&importedStru);
}

static
auto test_no_entry_category_not_allowed()
{
    static const char *cifText =
        "data_\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.segment_id\n"
        "_atom_site.cartn_x\n"
        "_atom_site.cartn_y\n"
        "_atom_site.cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.cartn_x_esd\n"
        "_atom_site.cartn_y_esd\n"
        "_atom_site.cartn_z_esd\n"
        "_atom_site.occupancy_esd\n"
        "_atom_site.B_iso_or_equiv_esd\n"
        "_atom_site.pdbx_formal_charge\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1  P P     . DA A 1 99 . . 8.258  14.006  12.104 1 30 ? ? ? ? ? ? 1 DT A P     1\n"
        "";

    LLKA_ImportedStructure importedStru{};
    char *error;

    auto tRet = LLKA_cifTextToStructure(cifText, &importedStru, &error, 0);
    EFF_expect(tRet, LLKA_E_BAD_DATA, "unexpected return value");
    EFF_expect(error, "Category entry is not present in block ", "unexpected error message");

    LLKA_destroyString(error);
}

static
auto test_allow_broken_atomsite()
{
    static const char *cifText =
        "data_\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.segment_id\n"
        "_atom_site.cartn_x\n"
        "_atom_site.cartn_y\n"
        "_atom_site.cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.cartn_x_esd\n"
        "_atom_site.cartn_y_esd\n"
        "_atom_site.cartn_z_esd\n"
        "_atom_site.occupancy_esd\n"
        "_atom_site.B_iso_or_equiv_esd\n"
        "_atom_site.pdbx_formal_charge\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        " ATOM 8  N .     . DA . 1 20 . . -2.3707326 -4.507768  -0.94517957  1 20 ? ? ? ? ? ?  1 DA A N9    1\n"
        "";

    LLKA_ImportedStructure importedStru{};
    char *error;

    auto tRet = LLKA_cifTextToStructure(cifText, &importedStru, &error, LLKA_MINICIF_ALLOW_NO_ENTRY_CATEGORY | LLKA_MINICIF_ALLOW_BROKEN_ATOMSITE);
    EFF_expect(tRet, LLKA_OK, "unexpected return value");
    EFF_expect(error, nullptr, "error should have been set to nullptr");

    EFF_expect(importedStru.structure.nAtoms, 1ULL, "wrong number of atoms in structure");
    EFF_expect(importedStru.entry.id, "", "unexpected entry id");
    EFF_expect(importedStru.structure.atoms[0].label_seq_id, 20, "unexpected label_seq_id");
    EFF_expect(importedStru.structure.atoms[0].label_atom_id, "N9", "unexpected label_atom_id");
    EFF_expect(importedStru.structure.atoms[0].label_asym_id, "A", "unexpected label_asym_id");

    LLKA_destroyImportedStructure(&importedStru);
}

static
auto test_allow_broken_atomsite_not_allowed()
{
    static const char *cifText =
        "data_\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.segment_id\n"
        "_atom_site.cartn_x\n"
        "_atom_site.cartn_y\n"
        "_atom_site.cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.cartn_x_esd\n"
        "_atom_site.cartn_y_esd\n"
        "_atom_site.cartn_z_esd\n"
        "_atom_site.occupancy_esd\n"
        "_atom_site.B_iso_or_equiv_esd\n"
        "_atom_site.pdbx_formal_charge\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        " ATOM 8  N .     . DA . 1 20 . . -2.3707326 -4.507768  -0.94517957  1 20 ? ? ? ? ? ?  1 DA A N9    1\n"
        "";

    LLKA_ImportedStructure importedStru{};
    char *error;

    auto tRet = LLKA_cifTextToStructure(cifText, &importedStru, &error, LLKA_MINICIF_ALLOW_NO_ENTRY_CATEGORY);
    EFF_expect(tRet, LLKA_E_BAD_DATA, "unexpected return value");

    LLKA_destroyString(error);
}

auto main(int, char **) -> int
{
    test_ok();
    test_broken();

    test_manipulate_categories();
    test_manipulate_entries();
    test_loop_not_loop();

    test_quote_in_quotes();
    test_unterminated_quote();

    test_no_entry_category();
    test_no_entry_category_not_allowed();

    test_allow_broken_atomsite();
    test_allow_broken_atomsite_not_allowed();
}
