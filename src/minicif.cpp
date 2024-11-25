/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "minicif/minicif_p.h"

#include "minicif/parser.h"
#include "minicif/writer.h"
#include "minicif/categories/atom-site.hpp"
#include "minicif/categories/entry.hpp"

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
    #if defined(LLKA_PLATFORM_UNIX) || defined(LLKA_PLATFORM_EMSCRIPTEN)
        #define ENABLE_READ_UNIX_MMAP

        #include <fcntl.h>
        #include <unistd.h>
        #include <sys/mman.h>
    #elif defined LLKA_PLATFORM_WIN32
        #define ENABLE_WIN32_MMAP

        #include <Windows.h>
    #endif // LLKA_PLATFORM_
#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

#include <cassert>
#include <cstring>
#include <filesystem>
#include <ranges>
#include <set>
#include <string_view>

namespace LLKAInternal::MiniCif {

static
auto compareAtoms(const LLKA_Atom &lhs, const LLKA_Atom &rhs)
{
    if (lhs.pdbx_PDB_model_num < rhs.pdbx_PDB_model_num)
        return true;

    return lhs.label_seq_id < rhs.label_seq_id;
}

static
auto destroyCifDataValues(LLKA_CifDataValue *values, size_t nValues)
{
    if (values == nullptr)
        return;

    for (size_t valIdx = 0; valIdx < nValues; valIdx++)
        LLKAInternal::destroyString(values[valIdx].text);
}

static
auto destroyCifDataItem(LLKA_CifDataItem *item)
{
    destroyCifDataValues(item->values, item->nValues);
    delete [] item->values;
    LLKAInternal::destroyString(item->keyword);

    delete item->p;
}

static
auto destroyCifDataCategory(LLKA_CifDataCategory *cat)
{
    auto item = cat->firstItem;
    while (item) {
        auto next = item->p->next;
        destroyCifDataItem(item);
        delete item;

        item = next;
    }

    LLKAInternal::destroyString(cat->name);
    delete cat->p;
}

static
auto destroyCifDataBlock(LLKA_CifDataBlock *block)
{
    auto cat = block->firstCategory;
    while (cat != nullptr) {
        auto next = cat->p->next;
        destroyCifDataCategory(cat);
        delete cat;

        cat = next;
    }

    LLKAInternal::destroyString(block->name);
    delete block->p;
}

static
auto normalize(std::unique_ptr<LLKA_Atom[]> &atoms, const size_t nAtoms)
{
    std::sort(atoms.get(), atoms.get() + nAtoms, compareAtoms);
}

static
auto toCifData(const std::vector<Block> &blocks)
{
    auto cifData = LLKA_cifData_empty();
    cifData->blocks = new LLKA_CifDataBlock[blocks.size()];
    cifData->nBlocks = blocks.size();

    size_t blockIdx = 0;
    for (const auto &block : blocks) {
        auto &cdBlock = cifData->blocks[blockIdx];
        cdBlock.name = duplicateString(block.name);
        cdBlock.firstCategory = nullptr;
        cdBlock.p = new LLKA_CifDataBlockPrivate;
        cdBlock.p->root = cifData;

        if (block.categories.empty())
            break;

        LLKA_CifDataCategory *catNext = nullptr;
        for (const auto &[ catName, catLwrName, catItems ] : std::ranges::reverse_view(block.categories)) {
            auto cdCat = new LLKA_CifDataCategory;
            cdCat->name = duplicateString(catName);
            cdCat->p = new LLKA_CifDataCategoryPrivate;
            cdCat->p->prev = nullptr;
            cdCat->p->next = catNext;
            cdCat->p->root = cifData;
            cdCat->firstItem = nullptr;

            LLKA_CifDataItem *itemNext = nullptr;
            for (const auto &[ keyword, lwrKeyword, values ] : std::ranges::reverse_view(catItems)) {
                auto cdItem = new LLKA_CifDataItem;
                cdItem->keyword = LLKAInternal::duplicateString(keyword);
                cdItem->values = values.size() == 0 ? nullptr : new LLKA_CifDataValue[values.size()];
                cdItem->nValues = values.size();
                cdItem->p = new LLKA_CifDataItemPrivate;
                cdItem->p->prev = nullptr;
                cdItem->p->next = itemNext;
                cdItem->p->root = cifData;

                cdCat->p->isLoop = cdItem->nValues > 1;

                for (size_t valueIdx = 0; valueIdx < values.size(); valueIdx++) {
                    const auto &v = values[valueIdx];
                    auto &cdValue = cdItem->values[valueIdx];

                    cdValue.state = static_cast<LLKA_CifDataValueState>(v.state);
                    if (v.state == Value::State::VALUE)
                        cdValue.text = LLKAInternal::duplicateString(v.text);
                    else
                        cdValue.text = nullptr;
                }
                if (itemNext != nullptr)
                    itemNext->p->prev = cdItem;
                itemNext = cdItem;
            }
            cdCat->firstItem = itemNext;

            if (catNext != nullptr)
                catNext->p->prev = cdCat;
            catNext = cdCat;
        }
        cdBlock.firstCategory = catNext;

        blockIdx++;
    }

    return cifData;
}

static
auto toStructure(const std::string_view &view, LLKA_ImportedStructure *importedStru, char **error, int32_t options)
{
    try {
        auto blocks = parse(view);

        // TODO: We should look into the potential memory leaks here if we get unexpected data

        auto [ entires, nEntries ] = LLKAInternal::MiniCif::Applier<LLKAInternal::MiniCif::Categories::Entry>::apply(
            blocks[0],
            NoopFixup<LLKA_StructureEntry>,
            [](const LLKA_StructureEntry &e) { LLKAInternal::destroyString(e.id);
        });

        if (nEntries < 1)
            return LLKA_E_BAD_DATA;

        auto fixupAtom = [](LLKA_Atom &atom) {
            if (atom.auth_atom_id == nullptr)
                atom.auth_atom_id = LLKAInternal::duplicateString(atom.label_atom_id);
            if (atom.auth_comp_id == nullptr)
                atom.auth_comp_id = LLKAInternal::duplicateString(atom.label_comp_id);
            if (atom.auth_asym_id == nullptr)
                atom.auth_asym_id = LLKAInternal::duplicateString(atom.label_asym_id);
            if (atom.pdbx_PDB_ins_code == nullptr)
                atom.pdbx_PDB_ins_code = LLKAInternal::duplicateString(LLKA_NO_INSCODE);
        };
        auto [ atoms, nAtoms ] = LLKAInternal::MiniCif::Applier<LLKAInternal::MiniCif::Categories::AtomSite>::apply(
            blocks[0],
            fixupAtom,
            [](const LLKA_Atom &atom) { LLKA_destroyAtom(&atom); }
        );

        if (options & LLKA_MINICIF_NORMALIZE)
            LLKAInternal::MiniCif::normalize(atoms, nAtoms);

        importedStru->entry.id = entires[0].id;
        importedStru->structure.atoms = atoms.release();
        importedStru->structure.nAtoms = nAtoms;

        if (options & LLKA_MINICIF_GET_CIFDATA)
            importedStru->cifData = toCifData(blocks);
        else
            importedStru->cifData = nullptr;

        *error = nullptr;

        return LLKA_OK;
    } catch (const LLKAInternal::MiniCif::CifParseError &ex) {
        const auto len = std::strlen(ex.what());
        *error = new char[len + 1];
        std::strcpy(*error, ex.what());

        return LLKA_E_BAD_DATA;
    }
}

} // name LLKAInternal::MiniCif

void LLKA_CC LLKA_cifData_addBlock(LLKA_CifData *data, const char *name)
{
    auto newBlocks = new LLKA_CifDataBlock[data->nBlocks + 1];

    if (data->blocks != nullptr)
        std::memcpy(newBlocks, data->blocks, sizeof(LLKA_CifDataBlock) * data->nBlocks);

    auto nb = &newBlocks[data->nBlocks];
    nb->name = LLKAInternal::duplicateString(name);
    nb->firstCategory = nullptr;
    nb->p = new LLKA_CifDataBlockPrivate;
    nb->p->root = data;

    delete [] data->blocks;
    data->blocks = newBlocks;
    data->nBlocks++;
}

LLKA_RetCode LLKA_CC LLKA_cifData_deleteBlock(LLKA_CifData *data, size_t idx)
{
    if (idx >= data->nBlocks)
        return LLKA_E_INVALID_ARGUMENT;

    if (data->nBlocks == 1) {
        LLKAInternal::MiniCif::destroyCifDataBlock(&data->blocks[0]);
        delete [] data->blocks;
        data->blocks = nullptr;
        data->nBlocks = 0;
    } else {
        LLKAInternal::MiniCif::destroyCifDataBlock(&data->blocks[idx]);
        auto newBlocks = new LLKA_CifDataBlock[data->nBlocks + 1];

        std::memcpy(newBlocks, data->blocks, sizeof(LLKA_CifDataBlock) * idx);
        std::memcpy(newBlocks + idx, data->blocks + idx + 1, sizeof(LLKA_CifDataBlock) * (data->nBlocks - idx));

        delete [] data->blocks;
        data->blocks = newBlocks;
        data->nBlocks--;
    }

    return LLKA_OK;
}

LLKA_RetCode LLKA_CC LLKA_cifData_detaint(LLKA_CifData *cifData)
{
    if (!cifData->p->tainted)
        return LLKA_OK;

    for (size_t blockIdx = 0; blockIdx < cifData->nBlocks; blockIdx++) {
        auto &block = cifData->blocks[blockIdx];

        auto cat = block.firstCategory;
        while (cat != nullptr) {
            assert(cat->firstItem);

            auto item = cat->firstItem;
            while (item != nullptr) {
                if (cat->firstItem->nValues != item->nValues)
                    return LLKA_E_MISMATCHING_SIZES;
                item = item->p->next;
            }

            cat->p->isLoop = cat->firstItem->nValues > 1;
            cat = cat->p->next;
        }
    }

    cifData->p->tainted = false;

    return LLKA_OK;
}

LLKA_CifData * LLKA_CC LLKA_cifData_empty(void)
{
    auto data = new LLKA_CifData;
    data->blocks = nullptr;
    data->nBlocks = 0;
    data->p = new LLKA_CifDataPrivate;
    data->p->tainted = false;

    return data;
}

LLKA_Bool LLKA_CC LLKA_cifData_isTainted(const LLKA_CifData *data)
{
    return data->p->tainted ? LLKA_TRUE : LLKA_FALSE;
}

LLKA_CifDataCategory * LLKA_CC LLKA_cifDataBlock_addCategory(LLKA_CifDataBlock *block, const char *name)
{
    if (std::strlen(name) == 0)
        return nullptr;
    if (LLKA_cifDataBlock_findCategory(block, name) != nullptr)
        return nullptr;

    auto cat = block->firstCategory;
    while (cat != nullptr && cat->p->next != nullptr)
        cat = cat->p->next;

    auto newCat = new LLKA_CifDataCategory;
    newCat->name = LLKAInternal::duplicateString(name);
    newCat->firstItem = nullptr;
    newCat->p = new LLKA_CifDataCategoryPrivate;
    newCat->p->prev = cat;
    newCat->p->next = nullptr;
    newCat->p->root = block->p->root;
    newCat->p->isLoop = false;

    if (cat != nullptr)
        cat->p->next = newCat;
    if (block->firstCategory == nullptr)
        block->firstCategory = newCat;

    block->p->root->p->tainted = true;

    return newCat;
}

LLKA_RetCode LLKA_CC LLKA_cifDataBlock_deleteCategory(LLKA_CifDataBlock *block, const char *name)
{
    if (block->firstCategory == nullptr)
        return LLKA_E_INVALID_ARGUMENT;

    auto cat = LLKA_cifDataBlock_findCategory(block, name);
    if (cat == nullptr)
        return LLKA_E_INVALID_ARGUMENT;

    auto catPrev = cat->p->prev;
    auto catNext = cat->p->next;
    if (catPrev)
        catPrev->p->next = catNext;
    if (catNext)
        catNext->p->prev = catPrev;

    if (block->firstCategory == cat)
        block->firstCategory = catNext;

    LLKAInternal::MiniCif::destroyCifDataCategory(cat);
    delete cat;

    return LLKA_OK;
}

LLKA_CifDataCategory * LLKA_CC LLKA_cifDataBlock_findCategory(const LLKA_CifDataBlock *block, const char *name)
{
    auto cat = block->firstCategory;
    while (cat != nullptr) {
        if (std::strcmp(cat->name, name) == 0)
            return cat;
        cat = cat->p->next;
    }

    return nullptr;
}

LLKA_API LLKA_RetCode LLKA_CC LLKA_cifDataBlock_initialize(LLKA_CifDataBlock *block, const char *name, LLKA_CifData *data)
{
    if (block->name != nullptr || block->p != nullptr)
        return LLKA_E_INVALID_ARGUMENT;

    block->name = LLKAInternal::duplicateString(name);
    block->p = new LLKA_CifDataBlockPrivate;
    block->p->root = data;

    return LLKA_OK;
}

LLKA_CifDataCategory * LLKA_CC LLKA_cifDataBlock_nextCategory(const LLKA_CifDataCategory *cat)
{
    return cat->p->next;
}

LLKA_CifDataCategory * LLKA_CC LLKA_cifDataBlock_previousCategory(const LLKA_CifDataCategory *cat)
{
    return cat->p->prev;
}

LLKA_CifDataItem * LLKA_CC LLKA_cifDataCategory_addItem(LLKA_CifDataCategory *cat, const char *keyword)
{
    if (std::strlen(keyword) == 0)
        return nullptr;
    if (LLKA_cifDataCategory_findItem(cat, keyword) != nullptr)
        return nullptr;

    auto item = cat->firstItem;
    while (item!= nullptr && item->p->next != nullptr)
        item = item ->p->next;

    auto newItem = new LLKA_CifDataItem;
    newItem->keyword = LLKAInternal::duplicateString(keyword);
    newItem->values = nullptr;
    newItem->nValues = 0;
    newItem->p = new LLKA_CifDataItemPrivate;
    newItem->p->root = cat->p->root;
    newItem->p->prev = item;
    newItem->p->next = nullptr;

    if (item != nullptr)
        item->p->next = newItem;
    if (cat->firstItem == nullptr)
        cat->firstItem = newItem;

    cat->p->root->p->tainted = true;

    return newItem;
}

LLKA_RetCode LLKA_CC LLKA_cifDataCategory_deleteItem(LLKA_CifDataCategory *cat, const char *keyword)
{
    auto item = LLKA_cifDataCategory_findItem(cat, keyword);
    if (item == nullptr)
        return LLKA_E_INVALID_ARGUMENT;

    auto prevItem = item->p->prev;
    auto nextItem = item->p->next;
    if (prevItem)
        prevItem->p->next = nextItem;
    if (nextItem)
        nextItem->p->prev = prevItem;
    if (cat->firstItem == item)
        cat->firstItem = nextItem;

    LLKAInternal::MiniCif::destroyCifDataItem(item);
    delete item;

    return LLKA_OK;
}

LLKA_CifDataItem * LLKA_CC LLKA_cifDataCategory_findItem(const LLKA_CifDataCategory *cat, const char *keyword)
{
    auto e = cat->firstItem;
    while (e != nullptr) {
        if (std::strcmp(e->keyword, keyword) == 0)
            return e;
        e = e->p->next;
    }

    return nullptr;
}

LLKA_CifDataItem * LLKA_CC LLKA_cifDataCategory_nextItem(const LLKA_CifDataItem *item)
{
    return item->p->next;
}

LLKA_CifDataItem * LLKA_CC LLKA_cifDataCategory_previousItem(const LLKA_CifDataItem *item)
{
    return item->p->prev;
}

void LLKA_CC LLKA_cifDataItem_addValue(LLKA_CifDataItem *item, const LLKA_CifDataValue *value)
{
    auto newValues = new LLKA_CifDataValue[item->nValues + 1];
    if (item->nValues > 0)
        std::memcpy(newValues, item->values, sizeof(LLKA_CifDataValue) * item->nValues);

    auto nv = &newValues[item->nValues];
    nv->state = value->state;
    if (value->state == LLKA_MINICIF_VALUE_SET)
        nv->text = LLKAInternal::duplicateString(value->text);
    else
        nv->text = nullptr;

    delete [] item->values;
    item->values = newValues;
    item->nValues++;

    item->p->root->p->tainted = true;
}

void LLKA_CC LLKA_cifDataItem_addValues(LLKA_CifDataItem *item, const LLKA_CifDataValue *values, size_t nValues)
{
    if (nValues == 0)
        return;

    auto newValues = new LLKA_CifDataValue[item->nValues + nValues];
    if (item->nValues > 0)
        std::memcpy(newValues, item->values, sizeof(LLKA_CifDataValue) * item->nValues);

    for (size_t idx = 0; idx < nValues; idx++) {
        auto v = &values[idx];
        auto nv = &newValues[idx + item->nValues];
        nv->state = v->state;

        if (v->state == LLKA_MINICIF_VALUE_SET)
            nv->text = LLKAInternal::duplicateString(v->text);
        else
            nv->text = nullptr;
    }

    delete [] item->values;
    item->values = newValues;
    item->nValues += nValues;

    item->p->root->p->tainted = true;
}

LLKA_Bool LLKA_CC LLKA_cifDataCategory_isLoop(LLKA_CifDataCategory *cat)
{
    return cat->p->isLoop ? LLKA_TRUE : LLKA_FALSE;
}

void LLKA_CC LLKA_cifDataItem_setValues(LLKA_CifDataItem *item, const LLKA_CifDataValue *values, size_t nValues)
{
    LLKAInternal::MiniCif::destroyCifDataValues(item->values, item->nValues);
    delete [] item->values;

    item->values = new LLKA_CifDataValue[nValues];
    for (size_t idx = 0; idx < nValues; idx++) {
        auto v = &values[idx];
        auto nv = &item->values[idx];

        nv->state = v->state;
        if (v->state == LLKA_MINICIF_VALUE_SET)
            nv->text = LLKAInternal::duplicateString(v->text);
        else
            nv->text = nullptr;
    }

    item->nValues = nValues;

    item->p->root->p->tainted = true;
}

LLKA_RetCode LLKA_CC LLKA_cifDataToString(const LLKA_CifData *cifData, LLKA_Bool pretty, char **cifString)
{
    if (cifData->p->tainted)
        return LLKA_E_BAD_DATA;

    try {
        auto str = LLKAInternal::MiniCif::dataToString(*cifData, pretty == LLKA_TRUE);
        *cifString = new char[str.length() + 1];

        std::copy_n(str.c_str(), str.length(), *cifString);
        (*cifString)[str.length()] = '\0';

        return LLKA_OK;
    } catch (const LLKA_RetCode tRet) {
        return tRet;
    }
}

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
LLKA_RetCode LLKA_CC LLKA_cifFileToData(const LLKA_PathChar *path, LLKA_CifData **data, char **error)
{
    *error = nullptr;

#ifdef ENABLE_READ_UNIX_MMAP
    if (!std::filesystem::is_regular_file(path))
        return LLKA_E_NO_FILE;

    off_t len = std::filesystem::file_size(path);
    if (len < 1)
        return LLKA_E_NO_DATA;

    int fd = open(path, O_RDONLY);
    if (fd == -1)
        return LLKA_E_CANNOT_READ_FILE;

    void *mapped = mmap(0, len, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        return LLKA_E_CANNOT_READ_FILE;
    }
#elif defined ENABLE_WIN32_MMAP
    HANDLE hFile = CreateFileW(path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, 0, 0);
    if (hFile == INVALID_HANDLE_VALUE) {
        DWORD err = GetLastError();
        if (err == ERROR_FILE_NOT_FOUND)
            return LLKA_E_NO_FILE;
        return LLKA_E_CANNOT_READ_FILE;
    }

    LARGE_INTEGER fileSize;
    if (!GetFileSizeEx(hFile, &fileSize)) {
        CloseHandle(hFile);
        return LLKA_E_CANNOT_READ_FILE;
    }

    // Windows will not map files with zero size
    if (fileSize.QuadPart == 0LL) {
        CloseHandle(hFile);
        return LLKA_E_CANNOT_READ_FILE;
    }

    HANDLE hMapping = CreateFileMappingW(hFile, NULL, PAGE_READONLY, fileSize.HighPart, fileSize.LowPart, NULL);
    if (hMapping == NULL) {
        CloseHandle(hFile);
        return LLKA_E_CANNOT_READ_FILE;
    }

    void *mapped = MapViewOfFile(hFile, FILE_MAP_READ, 0, 0, 0);
    if (mapped == NULL) {
        CloseHandle(hMapping);
        CloseHandle(hFile);
        return LLKA_E_CANNOT_READ_FILE;
    }

    size_t len = fileSize.QuadPart;
#else
    (void)path; (void)data; (void)error;
    return LLKA_E_NOT_IMPLEMENTED;
#endif // ENABLE_*_MMAP

    // A bit of common code to read the mapped file
    auto view = std::string_view(static_cast<const char *>(mapped), len);

    LLKA_RetCode tRet;
    try {
        auto blocks = LLKAInternal::MiniCif::parse(view);
        *data = LLKAInternal::MiniCif::toCifData(blocks);

        tRet = LLKA_OK;
    } catch (const LLKAInternal::MiniCif::CifParseError &ex) {
        const auto len = std::strlen(ex.what());
        *error = new char[len + 1];
        std::strcpy(*error, ex.what());

        tRet = LLKA_E_BAD_DATA;
    }

#ifdef ENABLE_READ_UNIX_MMAP
    munmap(mapped, len);
    close(fd);
#elif defined ENABLE_WIN32_MMAP
    UnmapViewOfFile(mapped);
    CloseHandle(hMapping);
    CloseHandle(hFile);
#endif // ENABLE_*_MMAP

    return tRet;
}
#endif // LLKA_FILESYSTEM_ACCESS_DISABLE

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
LLKA_RetCode LLKA_CC LLKA_cifFileToStructure(const LLKA_PathChar *path, LLKA_ImportedStructure *importedStru, char **error, int32_t options)
{
    *error = nullptr;

#ifdef ENABLE_READ_UNIX_MMAP
    if (!std::filesystem::is_regular_file(path))
        return LLKA_E_NO_FILE;

    off_t len = std::filesystem::file_size(path);
    if (len < 1)
        return LLKA_E_NO_DATA;

    int fd = open(path, O_RDONLY);
    if (fd == -1)
        return LLKA_E_CANNOT_READ_FILE;

    void *mapped = mmap(0, len, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapped == MAP_FAILED) {
        close(fd);
        return LLKA_E_CANNOT_READ_FILE;
    }
#elif defined ENABLE_WIN32_MMAP
    HANDLE hFile = CreateFileW(path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, 0, 0);
    if (hFile == INVALID_HANDLE_VALUE) {
        DWORD err = GetLastError();
        if (err == ERROR_FILE_NOT_FOUND)
            return LLKA_E_NO_FILE;
        return LLKA_E_CANNOT_READ_FILE;
    }

    LARGE_INTEGER fileSize;
    if (!GetFileSizeEx(hFile, &fileSize)) {
        CloseHandle(hFile);
        return LLKA_E_CANNOT_READ_FILE;
    }

    // Windows will not map files with zero size
    if (fileSize.QuadPart == 0LL) {
        CloseHandle(hFile);
        return LLKA_E_CANNOT_READ_FILE;
    }

    HANDLE hMapping = CreateFileMappingW(hFile, NULL, PAGE_READONLY, fileSize.HighPart, fileSize.LowPart, NULL);
    if (hMapping == NULL) {
        CloseHandle(hFile);
        return LLKA_E_CANNOT_READ_FILE;
    }

    void *mapped = MapViewOfFile(hMapping, FILE_MAP_READ, 0, 0, 0);
    if (mapped == NULL) {
        CloseHandle(hMapping);
        CloseHandle(hFile);
        return LLKA_E_CANNOT_READ_FILE;
    }

    size_t len = fileSize.QuadPart;
#else
    (void)path; (void)importedStru; (void)error; (void)options;
    return LLKA_E_NOT_IMPLEMENTED;
#endif // ENABLE_READ_UNIX_MMAP

    // A bit of common code
    auto view = std::string_view(static_cast<const char *>(mapped), len);
    auto tRet = LLKAInternal::MiniCif::toStructure(view, importedStru, error, options);

#ifdef ENABLE_READ_UNIX_MMAP
    munmap(mapped, len);
    close(fd);
#elif defined ENABLE_WIN32_MMAP
    UnmapViewOfFile(mapped);
    CloseHandle(hMapping);
    CloseHandle(hFile);
#endif // ENABLE_READ_UNIX_MMAP

    return tRet;
}
#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

LLKA_RetCode LLKA_CC LLKA_cifTextToData(const char *text, LLKA_CifData **data, char **error)
{
    auto view = std::string_view(text);

    LLKA_RetCode tRet;
    try {
        auto blocks = LLKAInternal::MiniCif::parse(view);
        *data = LLKAInternal::MiniCif::toCifData(blocks);

        tRet = LLKA_OK;
    } catch (const LLKAInternal::MiniCif::CifParseError &ex) {
        const auto len = std::strlen(ex.what());
        *error = new char[len + 1];
        std::strcpy(*error, ex.what());

        tRet = LLKA_E_BAD_DATA;
    }

    return tRet;
}

LLKA_RetCode LLKA_CC LLKA_cifTextToStructure(const char *text, LLKA_ImportedStructure *importedStru, char **error, int32_t options)
{
    auto view = std::string_view(text);

    return LLKAInternal::MiniCif::toStructure(view, importedStru, error, options);
}

void LLKA_CC LLKA_destroyImportedStructure(LLKA_ImportedStructure *importedStru)
{
    LLKAInternal::destroyString(importedStru->entry.id);
    LLKA_destroyStructure(&importedStru->structure);
    LLKA_destroyCifData(importedStru->cifData);
}

void LLKA_CC LLKA_destroyCifData(LLKA_CifData *cifData)
{
    if (cifData == nullptr)
        return;

    for (size_t blockIdx = 0; blockIdx < cifData->nBlocks; blockIdx++) {
        auto block = &cifData->blocks[blockIdx];
        LLKAInternal::MiniCif::destroyCifDataBlock(block);
    }

    delete [] cifData->blocks;
    delete cifData->p;
    delete cifData;
}
