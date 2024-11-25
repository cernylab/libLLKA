/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "writer.h"

#include "minicif_p.h"
#include "../fast_float/fast_float.h"

#include <array>
#include <cstring>
#include <sstream>
#include <string_view>
#include <tuple>
#include <vector>

namespace LLKAInternal::MiniCif {

inline constinit std::array<char, 6> WHITESPACES{ '\x20', '\x09', '\xa0', '\x0b', '\x0c', '\x0d' };

static
auto valueToCifFormat(const LLKA_CifDataValue &value) -> std::tuple<std::string, bool>
{
    if (value.state == LLKA_MINICIF_VALUE_NONE)
        return { ".", false };
    else if (value.state == LLKA_MINICIF_VALUE_UNKW)
        return { "?", false };

    std::string_view sv{value.text};
    // Empty value
    if (sv.empty())
        return { "\"\"", false };  // Rather odd case of an empty string that is not an empty value.
                                   // We need to represent this as a quoted string with no content.

    if (sv.find_first_of('\n') != std::string_view::npos) {
        /* Value contains newline characters - treat it as a multiline value */
        return {
            ";" + std::string{sv} + "\n;",
            true
        };
    }

    if (sv.find_first_of('"') != std::string_view::npos) {
        // Value contains double quotes, we may need to quote it
        if (sv.front() == '"' && sv.back() == '"' && sv.length() > 1) // Already quoted, just copy-return
            return { std::string{sv}, false };
        if (sv.find_first_of('\'') != std::string_view::npos) {
            // Value contains both single and double quotes, convert it to a multiline value
            return {
                ";" + std::string{sv} + "\n;",
                true
            };
        }
        return { "'" + std::string{sv} + "'", false }; // Use single quotes to quote
    }

    if (sv.find_first_of('\'') != std::string_view::npos) {
        // Value contains single quotes, we may need to quote it
        if (sv.front() == '\'' && sv.back() == '\'' && sv.length() > 1) // Already quoted, just copy-return
            return { std::string{sv}, false };
        if (sv.find_first_of('"') != std::string_view::npos) {
            // Value contains both single and double quotes, convert it to a multiline value
            return {
                ";" + std::string{sv} + "\n;",
                true
            };
        }
        return { "\"" + std::string{sv} + "\"", false }; // Use double quotes to quote
    }

    bool whitespaces = false;
    for (const auto &ch : sv) {
        for (const auto ws : WHITESPACES) {
            whitespaces = ws == ch;
            if (whitespaces)
                break;
        }
        if (whitespaces)
            break;
    }
    if (whitespaces) {
        /* Value contains whitespace, quote it with double quotes*/
        return { "\"" + std::string{sv} + "\"", false };
    }

    return { std::string{sv}, false };
}

static
auto writePrettyLoop(const LLKA_CifDataCategory *cat, std::string &cifStr)
{
    struct ColumnFormat {
        size_t maxLength;
        bool padOnLeft;
    };
    std::vector<ColumnFormat> columnFormats{};
    std::vector<std::vector<std::tuple<std::string, bool>>> columnsValues{};

    auto item = cat->firstItem;
    if (item == nullptr)
        throw LLKA_E_BAD_DATA;

    // Open the loop
    cifStr += "loop_\n";

    size_t nRows = item->nValues;
    size_t textLen = 0;
    while (item != nullptr) {
        // Write out the loop header
        cifStr += "_" + std::string{cat->name} + "." + std::string{item->keyword} + "\n";

        // Create Cif values and calculate formatting
        size_t maxLength = 0;
        bool canBeNumeric = true;
        std::vector<std::tuple<std::string, bool>> colValues(item->nValues);
        for (size_t idx = 0; idx < item->nValues; idx++) {
            auto cifValue = valueToCifFormat(item->values[idx]);
            auto len = std::get<0>(cifValue).length();
            if (len > maxLength)
                maxLength = len;

            if (canBeNumeric && item->values[idx].state == LLKA_MINICIF_VALUE_SET) {
                float f;
                const auto &text = std::get<0>(cifValue);
		        auto ret = fast_float::from_chars(text.data(), text.data() + text.size(), f);
		        canBeNumeric = ret.ec == std::errc();
            }

            textLen += std::get<0>(cifValue).length();
            colValues[idx] = std::move(cifValue);
        }

        columnFormats.push_back({.maxLength = maxLength, .padOnLeft = canBeNumeric});
        columnsValues.push_back(std::move(colValues));

        item = item->p->next;
    }

    // Ensure that we have enough space to store the text
    size_t toReserve = (125 * textLen) / 100;
    if (cifStr.capacity() < cifStr.length() + toReserve)
        cifStr.reserve(cifStr.length() + toReserve);

    for (size_t rowIdx = 0; rowIdx < nRows; rowIdx++) {
        bool previousValueWasMultiline = false;

        for (size_t colIdx = 0; colIdx < columnFormats.size(); colIdx++) {
            const auto &format = columnFormats[colIdx];
            const auto &[ text, isMultiline ] = columnsValues[colIdx][rowIdx];

            if (isMultiline) {
                cifStr
                    += (previousValueWasMultiline ? "" : "\n")
                    +  text + "\n";
                previousValueWasMultiline = true;
            } else {
                std::string padding(format.maxLength - text.length(), ' ');

                if (format.padOnLeft)
                    cifStr += padding + text + "  ";
                else
                    cifStr += text + padding + "  ";

                previousValueWasMultiline = false;
            }
        }
        if (!previousValueWasMultiline)
            cifStr += "\n";
    }
    cifStr += "##\n";
}

static
auto writeLoop(const LLKA_CifDataCategory *cat, std::string &cifStr)
{
    auto entry = cat->firstItem;
    if (entry == nullptr)
        throw LLKA_E_BAD_DATA;

    cifStr += "loop_\n";

    // Write out the loop header first
    while (entry != nullptr) {
        cifStr += "_" + std::string{cat->name} + "." + std::string{entry->keyword} + "\n";
        entry = entry->p->next;
    }

    const size_t nRows = cat->firstItem->nValues;

    if (cifStr.capacity() < cifStr.length() + 100)
        cifStr.reserve(cifStr.length() + 100); // Reserve space for additional 100 characters

    // Transpose contents from columns to rows
    for (size_t rowIdx = 0; rowIdx < nRows; rowIdx++) {
        bool previousValueWasMultiline = false;

        entry = cat->firstItem;
        while (entry != nullptr) {
            if (entry->nValues < nRows)
                throw LLKA_E_BAD_DATA;

            auto [ text, isMultiline ] = valueToCifFormat(entry->values[rowIdx]);
            if (isMultiline) {
                cifStr
                    += (previousValueWasMultiline ? "" : "\n")
                    +  text + "\n";
                previousValueWasMultiline = true;
            } else {
                cifStr += text + " ";
                previousValueWasMultiline = false;
            }

            entry = entry->p->next;
            if (cifStr.capacity() < cifStr.length() + 100)
                cifStr.reserve(cifStr.capacity() + 100); // Reserve space for additional 100 characters
        }
        if (!previousValueWasMultiline)
            cifStr += "\n";
    }
    cifStr += "##\n";
}

static
auto writePrettySingles(const LLKA_CifDataCategory *cat, std::string &cifStr)
{
    std::vector<std::tuple<std::string, std::tuple<std::string, bool>>> items{};
    size_t maxNameLength = 0;
    size_t entriesTextLen = 0;
    size_t nEntries = 0;

    auto entry = cat->firstItem;
    while (entry != nullptr) {
        if (entry->nValues == 0)
            throw LLKA_E_BAD_DATA;

        auto name = std::string{"_"} + cat->name + "." + entry->keyword;
        if (name.length() > maxNameLength)
            maxNameLength = name.length();

        auto cifValue = valueToCifFormat(entry->values[0]);
        items.emplace_back(std::move(name), std::move(cifValue));
        entriesTextLen += std::get<0>(cifValue).length();

        entry = entry->p->next;
        nEntries++;
    }

    size_t toReserve = (maxNameLength + 1) * nEntries + entriesTextLen + 3;
    if (cifStr.capacity() < cifStr.length() + toReserve)
        cifStr.reserve(cifStr.length() + toReserve);

    for (const auto &[ name, cifValue ] : items) {
        std::string padding(maxNameLength - name.length() + 2, ' ');

        cifStr += name + padding;
        if (std::get<1>(cifValue))
            cifStr += "\n" + std::get<0>(cifValue);
        else
            cifStr += std::get<0>(cifValue);
        cifStr += "\n";
    }
    cifStr += "##\n";
}

static
auto writeSingles(const LLKA_CifDataCategory *cat, std::string &cifStr)
{
    auto catName = std::string{cat->name};

    auto entry = cat->firstItem;
    while (entry != nullptr) {
        if (entry->nValues == 0)
            throw LLKA_E_BAD_DATA;

        auto [ text, isMultiline ] = valueToCifFormat(entry->values[0]);
        cifStr += "_" + std::string{cat->name} + "." + std::string{entry->keyword} + " ";
        if (isMultiline)
            cifStr += "\n" + text;
        else
            cifStr += text;
        cifStr += "\n";

        entry = entry->p->next;
    }
    cifStr += "##\n";
}

auto dataToString(const LLKA_CifData &cifData, const bool pretty) -> std::string
{
    std::string cifStr{};

    for (size_t blockIdx = 0; blockIdx < cifData.nBlocks; blockIdx++) {
        const auto &block = cifData.blocks[blockIdx];

        // Open data block
        cifStr += "data_" + std::string{block.name} + "\n#\n";

        auto cat = block.firstCategory;
        while (cat != nullptr) {
            if (!cat->p->isLoop)
                pretty ? writePrettySingles(cat, cifStr) : writeSingles(cat, cifStr);
            else
                pretty ? writePrettyLoop(cat, cifStr) : writeLoop(cat, cifStr);

            cat = cat->p->next;
        }
        cifStr += "#\n";
    }

    return cifStr;
}

} // namespace LLKAInternal::MiniCif
