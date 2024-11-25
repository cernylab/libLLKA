/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_UTIL_CSV
#define _LLKA_UTIL_CSV

#include <llka_main.h>

#include "elementaries.h"
#include "templates.hpp"
#include "../fast_float/fast_float.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

/*
  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  v                                            v
  v    The girls just wanted to have fun...    v
  v                                            v
  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
*/

namespace LLKAInternal::Csv {

class TokenParseError : public std::runtime_error {
public:
    TokenParseError(const std::string &tag) : std::runtime_error{"Cannot parse token " + tag}
    {
    }
};

inline
auto tagStr(const std::string &tag) -> const std::string & { return tag; }

inline
auto tagStr(size_t tag) { return std::to_string(tag); }

template <typename T>
struct Convert {
    template <typename Tag>
    static auto call(const std::string &s, const Tag &) { return s; }
};

template <>
struct Convert<const char *> {
    template <typename Tag>
    static auto call(const std::string &s, const Tag &)
    {
        return duplicateString(s);
    }
};

template <>
struct Convert<double> {
    template <typename Tag>
    static auto call(const std::string &s, const Tag &tag)
    {
		double d;
        auto res = fast_float::from_chars(s.data(), s.data() + s.size(), d);
		if (res.ec != std::errc()) [[ unlikely ]]
			throw TokenParseError{tagStr(tag)};

		return d;
    }
};

template <>
struct Convert<float> {
    template <typename Tag>
    static auto call(const std::string &s, const Tag &tag)
    {
		float f;
        auto res = fast_float::from_chars(s.data(), s.data() + s.size(), f);
		if (res.ec != std::errc()) [[ unlikely ]]
			throw TokenParseError{tagStr(tag)};

		return f;
    }
};

template <>
struct Convert<int32_t> {
    template <typename Tag>
    static auto call(const std::string &s, const Tag &tag)
    {
        try {
            return std::stoi(s);
        } catch (const std::invalid_argument &) {
            throw TokenParseError{tagStr(tag)};
        } catch (const std::out_of_range &) {
            throw TokenParseError{tagStr(tag)};
        }
    }
};

template <StringLiteral Tag, typename ValueSetter>
struct NamedField {
    using Type = typename ValueSetter::Type;
    using Setter = ValueSetter;
    static constexpr auto tag = Tag.value;
};

template <size_t _Index, typename ValueSetter>
struct NumberedField {
    using Type = typename ValueSetter::Type;
    using Setter = ValueSetter;
    static constexpr auto Index = _Index;
    static constexpr auto tag = _Index;
};

enum class SchemaKind {
    Named,
    Numbered
};

template <SchemaKind _Kind, typename _Line, typename ..._Fields>
struct Schema {
    using Line = _Line;
    using Fields = std::tuple<_Fields...>;
    static constexpr size_t nFields = sizeof...(_Fields);
    static constexpr SchemaKind Kind = _Kind;
};

template <typename Line, typename Fields>
struct FieldsMapper {
private:
    static constexpr auto NFields = std::tuple_size<Fields>();

    template <size_t Index>
    static auto mapHeaderNamed(const std::vector<std::string> &tags, std::array<size_t, NFields> &mapping)
    {
        using E = std::tuple_element_t<Index, Fields>;

        for (size_t idx = 0; idx < tags.size(); idx++) {
            const auto &tag = tags[idx];
            if (tag == E::tag) {
                mapping[Index] = idx;
                if constexpr (Index == 0)
                    return mapping;
                else
                    return FieldsMapper::mapHeaderNamed<Index - 1>(tags, mapping);
            }
        }

        throw std::runtime_error{std::string{"Tag "} + E::tag + " not found"};
    }

    template <size_t Index>
    static auto mapHeaderNumbered(std::array<size_t, NFields> &mapping)
    {
        using E = std::tuple_element_t<Index, Fields>;
        mapping[Index] = E::Index;

        if constexpr (Index == 0)
            return mapping;
        else
            return mapHeaderNumbered<Index - 1>(mapping);
    }

    template <size_t Index>
    static auto mapLine(const std::vector<std::string> &toks, const std::array<size_t, NFields> &mapping, Line &line)
    {
        using E = std::tuple_element_t<Index, Fields>;
        const auto MappedIndex = mapping[Index];

        if (toks.size() < MappedIndex)
            throw LLKA_E_MISMATCHING_SIZES;

        auto v = Convert<typename E::Type>::call(toks[mapping[Index]], E::tag);
        E::Setter::set(line, std::move(v));
        if constexpr (Index == 0)
            return line;
        else
            return mapLine<Index - 1>(toks, mapping, line);
    }

public:
    template <SchemaKind Kind>
    static auto mapHeader(const std::vector<std::string> &tags = {})
    {
        std::array<size_t, NFields> mapping{};
        if constexpr (Kind == SchemaKind::Named)
            return mapHeaderNamed<NFields - 1>(tags, mapping);
        else
            return mapHeaderNumbered<NFields - 1>(mapping);
    }

    static auto mapLine(const std::vector<std::string> &toks, const std::array<size_t, NFields> &mapping)
    {
        Line line{};
        return mapLine<NFields - 1>(toks, mapping, line);
    }
};

template <typename Schema, typename Input>
static
auto getMapping(Input &&strm, const char delim)
{
    if constexpr (Schema::Kind == SchemaKind::Named) {
        std::string buf{};

        std::getline(strm, buf);
        auto tags = split(buf, delim, Schema::nFields);
        if (tags.size() < Schema::nFields)
            throw LLKA_E_MISMATCHING_SIZES;

        return FieldsMapper<typename Schema::Line, typename Schema::Fields>::template mapHeader<Schema::Kind>(tags);
    } else if (Schema::Kind == SchemaKind::Numbered) {
        return FieldsMapper<typename Schema::Line, typename Schema::Fields>::template mapHeader<Schema::Kind>();
    }
}

inline
auto splitLine(const std::string &line, const size_t nToks, const char delim, const char quote)
{
    if (quote == 0)
        return split(line, delim, nToks);

    std::vector<std::string> toks{};
    toks.reserve(nToks);

    size_t idx = 0;
    while (idx < line.length()) {
        const auto ch = line[idx];
        if (ch == quote) {
            auto next = line.find_first_of(quote, idx + 1);
            if (next == std::string::npos)
                throw LLKA_E_BAD_DATA;    // Unterminated quote
            else if ((next < line.length() - 1) && (line[next + 1] != delim))
                throw LLKA_E_BAD_DATA;    // Field separator does not follow the end quote

            toks.emplace_back(line.substr(idx + 1, next - idx - 1));
            idx = next + 2;
        } else {
            const auto next = line.find_first_of(delim, idx);

            toks.emplace_back(line.substr(idx, next - idx));

            if (next == std::string::npos)
                break;
            idx = next + 1;
        }
    }

    return toks;
}

template <typename Schema, typename Input>
static
auto parse(Input &&strm, const char delim, const char quote)
{
    using FM = FieldsMapper<typename Schema::Line, typename Schema::Fields>;
    using Lines = std::vector<typename Schema::Line>;
    constexpr auto NToks = Schema::nFields;
    std::vector<std::string> warnings{};

    Lines lines{};

    try {
        std::string buf{};

        const auto mapping = getMapping<Schema>(strm, delim);

        size_t lineCounter = 2;
        while (strm) {
            std::getline(strm, buf);

            if (buf.empty())
                break;

            auto toks = splitLine(buf, NToks, delim, quote);

            try {
                lines.emplace_back(FM::mapLine(toks, mapping));
            } catch (const TokenParseError &ex) {
                std::string warn = "Warning on line " + std::to_string(lineCounter) + ": " + ex.what() + ", skipping that line";
                warnings.push_back(std::move(warn));
            }

            lineCounter++;
        }
    } catch (const std::runtime_error &) {
        throw LLKA_E_MISMATCHING_DATA;
    }

    return std::make_tuple(lines, warnings);
}

template <typename Schema>
static
auto parse(const std::filesystem::path &path, const char delim, const char quote)
{
    if (!std::filesystem::is_regular_file(path))
        throw LLKA_E_NO_FILE;

    std::ifstream ifs{path};
    if (!ifs)
        throw LLKA_E_CANNOT_READ_FILE;

    return parse<Schema>(ifs, delim, quote);
}

template <typename Schema>
static
auto parse(const std::string &text, const char delim, const char quote)
{
    return parse<Schema>(std::istringstream{text}, delim, quote);
}

} // namespace LLKAInternal::Csv

#endif // _LLKA_UTIL_CSV
