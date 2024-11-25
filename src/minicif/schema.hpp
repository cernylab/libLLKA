/* vim: set sw=4 ts=4 sts=4 : */

#ifndef _LLKA_MINICIF_SCHEMA_H
#define _LLKA_MINICIF_SCHEMA_H

#include "parser.h"
#include "../fast_float/fast_float.h"

#include "../util/elementaries.h"
#include "../util/templates.hpp"

#include <array>
#include <clocale>
#include <memory>
#include <sstream>

namespace LLKAInternal::MiniCif {

class CifSchemaError : public std::runtime_error {
	using std::runtime_error::runtime_error;
};

class InvalidValueError : public CifSchemaError {
public:
	InvalidValueError() : CifSchemaError{"Invalid value"}
	{
	}
};

template <typename T>
struct Convert;

template <>
struct Convert<char> {
	static auto call(const std::string &s)
	{
		if (s.empty()) [[unlikely]]
			throw InvalidValueError{};
		return s[0];
	}
};

template <>
struct Convert<const char *> {
	static auto call(const std::string &s)
	{
		const auto ds = dequote(s);
		return duplicateString(ds);
	}
};

template <>
struct Convert<double> {
	static auto call(const std::string &s)
	{
		double d;
		auto ret = fast_float::from_chars(s.data(), s.data() + s.size(), d);
		if (ret.ec != std::errc()) [[ unlikely ]]
			throw InvalidValueError{};

		return d;
	}
};

template <>
struct Convert<float> {
	static auto call(const std::string &s)
	{
		double f;
		auto ret = fast_float::from_chars(s.data(), s.data() + s.size(), f);
		if (ret.ec != std::errc()) [[ unlikely ]]
			throw InvalidValueError{};

		return f;
	}
};

template <>
struct Convert<int32_t> {
	static auto call(const std::string &s)
	{
		try {
			return std::stoi(s);
		} catch (const std::invalid_argument &) {
			throw InvalidValueError{};
		} catch (const std::out_of_range &) {
			throw InvalidValueError{};
		}
	}
};

template <>
struct Convert<uint32_t> {
	static auto call(const std::string &s)
	{
		try {
			return std::stoul(s);
		} catch (const std::invalid_argument &) {
			throw InvalidValueError{};
		} catch (const std::out_of_range &) {
			throw InvalidValueError{};
		}
	}
};

struct NoDefaultValue {
	using Type = void;
	static constexpr bool hasDefault = false;
};

template <typename T, T Value>
struct DefaultValue {
	using Type = T;
	static constexpr bool hasDefault = true;
	static constexpr auto value = Value;
};

template <StringLiteral Tag, typename ValueSetter, typename _Default = NoDefaultValue>
struct TagValue {
	using Type = typename ValueSetter::Type;
	using Setter  = ValueSetter;
	using Default = _Default;
	static constexpr auto tag = Tag.value;
};

template <StringLiteral Name, typename _Image, typename ..._Fields>
struct Schema {
	using Fields = std::tuple<_Fields...>;
	using Image = _Image;
	static constexpr auto name = Name.value;
	static constexpr auto nFields = sizeof...(_Fields);
};

template <typename T>
inline
auto NoopFixup(const T &) { /* NOOP */ }

template <typename T>
inline
auto NoopRelease(const T &) { /* NOOP */ }

template <typename Schema>
struct Applier {
private:
	using Mapping = std::array<const Values *, Schema::nFields>;

	template <size_t Index>
	static auto apply(const Mapping &mapping, const size_t row, typename Schema::Image &image)
	{
		using E = std::tuple_element_t<Index, typename Schema::Fields>;
		const auto valuesPtr = mapping[Index];

		if (valuesPtr == nullptr) {
			if constexpr (E::Default::hasDefault) {
				E::Setter::set(image, E::Default::value);
			} else {
				// We cannot hit this case here but we need to check anyway becase the template code does not know that this cannot happen
				throw CifSchemaError{"Category does not contain field " + std::string{E::tag}};
			}
		} else {
			const auto &value = (*valuesPtr)[row];

			if (value.state != Value::State::VALUE) {
				if constexpr (E::Default::hasDefault) {
					E::Setter::set(image, E::Default::value);
				} else {
					throw CifSchemaError{"Item in field " + std::string{E::tag} + " has empty value but does not have a default value to use instead"};
				}
			} else {
				try {
					E::Setter::set(image, Convert<typename E::Type>::call(value.text));
				} catch (const InvalidValueError &ex) {
					throw CifSchemaError{"Cannot set value for field " + std::string{E::tag} + ", " + ex.what() + " " + value.text};
				}
			}
		}

		if constexpr (Index > 0) {
			apply<Index - 1>(mapping, row, image);
		}
	}

	static auto getNumerOfRows(const Mapping &mapping) -> size_t
	{
		for (const auto &valuesPtr : mapping) {
			if (valuesPtr == nullptr)
				continue;
			return valuesPtr->size();
		}

		return 0UL;
	}

	template <size_t Index>
	static auto mapColumn(const Category &cat, Mapping &mapping)
	{
		using E = std::tuple_element_t<Index, typename Schema::Fields>;
		const auto fieldIt = std::find_if(cat.items.cbegin(), cat.items.cend(), [](const Item &im) { return im.lowercaseKeyword == E::tag; });

		if (fieldIt == cat.items.cend()) {
			if constexpr (E::Default::hasDefault) {
				mapping[Index] = nullptr;
			} else {
				throw CifSchemaError{"Category does not contain field " + std::string{E::tag}};
			}
		} else
			mapping[Index] = &fieldIt->values;

		if constexpr (Index > 0) {
			return mapColumn<Index - 1>(cat, mapping);
		} else {
			return mapping;
		}
	}

	static auto makeMapping(const Category &cat)
	{
		Mapping mapping{};

		return mapColumn<Schema::nFields - 1>(cat, mapping);
	}

public:
	template <typename FixerUpper, typename Releaser>
	static auto apply(const Block &block, FixerUpper fixerUpper = NoopFixup<typename Schema::Image>, Releaser releaseItem = NoopRelease<typename Schema::Image>)
	{
		auto catIt = std::find_if(block.categories.cbegin(), block.categories.cend(), [](const Category &c) { return c.lowecaseName == Schema::name; });
		if (catIt == block.categories.cend())
			throw CifSchemaError{"Category " + std::string{Schema::name} + " is not present in block " + block.name};

		const auto &cat = *catIt;

		if (cat.items.empty())
			throw CifSchemaError{"Empty category"};

		const auto mapping = makeMapping(cat);

		const auto NRows = getNumerOfRows(mapping);
		if (NRows == 0)
			throw CifSchemaError{"Category contains no values"};

		auto items = std::unique_ptr<typename Schema::Image[]>(new typename Schema::Image[NRows]);

		size_t row = 0;
		try {
			for (; row < NRows; row++) {
				apply<Schema::nFields - 1>(mapping, row, items[row]);
				fixerUpper(items[row]);
			}

			return std::make_tuple(std::move(items), NRows);
		} catch (const CifSchemaError &ex) {
			// Unwind
			for (size_t uwdx = 0; uwdx < row; uwdx++)
				releaseItem(items[uwdx]);

			throw ex;
		}
	}
};

} // namespace LLKAInternal::MiniCif

#endif // _LLKA_MINICIF_SCHEMA_H
