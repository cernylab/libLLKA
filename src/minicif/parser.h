/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_MINICIF_PARSER_H
#define _LLKA_MINICIF_PARSER_H

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace LLKAInternal::MiniCif {

class CifParseError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

class Value {
public:
    // Set the values explicitly so that we can safely cast them to LLKA_CifDataValueState
    enum class State : int32_t {
        VALUE = 0,
        NONE  = 1,
        UNKW  = 2
    };

    Value() noexcept;
    Value(const std::string_view &sv) noexcept;
    Value(std::string value) noexcept;
    Value(const Value &other);
    Value(Value &&other) noexcept;

    Value & operator=(const Value &other) noexcept;
    Value & operator=(Value &&other) noexcept;

    std::string text;
    State state; // This is the same thing as LLKA_CifDataValueState but
                 // we use a different enum here to avoid having to pull in
                 // the entire public API header.

    static Value NoneValue() noexcept;
    static Value UnkwValue() noexcept;

private:
    Value(State state) noexcept;
};
using Values = std::vector<Value>;

class Item {
public:
    Item(std::string keyword, std::string lowecaseKeyword, Values values) noexcept;

    std::string keyword;
    std::string lowercaseKeyword;
    Values values;
};
using Items = std::vector<Item>;

class Category {
public:
    Category(std::string name, std::string lowecaseName, Items items) noexcept;

    std::string name;
    std::string lowecaseName;
    Items items;
};

class Block {
public:
    explicit Block(std::string _name);

    auto add(std::string category, std::string keyword, Value value) -> void;
    auto addMultiple(std::string category, std::string keyword, Values values) -> void;
    auto nextAnonymousCategoryName() -> std::string;

    std::vector<Category> categories;
    std::string name;

private:
    size_t m_anonymousCategoriesCount;
};

auto parse(const std::string_view &data) -> std::vector<Block>;

} // namespace LLKAInternal::MiniCif

#endif // _LLKA_MINICIF_PARSER_H
