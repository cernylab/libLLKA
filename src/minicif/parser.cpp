/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "parser.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

namespace LLKAInternal::MiniCif {

Value::Value() noexcept :
    state{State::NONE}
{
}

Value::Value(State state) noexcept :
    state{state}
{
}

Value::Value(const std::string_view &sv) noexcept :
    text{sv},
    state{State::VALUE}
{
}

Value::Value(std::string text) noexcept :
    text{std::move(text)},
    state{State::VALUE}
{
}

Value::Value(const Value &other) :
    text{other.text},
    state{other.state}
{
}

Value::Value(Value &&other) noexcept :
    text{std::move(other.text)},
    state{other.state}
{
}

Value & Value::operator=(const Value &other) noexcept
{
    text = other.text;
    state = other.state;

    return *this;
}

Value & Value::operator=(Value &&other) noexcept
{
    text = std::move(other.text);
    state = other.state;

    return *this;
}

Value Value::NoneValue() noexcept
{
    return Value{State::NONE};
}

Value Value::UnkwValue() noexcept
{
    return Value{State::UNKW};
}

Item::Item(std::string keyword, std::string lowercaseKeyword, Values values) noexcept :
    keyword{std::move(keyword)},
    lowercaseKeyword{std::move(lowercaseKeyword)},
    values{std::move(values)}
{
}

Category::Category(std::string name, std::string lowercaseName, Items items) noexcept :
    name{std::move(name)},
    lowecaseName{std::move(lowercaseName)},
    items{std::move(items)}
{
}

static
auto toLowerCase(std::string &s)
{
    std::string lwr(s.length(), ' ');
    // PERF: Use SIMD here
    std::transform(s.begin(), s.end(), lwr.begin(), [](char ch) { return ::tolower(ch); });

    return lwr;
}

static
auto toValue(const std::string_view &sv) -> Value
{
    if (sv.length() == 1) {
        const auto ch = sv[0];
        if (ch == '.')
            return Value::NoneValue();
        else if (ch == '?')
            return Value::UnkwValue();
    }
    return Value{sv};
}

Block::Block(std::string _name) :
    name{std::move(_name)},
    m_anonymousCategoriesCount{0UL}
{
}

auto Block::add(std::string category, std::string keyword, Value value) -> void
{
    addMultiple(std::move(category), std::move(keyword), { std::move(value) });
}

auto Block::addMultiple(std::string category, std::string keyword, Values values) -> void
{
    if (category.empty())
        category = nextAnonymousCategoryName();

    auto lowercaseCategory = toLowerCase(category);
    auto lowercaseKeyword = toLowerCase(keyword);

    auto cat = std::find_if(categories.begin(), categories.end(), [&lowercaseCategory](const Category &c) { return c.lowecaseName == lowercaseCategory; });
    if (cat == categories.end()) {
        categories.emplace_back(std::move(category), std::move(lowercaseCategory), Items{ { std::move(keyword), std::move(lowercaseKeyword), std::move(values) } });
    } else {
        auto item = std::find_if(cat->items.begin(), cat->items.end(), [&lowercaseKeyword](const Item &im) { return im.lowercaseKeyword == lowercaseKeyword; });
        if (item != cat->items.end()) [[ unlikely ]]
            throw CifParseError{"Duplicit item " + keyword};

        cat->items.emplace_back(std::move(keyword), std::move(lowercaseKeyword), std::move(values));
    }
}

auto Block::nextAnonymousCategoryName() -> std::string
{
    return "anonymous_" + std::to_string(m_anonymousCategoriesCount++);
}

class PrimingState {
public:
    PrimingState() :
        primed{false},
        quotesIgnored{false}
    {
    }

    auto isPrimed(bool ignoreQuotes)
    {
        return primed && quotesIgnored == ignoreQuotes;
    }

    bool primed;
    bool quotesIgnored;
};

class Token {
public:
    enum class Kind : int {
        VALUE        = 0,
        TAG          = (1 << 0),
        COMMENT      = (1 << 1),
        MULTILINE    = (1 << 2),
        LOOP         = (1 << 3),
        DATA_BLOCK   = (1 << 4),
        SAVE_BLOCK   = (1 << 5),
        STOP         = (1 << 6),
        GLOBAL_BLOCK = (1 << 7),
        EMPTY        = (1 << 8)
    };

    Token(std::string_view _text, Kind _kind) :
        text{std::move(_text)},
        kind{_kind}
    {
    }

    std::string_view text;
    Kind kind;
};

class Stream {
public:
    static constexpr char DOUBLE_QUOTE = '"';
    static constexpr char HASH         = '#';
    static constexpr char NEW_LINE     = '\n';
    static constexpr char CARRIAGE_RET = '\r';
    static constexpr char SINGLE_QUOTE = '\'';
    static constexpr char SEMICOLON    = ';';
    static constexpr char SPACE        = ' ';
    static constexpr char TAB          = '\t';
    static constexpr char UNDERSCORE   = '_';

    Stream(const std::string_view &stream) :
        lineCounter{1},
        m_stream{stream},
        m_cursor{0UL},
        m_length{stream.length()}
    {
    }

    auto eat(bool ignoreQuotes = false)
    {
        auto token = [this, ignoreQuotes]() {
            if (m_priming.isPrimed(ignoreQuotes)) {
                m_priming.primed = false;
                return Token{m_primedText, m_primedKind};
            } else {
                auto text = getCifToken(ignoreQuotes);
                if (text.empty()) [[unlikely]]
                    return Token{"", Token::Kind::EMPTY};
                auto kind = tokenKind(text);
                return Token{text, kind};
            }
        }();

        m_cursor += token.text.length();
        skipWhitespaces();

        return token;
    }

    auto eatLine() -> std::string_view
    {
        m_priming.primed = false;

        auto end = std::string::npos;
        for (size_t idx = m_cursor; idx < m_length; idx++) {
            if (m_stream[idx] == NEW_LINE) {
                end = idx;
                lineCounter++;
                break;
            }
        }

        // If the file has \r\n line endings, this function will not remove the trailing \r.
        // This is okay because most callers of this function discard its output and trying
        // to fix this up would just slow things down.
        // If the caller wants to use the output, they are expected to handle the trailing \r
        // by themselves.

        if (end == std::string::npos) [[unlikely]] {
            auto line = m_stream.substr(m_cursor);
            m_cursor = m_length;
            return line;
        } else {
            auto line = m_stream.substr(m_cursor, end - m_cursor);
            m_cursor = end + 1;
            return line;
        }
    }

    auto exhausted() const { return m_cursor == m_length; }

    auto peekKind(bool ignoreQuotes = false)
    {
        if (m_priming.isPrimed(ignoreQuotes))
            return m_primedKind;

        skipWhitespaces();
        if (exhausted()) [[ unlikely ]]
            return Token::Kind::EMPTY;

        m_primedText = getCifToken(ignoreQuotes);
        m_primedKind = tokenKind(m_primedText);

        m_priming.primed = true;
        m_priming.quotesIgnored = ignoreQuotes;

        return m_primedKind;
    }

    size_t lineCounter;

private:
    auto getCifToken(bool ignoreQuotes) const -> std::string_view
    {
        auto ch = m_stream[m_cursor];
        const auto quote = (!ignoreQuotes && (ch == SINGLE_QUOTE || ch == DOUBLE_QUOTE)) * ch;    // Get quote char or set to zero if current char is not a quote
        bool unterminatedQuote = quote > 0;

        const auto from = m_cursor;
        auto idx = from + 1;
        for (; idx < m_length; idx++) {
            ch = m_stream[idx];
            if (quote) {
                if (ch == NEW_LINE || ch == CARRIAGE_RET)
                    throw CifParseError{"Quoted token that begins on line " + std::to_string(lineCounter) + " contains a new line"};
                else if (ch == quote) {
                    if (idx + 1 < m_length) [[ likely ]] {
                        auto nextCh = m_stream[idx + 1];

                        if ((nextCh == TAB) || (nextCh == NEW_LINE) || (nextCh == CARRIAGE_RET) || (nextCh == SPACE)) {
                            // This should be ending quote, jump one character ahead to get past the ending quote and break
                            // Note that we cannot dequote the string here because CIF names are allowed to be expressed as "Value"
                            // if they are quoted. Dequoting the string would therefore confuse the token kind detection.
                            unterminatedQuote = false;
                            idx++;
                            break;
                        } else {
                            // Continue parsing if the next character is not whitespace
                            // see points 15 and/or 56 on https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax
                            continue;
                        }
                    } else {
                        unterminatedQuote = false;
                        idx++;
                        break;
                    }
                }
            } else if ((ch == TAB) || (ch == NEW_LINE) || (ch == CARRIAGE_RET) || (ch == SPACE))
                break;
        }

        if (unterminatedQuote && idx == m_length)
            throw CifParseError{"Unterminated quoted token that begins on line " + std::to_string(lineCounter)};
        return m_stream.substr(from, idx - from);
    }

    auto isEmptyToken(const std::string_view &text) const
    {
        for (const char ch : text) {
            bool isWhiteSpace =
                (ch == TAB) || (ch == NEW_LINE) || (ch == CARRIAGE_RET) || (ch == SPACE);
            if (!isWhiteSpace)
                return false;
        }

        return true;
    }

    auto skipWhitespaces() -> void
    {
        while (m_cursor < m_length) {
            const auto ch = m_stream[m_cursor];

            const size_t isNewLine = ch == NEW_LINE;
            const bool isWhiteSpace = isNewLine || (ch == TAB) || (ch == CARRIAGE_RET) || (ch == SPACE);
            if (!isWhiteSpace)
                break;

            lineCounter += isNewLine;
            m_cursor++;
        }
    }

    auto tokenKind(const std::string_view &text) const -> Token::Kind
    {
        assert(!text.empty());

        const auto ch = text[0];
        const auto length = text.length();
        const auto data = text.data();

        int kind =
            ((ch == UNDERSCORE)) | // Tag?
            ((ch == HASH) << 1)  | // Comment?
            ((ch == SEMICOLON && (m_cursor == 0ULL || m_stream[m_cursor - 1] == NEW_LINE || m_stream[m_cursor - 1] == CARRIAGE_RET)) << 2); // Multiline?

        if (kind)
            return static_cast<Token::Kind>(kind);
        if (length < 5)
            return isEmptyToken(text) ? Token::Kind::EMPTY : Token::Kind::VALUE;

        if (std::memcmp(data, LOOP_.data(), 5) == 0)
            return Token::Kind::LOOP;
        else if (std::memcmp(data, DATA_.data(), 5) == 0)
            return Token::Kind::DATA_BLOCK;
        else if (std::memcmp(data, SAVE_.data(), 5) == 0)
            return Token::Kind::SAVE_BLOCK;
        else if (std::memcmp(data, STOP_.data(), 5) == 0)
            return Token::Kind::STOP;
        else if (length >= 7 && std::memcmp(data, GLOBAL_.data(), 7) == 0)
            return Token::Kind::GLOBAL_BLOCK;

        return isEmptyToken(text) ? Token::Kind::EMPTY : Token::Kind::VALUE;
    }

    const std::string_view &m_stream;
    size_t m_cursor;
    const size_t m_length;

    std::string_view m_primedText;
    Token::Kind m_primedKind;
    PrimingState m_priming;

    static constexpr std::string_view DATA_  {"data_"};
    static constexpr std::string_view GLOBAL_{"global_"};
    static constexpr std::string_view LOOP_  {"loop_"};
    static constexpr std::string_view SAVE_  {"save_"};
    static constexpr std::string_view STOP_  {"stop_"};
};

static
auto splitOnFirst(const std::string_view &str, const char delim) -> std::tuple<std::string_view, std::string_view>
{
    auto pos = str.find_first_of(delim);
    if (pos == std::string::npos)
        return std::make_tuple(std::string_view(str.cbegin(), str.cend()), "");
    return std::make_tuple(
        std::string_view(str.cbegin(), str.cbegin() + pos),
        std::string_view(str.cbegin() + pos + 1, str.cend())
    );
}

static
auto nextDataBlock(Stream &stream) -> std::optional<Block>
{
    while (!stream.exhausted()) {
        const auto kind = stream.peekKind();
        if (kind == Token::Kind::COMMENT)
            stream.eatLine();
        else if (kind == Token::Kind::DATA_BLOCK) {
            const auto [ _unused, name ] = splitOnFirst(stream.eat().text, '_');
            return std::make_optional<Block>(std::string{name});
        } else
            stream.eat();
    }

    return std::nullopt;
}

static
auto tagToCategoryKeyword(const std::string &key, size_t lineNo)
{
    // CONFORMANCE: Check that there is only one dot
    auto [ category, keyword ] = splitOnFirst(key, '.');

    if (category.length() < 2)
        throw CifParseError{"Invalid category name token on line " + std::to_string(lineNo)};

    if (keyword.empty())
        return std::make_tuple(std::string{}, std::string{category.substr(1)});    // "Swap" keyword for category because we need to have anonymous categories to deal with non-mmCIF data

    return std::make_tuple(std::string{category.substr(1)}, std::string{keyword});
}

static
auto doMultiline(const std::string &text, Stream &stream)
{
    if (stream.exhausted())
        throw CifParseError{"Unexpected end of line in the middle of a multiline entry on line " + std::to_string(stream.lineCounter)};

    std::string multiline = text.substr(1);
    while (!stream.exhausted()) {
        const auto kind = stream.peekKind(true);

        if (kind == Token::Kind::MULTILINE) {
            stream.eat(true);
            return multiline;
        } else if (kind == Token::Kind::COMMENT)
            stream.eatLine();
        else {
            auto line = stream.eatLine();
            if (line.length() > 0) [[ likely ]] {
                auto end = line.back() == Stream::CARRIAGE_RET ? line.length() - 1 : line.length();
                multiline += line.substr(0, end);
            }
        }
    }

    throw CifParseError{"Unterminated multiline entry"};
}

static
auto doLoop(Block &block, Stream &stream)
{
    if (stream.exhausted())
        throw CifParseError{"Unexpected end of file on line " + std::to_string(stream.lineCounter)};

    const auto token = stream.eat();
    if (token.kind != Token::Kind::TAG)
        throw CifParseError{"Loop on line " + std::to_string(stream.lineCounter) + " does not define any columns"};

    auto [ loopCategory, keyword ] = tagToCategoryKeyword(std::string{token.text}, stream.lineCounter);
    std::vector<std::string> columns = { std::move(keyword) };

    std::vector<Values> columnData{};
    size_t columnIndex = 0;

    while (!stream.exhausted()) {
        const auto kind = stream.peekKind();

        if (kind == Token::Kind::COMMENT)
            stream.eatLine();
        else if (kind == Token::Kind::EMPTY)
            stream.eat();
        else if (kind == Token::Kind::VALUE) {
            if (columns.size() == 0)
                throw CifParseError{"Loop on line " + std::to_string(stream.lineCounter) + " does not define any columns"};
            else if (columnData.size() == 0)
                columnData.resize(columns.size());

            auto value = stream.eat().text;
            columnData[columnIndex].emplace_back(toValue(value));
            columnIndex = (columnIndex + 1) % columns.size();
        } else if (kind == Token::Kind::MULTILINE) {
            if (columns.size() == 0)
                throw CifParseError{"Loop on line " + std::to_string(stream.lineCounter) + " does not define any columns"};
            else if (columnData.size() == 0)
                columnData.resize(columns.size());

            auto value = doMultiline(std::string{stream.eat().text}, stream);
            columnData[columnIndex].push_back(std::move(value));
            columnIndex = (columnIndex + 1) % columns.size();
        } else if (kind == Token::Kind::TAG) {
            // If we have data that can make up a loop, assume that that loop ends here
            if (columnData.size() > 0 && columnIndex == 0)
                break;

            const auto token = stream.eat();
            auto [ category, keyword ] = tagToCategoryKeyword(std::string{token.text}, stream.lineCounter);
            if (loopCategory != category) {
                if (columnData.size() == 0)
                    throw CifParseError{"Mismatching categories " + category + " vs. " + loopCategory + " in loop on line " + std::to_string(stream.lineCounter)};
                else
                    throw CifParseError{"Malformed loop on line " + std::to_string(stream.lineCounter) + ", unterminated data"};
            }
            columns.push_back(std::move(keyword));
        } else {
            // If we have data that can make up a loop, assume that that loop ends here
            if (columnData.size() > 0 && columnIndex == 0)
                break;

            throw CifParseError{"Malformed loop on line " + std::to_string(stream.lineCounter) + ", unexpected token kind " + std::to_string(static_cast<int>(kind))};
        }
    }

    if (!(columnData.size() > 0 && columnIndex == 0))
        throw CifParseError{"File ended in the middle of a loop"};

    auto actualCategory = loopCategory.empty() ? block.nextAnonymousCategoryName() : loopCategory;
    for (size_t colIdx = 0; colIdx < columnData.size(); colIdx++)
        block.addMultiple(actualCategory, std::move(columns[colIdx]), std::move(columnData[colIdx]));
}

static
auto doTagValue(const std::string &text, Block &block, Stream &stream)
{
    if (stream.exhausted())
        throw CifParseError{"Unexpected end of file on line " + std::to_string(stream.lineCounter)};

    const auto [ category, keyword ] = tagToCategoryKeyword(text, stream.lineCounter);

    while (!stream.exhausted()) {
        const auto kind = stream.peekKind();

        if (kind == Token::Kind::VALUE) {
            block.add(category, keyword, toValue(stream.eat().text));
            return;
        } else if (kind == Token::Kind::COMMENT)
            stream.eatLine();
        else if (kind == Token::Kind::MULTILINE) {
            block.add(category, keyword, doMultiline(std::string{stream.eat().text}, stream));
            return;
        } else
            throw CifParseError{"Unexpected token kind " + std::to_string(std::underlying_type_t<Token::Kind>(kind)) + " on line " + std::to_string(stream.lineCounter)};
    }

    throw CifParseError{"Unexpected end of file on line " + std::to_string(stream.lineCounter)};
}

auto parse(const std::string_view &data) -> std::vector<Block>
{
    auto stream = Stream(data);

    std::vector<Block> blocks{};

    auto maybeBlock = nextDataBlock(stream);
    if (!maybeBlock.has_value())
        throw CifParseError{"Cif does not contain any data blocks"};
    auto currentBlock = maybeBlock.value();

    while (!stream.exhausted()) {
        const auto kind = stream.peekKind();

        if (kind == Token::Kind::DATA_BLOCK) {
            auto [ _unused, name ] = splitOnFirst(stream.eat().text, '_');
            blocks.push_back(std::move(currentBlock));
            currentBlock = Block{std::string{name}};
        } else if (kind == Token::Kind::TAG)
            doTagValue(std::string{stream.eat().text}, currentBlock, stream);
        else if (kind == Token::Kind::LOOP) {
            stream.eat();
            doLoop(currentBlock, stream);
        } else if (kind == Token::Kind::COMMENT)
            stream.eatLine();
        else if (kind == Token::Kind::MULTILINE)
            throw CifParseError{"Unexpected multiline entry marker on line " + std::to_string(stream.lineCounter)};
        else if (kind == Token::Kind::VALUE)
            throw CifParseError{"Unexpected value without name on line " + std::to_string(stream.lineCounter)};
        else if (kind == Token::Kind::EMPTY)
            stream.eat();
        else if (kind == Token::Kind::SAVE_BLOCK || kind == Token::Kind::STOP || kind == Token::Kind::GLOBAL_BLOCK) {
            stream.eat();

            blocks.push_back(std::move(currentBlock));
            auto maybeBlock = nextDataBlock(stream);
            if (!maybeBlock.has_value())
                return blocks;

            currentBlock = maybeBlock.value();
        } else
            throw CifParseError{"Unknown or unhandled token on line " + std::to_string(stream.lineCounter)};
    }

    blocks.push_back(std::move(currentBlock));

    return blocks;
}

} // namespace LLKAInternal::MiniCif
