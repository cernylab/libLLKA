/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _ASSIGNMENT_H
#define _ASSIGNMENT_H

#include <llka_cpp.h>

#include <filesystem>
#include <optional>
#include <string>

class AssignmentError {
public:
    LLKA_RetCode tRet;
    std::string extraMessage;
};

class Step {
public:
    enum ClassificationState {
        NotYetClassified,
        Classified,
        ClassificationFailed
    };

    Step(const bool multipleModels, const std::string &id, const LLKA::Structure &structure) noexcept;
    Step(const bool multipleModels, const std::string &id, const LLKA::Structure &structure, const LLKA::ClassifiedStep &classification) noexcept;

    std::string name;
    ClassificationState state;
    LLKA::Structure structure;
    LLKA::ClassifiedStep classification;

    static auto DNATCOName(const bool multipleModels, const std::string &id, const LLKA::Structure &structure) -> std::string;
};

class LoadedStructure {
public:
    LoadedStructure() = default;
    LoadedStructure(const bool hasMultipleModels, std::string id, std::vector<Step> steps) noexcept;

    auto isValid() { return !id.empty(); }

    bool hasMultipleModels;
    std::string id;
    std::vector<Step> steps;
};

auto loadStructure(const std::filesystem::path &path) -> LoadedStructure;

#endif // _ASSIGNMENT_H
