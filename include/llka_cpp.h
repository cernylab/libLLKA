// vim: set sw=4 ts=4 sts=4 expandtab :

/*!
 * @file llka_cpp.h
 * This file provides C++ wrappers around the C API of libLLKA.
 * These bindings are provided primarily to allow for easy binding of JavaScript to libLLKA.
 * Some of the wrappers provided here come with a minor overhead compared to their C counterparts.
 * The C API is the preferred way how to call libLLKA
 */

#ifdef __cplusplus

#ifndef _LLKA_CPP_H
#define _LLKA_CPP_H

#include "llka_main.h"
#include "llka_connectivity_similarity.h"
#include "llka_classification.h"
#include "llka_minicif.h"
#include "llka_ntc.h"
#include "llka_structure.h"

#include <cassert>
#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
#include <filesystem>
#endif // LLKA_FILESYSTEM_ACCESS_DISABLED
#include <map>
#include <string>
#include <vector>

#ifdef LLKA_PLATFORM_EMSCRIPTEN
    #include <emscripten/bind.h>

    #define _EMX_GET(type, var) auto _emsGet_##var() const -> const type &;
    #define _EMX_SET(type, var) auto _emsSet_##var(const type &var) -> void;
    #define _EMX_GET_SET(type, var) _EMX_GET(type, var) _EMX_SET(type, var)
#else
    #define _EMX_GET(type, var)
    #define _EMX_SET(type, var)
    #define _EMX_GET_SET(type, var)
#endif // LLKA_PLATFORM_EMSCRIPTEN

namespace LLKA {

// Forward declarations (but I still love you, my C++17 sweetie...)
class Atom;
using Structure = std::vector<Atom>;
using Structures = std::vector<Structure>;
using StructureView = std::vector<const Atom *>;
class CifData;

template <typename T, bool>
struct ResultSuccessReturnType;

template <typename S>
struct ResultSuccessReturnType<S, true> {
    using RT = S &;
};

template <typename S>
struct ResultSuccessReturnType<S, false> {
    using RT = S &&;
};

// C++-only generic return type that we return from all function that return a value on success but may fail with an error code
template <typename S, typename F>
class Result {
public:
    template <typename ST = S>
    Result(ST success) :
        m_u{std::move(success)},
        m_isSuccess{true},
        m_movedAway{false}
    {
    }

    Result(F failure, bool) : // The dummy parameter is needed for Embind to be able to differentiate between S and F c-tors
        m_u{std::move(failure)},
        m_isSuccess{false},
        m_movedAway{false}
    {
    }

    Result(const Result &other) :
        m_u{S{}}  // Emscripten needs us to default-init this to something.
    {
        m_u.success.~S(); // Needed because m_u must be default-inited to something

        if (other.m_isSuccess)
            new (&this->m_u.success) S(other.m_u.success);
        else
            new (&this->m_u.failure) F(other.m_u.failure);
        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;
    }

    Result(Result &&other) noexcept :
        m_u{S{}}  // Emscripten needs us to default-init this to something.

    {
        m_u.success.~S(); // Needed because m_u must be default-inited to something

        if (other.m_isSuccess)
            new (&this->m_u.success) S(std::move(other.m_u.success));
        else
            new (&this->m_u.failure) F(std::move(other.m_u.failure));
        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;

        other.m_movedAway = true;
    }

    ~Result()
    {
        if (m_movedAway)
            return;
        if (m_isSuccess)
            m_u.success.~S();
        else
            m_u.failure.~F();
    }

    Result & operator=(const Result &other)
    {
        if (other.m_isSuccess)
            this->m_u.success = other.m_u.success;
        else
            this->m_u.failure = other.m_u.failure;
        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;

        return *this;
    }

    Result & operator=(Result &&other) noexcept
    {
        if (other.m_isSuccess)
            this->m_u.success = std::move(other.m_u.success);
        else
            this->m_u.failure = std::move(other.m_u.failure);
        this->m_isSuccess = other.m_isSuccess;
        other.m_movedAway = true;
        this->m_movedAway = false;

        return *this;
    }

    const F & failure() const
    {
        assert(!m_movedAway);

        if (m_isSuccess)
            throw std::runtime_error{"Cannot get failed value for a succesful result"};
        return m_u.failure;
    }

    typename ResultSuccessReturnType<S, std::is_copy_constructible_v<S>>::RT success()
    {
        assert(!m_movedAway);

        if (!m_isSuccess)
            throw std::runtime_error{"Cannot get success value for a failed result"};

        if constexpr (std::is_copy_constructible_v<S>) {
            return m_u.success;
        } else {
            return std::move(m_u.success);
        }
    }

#ifdef LLKA_PLATFORM_EMSCRIPTEN
    // const-overloaded methods make Emscripten cry for momma.
    // We may need to provide some kind of a macro to make the code build seamlessly
    // in C++ and JS...
    const S & const_success() const
#else
    template <typename ST = S>
    const ST & success() const
#endif // LLKA_PLATFORM_EMSCRIPTEN
    {
        assert(!m_movedAway);

        if (!m_isSuccess)
            throw std::runtime_error{"Cannot get success value for a failed result"};
        return m_u.success;
    }

    bool isSuccess() const
    {
        assert(!m_movedAway);

        return m_isSuccess;
    }

    template <typename ...Args>
    static Result fail(Args &&...args)
    {
        return Result{F{std::forward<Args>(args)...}, false};
    }

    template <typename ...Args>
    static Result succeed(Args &&...args)
    {
        return Result{S{std::forward<Args>(args)...}};
    }

private:
    union U {
        F failure;
        S success;
        U(S s) noexcept : success{std::move(s)} {}
        template <typename FT = std::enable_if_t<!std::is_same_v<S, F>, F>>
        U(FT f) noexcept : failure{std::move(f)} {}
        ~U() {}
    } m_u;

    bool m_isSuccess;
    bool m_movedAway;
};
template <typename S> using RCResult = Result<S, LLKA_RetCode>;

template <typename F>
class Result<void, F> {
public:
    Result() noexcept :
        m_isSuccess{true},
        m_movedAway{false}
    {
    }

    Result(F f, bool) noexcept :
        m_failure{std::move(f)},
        m_isSuccess{false},
        m_movedAway{false}
    {
    }

    Result(const Result &other)
    {
        if (!other.m_isSuccess)
            this->m_failure = other.m_failure;

        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;
    }

    Result(Result &&other) noexcept
    {
        if (!other.m_isSuccess)
            this->m_failure = std::move(other.m_failure);

        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;

        other.m_movedAway = true;
    }

    Result & operator=(const Result &other)
    {
        if (!other.m_isSuccess)
            this->m_failure = other.m_failure;

        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;

        return *this;
    }

    Result & operator=(Result &&other) noexcept
    {
        if (!other.m_isSuccess)
            this->m_failure = std::move(m_failure);

        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;

        other.m_movedAway = true;

        return *this;
    }

    const F & failure() const
    {
        assert(!m_movedAway);

        if (m_isSuccess)
            throw std::runtime_error{"Cannot get failed value for a succesful result"};
        return failure;
    }

    bool isSuccess() const
    {
        assert(!m_movedAway);

        return m_isSuccess;
    }

    template <typename ...Args>
    static Result fail(Args &&...args)
    {
        return Result{F{std::forward<Args>(args)...}, false};
    }

    static Result succeed()
    {
        return Result{};
    }

private:
    F m_failure;

    bool m_isSuccess;
    bool m_movedAway;
};

template <typename S>
class Result<S, void> {
public:
    Result(S s) noexcept :
        m_success{std::move(s)},
        m_isSuccess{true},
        m_movedAway{false}
    {
    }

    Result() noexcept :
        m_isSuccess{false},
        m_movedAway{false}
    {
    }

    Result(const Result &other)
    {
        if (other.m_isSuccess)
            this->m_success = other.m_success;

        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;
    }

    Result(Result &&other) noexcept
    {
        if (other.m_isSuccess)
            this->m_success = std::move(other.m_success);

        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;

        other.m_movedAway = true;
    }

    Result & operator=(const Result &other)
    {
        if (other.m_isSuccess)
            this->m_success = other.m_success;

        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;

        return *this;
    }

    Result & operator=(Result &&other) noexcept
    {
        if (other.m_isSuccess)
            this->m_success = std::move(m_success);

        this->m_isSuccess = other.m_isSuccess;
        this->m_movedAway = false;

        other.m_movedAway = true;

        return *this;
    }

    const S & success() const
    {
        assert(!m_movedAway);

        if (!m_isSuccess)
            throw std::runtime_error{"Cannot get success value for a failed result"};
        return m_success;
    }

    bool isSuccess() const
    {
        assert(!m_movedAway);

        return m_isSuccess;
    }

    static Result fail()
    {
        return Result{};
    }

    template <typename ...Args>
    static Result succeed(Args &&...args)
    {
        return Result{S{std::forward<Args>(args)...}};
    }

private:
    S m_success;

    bool m_isSuccess;
    bool m_movedAway;
};

//
// Main
//

using Points = std::vector<LLKA_Point>; // TODO: Revise this when we figure out how to do memory alignment of point arrays

LLKA_CPP_API
auto operator<<(std::ostream &os, const LLKA_Point &pt) -> std::ostream &;

LLKA_CPP_API
auto errorToString(LLKA_RetCode tRet) -> std::string;

//
// Structure
//

// We cannot take advantage of any "modern" features like constructors or class methods
// because we would have to export this object as a class instead of value_object (== POD)
// in Emscripten builds. Emscripten seems to have major memory management issues when
// complex objects are exported as classes and these issues disappear when we replace them
// with dumb value_objects. #software_is_in_decline
class LLKA_CPP_API Atom {
public:
    std::string type_symbol;
    std::string label_atom_id;
    std::string label_entity_id;
    std::string label_comp_id;
    std::string label_asym_id;
    std::string auth_atom_id;
    std::string auth_comp_id;
    std::string auth_asym_id;
    std::string pdbx_PDB_ins_code;
    LLKA_Point coords;
    uint32_t id;
    int32_t label_seq_id;
    int32_t auth_seq_id;
    int32_t pdbx_PDB_model_num;
    char label_alt_id;
};

LLKA_CPP_API
auto makeEmptyAtom() -> Atom;

LLKA_CPP_API
auto makeAtom(
    std::string type_symbol,
    std::string label_atom_id,
    std::string label_entity_id,
    std::string label_comp_id,
    std::string label_asym_id,
    std::string auth_atom_id,
    std::string auth_comp_id,
    std::string auth_asym_id,
    LLKA_Point coords,
    uint32_t id,
    int32_t label_seq_id,
    int32_t auth_seq_id,
    int32_t pdbx_PDB_model_num,
    std::string pdbx_PDB_ins_code,
    char label_alt_id
) -> Atom;

LLKA_CPP_API
auto atomsEqual(const Atom &a, const Atom &b) noexcept -> bool;

class LLKA_CPP_API AltIdSplit {
public:
    AltIdSplit() = default;

    Structure structure;
    char altId;
};

LLKA_CPP_API
auto operator<<(std::ostream &os, const Atom &atom) -> std::ostream &;

LLKA_CPP_API
auto operator<<(std::ostream &os, const Structure &stru) -> std::ostream &;

inline
auto atomMatches(
    const Atom &atom,
    const std::string &label_atom_id,
    const std::string &label_comp_id,
    const std::string &label_asym_id,
    int32_t label_seq_id,
    char label_alt_id = LLKA_NO_ALTID,
    const std::string &pdbx_PDB_ins_code = LLKA_NO_INSCODE,
    int32_t pdbx_PDB_model_num = 1
) {
    // This function must be kept in sync with atomMatchesCriteria() from structure.cpp
    // Conversion to C API objects and a function call would be too expensive to just
    // call that function.
    if (!label_atom_id.empty() && label_atom_id != atom.label_atom_id)
        return false;

    if (!label_comp_id.empty() && label_comp_id != atom.label_comp_id)
        return false;

    if (!label_asym_id.empty() && label_asym_id != atom.label_asym_id)
        return false;

    if (label_seq_id >= 0 && atom.label_seq_id != label_seq_id)
        return false;

    if (label_alt_id != LLKA_NO_ALTID && atom.label_alt_id != label_alt_id)
        return false;

    if (pdbx_PDB_ins_code != LLKA_NO_INSCODE && pdbx_PDB_ins_code != atom.pdbx_PDB_ins_code)
        return false;

    return atom.pdbx_PDB_model_num == pdbx_PDB_model_num;
}

LLKA_CPP_API
auto compareAtoms(const Atom &a, const Atom &b, bool ignoreId) noexcept -> bool;

LLKA_CPP_API
auto makeStructure(const LLKA_Atom *atoms, size_t nAtoms) -> Structure;

LLKA_CPP_API
auto splitByAltIds(const Structure &stru) -> std::vector<AltIdSplit>;

LLKA_CPP_API
auto splitStructureToDinucleotideSteps(const Structure &stru) -> RCResult<Structures>;

//
// Measurements
//

template <typename T> requires std::is_floating_point_v<T>
LLKA_CPP_API
auto measureAngle(const Atom &a, const Atom &b, const Atom &c) -> T;

template <typename T> requires std::is_floating_point_v<T>
LLKA_CPP_API
auto measureDihedral(const Atom &a, const Atom &b, const Atom &c, const Atom &d) -> T;

template <typename T> requires std::is_floating_point_v<T>
LLKA_CPP_API
auto measureDistance(const Atom &a, const Atom &b) -> T;

//
// Superposition
//

LLKA_CPP_API
auto centroid(const Points &points) noexcept -> LLKA_Point;

LLKA_CPP_API
auto centroid(const Structure &stru) noexcept -> LLKA_Point;

LLKA_CPP_API
auto rmsd(const Points &a, const Points &b) noexcept -> RCResult<double>;

LLKA_CPP_API
auto rmsd(const Structure &a, const Structure &b) noexcept -> RCResult<double>;

LLKA_CPP_API
auto superpose(Points &what, const Points &onto) noexcept -> RCResult<double>;

LLKA_CPP_API
auto superpose(Structure &what, const Structure &onto) noexcept -> RCResult<double>;

//
// Segmentation
//

class LLKA_CPP_API StructureSegments {
/*
 * To counter memory management issues that cause memory leaks with Emscripten,
 * we need to export LLKA::Atom class to JavaScript as value_object instead of a class.
 * This, apparently, has the unfortunate consequence of not being able to expose
 * a _pointer_ to LLKA::Atom in the JavaScript bindings. This breaks the originally
 * intended layout of Segmentation mapping.
 * To make Segmentation available in the JavaScript bindings, we need to provide
 * a simplified version that stores copies of atoms instead of pointers.
 * It is dumb but at the moment there does not appear to be a better solution.
 * This is why you will have to deal with the barrage of ifdefs.
 */
public:
#ifdef LLKA_PLATFORM_EMSCRIPTEN
    using Atoms = std::vector<Atom>;
#else
    using Atoms = std::vector<Atom *>;
#endif // LLKA_PLATFORM_EMSCRIPTEN

    struct Residue {
        Atoms atoms;
    };
    using Residues = std::map<int32_t, Residue>;

    struct Chain {
        Residues residues;
#ifndef LLKA_PLATFORM_EMSCRIPTEN
        Atoms atoms;
#endif // LLKA_PLATFORM_EMSCRIPTEN
    };
    using Chains = std::map<std::string, Chain>;

    struct Model {
        Chains chains;
#ifndef LLKA_PLATFORM_EMSCRIPTEN
        Atoms atoms;
#endif // LLKA_PLATFORM_EMSCRIPTEN
    };
    using Models = std::map<int32_t, Model>;

    // Default c-tor makes absolutely no sense but Emscripten requires it
#ifdef LLKA_PLATFORM_EMSCRIPTEN
    StructureSegments() noexcept;
#endif // LLKA_PLATFORM_EMSCRIPTEN
    StructureSegments(Structure &structure) noexcept;

    Models models;
    // We cannot use reference because Emscripten doesn't take it
    Structure *structure;

    _EMX_GET_SET(Models, models)

#ifdef LLKA_PLATFORM_EMSCRIPTEN
    auto _emsGet_structure() const -> const Structure &;
#endif // LLKA_PLATFORM_EMSCRIPTEN
};

//
// NtC
//

class LLKA_CPP_API AtomNameQuad {
public:
    AtomNameQuad()
    {
        a.reserve(sizeof(LLKA_AtomNameQuad::a));
        b.reserve(sizeof(LLKA_AtomNameQuad::b));
        c.reserve(sizeof(LLKA_AtomNameQuad::c));
        d.reserve(sizeof(LLKA_AtomNameQuad::d));
    }

    AtomNameQuad(const AtomNameQuad &) = default;

    AtomNameQuad(AtomNameQuad &&other) noexcept :
        a{std::move(other.a)},
        b{std::move(other.b)},
        c{std::move(other.c)},
        d{std::move(other.d)}
    {}

    AtomNameQuad(std::string _a, std::string _b, std::string _c, std::string _d) noexcept :
        a{std::move(_a)}, b{std::move(_b)}, c{std::move(_c)}, d{std::move(_d)}
    {}

    AtomNameQuad & operator=(const AtomNameQuad &) = default;
    AtomNameQuad & operator=(AtomNameQuad &&other) noexcept
    {
        this->a = std::move(other.a);
        this->b = std::move(other.b);
        this->c = std::move(other.c);
        this->d = std::move(other.d);

        return *this;
    }

    std::string a;
    std::string b;
    std::string c;
    std::string d;

    _EMX_GET_SET(std::string, a)
    _EMX_GET_SET(std::string, b)
    _EMX_GET_SET(std::string, c)
    _EMX_GET_SET(std::string, d)
};

class LLKA_CPP_API BackboneAtom {
public:
    BackboneAtom() = default;

    LLKA_ResidueInStep residue;
    std::string name;
};

LLKA_CPP_API
auto backboneAtomIndex(const BackboneAtom &bkbnAtom, const Structure &backbone) noexcept -> size_t;

LLKA_CPP_API
auto calculateStepMetrics(const Structure &stru) noexcept -> RCResult<LLKA_StepMetrics>;

LLKA_CPP_API
auto calculateStepMetricsDifferenceAgainstReference(const Structure &stru, LLKA_NtC ntc) noexcept -> RCResult<LLKA_StepMetrics>;

LLKA_CPP_API
auto crossResidueMetric(LLKA_CrossResidueMetric metric, const Structure &stru) noexcept -> RCResult<Structure>;

LLKA_CPP_API
auto crossResidueMetricAtoms(const std::string &firstBase, const std::string &secondBase, LLKA_CrossResidueMetric metric) -> RCResult<AtomNameQuad>;

LLKA_CPP_API
auto crossResidueMetricAtoms(const Structure &stru, LLKA_CrossResidueMetric metric) -> RCResult<AtomNameQuad>;

LLKA_CPP_API
auto crossResidueMetricName(LLKA_CrossResidueMetric metric, bool greek) noexcept -> std::string;

LLKA_CPP_API
auto dinucleotideTorsion(LLKA_DinucleotideTorsion torsion, const Structure &stru) -> RCResult<Structure>;

LLKA_CPP_API
auto dinucleotideTorsionAtoms(const std::string &firstBase, const std::string &secondBase, LLKA_DinucleotideTorsion torsion) -> RCResult<AtomNameQuad>;

LLKA_CPP_API
auto dinucleotideTorsionAtoms(const Structure &stru, LLKA_DinucleotideTorsion torsion) -> RCResult<AtomNameQuad>;

LLKA_CPP_API
auto dinucleotideTorsionName(LLKA_DinucleotideTorsion torsion, bool greek) noexcept -> std::string;

LLKA_CPP_API
auto extractBackbone(const Structure &stru) -> RCResult<Structure>;

LLKA_CPP_API
auto extractExtendedBackbone(const Structure &stru) noexcept -> RCResult<Structure>;

LLKA_CPP_API
auto extractMetricsStructure(const Structure &stru) noexcept -> RCResult<Structure>;

LLKA_CPP_API
auto nameToCANA(const std::string &name) noexcept -> LLKA_CANA;

LLKA_CPP_API
auto CANAToName(LLKA_CANA cana) noexcept -> std::string;

LLKA_CPP_API
auto nameToNtC(const std::string &name) noexcept -> LLKA_NtC;

LLKA_CPP_API
auto NtCToName(LLKA_NtC ntc) noexcept -> std::string;

using StepInfo = LLKA_StepInfo;
LLKA_CPP_API
auto structureIsStep(const Structure &stru) noexcept -> RCResult<LLKA_StepInfo>;

using StepMetrics = LLKA_StepMetrics;
auto operator<<(std::ostream &os, const StepMetrics &metrics) -> std::ostream &;

//
// Connectivity, similarity
//

using Connectivity = LLKA_Connectivity;
using Connectivities = std::vector<Connectivity>;
using Similarity = LLKA_Similarity;
using Similarities = std::vector<LLKA_Similarity>;

LLKA_CPP_API
auto measureStepConnectivity(const Structure &positionFirst, LLKA_NtC ntcFirst, const Structure &positionSecond, LLKA_NtC ntcSecond) noexcept -> RCResult<Connectivity>;

LLKA_CPP_API
auto measureStepConnectivity(const Structure &positionFirst, std::vector<LLKA_NtC> ntcsFirst, const Structure &positionSecond, LLKA_NtC ntcSecond) noexcept -> RCResult<Connectivities>;

LLKA_CPP_API
auto measureStepConnectivity(const Structure &positionFirst, LLKA_NtC ntcFirst, const Structure &positionSecond, std::vector<LLKA_NtC> ntcsSecond) noexcept -> RCResult<Connectivities>;

LLKA_CPP_API
auto measureStepConnectivity(const Structure &positionFirst, const Structure &dinuFirst, const Structure &positionSecond, const Structure &dinuSecond) noexcept -> RCResult<Connectivity>;

LLKA_CPP_API
auto measureStepConnectivity(const Structure &positionFirst, const Structure &dinuFirst, const Structure &positionSecond, const Structures &dinusSecond) noexcept -> RCResult<Connectivities>;

LLKA_CPP_API
auto measureStepSimilarity(const Structure &stepStru, LLKA_NtC ntc) noexcept -> RCResult<Similarity>;

LLKA_CPP_API
auto measureStepSimilarity(const Structure &stepStru, std::vector<LLKA_NtC> ntcs) noexcept -> RCResult<Similarities>;

LLKA_CPP_API
auto measureStepSimilarity(const Structure &stepStru, const Structure &refStru) noexcept -> RCResult<Similarity>;

//
// MiniCif
//

class CifData {
public:
    class Value {
    public:
        std::string text;
        LLKA_CifDataValueState state;
    };
    using Values = std::vector<Value>;

    class Item {
    public:
        std::string keyword;
        Values values;
    };
    using Items = std::vector<Item>;

    class Category {
    public:
        Category() = default;
        Category(std::string name, Items items);

        auto isLoop() -> bool;

        std::string name;
        Items items;

        _EMX_GET_SET(std::string, name)
        _EMX_GET_SET(Items, items)

    };
    using Categories = std::vector<Category>;

    class Block {
    public:
        std::string name;
        Categories categories;
    };
    using Blocks = std::vector<Block>;

    Blocks blocks;
};

class LLKA_CPP_API ImportedStructure {
public:
    ImportedStructure() = default;

    std::string id;
    Structure structure;
    CifData cifData;
};

class LLKA_CPP_API CifError {
public:
    CifError() = default;

    LLKA_RetCode tRet;
    std::string error;
};
using CifResult = Result<ImportedStructure, CifError>;

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED

LLKA_CPP_API
auto cifToStructure(const std::filesystem::path &path, int32_t options = 0) -> CifResult;

#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

LLKA_CPP_API
auto cifToStructure(const std::string &text, int32_t options = 0) -> CifResult;

LLKA_CPP_API
auto cifDataToString(const CifData &cifData, bool pretty) -> RCResult<std::string>;

//
// Nucleotide
//

LLKA_CPP_API
auto extractNucleotide(const Structure &stru, int32_t pdbx_PDB_model_num, const std::string &label_asym_id, int32_t label_seq_id) -> Structure;

LLKA_CPP_API
auto extractNucleotideView(const Structure &stru, int32_t pdbx_PDB_model_num, const std::string &label_asym_id, int32_t label_seq_id) -> StructureView;

LLKA_CPP_API
auto extractRibose(const Structure &stru) -> RCResult<Structure>;

LLKA_CPP_API
auto extractRiboseView(const Structure &stru) -> RCResult<StructureView>;

LLKA_CPP_API
auto isNucleotideCompound(const std::string &compId) -> bool;

LLKA_CPP_API
auto nameToSugarPucker(const std::string &name) -> LLKA_SugarPucker;

LLKA_CPP_API
auto riboseMetrics(const Structure &stru) -> RCResult<LLKA_RiboseMetrics>;

LLKA_CPP_API
auto riboseMetrics(const StructureView &view) -> RCResult<LLKA_RiboseMetrics>;

LLKA_CPP_API
auto sugarPucker(const Structure &stru) -> RCResult<LLKA_SugarPucker>;

LLKA_CPP_API
auto sugarPucker(const StructureView &stru) -> RCResult<LLKA_SugarPucker>;

LLKA_CPP_API
auto sugarPuckerToName(LLKA_SugarPucker pucker, LLKA_SugarPuckerNameBrevity brevity) -> std::string;

// Classification - declaration of value objects that cannot be mapped
// directly to the C structs

class LLKA_CPP_API GoldenStep {
public:
    GoldenStep() = default;
    GoldenStep(const LLKA_GoldenStep &gs);

    LLKA_SugarPucker pucker_1;
    LLKA_SugarPucker pucker_2;
    LLKA_NuAngles nuAngles_1;
    LLKA_NuAngles nuAngles_2;
    LLKA_StepMetrics metrics;
    std::string name;
    size_t clusterIdx;
    int32_t clusterNumber;
};

//
// Resource loaders
//

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED

LLKA_CPP_API
auto loadClusterNuAngles(const std::filesystem::path &path) -> RCResult<std::vector<LLKA_ClusterNuAngles>>;

LLKA_CPP_API
auto loadClusters(const std::filesystem::path &path) -> RCResult<std::vector<LLKA_ClassificationCluster>>;

LLKA_CPP_API
auto loadConfals(const std::filesystem::path &path) -> RCResult<std::vector<LLKA_Confal>>;

LLKA_CPP_API
auto loadConfalPercentiles(const std::filesystem::path &path) -> RCResult<std::vector<LLKA_ConfalPercentile>>;

LLKA_CPP_API
auto loadGoldenSteps(const std::filesystem::path &path) -> RCResult<std::vector<GoldenStep>>;

#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

LLKA_CPP_API
auto loadClusterNuAngles(const std::string &text) -> RCResult<std::vector<LLKA_ClusterNuAngles>>;

LLKA_CPP_API
auto loadClusters(const std::string &text) -> RCResult<std::vector<LLKA_ClassificationCluster>>;

LLKA_CPP_API
auto loadConfals(const std::string &text) -> RCResult<std::vector<LLKA_Confal>>;

LLKA_CPP_API
auto loadConfalPercentiles(const std::string &path) -> RCResult<std::vector<LLKA_ConfalPercentile>>;

LLKA_CPP_API
auto loadGoldenSteps(const std::string &text) -> RCResult<std::vector<GoldenStep>>;

//
// Tracing
//

class TracepointInfo {
public:
    int32_t TPID;
    std::string description;
};

LLKA_CPP_API
auto toggleAllTracepoints(bool enable) -> void;

LLKA_CPP_API
auto toggleTracepoint(int32_t TPID, bool enable) -> void;

LLKA_CPP_API
auto trace(bool dontClear = false) -> std::string;

LLKA_CPP_API
auto tracepointInfo() -> std::vector<TracepointInfo>;

LLKA_CPP_API
auto tracepointState(int32_t TPID) -> bool;

//
// Classification
//

class LLKA_CPP_API ClassificationContext {
public:
    ClassificationContext();
    ClassificationContext(LLKA_ClassificationContext *ctx);

    // We need to cater to Emscripten requiring copy c-tor
#ifdef LLKA_PLATFORM_EMSCRIPTEN
    ClassificationContext(const ClassificationContext &other);
#else
    ClassificationContext(const ClassificationContext &other) = delete;
#endif // LLKA_PLATFORM_EMSCRIPTEN

    ClassificationContext(ClassificationContext &&other) noexcept;

    ~ClassificationContext();

    ClassificationContext & operator=(const ClassificationContext &) = delete;
    ClassificationContext & operator=(ClassificationContext &&other) noexcept;

    const LLKA_ClassificationContext * get() const;

    auto isValid() const -> bool;

private:
    LLKA_ClassificationContext *m_ctx;

#ifdef LLKA_PLATFORM_EMSCRIPTEN
    mutable bool m_movedAway;
#else
    bool m_movedAway;
#endif// LLKA_PLATFORM_EMSCRIPTEN
};

class LLKA_CPP_API ClassifiedStep {
public:
    LLKA_NtC assignedNtC;
    LLKA_CANA assignedCANA;
    LLKA_NtC closestNtC;
    LLKA_CANA closestCANA;
    LLKA_ConfalScore confalScore;
    double euclideanDistanceNtCIdeal;
    LLKA_StepMetrics metrics;
    LLKA_StepMetrics differencesFromNtCAverages;
    LLKA_NuAngles nuAngles_1;
    LLKA_NuAngles nuAngles_2;
    double ribosePseudorotation_1;
    double ribosePseudorotation_2;
    double tau_1;
    double tau_2;
    LLKA_SugarPucker sugarPucker_1;
    LLKA_SugarPucker sugarPucker_2;
    LLKA_NuAngles nuAngleDifferences_1;
    LLKA_NuAngles nuAngleDifferences_2;
    double rmsdToClosestNtC;

    std::string closestGoldenStep;

    int32_t violations;

    int16_t violatingTorsionsAverage;
    int16_t violatingTorsionsNearest;

    ClassifiedStep();
    ClassifiedStep(const LLKA_ClassifiedStep &cStep);

    auto hasViolations() const -> bool;
    auto namedViolations() const -> std::vector<std::string>;

    _EMX_GET_SET(LLKA_NtC, assignedNtC)
    _EMX_GET_SET(LLKA_CANA, assignedCANA)
    _EMX_GET_SET(LLKA_NtC, closestNtC)
    _EMX_GET_SET(LLKA_CANA, closestCANA)
    _EMX_GET_SET(LLKA_ConfalScore, confalScore)
    _EMX_GET_SET(double, euclideanDistanceNtCIdeal)
    _EMX_GET_SET(LLKA_StepMetrics, metrics)
    _EMX_GET_SET(LLKA_StepMetrics, differencesFromNtCAverages)
    _EMX_GET_SET(LLKA_NuAngles, nuAngles_1)
    _EMX_GET_SET(LLKA_NuAngles, nuAngles_2)
    _EMX_GET_SET(double, ribosePseudorotation_1)
    _EMX_GET_SET(double, ribosePseudorotation_2)
    _EMX_GET_SET(double, tau_1)
    _EMX_GET_SET(double, tau_2)
    _EMX_GET_SET(LLKA_SugarPucker, sugarPucker_1)
    _EMX_GET_SET(LLKA_SugarPucker, sugarPucker_2)
    _EMX_GET_SET(LLKA_NuAngles, nuAngleDifferences_1)
    _EMX_GET_SET(LLKA_NuAngles, nuAngleDifferences_2)
    _EMX_GET_SET(double, rmsdToClosestNtC)
    _EMX_GET_SET(std::string, closestGoldenStep)
    _EMX_GET_SET(int32_t, violations)
    _EMX_GET_SET(int16_t, violatingTorsionsAverage)
    _EMX_GET_SET(int16_t, violatingTorsionsNearest)
};

class LLKA_CPP_API AttemptedClassifiedStep {
public:
    AttemptedClassifiedStep() = default;

    LLKA_RetCode status;
    ClassifiedStep step;
};
using AttemptedClassifiedSteps = std::vector<AttemptedClassifiedStep>;

LLKA_CPP_API
auto averageConfal(const std::vector<ClassifiedStep> &steps, const ClassificationContext &ctx) noexcept -> LLKA_AverageConfal;

LLKA_CPP_API
auto averageConfal(const std::vector<AttemptedClassifiedStep> &steps, const ClassificationContext &ctx) noexcept -> LLKA_AverageConfal;

LLKA_CPP_API
auto classificationClusterForNtC(LLKA_NtC ntc, const ClassificationContext &ctx) noexcept -> RCResult<LLKA_ClassificationCluster>;

LLKA_CPP_API
auto confalForNtC(LLKA_NtC ntc, const ClassificationContext &ctx) noexcept -> RCResult<LLKA_Confal>;

LLKA_CPP_API
auto confalPercentile(double confalScore, const ClassificationContext &ctx) noexcept -> double;

LLKA_CPP_API
auto classifyStep(const Structure &stru, const ClassificationContext &ctx) -> RCResult<ClassifiedStep>;

LLKA_CPP_API
auto classifySteps(const Structures &strus, const ClassificationContext &ctx) -> RCResult<AttemptedClassifiedSteps>;

LLKA_CPP_API
auto initializeClassificationContext(
    const std::vector<LLKA_ClassificationCluster> &clusters,
    const std::vector<GoldenStep> &goldenSteps,
    const std::vector<LLKA_Confal> &confals,
    const std::vector<LLKA_ClusterNuAngles> &clusterNuAngles,
    const std::vector<LLKA_ConfalPercentile> &confalPercentiles,
    const LLKA_ClassificationLimits &limits,
    double maxCloseEnoughRmsd
) -> RCResult<ClassificationContext>;

} // namespace LLKA


#ifdef LLKA_PLATFORM_EMSCRIPTEN

namespace LLKA {

template <typename T>
LLKA_CPP_API
auto makeStdVector()
{
    return std::vector<T>{};
}

template LLKA_CPP_API auto makeStdVector<double>();
template LLKA_CPP_API auto makeStdVector<Atom>();
template LLKA_CPP_API auto makeStdVector<Atom *>();
template LLKA_CPP_API auto makeStdVector<Structure>();
template LLKA_CPP_API auto makeStdVector<LLKA_NtC>();
template LLKA_CPP_API auto makeStdVector<CifData::Value>();
template LLKA_CPP_API auto makeStdVector<CifData::Item>();
template LLKA_CPP_API auto makeStdVector<CifData::Category>();
template LLKA_CPP_API auto makeStdVector<std::string>();
template LLKA_CPP_API auto makeStdVector<int32_t>();

} // LLKA

#endif // LLKA_PLATFORM_EMSCRIPTEN

#ifdef LLKA_GENERATE_EMSCRIPTEN_BINDINGS

#define _EMX_CLS_MTH(type, cls, mth) .function(#mth, &cls::mth)
#define _EMX_CLS_PROP(prop, cls) .property(#prop, &cls::_emsGet_##prop, &cls::_emsSet_##prop)
#define _EMX_CLS_PROP_READONLY(prop, cls) .property(#prop, &cls::_emsGet_##prop)
#define _EMX_ENUM_VAL(ev) .value(#ev, ev)
#define _EMX_VOF(f, vo) .field(#f, &vo::f)
#define _EMX_MK_RCRESULT(S) \
    emscripten::class_<LLKA::RCResult<S>>("RCResult_"#S) \
        .constructor<S>() \
        .constructor<LLKA_RetCode, bool>() \
        .function("failure", &LLKA::RCResult<S>::failure) \
        .function("success", std::function<typename LLKA::ResultSuccessReturnType<S, std::is_copy_constructible_v<S>>::RT (LLKA::RCResult<S>&)>(&LLKA::RCResult<S>::success)) \
        .function("const_success", std::function<const S &(const LLKA::RCResult<S>&)>(&LLKA::RCResult<S>::const_success)) \
        .function("isSuccess", &LLKA::RCResult<S>::isSuccess)

EMSCRIPTEN_BINDINGS(LLKA)
{
    //
    // Main
    //

    emscripten::register_vector<std::string>("StringVector");
    emscripten::register_vector<int32_t>("Int32Vector");

    emscripten::enum_<LLKA_RetCode>("RetCode")
        _EMX_ENUM_VAL(LLKA_OK)
        _EMX_ENUM_VAL(LLKA_E_INVALID_ARGUMENT)
        _EMX_ENUM_VAL(LLKA_E_MISMATCHING_SIZES)
        _EMX_ENUM_VAL(LLKA_E_NOT_IMPLEMENTED)
        _EMX_ENUM_VAL(LLKA_E_MULTIPLE_ALT_IDS)
        _EMX_ENUM_VAL(LLKA_E_MISSING_ATOMS)
        _EMX_ENUM_VAL(LLKA_E_MISMATCHING_DATA)
        _EMX_ENUM_VAL(LLKA_E_BAD_CLASSIFICATION_CLUSTERS)
        _EMX_ENUM_VAL(LLKA_E_BAD_GOLDEN_STEPS)
        _EMX_ENUM_VAL(LLKA_E_BAD_CONFALS)
        _EMX_ENUM_VAL(LLKA_E_BAD_AVERAGE_NU_ANGLES)
        _EMX_ENUM_VAL(LLKA_E_BAD_CLASSIFICATION_LIMITS)
        _EMX_ENUM_VAL(LLKA_E_NO_FILE)
        _EMX_ENUM_VAL(LLKA_E_CANNOT_READ_FILE)
        _EMX_ENUM_VAL(LLKA_E_BAD_DATA)
        _EMX_ENUM_VAL(LLKA_E_NO_DATA)
        _EMX_ENUM_VAL(LLKA_E_NOTHING_TO_CLASSIFY)
    ;

    emscripten::function("errorToString", &LLKA::errorToString);

    //
    // Structure
    //

    emscripten::value_object<LLKA_Point>("Point")
        _EMX_VOF(x, LLKA_Point)
        _EMX_VOF(y, LLKA_Point)
        _EMX_VOF(z, LLKA_Point)
    ;
    emscripten::register_vector<LLKA_Point>("Points");

    emscripten::value_object<LLKA::Atom>("Atom")
        _EMX_VOF(type_symbol, LLKA::Atom)
        _EMX_VOF(label_atom_id, LLKA::Atom)
        _EMX_VOF(label_entity_id, LLKA::Atom)
        _EMX_VOF(label_comp_id, LLKA::Atom)
        _EMX_VOF(label_asym_id, LLKA::Atom)
        _EMX_VOF(auth_atom_id, LLKA::Atom)
        _EMX_VOF(auth_comp_id, LLKA::Atom)
        _EMX_VOF(auth_asym_id, LLKA::Atom)
        _EMX_VOF(coords, LLKA::Atom)
        _EMX_VOF(id, LLKA::Atom)
        _EMX_VOF(label_seq_id, LLKA::Atom)
        _EMX_VOF(auth_seq_id, LLKA::Atom)
        _EMX_VOF(pdbx_PDB_model_num, LLKA::Atom)
        _EMX_VOF(pdbx_PDB_ins_code, LLKA::Atom)
        _EMX_VOF(label_alt_id, LLKA::Atom)
    ;

    emscripten::register_vector<LLKA::Atom>("Structure");
    emscripten::register_vector<LLKA::Atom *>("StructureView"); // LLKA::Atom * is supposed to be const-qualified here but some versions of Emscripten have a problem with that
    emscripten::register_vector<LLKA::Structure>("Structures");

    _EMX_MK_RCRESULT(LLKA::Structure);
    _EMX_MK_RCRESULT(LLKA::Structures);

    emscripten::function("makeEmptyAtom", &LLKA::makeEmptyAtom);
    emscripten::function("makeAtom", &LLKA::makeAtom);
    emscripten::function("atomsEqual", &LLKA::atomsEqual);
    emscripten::function("atomMatches", &LLKA::atomMatches);
    emscripten::function("compareAtoms", &LLKA::compareAtoms);

    //
    // Measurements
    //

    emscripten::function("measureAnglef", &LLKA::measureAngle<float>);
    emscripten::function("measureAngle", &LLKA::measureAngle<double>);
    emscripten::function("measureAngleld", &LLKA::measureAngle<long double>);

    emscripten::function("measureDihedralf", &LLKA::measureDihedral<float>);
    emscripten::function("measureDihedral", &LLKA::measureDihedral<double>);
    emscripten::function("measureDihedralld", &LLKA::measureDihedral<long double>);

    emscripten::function("measureDistancef", &LLKA::measureDistance<float>);
    emscripten::function("measureDistance", &LLKA::measureDistance<double>);
    emscripten::function("measureDistanceld", &LLKA::measureDistance<long double>);

    //
    // Superposition
    //

    _EMX_MK_RCRESULT(double);

    emscripten::function("centroidPoints", emscripten::select_overload<LLKA_Point(const LLKA::Points&)>(&LLKA::centroid));
    emscripten::function("centroidStructure", emscripten::select_overload<LLKA_Point(const LLKA::Structure&)>(&LLKA::centroid));
    emscripten::function("rmsdPoints", emscripten::select_overload<LLKA::RCResult<double>(const LLKA::Points&, const LLKA::Points&)>(&LLKA::rmsd));
    emscripten::function("rmsd", emscripten::select_overload<LLKA::RCResult<double>(const LLKA::Structure&, const LLKA::Structure&)>(&LLKA::rmsd));
    emscripten::function("superposePoints", emscripten::select_overload<LLKA::RCResult<double>(LLKA::Points&, const LLKA::Points &)>(&LLKA::superpose));
    emscripten::function("superposeStructures", emscripten::select_overload<LLKA::RCResult<double>(LLKA::Structure &, const LLKA::Structure &)>(&LLKA::superpose));
    emscripten::function("splitByAltIds", &LLKA::splitByAltIds);
    emscripten::function("splitStructureToDinucleotideSteps", &LLKA::splitStructureToDinucleotideSteps);

    //
    // Segmentation
    //

    emscripten::register_map<int32_t, LLKA::StructureSegments::Residue>("StructureSegmentsResidueMap");
    emscripten::register_map<std::string, LLKA::StructureSegments::Chain>("StructureSegmentsChainMap");
    emscripten::register_map<int32_t, LLKA::StructureSegments::Model>("StructureSegmentsModelMap");

    emscripten::value_object<LLKA::StructureSegments::Residue>("StructureSegmentsResidue")
        _EMX_VOF(atoms, LLKA::StructureSegments::Residue)
    ;

    emscripten::value_object<LLKA::StructureSegments::Chain>("StructureSegmentsChain")
        _EMX_VOF(residues, LLKA::StructureSegments::Chain)
    ;

    emscripten::value_object<LLKA::StructureSegments::Model>("StructureSegmentsModel")
        _EMX_VOF(chains, LLKA::StructureSegments::Model)
    ;

    emscripten::class_<LLKA::StructureSegments>("StructureSegments")
        .constructor<>()
        .constructor<LLKA::Structure &>()
        _EMX_CLS_PROP(models, LLKA::StructureSegments)
        _EMX_CLS_PROP_READONLY(structure, LLKA::StructureSegments)
    ;

    //
    // NtC
    //

    emscripten::enum_<LLKA_CrossResidueMetric>("CrossResidueMetric")
        _EMX_ENUM_VAL(LLKA_XR_DIST_CC)
        _EMX_ENUM_VAL(LLKA_XR_DIST_NN)
        _EMX_ENUM_VAL(LLKA_XR_TOR_MU)
    ;

    emscripten::enum_<LLKA_DinucleotideTorsion>("DinucleotideTorsion")
        _EMX_ENUM_VAL(LLKA_TOR_DELTA_1)
        _EMX_ENUM_VAL(LLKA_TOR_EPSILON_1)
        _EMX_ENUM_VAL(LLKA_TOR_ZETA_1)
        _EMX_ENUM_VAL(LLKA_TOR_ALPHA_2)
        _EMX_ENUM_VAL(LLKA_TOR_BETA_2)
        _EMX_ENUM_VAL(LLKA_TOR_GAMMA_2)
        _EMX_ENUM_VAL(LLKA_TOR_DELTA_2)
        _EMX_ENUM_VAL(LLKA_TOR_CHI_1)
        _EMX_ENUM_VAL(LLKA_TOR_CHI_2)
    ;

    emscripten::enum_<LLKA_CANA>("CANA")
        _EMX_ENUM_VAL(LLKA_AAA)
        _EMX_ENUM_VAL(LLKA_AAw)
        _EMX_ENUM_VAL(LLKA_AAu)
        _EMX_ENUM_VAL(LLKA_AB)
        _EMX_ENUM_VAL(LLKA_BA)
        _EMX_ENUM_VAL(LLKA_BBB)
        _EMX_ENUM_VAL(LLKA_BBw)
        _EMX_ENUM_VAL(LLKA_B12)
        _EMX_ENUM_VAL(LLKA_BB2)
        _EMX_ENUM_VAL(LLKA_miB)
        _EMX_ENUM_VAL(LLKA_ICL)
        _EMX_ENUM_VAL(LLKA_OPN)
        _EMX_ENUM_VAL(LLKA_SYN)
        _EMX_ENUM_VAL(LLKA_ZZZ)
        _EMX_ENUM_VAL(LLKA_INVALID_CANA)
    ;

    emscripten::enum_<LLKA_NtC>("NtC")
        _EMX_ENUM_VAL(LLKA_AA00)
        _EMX_ENUM_VAL(LLKA_AA02)
        _EMX_ENUM_VAL(LLKA_AA03)
        _EMX_ENUM_VAL(LLKA_AA04)
        _EMX_ENUM_VAL(LLKA_AA08)
        _EMX_ENUM_VAL(LLKA_AA09)
        _EMX_ENUM_VAL(LLKA_AA01)
        _EMX_ENUM_VAL(LLKA_AA05)
        _EMX_ENUM_VAL(LLKA_AA06)
        _EMX_ENUM_VAL(LLKA_AA10)
        _EMX_ENUM_VAL(LLKA_AA11)
        _EMX_ENUM_VAL(LLKA_AA07)
        _EMX_ENUM_VAL(LLKA_AA12)
        _EMX_ENUM_VAL(LLKA_AA13)
        _EMX_ENUM_VAL(LLKA_AB01)
        _EMX_ENUM_VAL(LLKA_AB02)
        _EMX_ENUM_VAL(LLKA_AB03)
        _EMX_ENUM_VAL(LLKA_AB04)
        _EMX_ENUM_VAL(LLKA_AB05)
        _EMX_ENUM_VAL(LLKA_BA01)
        _EMX_ENUM_VAL(LLKA_BA05)
        _EMX_ENUM_VAL(LLKA_BA09)
        _EMX_ENUM_VAL(LLKA_BA08)
        _EMX_ENUM_VAL(LLKA_BA10)
        _EMX_ENUM_VAL(LLKA_BA13)
        _EMX_ENUM_VAL(LLKA_BA16)
        _EMX_ENUM_VAL(LLKA_BA17)
        _EMX_ENUM_VAL(LLKA_BB00)
        _EMX_ENUM_VAL(LLKA_BB01)
        _EMX_ENUM_VAL(LLKA_BB17)
        _EMX_ENUM_VAL(LLKA_BB02)
        _EMX_ENUM_VAL(LLKA_BB03)
        _EMX_ENUM_VAL(LLKA_BB11)
        _EMX_ENUM_VAL(LLKA_BB16)
        _EMX_ENUM_VAL(LLKA_BB04)
        _EMX_ENUM_VAL(LLKA_BB05)
        _EMX_ENUM_VAL(LLKA_BB07)
        _EMX_ENUM_VAL(LLKA_BB08)
        _EMX_ENUM_VAL(LLKA_BB10)
        _EMX_ENUM_VAL(LLKA_BB12)
        _EMX_ENUM_VAL(LLKA_BB13)
        _EMX_ENUM_VAL(LLKA_BB14)
        _EMX_ENUM_VAL(LLKA_BB15)
        _EMX_ENUM_VAL(LLKA_BB20)
        _EMX_ENUM_VAL(LLKA_IC01)
        _EMX_ENUM_VAL(LLKA_IC02)
        _EMX_ENUM_VAL(LLKA_IC03)
        _EMX_ENUM_VAL(LLKA_IC04)
        _EMX_ENUM_VAL(LLKA_IC05)
        _EMX_ENUM_VAL(LLKA_IC06)
        _EMX_ENUM_VAL(LLKA_IC07)
        _EMX_ENUM_VAL(LLKA_OP01)
        _EMX_ENUM_VAL(LLKA_OP02)
        _EMX_ENUM_VAL(LLKA_OP03)
        _EMX_ENUM_VAL(LLKA_OP04)
        _EMX_ENUM_VAL(LLKA_OP05)
        _EMX_ENUM_VAL(LLKA_OP06)
        _EMX_ENUM_VAL(LLKA_OP07)
        _EMX_ENUM_VAL(LLKA_OP08)
        _EMX_ENUM_VAL(LLKA_OP09)
        _EMX_ENUM_VAL(LLKA_OP10)
        _EMX_ENUM_VAL(LLKA_OP11)
        _EMX_ENUM_VAL(LLKA_OP12)
        _EMX_ENUM_VAL(LLKA_OP13)
        _EMX_ENUM_VAL(LLKA_OP14)
        _EMX_ENUM_VAL(LLKA_OP15)
        _EMX_ENUM_VAL(LLKA_OP16)
        _EMX_ENUM_VAL(LLKA_OP17)
        _EMX_ENUM_VAL(LLKA_OP18)
        _EMX_ENUM_VAL(LLKA_OP19)
        _EMX_ENUM_VAL(LLKA_OP20)
        _EMX_ENUM_VAL(LLKA_OP21)
        _EMX_ENUM_VAL(LLKA_OP22)
        _EMX_ENUM_VAL(LLKA_OP23)
        _EMX_ENUM_VAL(LLKA_OP24)
        _EMX_ENUM_VAL(LLKA_OP25)
        _EMX_ENUM_VAL(LLKA_OP26)
        _EMX_ENUM_VAL(LLKA_OP27)
        _EMX_ENUM_VAL(LLKA_OP28)
        _EMX_ENUM_VAL(LLKA_OP29)
        _EMX_ENUM_VAL(LLKA_OP30)
        _EMX_ENUM_VAL(LLKA_OP31)
        _EMX_ENUM_VAL(LLKA_OPS1)
        _EMX_ENUM_VAL(LLKA_OP1S)
        _EMX_ENUM_VAL(LLKA_AAS1)
        _EMX_ENUM_VAL(LLKA_AB1S)
        _EMX_ENUM_VAL(LLKA_AB2S)
        _EMX_ENUM_VAL(LLKA_BB1S)
        _EMX_ENUM_VAL(LLKA_BB2S)
        _EMX_ENUM_VAL(LLKA_BBS1)
        _EMX_ENUM_VAL(LLKA_ZZ01)
        _EMX_ENUM_VAL(LLKA_ZZ02)
        _EMX_ENUM_VAL(LLKA_ZZ1S)
        _EMX_ENUM_VAL(LLKA_ZZ2S)
        _EMX_ENUM_VAL(LLKA_ZZS1)
        _EMX_ENUM_VAL(LLKA_ZZS2)
        _EMX_ENUM_VAL(LLKA_INVALID_NTC)
    ;

    _EMX_MK_RCRESULT(LLKA::AtomNameQuad);

    emscripten::class_<LLKA::AtomNameQuad>("AtomNameQuad")
        .constructor<>()
        .constructor<
            std::string, std::string, std::string, std::string
        >()
        _EMX_CLS_PROP(a, LLKA::AtomNameQuad)
        _EMX_CLS_PROP(b, LLKA::AtomNameQuad)
        _EMX_CLS_PROP(c, LLKA::AtomNameQuad)
        _EMX_CLS_PROP(d, LLKA::AtomNameQuad)
    ;

    emscripten::value_object<LLKA::BackboneAtom>("BackboneAtom")
        _EMX_VOF(residue, LLKA::BackboneAtom)
        _EMX_VOF(name, LLKA::BackboneAtom)
    ;

    emscripten::enum_<LLKA_BaseKind>("BaseKind")
        _EMX_ENUM_VAL(LLKA_PURINE)
        _EMX_ENUM_VAL(LLKA_PYRIMIDINE)
        _EMX_ENUM_VAL(LLKA_NON_STANDARD_BASE)
    ;

    emscripten::value_object<LLKA_StepInfo>("StepInfo")
        _EMX_VOF(firstSeqId, LLKA_StepInfo)
        _EMX_VOF(firstSeqIdAuth, LLKA_StepInfo)
        _EMX_VOF(secondSeqId, LLKA_StepInfo)
        _EMX_VOF(secondSeqIdAuth, LLKA_StepInfo)
        _EMX_VOF(firstBaseKind, LLKA_StepInfo)
        _EMX_VOF(secondBaseKind, LLKA_StepInfo)
    ;

    emscripten::value_object<LLKA_StepMetrics>("StepMetrics")
        _EMX_VOF(delta_1, LLKA_StepMetrics)
        _EMX_VOF(epsilon_1, LLKA_StepMetrics)
        _EMX_VOF(zeta_1, LLKA_StepMetrics)
        _EMX_VOF(alpha_2, LLKA_StepMetrics)
        _EMX_VOF(beta_2, LLKA_StepMetrics)
        _EMX_VOF(gamma_2, LLKA_StepMetrics)
        _EMX_VOF(delta_2, LLKA_StepMetrics)
        _EMX_VOF(chi_1, LLKA_StepMetrics)
        _EMX_VOF(chi_2, LLKA_StepMetrics)
        _EMX_VOF(CC, LLKA_StepMetrics)
        _EMX_VOF(NN, LLKA_StepMetrics)
        _EMX_VOF(mu, LLKA_StepMetrics)
    ;

    _EMX_MK_RCRESULT(LLKA_StepInfo);

    emscripten::function("backboneAtomIndex", &LLKA::backboneAtomIndex);
    emscripten::function("calculateStepMetrics", &LLKA::calculateStepMetrics);
    emscripten::function("crossResidueMetric", &LLKA::crossResidueMetric);
    emscripten::function(
        "crossResidueMetricAtomsFromBases",
        emscripten::select_overload<LLKA::RCResult<LLKA::AtomNameQuad>(const std::string &, const std::string &, LLKA_CrossResidueMetric)>(&LLKA::crossResidueMetricAtoms)
    );
    emscripten::function(
        "crossResidueMetricAtomsFromStructure",
        emscripten::select_overload<LLKA::RCResult<LLKA::AtomNameQuad>(const LLKA::Structure &, LLKA_CrossResidueMetric)>(&LLKA::crossResidueMetricAtoms)
    );
    emscripten::function("crossResidueMetricName", &LLKA::crossResidueMetricName);
    emscripten::function("dinucleotideTorsion", &LLKA::dinucleotideTorsion);
    emscripten::function(
        "dinucleotideTorsionAtomsFromBases",
        emscripten::select_overload<LLKA::RCResult<LLKA::AtomNameQuad>(const std::string &, const std::string &, LLKA_DinucleotideTorsion)>(&LLKA::dinucleotideTorsionAtoms)
    );
    emscripten::function(
        "dinucleotideTorsionAtomsFromStructure",
        emscripten::select_overload<LLKA::RCResult<LLKA::AtomNameQuad>(const LLKA::Structure &, LLKA_DinucleotideTorsion)>(&LLKA::dinucleotideTorsionAtoms)
    );
    emscripten::function("dinucleotideTorsionName", &LLKA::dinucleotideTorsionName);
    emscripten::function("extractBackbone", &LLKA::extractBackbone);
    emscripten::function("extractExtendedBackbone", &LLKA::extractExtendedBackbone);
    emscripten::function("extractMetricsStructure", &LLKA::extractMetricsStructure);
    emscripten::function("nameToCANA", &LLKA::nameToCANA);
    emscripten::function("CANAToName", &LLKA::CANAToName);
    emscripten::function("nameToNtC", &LLKA::nameToNtC);
    emscripten::function("NtCToName", &LLKA::NtCToName);
    emscripten::function("structureIsStep", &LLKA::structureIsStep);

    //
    // Connectivity, similarity
    //

    emscripten::register_vector<LLKA_NtC>("NtCs");

    emscripten::value_object<LLKA_Connectivity>("Connectivity")
        _EMX_VOF(C5PrimeDistance, LLKA_Connectivity)
        _EMX_VOF(O3PrimeDistance, LLKA_Connectivity)
    ;
    emscripten::register_vector<LLKA_Connectivity>("Connectivities");

    emscripten::value_object<LLKA_Similarity>("Similarity")
        _EMX_VOF(rmsd, LLKA_Similarity)
        _EMX_VOF(euclideanDistance, LLKA_Similarity)
    ;
    emscripten::register_vector<LLKA_Similarity>("Similarities");

    _EMX_MK_RCRESULT(LLKA_Connectivity);
    _EMX_MK_RCRESULT(LLKA::Connectivities);
    _EMX_MK_RCRESULT(LLKA_Similarity);
    _EMX_MK_RCRESULT(LLKA::Similarities);

    emscripten::function(
        "measureStepConnectivityNtCs",
        emscripten::select_overload<LLKA::RCResult<LLKA_Connectivity>(const LLKA::Structure &, LLKA_NtC, const LLKA::Structure &, LLKA_NtC)>(&LLKA::measureStepConnectivity)
    );
    emscripten::function(
        "measureStepConnectivityNtCsMultipleFirst",
        emscripten::select_overload<LLKA::RCResult<LLKA::Connectivities>(const LLKA::Structure &, std::vector<LLKA_NtC>, const LLKA::Structure &, LLKA_NtC)>(&LLKA::measureStepConnectivity)
    );
    emscripten::function(
        "measureStepConnectivityNtCsMultipleSecond",
        emscripten::select_overload<LLKA::RCResult<LLKA::Connectivities>(const LLKA::Structure &, LLKA_NtC, const LLKA::Structure &, std::vector<LLKA_NtC>)>(&LLKA::measureStepConnectivity)
    );
    emscripten::function(
        "measureStepConnectivityStructures",
        emscripten::select_overload<LLKA::RCResult<LLKA_Connectivity>(const LLKA::Structure &, const LLKA::Structure &, const LLKA::Structure &, const LLKA::Structure &)>(&LLKA::measureStepConnectivity)
    );
    emscripten::function(
        "measureStepConnectivityStructuresMultiple",
        emscripten::select_overload<LLKA::RCResult<LLKA::Connectivities>(const LLKA::Structure &, const LLKA::Structure &, const LLKA::Structure &, const LLKA::Structures &)>(&LLKA::measureStepConnectivity)
    );
    emscripten::function(
        "measureStepSimilarityNtC",
        emscripten::select_overload<LLKA::RCResult<LLKA_Similarity>(const LLKA::Structure &, LLKA_NtC)>(&LLKA::measureStepSimilarity)
    );
    emscripten::function(
        "measureStepSimilarityNtCMultiple",
        emscripten::select_overload<LLKA::RCResult<LLKA::Similarities>(const LLKA::Structure &, std::vector<LLKA_NtC>)>(&LLKA::measureStepSimilarity)
    );
    emscripten::function(
        "measureStepSimilarityNtCStructure",
        emscripten::select_overload<LLKA::RCResult<LLKA_Similarity>(const LLKA::Structure &, const LLKA::Structure &)>(&LLKA::measureStepSimilarity)
    );

    emscripten::function("makeStdVectorDouble", &LLKA::makeStdVector<double>);
    emscripten::function("makeStdVectorAtom", &LLKA::makeStdVector<LLKA::Atom>);
    emscripten::function("makeStdVectorAtomPtr", &LLKA::makeStdVector<LLKA::Atom *>);
    emscripten::function("makeStdVectorNtC", &LLKA::makeStdVector<LLKA_NtC>);
    emscripten::function("makeStdVectorStructure", &LLKA::makeStdVector<LLKA::Atom>);
    emscripten::function("makeStdVectorStructures", &LLKA::makeStdVector<LLKA::Structure>);
    emscripten::function("makeStdVectorString", &LLKA::makeStdVector<std::string>);
    emscripten::function("makeStdVectorInt32", &LLKA::makeStdVector<int32_t>);

    //
    // MiniCif
    //

    emscripten::function("makeStdVectorCifDataValue", &LLKA::makeStdVector<LLKA::CifData::Value>);
    emscripten::function("makeStdVectorCifDataItem", &LLKA::makeStdVector<LLKA::CifData::Item>);
    emscripten::function("makeStdVectorCifDataCategory", &LLKA::makeStdVector<LLKA::CifData::Category>);
    emscripten::function("makeStdVectorCifDataBlock", &LLKA::makeStdVector<LLKA::CifData::Block>);

    emscripten::enum_<LLKA_CifDataValueState>("CifDataValueState")
        _EMX_ENUM_VAL(LLKA_MINICIF_VALUE_SET)
        _EMX_ENUM_VAL(LLKA_MINICIF_VALUE_NONE)
        _EMX_ENUM_VAL(LLKA_MINICIF_VALUE_UNKW)
    ;

    emscripten::value_object<LLKA::CifData::Value>("CifDataValue")
        _EMX_VOF(text, LLKA::CifData::Value)
        _EMX_VOF(state, LLKA::CifData::Value)
    ;
    emscripten::register_vector<LLKA::CifData::Value>("CifDataValues");

    emscripten::value_object<LLKA::CifData::Item>("CifDataItem")
        _EMX_VOF(keyword, LLKA::CifData::Item)
        _EMX_VOF(values, LLKA::CifData::Item)
    ;
    emscripten::register_vector<LLKA::CifData::Item>("CifDataItems");

    emscripten::class_<LLKA::CifData::Category>("CifDataCategory")
        .constructor<>()
        .constructor<std::string, LLKA::CifData::Items>()
        _EMX_CLS_PROP(name, LLKA::CifData::Category)
        _EMX_CLS_PROP(items, LLKA::CifData::Category)
        _EMX_CLS_MTH(bool, LLKA::CifData::Category, isLoop)
    ;
    emscripten::register_vector<LLKA::CifData::Category>("CifDataCategories");

    emscripten::value_object<LLKA::CifData::Block>("CifDataBlock")
        _EMX_VOF(name, LLKA::CifData::Block)
        _EMX_VOF(categories, LLKA::CifData::Block)
    ;
    emscripten::register_vector<LLKA::CifData::Block>("CifDataBlocks");

    emscripten::value_object<LLKA::CifData>("CifData")
        _EMX_VOF(blocks, LLKA::CifData)
    ;

    emscripten::value_object<LLKA::ImportedStructure>("ImportedStructure")
        _EMX_VOF(id, LLKA::ImportedStructure)
        _EMX_VOF(structure, LLKA::ImportedStructure)
        _EMX_VOF(cifData, LLKA::ImportedStructure)
    ;

    emscripten::value_object<LLKA::CifError>("CifError")
        _EMX_VOF(tRet, LLKA::CifError)
        _EMX_VOF(error, LLKA::CifError)
    ;

    emscripten::class_<LLKA::CifResult>("CifResult")
        .constructor<LLKA::ImportedStructure>()
        .constructor<LLKA::CifError, bool>() \
        .function("failure", &LLKA::CifResult::failure) \
        .function("success", std::function<typename LLKA::ResultSuccessReturnType<LLKA::ImportedStructure, std::is_copy_constructible_v<LLKA::ImportedStructure>>::RT (LLKA::CifResult&)>(&LLKA::CifResult::success))
        .function("const_success", std::function<const LLKA::ImportedStructure &(const LLKA::CifResult&)>(&LLKA::CifResult::const_success))
        .function("isSuccess", &LLKA::CifResult::isSuccess)
    ;

    // Overload that takes path to file as a parameter is disabled in Emscripten builds
    // so we do not need to use select_overload here

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
    // Read from file
    emscripten::function(
        "cifToStructueFromFile",
        emscripten::select_overload<LLKA::CifResult(const std::filesystem::path &, int32_t options)>(&LLKA::cifToStructure)
    );
#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

    // Read from text
    emscripten::function(
        "cifToStructure",
        emscripten::select_overload<LLKA::CifResult(const std::string &, int32_t options)>(&LLKA::cifToStructure)
    );

    emscripten::function("cifDataToString", &LLKA::cifDataToString);

    _EMX_MK_RCRESULT(std::string);

    //
    // Nucleotide
    //

    emscripten::enum_<LLKA_SugarPucker>("SugarPucker")
        _EMX_ENUM_VAL(LLKA_C3_ENDO)
        _EMX_ENUM_VAL(LLKA_C4_EXO)
        _EMX_ENUM_VAL(LLKA_O4_ENDO)
        _EMX_ENUM_VAL(LLKA_C1_EXO)
        _EMX_ENUM_VAL(LLKA_C2_ENDO)
        _EMX_ENUM_VAL(LLKA_C3_EXO)
        _EMX_ENUM_VAL(LLKA_C4_ENDO)
        _EMX_ENUM_VAL(LLKA_O4_EXO)
        _EMX_ENUM_VAL(LLKA_C1_ENDO)
        _EMX_ENUM_VAL(LLKA_C2_EXO)
        _EMX_ENUM_VAL(LLKA_INVALID_SUGAR_PUCKER)
    ;

    emscripten::enum_<LLKA_SugarPuckerNameBrevity>("SugarPuckerNameBrevity")
        _EMX_ENUM_VAL(LLKA_SPN_VERY_TERSE)
        _EMX_ENUM_VAL(LLKA_SPN_TERSE)
        _EMX_ENUM_VAL(LLKA_SPN_FANCY)
    ;

    emscripten::value_object<LLKA_RiboseMetrics>("RiboseMetrics")
        _EMX_VOF(nus, LLKA_RiboseMetrics)
        _EMX_VOF(P, LLKA_RiboseMetrics)
        _EMX_VOF(tMax, LLKA_RiboseMetrics)
        _EMX_VOF(pucker, LLKA_RiboseMetrics)
    ;

    _EMX_MK_RCRESULT(LLKA_RiboseMetrics);
    _EMX_MK_RCRESULT(LLKA_SugarPucker);

    emscripten::function("extractNucleotide", &LLKA::extractNucleotide);
    emscripten::function("extractNucleotideView", &LLKA::extractNucleotideView);
    emscripten::function("extractRibose", &LLKA::extractRibose);
    emscripten::function("extractRiboseView", &LLKA::extractRiboseView);
    emscripten::function("isNucleotideCompound", &LLKA::isNucleotideCompound);
    emscripten::function("nameToSugarPucker", &LLKA::nameToSugarPucker);
    emscripten::function(
        "riboseMetrics",
        emscripten::select_overload<LLKA::RCResult<LLKA_RiboseMetrics>(const LLKA::Structure &)>(&LLKA::riboseMetrics)
    );
    emscripten::function(
        "riboseMetricsView",
        emscripten::select_overload<LLKA::RCResult<LLKA_RiboseMetrics>(const LLKA::StructureView &)>(&LLKA::riboseMetrics)
    );
    emscripten::function(
        "sugarPucker",
        emscripten::select_overload<LLKA::RCResult<LLKA_SugarPucker>(const LLKA::Structure &)>(&LLKA::sugarPucker)
    );
    emscripten::function(
        "sugarPuckerView",
        emscripten::select_overload<LLKA::RCResult<LLKA_SugarPucker>(const LLKA::StructureView &)>(&LLKA::sugarPucker)
    );
    emscripten::function("sugarPuckerToName", &LLKA::sugarPuckerToName);

    //
    // Classification
    //

    // We try to use raw C structs directly where possible
    emscripten::value_object<LLKA_NuAngles>("NuAngles")
        _EMX_VOF(nu_0, LLKA_NuAngles)
        _EMX_VOF(nu_1, LLKA_NuAngles)
        _EMX_VOF(nu_2, LLKA_NuAngles)
        _EMX_VOF(nu_3, LLKA_NuAngles)
        _EMX_VOF(nu_4, LLKA_NuAngles)
    ;

    emscripten::value_object<LLKA_NuAnglesMetrics>("NuAnglesMetrics")
        _EMX_VOF(nu_0, LLKA_NuAnglesMetrics)
        _EMX_VOF(nu_1, LLKA_NuAnglesMetrics)
        _EMX_VOF(nu_2, LLKA_NuAnglesMetrics)
        _EMX_VOF(nu_3, LLKA_NuAnglesMetrics)
        _EMX_VOF(nu_4, LLKA_NuAnglesMetrics)
    ;

    emscripten::value_object<LLKA_ClusterNuAngles>("ClusterNuAngles")
        _EMX_VOF(firstNucleotide, LLKA_ClusterNuAngles)
        _EMX_VOF(secondNucleotide, LLKA_ClusterNuAngles)
        _EMX_VOF(clusterNumber, LLKA_ClusterNuAngles)
    ;

    emscripten::value_object<LLKA_ClassificationLimits>("ClassificationLimits")
        _EMX_VOF(minimumNearestNeighbors, LLKA_ClassificationLimits)
        _EMX_VOF(numberOfUsedNearestNeighbors, LLKA_ClassificationLimits)
        _EMX_VOF(averageNeighborsTorsionCutoff, LLKA_ClassificationLimits)
        _EMX_VOF(nearestNeighborTorsionsCutoff, LLKA_ClassificationLimits)
        _EMX_VOF(totalDistanceCutoff, LLKA_ClassificationLimits)
        _EMX_VOF(pseudorotationCutoff, LLKA_ClassificationLimits)
        _EMX_VOF(minimumClusterVotes, LLKA_ClassificationLimits)
    ;

    emscripten::value_object<LLKA_ClassificationMetric>("ClassificationMetric")
        _EMX_VOF(deviation, LLKA_ClassificationMetric)
        _EMX_VOF(minValue, LLKA_ClassificationMetric)
        _EMX_VOF(meanValue, LLKA_ClassificationMetric)
        _EMX_VOF(maxValue, LLKA_ClassificationMetric)
    ;

    emscripten::value_object<LLKA_ClassificationCluster>("ClassificationCluster")
        _EMX_VOF(delta_1, LLKA_ClassificationCluster)
        _EMX_VOF(epsilon_1,LLKA_ClassificationCluster)
        _EMX_VOF(zeta_1,LLKA_ClassificationCluster)
        _EMX_VOF(alpha_2, LLKA_ClassificationCluster)
        _EMX_VOF(beta_2, LLKA_ClassificationCluster)
        _EMX_VOF(gamma_2, LLKA_ClassificationCluster)
        _EMX_VOF(delta_2, LLKA_ClassificationCluster)
        _EMX_VOF(chi_1, LLKA_ClassificationCluster)
        _EMX_VOF(chi_2, LLKA_ClassificationCluster)
        _EMX_VOF(CC, LLKA_ClassificationCluster)
        _EMX_VOF(NN, LLKA_ClassificationCluster)
        _EMX_VOF(mu, LLKA_ClassificationCluster)
        _EMX_VOF(nusFirst, LLKA_ClassificationCluster)
        _EMX_VOF(nusSecond, LLKA_ClassificationCluster)
        _EMX_VOF(ribosePseudorotation_1, LLKA_ClassificationCluster)
        _EMX_VOF(ribosePseudorotation_2, LLKA_ClassificationCluster)
        _EMX_VOF(number, LLKA_ClassificationCluster)
        _EMX_VOF(NtC, LLKA_ClassificationCluster)
        _EMX_VOF(CANA, LLKA_ClassificationCluster)
    ;

    emscripten::value_object<LLKA_ConfalScore>("ConfalScore")
        _EMX_VOF(delta_1, LLKA_ConfalScore)
        _EMX_VOF(epsilon_1,  LLKA_ConfalScore)
        _EMX_VOF(zeta_1,  LLKA_ConfalScore)
        _EMX_VOF(alpha_2, LLKA_ConfalScore)
        _EMX_VOF(beta_2, LLKA_ConfalScore)
        _EMX_VOF(gamma_2, LLKA_ConfalScore)
        _EMX_VOF(delta_2, LLKA_ConfalScore)
        _EMX_VOF(chi_1, LLKA_ConfalScore)
        _EMX_VOF(chi_2, LLKA_ConfalScore)
        _EMX_VOF(CC, LLKA_ConfalScore)
        _EMX_VOF(NN, LLKA_ConfalScore)
        _EMX_VOF(mu, LLKA_ConfalScore)
        _EMX_VOF(total, LLKA_ConfalScore)
    ;

    emscripten::value_object<LLKA_Confal>("Confal")
        _EMX_VOF(delta_1, LLKA_Confal)
        _EMX_VOF(epsilon_1, LLKA_Confal)
        _EMX_VOF(zeta_1, LLKA_Confal)
        _EMX_VOF(alpha_2, LLKA_Confal)
        _EMX_VOF(beta_2, LLKA_Confal)
        _EMX_VOF(gamma_2, LLKA_Confal)
        _EMX_VOF(delta_2, LLKA_Confal)
        _EMX_VOF(chi_1, LLKA_Confal)
        _EMX_VOF(chi_2, LLKA_Confal)
        _EMX_VOF(CC, LLKA_Confal)
        _EMX_VOF(NN, LLKA_Confal)
        _EMX_VOF(mu, LLKA_Confal)
        _EMX_VOF(nusFirst, LLKA_Confal)
        _EMX_VOF(nusSecond, LLKA_Confal)
        _EMX_VOF(clusterNumber, LLKA_Confal)
    ;

    emscripten::value_object<LLKA_ConfalPercentile>("ConfalPercentile")
        _EMX_VOF(value, LLKA_ConfalPercentile)
    ;

    emscripten::value_object<LLKA::GoldenStep>("GoldenStep")
        _EMX_VOF(pucker_1, LLKA::GoldenStep)
        _EMX_VOF(pucker_2, LLKA::GoldenStep)
        _EMX_VOF(nuAngles_1, LLKA::GoldenStep)
        _EMX_VOF(nuAngles_2, LLKA::GoldenStep)
        _EMX_VOF(metrics, LLKA::GoldenStep)
        _EMX_VOF(name, LLKA::GoldenStep)
        _EMX_VOF(clusterIdx, LLKA::GoldenStep)
    ;

    //
    // Resource loaders
    //

    emscripten::register_vector<LLKA_ClusterNuAngles>("LLKA_ClusterNuAnglesVec");
    emscripten::register_vector<LLKA_ClassificationCluster>("LLKA_ClassificationClusterVec");
    emscripten::register_vector<LLKA_Confal>("LLKA_ConfalVec");
    emscripten::register_vector<LLKA::GoldenStep>("GoldenStepVec");

    _EMX_MK_RCRESULT(std::vector<LLKA_ClusterNuAngles>);
    _EMX_MK_RCRESULT(std::vector<LLKA_ClassificationCluster>);
    _EMX_MK_RCRESULT(std::vector<LLKA_Confal>);
    _EMX_MK_RCRESULT(std::vector<LLKA_ConfalPercentile>);
    _EMX_MK_RCRESULT(std::vector<LLKA::GoldenStep>);

    // Read data from files
#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
    emscripten::function(
        "loadClusterNuAnglesFromFile",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA_ClusterNuAngles>>(const std::filesystem::path &)>(&LLKA::loadClusterNuAngles)
    );
    emscripten::function(
        "loadClustersFromFile",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA_ClassificationCluster>>(const std::filesystem::path &)>(&LLKA::loadClusters)
    );
    emscripten::function(
        "loadConfalsFromFile",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA_Confal>>(const std::filesystem::path &)>(&LLKA::loadConfals)
    );
    emscripten::function(
        "loadConfalPercentilesFromFile",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA_ConfalPercentile>>(const std::filesystem::path &)>(&LLKA::loadConfalPercentiles)
    );
    emscripten::function(
        "loadGoldenStepsFromFile",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA::GoldenStep>>(const std::filesystem::path &)>(&LLKA::loadGoldenSteps)
    );
#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

    // Read data from text
    emscripten::function(
        "loadClusterNuAngles",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA_ClusterNuAngles>>(const std::string &)>(&LLKA::loadClusterNuAngles)
    );
    emscripten::function(
        "loadClusters",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA_ClassificationCluster>>(const std::string &)>(&LLKA::loadClusters)
    );
    emscripten::function(
        "loadConfals",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA_Confal>>(const std::string &)>(&LLKA::loadConfals)
    );
    emscripten::function(
        "loadConfalPercentiles",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA_ConfalPercentile>>(const std::string &)>(&LLKA::loadConfalPercentiles)
    );
    emscripten::function(
        "loadGoldenSteps",
        emscripten::select_overload<LLKA::RCResult<std::vector<LLKA::GoldenStep>>(const std::string &)>(&LLKA::loadGoldenSteps)
    );

    //
    // Tracing
    //

    emscripten::value_object<LLKA::TracepointInfo>("TracepointInfo")
        _EMX_VOF(TPID, LLKA::TracepointInfo)
        _EMX_VOF(description, LLKA::TracepointInfo)
    ;

    emscripten::function("toggleAllTracepoints", &LLKA::toggleAllTracepoints);
    emscripten::function("toggleTracepoint", &LLKA::toggleTracepoint);
    emscripten::function("trace", &LLKA::trace);
    emscripten::function("tracepointInfo", &LLKA::tracepointInfo);
    emscripten::function("tracepointState", &LLKA::tracepointState);

    //
    // Classification
    //

    _EMX_MK_RCRESULT(LLKA_ClassificationCluster);
    _EMX_MK_RCRESULT(LLKA_Confal);

    emscripten::value_object<LLKA_AverageConfal>("AverageConfal")
        _EMX_VOF(score, LLKA_AverageConfal)
        _EMX_VOF(percentile, LLKA_AverageConfal)
    ;

    emscripten::class_<LLKA::ClassificationContext>("ClassificationContext")
        .constructor<>()
        .constructor<LLKA_ClassificationContext *>()
        _EMX_CLS_MTH(boo, LLKA::ClassificationContext, isValid)
    ;

    _EMX_MK_RCRESULT(LLKA::ClassificationContext);

    emscripten::class_<LLKA::ClassifiedStep>("ClassifiedStep")
        .constructor<>()
        .constructor<LLKA_ClassifiedStep &>()
        _EMX_CLS_PROP(assignedNtC, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(assignedCANA, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(closestNtC, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(closestCANA, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(confalScore, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(euclideanDistanceNtCIdeal, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(metrics, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(differencesFromNtCAverages, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(nuAngles_1, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(nuAngles_2, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(ribosePseudorotation_1, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(ribosePseudorotation_2, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(tau_1, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(tau_2, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(sugarPucker_1, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(sugarPucker_2, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(nuAngleDifferences_1, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(nuAngleDifferences_2, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(rmsdToClosestNtC, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(closestGoldenStep, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(violations, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(violatingTorsionsAverage, LLKA::ClassifiedStep)
        _EMX_CLS_PROP(violatingTorsionsNearest, LLKA::ClassifiedStep)
        _EMX_CLS_MTH(bool, LLKA::ClassifiedStep, hasViolations)
        _EMX_CLS_MTH(std::vector<std::string>, LLKA::ClassifiedStep, namedViolations)
    ;

    emscripten::value_object<LLKA::AttemptedClassifiedStep>("AttemptedClassifiedStep")
        _EMX_VOF(status, LLKA::AttemptedClassifiedStep)
        _EMX_VOF(step, LLKA::AttemptedClassifiedStep)
    ;

    emscripten::register_vector<LLKA::ClassifiedStep>("ClassifiedSteps");
    emscripten::register_vector<LLKA_ConfalPercentile>("ConfalPercentiles");
    emscripten::register_vector<LLKA::AttemptedClassifiedStep>("AttemptedClassifiedSteps");

    _EMX_MK_RCRESULT(LLKA::AttemptedClassifiedSteps);

    emscripten::function(
        "averageConfal",
        emscripten::select_overload<LLKA_AverageConfal(const std::vector<LLKA::ClassifiedStep> &, const LLKA::ClassificationContext &)>(&LLKA::averageConfal)
    );
    emscripten::function(
        "averageConfalAttempted",
        emscripten::select_overload<LLKA_AverageConfal(const std::vector<LLKA::AttemptedClassifiedStep> &, const LLKA::ClassificationContext &)>(&LLKA::averageConfal)
    );
    emscripten::function("classificationClusterForNtC", &LLKA::classificationClusterForNtC);
    emscripten::function("confalForNtC", &LLKA::confalForNtC);
    emscripten::function("confalPercentile", &LLKA::confalPercentile);
    emscripten::function("classifyStep", &LLKA::classifyStep);
    emscripten::function("classifySteps", &LLKA::classifySteps);
    emscripten::function("initializeClassificationContext", &LLKA::initializeClassificationContext);
}

#endif // LLKA_GENERATE_EMSCRIPTEN_BINDINGS

#endif // _LLKA_CPP_H

#endif // __cplusplus
