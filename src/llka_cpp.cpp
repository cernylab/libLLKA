/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifdef __cplusplus

#include <llka_config.h>
#ifdef LLKA_PLATFORM_EMSCRIPTEN
#define LLKA_GENERATE_EMSCRIPTEN_BINDINGS
#endif // LLKA_PLATFORM_EMSCRIPTEN

#include <llka_cpp.h>
#include <llka_nucleotide.h>
#include <llka_resource_loaders.h>
#include <llka_superposition.h>
#include <llka_tracing.h>

#include <segmentation.h>
#include <util/geometry.h>

#include <algorithm>
#include <cstring>
#include <memory>
#include <ostream>
#include <span>

#ifdef LLKA_PLATFORM_EMSCRIPTEN
    #define _EMX_GET_DEF(type, cls, var) auto cls::_emsGet_##var() const -> const type & { return var; }
    #define _EMX_SET_DEF(type, cls, var) auto cls::_emsSet_##var(const type &var) -> void { this->var = var; }
    #define _EMX_GET_SET_DEF(type, cls, var) _EMX_GET_DEF(type, cls, var) _EMX_SET_DEF(type, cls, var)
    #define _EMX_GET_SET_DEF_NOCONST(type, cls, var) _EMX_GET_DEF(type, cls, var) _EMX_SET_DEF_NOCONST(type, cls, var)
#else
    #define _EMX_GET_DEF(type, cls, var)
    #define _EMX_SET_DEF(type, cls, var)
    #define _EMX_GET_SET_DEF(type, cls, var)
#endif // LLKA_PLATFORM_EMSCRIPTEN


namespace LLKA {

namespace helpers {
    class WrappedCStructure {
    public:
        LLKA_Structure cStru;
        std::vector<LLKA_Atom> cAtoms;
    };

    class WrappedCStructures {
    public:
        explicit WrappedCStructures(std::vector<LLKA_Structure> cStrus) :
            cStrus{std::move(cStrus)}
        {
        }

        ~WrappedCStructures() {
            // We do not own the data of the atoms but we still need to delete
            // the arrays that we copied the source atoms into
            for (auto &cS : cStrus)
                delete [] cS.atoms;
        }

        auto get() -> LLKA_Structures {
            return { cStrus.data(), cStrus.size() };
        }
        auto get() const -> const LLKA_Structures {
            return { const_cast<LLKA_Structure *>(cStrus.data()), cStrus.size() };
        }

        std::vector<LLKA_Structure> cStrus;
    };

    inline auto atomToCAtom(const Atom &atom, LLKA_Atom &cAtom) -> void;
    inline auto stepToCStep(const ClassifiedStep &step) -> LLKA_ClassifiedStep;
    inline auto struToCStru(const Structure &stru, LLKA_Structure &cStru, std::vector<LLKA_Atom> &cAtoms) -> void;
    inline auto struToWrappedCStru(const Structure &stru) -> WrappedCStructure;
    inline auto struViewToWrappedCStru(const StructureView &view) -> WrappedCStructure;
    inline auto strusToWrappedCStrus(const Structures &strus) -> WrappedCStructures;
    inline auto cStruToStru(const LLKA_Structure &cStru) -> Structure;
    inline auto cStruToStru(const LLKA_Structure &cStru, Structure &stru) -> void;
    inline auto cStruViewToStruView(const LLKA_StructureView &cView, const WrappedCStructure &wStru, StructureView &view, const Structure &stru) -> void;
    inline auto cCifDataToCifData(const LLKA_CifData *cCifData) -> CifData;
} // namespace helpers

//
// Main
//

using Points = std::vector<LLKA_Point>; // TODO: Revise this when we figure out how to do memory alignment of point arrays
auto operator<<(std::ostream &os, const LLKA_Point &pt) -> std::ostream &
{
    os << "[ " << pt.x << "; " << pt.y << "; " << pt.z << " ]";
    return os;
}

auto errorToString(LLKA_RetCode tRet) -> std::string
{
    std::string s = LLKA_errorToString(tRet);
    return s;
}

//
// Structure
//

auto makeEmptyAtom() -> Atom
{
    Atom atom;

    atom.pdbx_PDB_ins_code = LLKA_NO_INSCODE;
    atom.label_atom_id = LLKA_NO_ALTID;

    return atom;
}

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
) -> Atom
{
    Atom atom;

    atom.type_symbol = std::move(type_symbol);
    atom.label_atom_id = std::move(label_atom_id);
    atom.label_entity_id = std::move(label_entity_id);
    atom.label_comp_id = std::move(label_comp_id);
    atom.label_asym_id = std::move(label_asym_id);
    atom.pdbx_PDB_ins_code = std::move(pdbx_PDB_ins_code);
    atom.coords = coords;
    atom.id = id;
    atom.label_seq_id = label_seq_id;
    atom.auth_seq_id = auth_seq_id;
    atom.pdbx_PDB_model_num = pdbx_PDB_model_num;
    atom.label_alt_id = label_alt_id;

    atom.auth_atom_id = auth_atom_id.empty() ? atom.label_atom_id : std::move(auth_atom_id);
    atom.auth_comp_id = auth_comp_id.empty() ? atom.label_comp_id : std::move(auth_comp_id);
    atom.auth_asym_id = auth_asym_id.empty() ? atom.label_asym_id : std::move(auth_asym_id);

    return atom;
}

auto atomsEqual(const Atom &a, const Atom &b) noexcept -> bool
{
    LLKA_Atom cA{
        a.type_symbol.c_str(),
        a.label_atom_id.c_str(),
        a.label_entity_id.c_str(),
        a.label_comp_id.c_str(),
        a.label_asym_id.c_str(),
        a.auth_atom_id.c_str(),
        a.auth_comp_id.c_str(),
        a.auth_asym_id.c_str(),
        { a.coords.x, a.coords.y, a.coords.z },
        a.id,
        a.label_seq_id,
        a.auth_seq_id,
        a.pdbx_PDB_model_num,
        a.pdbx_PDB_ins_code.c_str(),
        a.label_alt_id
    };

    LLKA_Atom cB{
        b.type_symbol.c_str(),
        b.label_atom_id.c_str(),
        b.label_entity_id.c_str(),
        b.label_comp_id.c_str(),
        b.label_asym_id.c_str(),
        b.auth_atom_id.c_str(),
        b.auth_comp_id.c_str(),
        b.auth_asym_id.c_str(),
        { b.coords.x, b.coords.y, b.coords.z },
        b.id,
        b.label_seq_id,
        b.auth_seq_id,
        b.pdbx_PDB_model_num,
        b.pdbx_PDB_ins_code.c_str(),
        b.label_alt_id
    };

    return LLKA_compareAtoms(&cA, &cB, LLKA_FALSE) == LLKA_TRUE;
}

auto operator<<(std::ostream &os, const Atom &atom) -> std::ostream &
{
    os
        << atom.pdbx_PDB_model_num << " "
        << atom.label_asym_id << " "
        << atom.label_comp_id << " "
        << atom.label_seq_id << "(" << atom.auth_seq_id << ") "
        << (!atom.pdbx_PDB_ins_code.empty() ? std::string{"<"} + atom.pdbx_PDB_ins_code + std::string{"> "} : "")
        << (atom.label_alt_id ? std::string{"["} + atom.label_alt_id + std::string{"] "} : "")
        << atom.label_atom_id
        << atom.coords;
    return os;
}

auto operator<<(std::ostream &os, const Structure &stru) -> std::ostream &
{
    for (const auto &atom : stru)
        os << atom << "\n";
    return os;
}

auto compareAtoms(const Atom &a, const Atom &b, bool ignoreId) noexcept -> bool
{
    LLKA_Atom ca{
        a.type_symbol.c_str(),
        a.label_atom_id.c_str(),
        a.label_entity_id.c_str(),
        a.label_comp_id.c_str(),
        a.label_asym_id.c_str(),
        a.auth_atom_id.c_str(),
        a.auth_comp_id.c_str(),
        a.auth_asym_id.c_str(),
        a.coords,
        a.id,
        a.label_seq_id,
        a.auth_seq_id,
        a.pdbx_PDB_model_num,
        a.pdbx_PDB_ins_code.c_str(),
        a.label_alt_id
    };

    LLKA_Atom cb{
        b.type_symbol.c_str(),
        b.label_atom_id.c_str(),
        b.label_entity_id.c_str(),
        b.label_comp_id.c_str(),
        b.label_asym_id.c_str(),
        b.auth_atom_id.c_str(),
        b.auth_comp_id.c_str(),
        b.auth_asym_id.c_str(),
        b.coords,
        b.id,
        b.label_seq_id,
        b.auth_seq_id,
        b.pdbx_PDB_model_num,
        b.pdbx_PDB_ins_code.c_str(),
        b.label_alt_id
    };

    return LLKA_compareAtoms(&ca, &cb, ignoreId);
}

auto makeStructure(const LLKA_Atom *atoms, size_t nAtoms) -> Structure
{
    Structure stru;

    for (size_t idx = 0; idx < nAtoms; idx++) {
        const auto &cAtom = atoms[idx];
        stru.push_back(
            makeAtom(
                cAtom.type_symbol,
                cAtom.label_atom_id,
                cAtom.label_entity_id,
                cAtom.label_comp_id,
                cAtom.label_asym_id,
                cAtom.auth_atom_id,
                cAtom.auth_comp_id,
                cAtom.auth_asym_id,
                cAtom.coords,
                cAtom.id,
                cAtom.label_seq_id,
                cAtom.auth_seq_id,
                cAtom.pdbx_PDB_model_num,
                cAtom.pdbx_PDB_ins_code,
                cAtom.label_alt_id
            )
        );
    }

    return stru;
}

auto splitByAltIds(const Structure &stru) -> std::vector<AltIdSplit>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_AlternatePositions alts;
    LLKA_Structures strus = LLKA_splitByAltIds(&wStru.cStru, &alts);

    std::vector<AltIdSplit> splits{};

    if (alts.nPositions == 0) {
        splits.resize(1);

        splits[0].structure = helpers::cStruToStru(strus.strus[0]);
        splits[0].altId = LLKA_NO_ALTID;
    } else {
        splits.resize(strus.nStrus);

        for (size_t idx = 0; idx < strus.nStrus; idx++) {
            splits[idx].structure = helpers::cStruToStru(strus.strus[idx]);
            splits[idx].altId = alts.positions[idx];
        }
    }

    LLKA_destroyStructures(&strus);
    LLKA_destroyAlternatePositions(&alts);

    return splits;
}

auto splitStructureToDinucleotideSteps(const Structure &stru) -> RCResult<Structures>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_Structures steps;
    auto tRet = LLKA_splitStructureToDinucleotideSteps(&wStru.cStru, &steps);
    if (tRet != LLKA_OK)
        return RCResult<Structures>::fail(tRet);

    Structures strus(steps.nStrus);
    for (size_t idx = 0; idx < steps.nStrus; idx++)
        strus[idx] = helpers::cStruToStru(steps.strus[idx]);

    LLKA_destroyStructures(&steps);

    return RCResult<Structures>::succeed(std::move(strus));
}

//
// Measurements
//

template <typename T> requires std::is_floating_point_v<T>
LLKA_CPP_API
auto measureAngle(const Atom &a, const Atom &b, const Atom &c) -> T
{
    return LLKAInternal::angle<T>(a.coords, b.coords, c.coords);
}
template LLKA_CPP_API auto measureAngle<float>(const Atom &a, const Atom &b, const Atom &c) -> float;
template LLKA_CPP_API auto measureAngle<double>(const Atom &a, const Atom &b, const Atom &c) -> double;
template LLKA_CPP_API auto measureAngle<long double>(const Atom &a, const Atom &b, const Atom &c) -> long double;

template <typename T> requires std::is_floating_point_v<T>
LLKA_CPP_API
auto measureDihedral(const Atom &a, const Atom &b, const Atom &c, const Atom &d ) -> T
{
    return LLKAInternal::dihedralAngle<T>(a.coords, b.coords, c.coords, d.coords);
}
template LLKA_CPP_API auto measureDihedral<float>(const Atom &a, const Atom &b, const Atom &c, const Atom &d) -> float;
template LLKA_CPP_API auto measureDihedral<double>(const Atom &a, const Atom &b, const Atom &c, const Atom &d) -> double;
template LLKA_CPP_API auto measureDihedral<long double>(const Atom &a, const Atom &b, const Atom &c, const Atom &d) -> long double;

template <typename T> requires std::is_floating_point_v<T>
LLKA_CPP_API
auto measureDistance(const Atom &a, const Atom &b) -> T
{
    return LLKAInternal::spatialDistance<T>(a.coords, b.coords);
}
template LLKA_CPP_API auto measureDistance<float>(const Atom &a, const Atom &b) -> float;
template LLKA_CPP_API auto measureDistance<double>(const Atom &a, const Atom &b) -> double;
template LLKA_CPP_API auto measureDistance<long double>(const Atom &a, const Atom &b) -> long double;

//
// Superposition
//

auto centroid(const Points &points) noexcept -> LLKA_Point
{
    // TODO: We are mapping the points vector memory directly onto
    // the LLKA_Points structure. This will work okay until we start
    // to align memory. Once we start doing that, this will have to be
    // revised.

    const LLKA_Points pts {
        { const_cast<LLKA_Point *>(points.data()) },   // LLKA_Points expects mutable array
        points.size()
    };

    return LLKA_centroidPoints(&pts);
}

auto centroid(const Structure &stru) noexcept -> LLKA_Point
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    return LLKA_centroidStructure(&wStru.cStru);
}

auto rmsd(const Points &a, const Points &b) noexcept -> RCResult<double>
{
    const LLKA_Points cA{
        { const_cast<LLKA_Point *>(a.data()) },   // LLKA_Points expects mutable array
        a.size()
    };
    const LLKA_Points cB{
        { const_cast<LLKA_Point *>(b.data()) },
        b.size()
    };

    double rmsd;
    auto tRet = LLKA_rmsdPoints(&cA, &cB, &rmsd);
    if (tRet != LLKA_OK)
        return RCResult<double>::fail(tRet);
    return RCResult<double>::succeed(rmsd);
}

auto rmsd(const Structure &a, const Structure &b) noexcept -> RCResult<double>
{
    const auto cStruA = helpers::struToWrappedCStru(a);
    const auto cStruB = helpers::struToWrappedCStru(b);

    double rmsd;
    auto tRet = LLKA_rmsdStructures(&cStruA.cStru, &cStruB.cStru, &rmsd);
    if (tRet != LLKA_OK)
        return RCResult<double>::fail(tRet);
    return RCResult<double>::succeed(rmsd);
}

auto superpose(Points &what, const Points &onto) noexcept -> RCResult<double>
{
    LLKA_Points cWhat{
        { what.data() },
        what.size()
    };
    const LLKA_Points cOnto{
        { const_cast<LLKA_Point *>(onto.data()) },   // LLKA_Points expects mutable array
        onto.size()
    };

    double rmsd;
    auto tRet = LLKA_superposePoints(&cWhat, &cOnto, &rmsd);
    if (tRet != LLKA_OK)
        return RCResult<double>::fail(tRet);
    return RCResult<double>::succeed(rmsd);
}

auto superpose(Structure &what, const Structure &onto) noexcept -> RCResult<double>
{
    auto wWhat = helpers::struToWrappedCStru(what);
    const auto wOnto = helpers::struToWrappedCStru(onto);

    double rmsd;
    auto tRet = LLKA_superposeStructures(&wWhat.cStru, &wOnto.cStru, &rmsd);
    if (tRet != LLKA_OK)
        return RCResult<double>::fail(tRet);

    // Map the results back
    for (size_t idx = 0; idx < what.size(); idx++)
        what[idx].coords = wWhat.cAtoms[idx].coords;

    return RCResult<double>::succeed(rmsd);
}

#ifdef LLKA_PLATFORM_EMSCRIPTEN
StructureSegments::StructureSegments() noexcept :
    structure{nullptr}
{}
#endif // LLKA_PLATFORM_EMSCRIPTEN

//
// Segmentation
//

StructureSegments::StructureSegments(Structure &structure) noexcept :
    structure{&structure}
{
    auto counts = LLKAInternal::countAtomsForSegmentation(structure.data(), structure.size());

    // Prepare the mapping
    for (auto cntModelIt = counts.cbegin(); cntModelIt != counts.cend(); cntModelIt++) {
        models.emplace(cntModelIt->first, Model{});
        const auto &cntModel = cntModelIt->second;
        auto &model = models.at(cntModelIt->first);
#ifndef LLKA_PLATFORM_EMSCRIPTEN
        model.atoms.reserve(cntModel.nAtoms);
#endif // LLKA_PLATFORM_EMSCRIPTEN

        for (auto cntChainIt = cntModel.chains.cbegin(); cntChainIt != cntModel.chains.cend(); cntChainIt++) {
            model.chains.emplace(cntChainIt->first, Chain{});
            const auto &cntChain = cntChainIt->second;
            auto &chain = model.chains.at(cntChainIt->first);
#ifndef LLKA_PLATFORM_EMSCRIPTEN
            chain.atoms.reserve(cntChain.nAtoms);
#endif // LLKA_PLATFORM_EMSCRIPTEN

            for (auto cntResidueIt = cntChain.residues.cbegin(); cntResidueIt != cntChain.residues.cend(); cntResidueIt++) {
                chain.residues.emplace(cntResidueIt->first, Residue{});
                const auto &cntResidue = cntResidueIt->second; // cntResidue is just the number of atoms in the residue
                chain.residues.at(cntResidueIt->first).atoms.reserve(cntResidue);
            }
        }
    }

    // Do the mapping
    for (auto &atom : structure) {
#ifdef LLKA_PLATFORM_EMSCRIPTEN
        auto atomPtr = atom;
#else
        auto *atomPtr = &atom;
#endif // LLKA_PLATFORM_EMSCRIPTEN

        auto &model = models.at(atom.pdbx_PDB_model_num);
        auto &chain = model.chains.at(atom.label_asym_id);
        auto &residue = chain.residues.at(atom.label_seq_id);

#ifndef LLKA_PLATFORM_EMSCRIPTEN
        model.atoms.push_back(atomPtr);
        chain.atoms.push_back(atomPtr);
#endif // LLKA_PLATFORM_EMSCRIPTEN
        residue.atoms.push_back(atomPtr);
    }
}

_EMX_GET_SET_DEF(StructureSegments::Models, StructureSegments, models)
#ifdef LLKA_PLATFORM_EMSCRIPTEN
auto StructureSegments::_emsGet_structure() const -> const Structure &
{
    return *this->structure;
}
#endif // LLKA_PLATFORM_EMSCRIPTEN

//
// NtC
//

_EMX_GET_SET_DEF(std::string, AtomNameQuad, a)
_EMX_GET_SET_DEF(std::string, AtomNameQuad, b)
_EMX_GET_SET_DEF(std::string, AtomNameQuad, c)
_EMX_GET_SET_DEF(std::string, AtomNameQuad, d)

auto backboneAtomIndex(const BackboneAtom &bkbnAtom, const Structure &backbone) noexcept -> size_t
{
    const LLKA_BackboneAtom cBkbnAtom{
        bkbnAtom.residue,
        bkbnAtom.name.c_str()
    };
    const auto wBackbone = helpers::struToWrappedCStru(backbone);

    return LLKA_backboneAtomIndex(&cBkbnAtom, &wBackbone.cStru);
}

auto calculateStepMetrics(const Structure &stru) noexcept -> RCResult<LLKA_StepMetrics>
{
    auto wStru = helpers::struToWrappedCStru(stru);
    LLKA_StepMetrics metrics;

    auto tRet = LLKA_calculateStepMetrics(&wStru.cStru, &metrics);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_StepMetrics>::fail(tRet);
    return RCResult<LLKA_StepMetrics>::succeed(std::move(metrics));
}

auto calculateStepMetricsDifferenceAgainstReference(const Structure &stru, LLKA_NtC ntc) noexcept -> RCResult<LLKA_StepMetrics>
{
    auto wStru = helpers::struToWrappedCStru(stru);
    LLKA_StepMetrics metrics;

    auto tRet = LLKA_calculateStepMetricsDifferenceAgainstReference(&wStru.cStru, ntc, &metrics);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_StepMetrics>::fail(tRet);
    return RCResult<LLKA_StepMetrics>::succeed(std::move(metrics));
}

auto crossResidueMetric(LLKA_CrossResidueMetric metric, const Structure &stru) noexcept -> RCResult<Structure>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_Structure cMetricStru;
    auto tRet = LLKA_crossResidueMetric(metric, &wStru.cStru, &cMetricStru);
    if (tRet != LLKA_OK)
        return RCResult<Structure>::fail(tRet);

    Structure metricStru;
    helpers::cStruToStru(cMetricStru, metricStru);
    LLKA_destroyStructure(&cMetricStru);

    return RCResult<Structure>::succeed(std::move(metricStru));
}

auto crossResidueMetricAtoms(const std::string &firstBase, const std::string &secondBase, LLKA_CrossResidueMetric metric) -> RCResult<AtomNameQuad>
{
    LLKA_AtomNameQuad cQuad;
    auto tRet = LLKA_crossResidueMetricAtomsFromBases(firstBase.c_str(), secondBase.c_str(), metric, &cQuad);
    if (tRet != LLKA_OK)
        return RCResult<AtomNameQuad>::fail(tRet);

    AtomNameQuad quad{};
    quad.a = cQuad.a;
    quad.b = cQuad.b;
    if (metric == LLKA_XR_TOR_MU) {
        quad.c = cQuad.c;
        quad.d = cQuad.d;
    }

    return RCResult<AtomNameQuad>::succeed(std::move(quad));
}

auto crossResidueMetricAtoms(const Structure &stru, LLKA_CrossResidueMetric metric) -> RCResult<AtomNameQuad>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_AtomNameQuad cQuad;
    auto tRet = LLKA_crossResidueMetricAtomsFromStructure(&wStru.cStru, metric, &cQuad);
    if (tRet != LLKA_OK)
        return RCResult<AtomNameQuad>::fail(tRet);

    AtomNameQuad quad{};
    quad.a = cQuad.a;
    quad.b = cQuad.b;
    if (metric == LLKA_XR_TOR_MU) {
        quad.c = cQuad.c;
        quad.d = cQuad.d;
    }

    return RCResult<AtomNameQuad>::succeed(std::move(quad));
}

auto crossResidueMetricName(LLKA_CrossResidueMetric metric, bool greek) noexcept -> std::string
{
    return LLKA_crossResidueMetricName(metric, greek);
}

auto dinucleotideTorsion(LLKA_DinucleotideTorsion torsion, const Structure &stru) -> RCResult<Structure>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_Structure cTorsionStru;
    auto tRet = LLKA_dinucleotideTorsion(torsion, &wStru.cStru, &cTorsionStru);
    if (tRet != LLKA_OK)
        return RCResult<Structure>::fail(tRet);

    Structure torsionStru;
    helpers::cStruToStru(cTorsionStru, torsionStru);
    LLKA_destroyStructure(&cTorsionStru);

    return RCResult<Structure>::succeed(std::move(torsionStru));
}

auto dinucleotideTorsionAtoms(const std::string &firstBase, const std::string &secondBase, LLKA_DinucleotideTorsion torsion) -> RCResult<AtomNameQuad>
{
    LLKA_AtomNameQuad cQuad;
    auto tRet = LLKA_dinucleotideTorsionAtomsFromBases(firstBase.c_str(), secondBase.c_str(), torsion, &cQuad);
    if (tRet != LLKA_OK)
        return RCResult<AtomNameQuad>::fail(tRet);

    AtomNameQuad quad{};
    quad.a = cQuad.a;
    quad.b = cQuad.b;
    quad.c = cQuad.c;
    quad.d = cQuad.d;

    return RCResult<AtomNameQuad>::succeed(std::move(quad));
}

auto dinucleotideTorsionAtoms(const Structure &stru, LLKA_DinucleotideTorsion torsion) -> RCResult<AtomNameQuad>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_AtomNameQuad cQuad;
    auto tRet = LLKA_dinucleotideTorsionAtomsFromStructure(&wStru.cStru, torsion, &cQuad);
    if (tRet != LLKA_OK)
        return RCResult<AtomNameQuad>::fail(tRet);

    AtomNameQuad quad{};
    quad.a = cQuad.a;
    quad.b = cQuad.b;
    quad.c = cQuad.c;
    quad.d = cQuad.d;

    return RCResult<AtomNameQuad>::succeed(std::move(quad));
}

auto dinucleotideTorsionName(LLKA_DinucleotideTorsion torsion, bool greek) noexcept -> std::string
{
    return LLKA_dinucleotideTorsionName(torsion, greek);
}

auto extractBackbone(const Structure &stru) -> RCResult<Structure>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_Structure cBackbone;
    auto tRet = LLKA_extractBackbone(&wStru.cStru, &cBackbone);
    if (tRet != LLKA_OK)
        return RCResult<Structure>::fail(tRet);

    Structure backbone;
    helpers::cStruToStru(cBackbone, backbone);
    LLKA_destroyStructure(&cBackbone);

    return RCResult<Structure>::succeed(std::move(backbone));
}

auto extractExtendedBackbone(const Structure &stru) noexcept -> RCResult<Structure>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_Structure cBackbone;
    auto tRet = LLKA_extractExtendedBackbone(&wStru.cStru, &cBackbone);
    if (tRet != LLKA_OK)
        return RCResult<Structure>::fail(tRet);

    Structure backbone;
    helpers::cStruToStru(cBackbone, backbone);
    LLKA_destroyStructure(&cBackbone);

    return RCResult<Structure>::succeed(std::move(backbone));
}

auto extractMetricsStructure(const Structure &stru) noexcept -> RCResult<Structure>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_Structure cBackbone;
    auto tRet = LLKA_extractMetricsStructure(&wStru.cStru, &cBackbone);
    if (tRet != LLKA_OK)
        return RCResult<Structure>::fail(tRet);

    Structure backbone;
    helpers::cStruToStru(cBackbone, backbone);
    LLKA_destroyStructure(&cBackbone);

    return RCResult<Structure>::succeed(std::move(backbone));
}

auto nameToCANA(const std::string &name) noexcept -> LLKA_CANA
{
    return LLKA_nameToCANA(name.c_str());
}

auto CANAToName(LLKA_CANA cana) noexcept -> std::string
{
    return LLKA_CANAToName(cana);
}

auto nameToNtC(const std::string &name) noexcept -> LLKA_NtC
{
    return LLKA_nameToNtC(name.c_str());
}

auto NtCToName(LLKA_NtC ntc) noexcept -> std::string
{
    return LLKA_NtCToName(ntc);
}

auto structureIsStep(const Structure &stru) noexcept -> RCResult<LLKA_StepInfo>
{
    const auto wStru = helpers::struToWrappedCStru(stru);
    StepInfo info;

    auto tRet  = LLKA_structureIsStep(&wStru.cStru, &info);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_StepInfo>::fail(tRet);
    return RCResult<LLKA_StepInfo>::succeed(std::move(info));
}

auto operator<<(std::ostream &os, const StepMetrics &metrics) -> std::ostream &
{
    os
        << LLKA::dinucleotideTorsionName(LLKA_TOR_DELTA_1, LLKA_FALSE) << " " << metrics.delta_1 << "\n"
        << LLKA::dinucleotideTorsionName(LLKA_TOR_EPSILON_1, LLKA_FALSE) << " " << metrics.epsilon_1 << "\n"
        << LLKA::dinucleotideTorsionName(LLKA_TOR_ZETA_1, LLKA_FALSE) << " " << metrics.zeta_1 << "\n"
        << LLKA::dinucleotideTorsionName(LLKA_TOR_ALPHA_2, LLKA_FALSE) << " " << metrics.alpha_2 << "\n"
        << LLKA::dinucleotideTorsionName(LLKA_TOR_BETA_2, LLKA_FALSE) << " " << metrics.beta_2 << "\n"
        << LLKA::dinucleotideTorsionName(LLKA_TOR_GAMMA_2, LLKA_FALSE) << " " << metrics.gamma_2 << "\n"
        << LLKA::dinucleotideTorsionName(LLKA_TOR_DELTA_2, LLKA_FALSE) << " " << metrics.delta_2 << "\n"
        << LLKA::dinucleotideTorsionName(LLKA_TOR_CHI_1, LLKA_FALSE) << " " << metrics.chi_1 << "\n"
        << LLKA::dinucleotideTorsionName(LLKA_TOR_CHI_2, LLKA_FALSE) << " " << metrics.chi_2 << "\n"
        << LLKA::crossResidueMetricName(LLKA_XR_DIST_CC, LLKA_FALSE) << " " << metrics.CC << "\n"
        << LLKA::crossResidueMetricName(LLKA_XR_DIST_NN, LLKA_FALSE) << " " << metrics.NN << "\n"
        << LLKA::crossResidueMetricName(LLKA_XR_TOR_MU, LLKA_FALSE) << " " << metrics.mu << "\n";

    return os;
}

//
// Connectivity, similarity
//

auto measureStepConnectivity(const Structure &positionFirst, LLKA_NtC ntcFirst, const Structure &positionSecond, LLKA_NtC ntcSecond) noexcept -> RCResult<LLKA_Connectivity>
{
    const auto wPositionFirst = helpers::struToWrappedCStru(positionFirst);
    const auto wPositionSecond = helpers::struToWrappedCStru(positionSecond);

    LLKA_Connectivity conn;
    auto tRet = LLKA_measureStepConnectivityNtCs(&wPositionFirst.cStru, ntcFirst, &wPositionSecond.cStru, ntcSecond, &conn);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_Connectivity>::fail(tRet);
    return RCResult<LLKA_Connectivity>::succeed(std::move(conn));
}

auto measureStepConnectivity(const Structure &positionFirst, std::vector<LLKA_NtC> ntcsFirst, const Structure &positionSecond, LLKA_NtC ntcSecond) noexcept -> RCResult<Connectivities>
{
    const auto wPositionFirst = helpers::struToWrappedCStru(positionFirst);
    const auto wPositionSecond = helpers::struToWrappedCStru(positionSecond);

    Connectivities conns;
    conns.resize(ntcsFirst.size());

    LLKA_Connectivities cConns {
        conns.data(),
        conns.size()
    };

    ntcsFirst.push_back(LLKA_INVALID_NTC);
    auto tRet = LLKA_measureStepConnectivityNtCsMultipleFirst(&wPositionFirst.cStru, ntcsFirst.data(), &wPositionSecond.cStru, ntcSecond, &cConns);
    if (tRet != LLKA_OK)
        return RCResult<Connectivities>::fail(tRet);
    return RCResult<Connectivities>::succeed(std::move(conns));
}

auto measureStepConnectivity(const Structure &positionFirst, LLKA_NtC ntcFirst, const Structure &positionSecond, std::vector<LLKA_NtC> ntcsSecond) noexcept -> RCResult<Connectivities>
{
    const auto wPositionFirst = helpers::struToWrappedCStru(positionFirst);
    const auto wPositionSecond = helpers::struToWrappedCStru(positionSecond);

    Connectivities conns;
    conns.resize(ntcsSecond.size());

    LLKA_Connectivities cConns {
        conns.data(),
        conns.size()
    };

    ntcsSecond.push_back(LLKA_INVALID_NTC);
    auto tRet = LLKA_measureStepConnectivityNtCsMultipleSecond(&wPositionFirst.cStru, ntcFirst, &wPositionSecond.cStru, ntcsSecond.data(), &cConns);
    if (tRet != LLKA_OK)
        return RCResult<Connectivities>::fail(tRet);
    return RCResult<Connectivities>::succeed(std::move(conns));
}

auto measureStepConnectivity(const Structure &positionFirst, const Structure &dinuFirst, const Structure &positionSecond, const Structure &dinuSecond) noexcept -> RCResult<Connectivity>
{
    const auto wPositionFirst = helpers::struToWrappedCStru(positionFirst);
    const auto wDinuFirst = helpers::struToWrappedCStru(dinuFirst);
    const auto wPositionSecond = helpers::struToWrappedCStru(positionSecond);
    const auto wDinuSecond = helpers::struToWrappedCStru(dinuSecond);

    LLKA_Connectivity conn;
    auto tRet = LLKA_measureStepConnectivityStructures(&wPositionFirst.cStru, &wDinuFirst.cStru, &wPositionSecond.cStru, &wDinuSecond.cStru, &conn);

    if (tRet != LLKA_OK)
        return RCResult<Connectivity>::fail(tRet);
    return RCResult<Connectivity>::succeed(std::move(conn));
}

auto measureStepConnectivity(const Structure &positionFirst, const Structure &dinuFirst, const Structure &positionSecond, const Structures &dinusSecond) noexcept -> RCResult<Connectivities>
{
    const auto wPositionFirst = helpers::struToWrappedCStru(positionFirst);
    const auto wDinuFirst = helpers::struToWrappedCStru(dinuFirst);
    const auto wPositionSecond = helpers::struToWrappedCStru(positionSecond);
    const auto wDinusSecond = helpers::strusToWrappedCStrus(dinusSecond);
    auto ds = wDinusSecond.get();

    Connectivities conns;
    conns.resize(wDinusSecond.cStrus.size());

    LLKA_Connectivities cConns{
        conns.data(),
        conns.size()
    };

    auto tRet = LLKA_measureStepConnectivityStructuresMultiple(&wPositionFirst.cStru, &wDinuFirst.cStru, &wPositionSecond.cStru, &ds, &cConns);

    if (tRet != LLKA_OK)
        return RCResult<Connectivities>::fail(tRet);
    return RCResult<Connectivities>::succeed(std::move(conns));
}

auto measureStepSimilarity(const Structure &stepStru, LLKA_NtC ntc) noexcept -> RCResult<Similarity>
{
    const auto wStepStru = helpers::struToWrappedCStru(stepStru);
    LLKA_Similarity simil;

    auto tRet = LLKA_measureStepSimilarityNtC(&wStepStru.cStru, ntc, &simil);
    if (tRet != LLKA_OK)
        return RCResult<Similarity>::fail(tRet);
    return RCResult<Similarity>::succeed(std::move(simil));
}

auto measureStepSimilarity(const Structure &stepStru, std::vector<LLKA_NtC> ntcs) noexcept -> RCResult<Similarities>
{
    const auto wStepStru = helpers::struToWrappedCStru(stepStru);
    Similarities simils;
    simils.resize(ntcs.size());

    LLKA_Similarities cResults{
        simils.data(),
        simils.size()
    };

    ntcs.push_back(LLKA_INVALID_NTC);
    auto tRet = LLKA_measureStepSimilarityNtCMultiple(&wStepStru.cStru, ntcs.data(), &cResults);
    if (tRet != LLKA_OK)
        return RCResult<Similarities>::fail(tRet);
    return RCResult<Similarities>::succeed(std::move(simils));
}

auto measureStepSimilarity(const Structure &stepStru, const Structure &refStru) noexcept -> RCResult<LLKA_Similarity>
{
    const auto wStepStru = helpers::struToWrappedCStru(stepStru);
    const auto wRefStru = helpers::struToWrappedCStru(refStru);

    LLKA_Similarity simil;

    auto tRet = LLKA_measureStepSimilarityStructure(&wStepStru.cStru, &wRefStru.cStru, &simil);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_Similarity>::fail(tRet);
    return RCResult<LLKA_Similarity>::succeed(std::move(simil));
}

//
// MiniCif
//

CifData::Category::Category(std::string name, Items items) :
    name{std::move(name)},
    items{std::move(items)}
{
}

auto CifData::Category::isLoop() -> bool
{
    if (items.empty())
        return false;
    return items.front().values.size() > 1;
}

_EMX_GET_SET_DEF(std::string, CifData::Category, name)
_EMX_GET_SET_DEF(Items, CifData::Category, items)

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED

auto cifToStructure(const std::filesystem::path &path, int32_t options) -> CifResult
{
    LLKA_ImportedStructure cImportedStru;
    char *error;

    auto tRet = LLKA_cifFileToStructure(path.c_str(), &cImportedStru, &error, options);
    if (tRet != LLKA_OK) {
        CifError err;
        err.tRet = tRet;

        if (error) {
            err.error = std::string{error};
            LLKA_destroyString(error);
        }

        return CifResult::fail(std::move(err));
    } else {
        ImportedStructure importedStru;

        importedStru.id = cImportedStru.entry.id;
        importedStru.structure = helpers::cStruToStru(cImportedStru.structure);
        if (cImportedStru.cifData != nullptr)
            importedStru.cifData = helpers::cCifDataToCifData(cImportedStru.cifData);

        LLKA_destroyImportedStructure(&cImportedStru);

        return CifResult::succeed(std::move(importedStru));
    }
}

#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

auto cifToStructure(const std::string &text, int32_t options) -> CifResult
{
    LLKA_ImportedStructure cImportedStru;
    char *error;

    auto tRet = LLKA_cifTextToStructure(text.c_str(), &cImportedStru, &error, options);
    if (tRet != LLKA_OK) {
        CifError err;
        err.tRet = tRet;

        if (error) {
            err.error = std::string{error};
            LLKA_destroyString(error);
        }

        return CifResult::fail(std::move(err));
    } else {
        ImportedStructure importedStru;

        importedStru.id = cImportedStru.entry.id;
        importedStru.structure = helpers::cStruToStru(cImportedStru.structure);
        if (cImportedStru.cifData != nullptr)
            importedStru.cifData = helpers::cCifDataToCifData(cImportedStru.cifData);

        LLKA_destroyImportedStructure(&cImportedStru);

        return CifResult::succeed(std::move(importedStru));
    }
}

auto cifDataToString(const CifData &cifData, bool pretty) -> RCResult<std::string>
{
    LLKA_RetCode tRet = LLKA_OK;

    auto cCifData = LLKA_cifData_empty();

    // Map the cifData to the corresponding C representation
    for (size_t blockIdx = 0; blockIdx < cifData.blocks.size(); blockIdx++) {
        const auto &block = cifData.blocks[blockIdx];

        LLKA_cifData_addBlock(cCifData, block.name.c_str());
        auto &cBlock = cCifData->blocks[blockIdx];

        for (const auto &cat : block.categories) {
            auto cCat = LLKA_cifDataBlock_addCategory(&cBlock, cat.name.c_str());
            if (cCat == nullptr) {
                tRet = LLKA_E_BAD_DATA;
                goto mapping_done;
            }

            for (const auto &item : cat.items) {
                auto cItem = LLKA_cifDataCategory_addItem(cCat, item.keyword.c_str());
                if (cItem == nullptr) {
                    tRet = LLKA_E_BAD_DATA;
                    goto mapping_done;
                }

                LLKA_CifDataValue *cValues = new LLKA_CifDataValue[item.values.size()];
                for (size_t valIdx = 0; valIdx < item.values.size(); valIdx++) {
                    const auto &value = item.values[valIdx];
                    auto cValue = &cValues[valIdx];

                    cValue->state = value.state;
                    if (value.state == LLKA_MINICIF_VALUE_SET)
                        cValue->text = value.text.c_str();
                    else
                        cValue->text = nullptr;
                }

                LLKA_cifDataItem_setValues(cItem, cValues, item.values.size());
                delete [] cValues;
            }
        }
    }

mapping_done:
    if (tRet != LLKA_OK) {
        auto res = RCResult<std::string>::fail(tRet);
        LLKA_destroyCifData(cCifData);

        return res;
    }

    tRet = LLKA_cifData_detaint(cCifData);
    if (tRet != LLKA_OK) {
        auto res = RCResult<std::string>::fail(tRet);
        LLKA_destroyCifData(cCifData);

        return res;
    }

    char *cCifString;
    tRet = LLKA_cifDataToString(cCifData, pretty ? LLKA_TRUE : LLKA_FALSE, &cCifString);

    LLKA_destroyCifData(cCifData);

    if (tRet != LLKA_OK)
        return RCResult<std::string>::fail(tRet);
    auto res = RCResult<std::string>::succeed(cCifString);

    LLKA_destroyString(cCifString);

    return res;
}

//
// Nucleotide
//

auto extractNucleotide(const Structure &stru, int32_t pdbx_PDB_model_num, const std::string &label_asym_id, int32_t label_seq_id) -> Structure
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_Structure cNucl = LLKA_extractNucleotide(&wStru.cStru, pdbx_PDB_model_num, label_asym_id.c_str(), label_seq_id);

    Structure nucl;
    helpers::cStruToStru(cNucl, nucl);
    LLKA_destroyStructure(&cNucl);

    return nucl;
}

auto extractNucleotideView(const Structure &stru, int32_t pdbx_PDB_model_num, const std::string &label_asym_id, int32_t label_seq_id) -> StructureView
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_StructureView cView = LLKA_extractNucleotideView(&wStru.cStru, pdbx_PDB_model_num, label_asym_id.c_str(), label_seq_id);

    StructureView view;
    helpers::cStruViewToStruView(cView, wStru, view, stru);
    LLKA_destroyStructureView(&cView);

    return view;
}

auto extractRibose(const Structure &stru) -> RCResult<Structure>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_Structure cRibose;
    auto tRet = LLKA_extractRibose(&wStru.cStru, &cRibose);
    if (tRet != LLKA_OK)
        return RCResult<Structure>::fail(tRet);

    Structure ribose;
    helpers::cStruToStru(cRibose, ribose);
    LLKA_destroyStructure(&cRibose);

    return RCResult<Structure>::succeed(std::move(ribose));
}

auto extractRiboseView(const Structure &stru) -> RCResult<StructureView>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_StructureView cRibose;
    auto tRet = LLKA_extractRiboseView(&wStru.cStru, &cRibose);
    if (tRet != LLKA_OK)
        return RCResult<StructureView>::fail(tRet);

    StructureView ribose;
    helpers::cStruViewToStruView(cRibose, wStru, ribose, stru);
    LLKA_destroyStructureView(&cRibose);

    return RCResult<StructureView>::succeed(std::move(ribose));
}

auto isNucleotideCompound(const std::string &compId) -> bool
{
    return LLKA_isNucleotideCompound(compId.c_str()) == LLKA_TRUE;
}

auto riboseMetrics(const Structure &stru) -> RCResult<LLKA_RiboseMetrics>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_RiboseMetrics metrics;
    auto tRet = LLKA_riboseMetrics(&wStru.cStru, &metrics);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_RiboseMetrics>::fail(tRet);

    return RCResult<LLKA_RiboseMetrics>::succeed(metrics);
}

auto riboseMetrics(const StructureView &view) -> RCResult<LLKA_RiboseMetrics>
{
    const auto wStru = helpers::struViewToWrappedCStru(view);

    LLKA_RiboseMetrics metrics;
    auto tRet = LLKA_riboseMetrics(&wStru.cStru, &metrics);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_RiboseMetrics>::fail(tRet);

    return RCResult<LLKA_RiboseMetrics>::succeed(metrics);
}

auto sugarPucker(const Structure &stru) -> RCResult<LLKA_SugarPucker>
{
    const auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_SugarPucker pucker;
    auto tRet = LLKA_sugarPucker(&wStru.cStru, &pucker);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_SugarPucker>::fail(tRet);

    return RCResult<LLKA_SugarPucker>::succeed(pucker);
}

auto sugarPucker(const StructureView &view) -> RCResult<LLKA_SugarPucker>
{
    const auto wStru = helpers::struViewToWrappedCStru(view);

    LLKA_SugarPucker pucker;
    auto tRet = LLKA_sugarPucker(&wStru.cStru, &pucker);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_SugarPucker>::fail(tRet);

    return RCResult<LLKA_SugarPucker>::succeed(pucker);
}

// Classification - declaration of value objects that cannot be mapped
// directly to the C structs

GoldenStep::GoldenStep(const LLKA_GoldenStep &gs) :
    pucker_1{gs.pucker_1},
    pucker_2{gs.pucker_2},
    nuAngles_1{gs.nuAngles_1},
    nuAngles_2{gs.nuAngles_2},
    metrics{gs.metrics},
    name{gs.name},
    clusterIdx{gs.clusterIdx},
    clusterNumber{gs.clusterNumber}
{
}

//
// Resource loaders
//

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED

auto loadClusterNuAngles(const std::filesystem::path &path) -> RCResult<std::vector<LLKA_ClusterNuAngles>>
{
    using RT = RCResult<std::vector<LLKA_ClusterNuAngles>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_AVERAGE_NU_ANGLES;

    auto tRet = LLKA_loadResourceFile(path.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<LLKA_ClusterNuAngles> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = cResource.data.clusterNuAngles[idx];

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

auto loadClusters(const std::filesystem::path &path) -> RCResult<std::vector<LLKA_ClassificationCluster>>
{
    using RT = RCResult<std::vector<LLKA_ClassificationCluster>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_CLUSTERS;

    auto tRet = LLKA_loadResourceFile(path.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<LLKA_ClassificationCluster> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = cResource.data.clusters[idx];

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

auto loadConfals(const std::filesystem::path &path) -> RCResult<std::vector<LLKA_Confal>>
{
    using RT = RCResult<std::vector<LLKA_Confal>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_CONFALS;

    auto tRet = LLKA_loadResourceFile(path.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<LLKA_Confal> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = cResource.data.confals[idx];

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

auto loadConfalPercentiles(const std::filesystem::path &path) -> RCResult<std::vector<LLKA_ConfalPercentile>>
{
    using RT = RCResult<std::vector<LLKA_ConfalPercentile>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_CONFAL_PERCENTILES;

    auto tRet = LLKA_loadResourceFile(path.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<LLKA_ConfalPercentile> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = cResource.data.confalPercentiles[idx];

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

auto loadGoldenSteps(const std::filesystem::path &path) -> RCResult<std::vector<GoldenStep>>
{
    using RT = RCResult<std::vector<GoldenStep>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_GOLDEN_STEPS;

    auto tRet = LLKA_loadResourceFile(path.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<GoldenStep> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = GoldenStep{cResource.data.goldenSteps[idx]};

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

#endif // LLKA_FILESYSTEM_ACCESS_DISABLED

auto loadClusterNuAngles(const std::string &text) -> RCResult<std::vector<LLKA_ClusterNuAngles>>
{
    using RT = RCResult<std::vector<LLKA_ClusterNuAngles>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_AVERAGE_NU_ANGLES;

    auto tRet = LLKA_loadResourceText(text.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<LLKA_ClusterNuAngles> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = cResource.data.clusterNuAngles[idx];

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

auto loadClusters(const std::string &text) -> RCResult<std::vector<LLKA_ClassificationCluster>>
{
    using RT = RCResult<std::vector<LLKA_ClassificationCluster>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_CLUSTERS;

    auto tRet = LLKA_loadResourceText(text.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<LLKA_ClassificationCluster> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = cResource.data.clusters[idx];

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

auto loadConfals(const std::string &text) -> RCResult<std::vector<LLKA_Confal>>
{
    using RT = RCResult<std::vector<LLKA_Confal>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_CONFALS;

    auto tRet = LLKA_loadResourceText(text.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<LLKA_Confal> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = cResource.data.confals[idx];

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

auto loadConfalPercentiles(const std::string &text) -> RCResult<std::vector<LLKA_ConfalPercentile>>
{
    using RT = RCResult<std::vector<LLKA_ConfalPercentile>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_CONFAL_PERCENTILES;

    auto tRet = LLKA_loadResourceText(text.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<LLKA_ConfalPercentile> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = cResource.data.confalPercentiles[idx];

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

auto loadGoldenSteps(const std::string &text) -> RCResult<std::vector<GoldenStep>>
{
    using RT = RCResult<std::vector<GoldenStep>>;

    LLKA_Resource cResource = {};
    cResource.type = LLKA_RES_GOLDEN_STEPS;

    auto tRet = LLKA_loadResourceText(text.c_str(), &cResource);
    if (tRet != LLKA_OK)
        return RT::fail(tRet);

    std::vector<GoldenStep> resource{cResource.count};
    for (size_t idx = 0; idx < cResource.count; idx++)
        resource[idx] = GoldenStep{cResource.data.goldenSteps[idx]};

    LLKA_destroyResource(&cResource);

    return RT::succeed(std::move(resource));
}

//
// Tracing
//

auto toggleAllTracepoints(bool enable) -> void
{
    LLKA_toggleAllTracepoints(enable ? LLKA_TRUE : LLKA_FALSE);
}

auto toggleTracepoint(int32_t TPID, bool enable) -> void
{
    LLKA_toggleTracepoint(TPID, enable ? LLKA_TRUE : LLKA_FALSE);
}

auto trace(bool dontClear) -> std::string
{
    auto cTrace = LLKA_trace(dontClear ? LLKA_TRUE : LLKA_FALSE);
    std::string str{cTrace};

    LLKA_destroyTrace(cTrace);

    return str;
}

auto tracepointInfo() -> std::vector<TracepointInfo>
{
  auto cTpInfos = LLKA_tracepointInfo();

  std::vector<TracepointInfo> tpInfos{cTpInfos.nInfos};
  for (size_t idx = 0; idx < cTpInfos.nInfos; idx++) {
      tpInfos[idx].TPID = cTpInfos.infos[idx].TPID;
      tpInfos[idx].description = cTpInfos.infos[idx].description;
  }

  LLKA_destroyTracepointInfo(&cTpInfos);

  return tpInfos;
}

LLKA_CPP_API
auto tracepointState(int32_t TPID) -> bool
{
    return LLKA_tracepointState(TPID);
}

// Classification

ClassificationContext::ClassificationContext() :
    m_ctx{nullptr},
    m_movedAway{false}
{
    // Needed for Emscripten builds.
}

ClassificationContext::ClassificationContext(LLKA_ClassificationContext *ctx) :
    m_ctx{ctx},
    m_movedAway{false}
{
}

ClassificationContext::~ClassificationContext()
{
    if (!m_movedAway && m_ctx)
        LLKA_destroyClassificationContext(m_ctx);
}

// We need to cater to Emscripten requiring copy c-tor
#ifdef LLKA_PLATFORM_EMSCRIPTEN
ClassificationContext::ClassificationContext(const ClassificationContext &other)
{
    this->m_ctx = other.m_ctx;
    this->m_movedAway = false;

    other.m_movedAway = true;
}
#endif // LLKA_PLATFORM_EMSCRIPTEN

ClassificationContext::ClassificationContext(ClassificationContext &&other) noexcept
{
    this->m_ctx = other.m_ctx;
    this->m_movedAway = false;

    other.m_movedAway = true;
}

ClassificationContext & ClassificationContext::operator=(ClassificationContext &&other) noexcept
{
    this->m_ctx = other.m_ctx;
    this->m_movedAway = false;

    other.m_movedAway = true;

    return *this;
}

const LLKA_ClassificationContext * ClassificationContext::get() const {
    if (!isValid())
        throw std::runtime_error{"Attempted to get invalid ClassificationContext"};
    return m_ctx;
}

auto ClassificationContext::isValid() const -> bool
{
    return !m_movedAway && m_ctx;
}

ClassifiedStep::ClassifiedStep() :
    assignedNtC{LLKA_INVALID_NTC},
    assignedCANA{LLKA_INVALID_CANA},
    closestNtC{LLKA_INVALID_NTC},
    closestCANA{LLKA_INVALID_CANA},
    sugarPucker_1{LLKA_INVALID_SUGAR_PUCKER},
    sugarPucker_2{LLKA_INVALID_SUGAR_PUCKER},
    violations{0},
    violatingTorsionsAverage{0},
    violatingTorsionsNearest{0}
{
}

ClassifiedStep::ClassifiedStep(const LLKA_ClassifiedStep &cStep)
{
    assignedNtC = cStep.assignedNtC;
    assignedCANA = cStep.assignedCANA;
    closestNtC = cStep.closestNtC;
    closestCANA = cStep.closestCANA;
    confalScore = cStep.confalScore;
    euclideanDistanceNtCIdeal = cStep.euclideanDistanceNtCIdeal;
    metrics = cStep.metrics;
    differencesFromNtCAverages = cStep.differencesFromNtCAverages;
    nuAngles_1 = cStep.nuAngles_1;
    nuAngles_2 = cStep.nuAngles_2;
    ribosePseudorotation_1 = cStep.ribosePseudorotation_1;
    ribosePseudorotation_2 = cStep.ribosePseudorotation_2;
    tau_1 = cStep.tau_1;
    tau_2 = cStep.tau_2;
    sugarPucker_1 = cStep.sugarPucker_1;
    sugarPucker_2 = cStep.sugarPucker_2;
    nuAngleDifferences_1 = cStep.nuAngleDifferences_1;
    nuAngleDifferences_2 = cStep.nuAngleDifferences_2;
    rmsdToClosestNtC = cStep.rmsdToClosestNtC;
    closestGoldenStep = cStep.closestGoldenStep;
    violations = cStep.violations;
    violatingTorsionsAverage = cStep.violatingTorsionsAverage;
    violatingTorsionsNearest = cStep.violatingTorsionsNearest;
}

auto ClassifiedStep::hasViolations() const -> bool
{
    return violations != LLKA_CLASSIFICATION_OK;
}

auto ClassifiedStep::namedViolations() const -> std::vector<std::string>
{
    std::vector<std::string> names{};

    for (int32_t bit = 0; bit < 32; bit++) {
        auto viol = (1 << bit);
        if (violations & viol)
            names.push_back(LLKA_classificationViolationToName(viol));
    }

    return names;
}

_EMX_GET_SET_DEF(LLKA_NtC, ClassifiedStep, assignedNtC)
_EMX_GET_SET_DEF(LLKA_CANA, ClassifiedStep, assignedCANA)
_EMX_GET_SET_DEF(LLKA_NtC, ClassifiedStep, closestNtC)
_EMX_GET_SET_DEF(LLKA_CANA, ClassifiedStep, closestCANA)
_EMX_GET_SET_DEF(LLKA_ConfalScore, ClassifiedStep, confalScore)
_EMX_GET_SET_DEF(double, ClassifiedStep, euclideanDistanceNtCIdeal)
_EMX_GET_SET_DEF(LLKA_StepMetrics, ClassifiedStep, metrics)
_EMX_GET_SET_DEF(LLKA_StepMetrics, ClassifiedStep, differencesFromNtCAverages)
_EMX_GET_SET_DEF(LLKA_NuAngles, ClassifiedStep, nuAngles_1)
_EMX_GET_SET_DEF(LLKA_NuAngles, ClassifiedStep, nuAngles_2)
_EMX_GET_SET_DEF(double, ClassifiedStep, ribosePseudorotation_1)
_EMX_GET_SET_DEF(double, ClassifiedStep, ribosePseudorotation_2)
_EMX_GET_SET_DEF(double, ClassifiedStep, tau_1)
_EMX_GET_SET_DEF(double, ClassifiedStep, tau_2)
_EMX_GET_SET_DEF(LLKA_SugarPucker, ClassifiedStep, sugarPucker_1)
_EMX_GET_SET_DEF(LLKA_SugarPucker, ClassifiedStep, sugarPucker_2)
_EMX_GET_SET_DEF(LLKA_NuAngles, ClassifiedStep, nuAngleDifferences_1)
_EMX_GET_SET_DEF(LLKA_NuAngles, ClassifiedStep, nuAngleDifferences_2)
_EMX_GET_SET_DEF(double, ClassifiedStep, rmsdToClosestNtC)
_EMX_GET_SET_DEF(std::string, ClassifiedStep, closestGoldenStep)
_EMX_GET_SET_DEF(int32_t, ClassifiedStep, violations)
_EMX_GET_SET_DEF(int16_t, ClassifiedStep, violatingTorsionsAverage)
_EMX_GET_SET_DEF(int16_t, ClassifiedStep, violatingTorsionsNearest)

auto averageConfal(const std::vector<ClassifiedStep> &steps, const ClassificationContext &ctx) noexcept -> LLKA_AverageConfal
{
    std::vector<LLKA_ClassifiedStep> cClassifiedSteps(steps.size());
    for (size_t idx = 0; idx < steps.size(); idx++)
        cClassifiedSteps[idx] = helpers::stepToCStep(steps[idx]);

    return LLKA_averageConfal(cClassifiedSteps.data(), cClassifiedSteps.size(), ctx.get());
}

auto averageConfal(const std::vector<AttemptedClassifiedStep> &steps, const ClassificationContext &ctx) noexcept -> LLKA_AverageConfal
{
    std::vector<LLKA_AttemptedClassifiedStep> cAttemptedClassifiedSteps(steps.size());
    for (size_t idx = 0; idx < steps.size(); idx++) {
        cAttemptedClassifiedSteps[idx].step = helpers::stepToCStep(steps[idx].step);
        cAttemptedClassifiedSteps[idx].status = steps[idx].status;
    }

    LLKA_ClassifiedSteps cSteps {
        .attemptedSteps = cAttemptedClassifiedSteps.data(),
        .nAttemptedSteps = cAttemptedClassifiedSteps.size(),
    };

    return LLKA_averageConfalAttempted(&cSteps, ctx.get());
}

auto classificationClusterForNtC(LLKA_NtC ntc, const ClassificationContext &ctx) noexcept -> RCResult<LLKA_ClassificationCluster>
{
    LLKA_ClassificationCluster cluster;

    auto tRet = LLKA_classificationClusterForNtC(ntc, ctx.get(), &cluster);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_ClassificationCluster>::fail(tRet);
    return RCResult<LLKA_ClassificationCluster>::succeed(std::move(cluster));
}

auto classifyStep(const Structure &stru, const ClassificationContext &ctx) -> RCResult<ClassifiedStep>
{
    auto wStru = helpers::struToWrappedCStru(stru);

    LLKA_ClassifiedStep cClassifiedStep;
    auto tRet = LLKA_classifyStep(&wStru.cStru, ctx.get(), &cClassifiedStep);
    if (tRet != LLKA_OK)
        return RCResult<ClassifiedStep>::fail(tRet);
    return RCResult<ClassifiedStep>::succeed(cClassifiedStep);
}

auto classifySteps(const Structures &strus, const ClassificationContext &ctx) -> RCResult<AttemptedClassifiedSteps>
{
    auto wStrus = helpers::strusToWrappedCStrus(strus);
    const auto ds = wStrus.get();

    LLKA_ClassifiedSteps cClassifiedSteps;
    auto tRet = LLKA_classifyStepsMultiple(&ds, ctx.get(), &cClassifiedSteps);
    if (tRet != LLKA_OK)
        return RCResult<AttemptedClassifiedSteps>::fail(tRet);

    std::vector<AttemptedClassifiedStep> attemptedSteps{cClassifiedSteps.nAttemptedSteps};
    for (size_t idx = 0; idx < cClassifiedSteps.nAttemptedSteps; idx++) {
        attemptedSteps[idx].status = cClassifiedSteps.attemptedSteps[idx].status;
        attemptedSteps[idx].step = ClassifiedStep(cClassifiedSteps.attemptedSteps[idx].step);
    }

    LLKA_destroyClassifiedSteps(&cClassifiedSteps);

    return RCResult<AttemptedClassifiedSteps>::succeed(std::move(attemptedSteps));
}

LLKA_CPP_API
auto confalPercentile(double confalScore, const ClassificationContext &ctx) noexcept -> double
{
    return LLKA_confalPercentile(confalScore, ctx.get());
}

LLKA_CPP_API
auto confalForNtC(LLKA_NtC ntc, const ClassificationContext &ctx) noexcept -> RCResult<LLKA_Confal>
{
    LLKA_Confal confal;

    auto tRet = LLKA_confalForNtC(ntc, ctx.get(), &confal);
    if (tRet != LLKA_OK)
        return RCResult<LLKA_Confal>::fail(tRet);
    return RCResult<LLKA_Confal>::succeed(std::move(confal));
}

auto initializeClassificationContext(
    const std::vector<LLKA_ClassificationCluster> &clusters,
    const std::vector<GoldenStep> &goldenSteps,
    const std::vector<LLKA_Confal> &confals,
    const std::vector<LLKA_ClusterNuAngles> &clusterNuAngles,
    const std::vector<LLKA_ConfalPercentile> &confalPercentiles,
    const LLKA_ClassificationLimits &limits,
    double maxCloseEnoughRmsd
) -> RCResult<ClassificationContext>
{
    LLKA_ClassificationContext *ctx;

    // Do a shallow-ish copy of the GoldenSteps to match get to the expected LLKA_GoldenStep type
    std::vector<LLKA_GoldenStep> _goldenSteps{goldenSteps.size()};
    for (size_t idx = 0; idx < goldenSteps.size(); idx++) {
        const auto &gs = goldenSteps[idx];
        auto &lgs = _goldenSteps[idx];

        lgs.pucker_1 = gs.pucker_1;
        lgs.pucker_2 = gs.pucker_2;
        lgs.nuAngles_1 = gs.nuAngles_1;
        lgs.nuAngles_2 = gs.nuAngles_2;
        lgs.metrics = gs.metrics;
        lgs.name = gs.name.c_str();   // BEWARE of this pointer. goldenSteps vector must stay alive longer than _goldenSteps!
        lgs.clusterIdx = gs.clusterIdx;
        lgs.clusterNumber = gs.clusterNumber;
    }

    auto tRet = LLKA_initializeClassificationContext(
        clusters.data(), clusters.size(),
        _goldenSteps.data(), _goldenSteps.size(),
        confals.data(), confals.size(),
        clusterNuAngles.data(), clusterNuAngles.size(),
        confalPercentiles.data(), confalPercentiles.size(),
        &limits,
        maxCloseEnoughRmsd,
        &ctx
    );
    if (tRet != LLKA_OK)
        return RCResult<ClassificationContext>::fail(tRet);
    return RCResult<ClassificationContext>::succeed(ctx);
}

auto nameToSugarPucker(const std::string &name) -> LLKA_SugarPucker
{
    return LLKA_nameToSugarPucker(name.c_str());
}

auto sugarPuckerToName(LLKA_SugarPucker pucker, LLKA_SugarPuckerNameBrevity brevity) -> std::string
{
    return std::string{LLKA_sugarPuckerToName(pucker, brevity)};
}

namespace helpers {

    auto atomToCAtom(const Atom &atom, LLKA_Atom &cAtom) -> void
    {
        cAtom.type_symbol = atom.type_symbol.c_str();
        cAtom.label_atom_id = atom.label_atom_id.c_str();
        cAtom.label_entity_id = atom.label_entity_id.c_str();
        cAtom.label_comp_id = atom.label_comp_id.c_str();
        cAtom.label_asym_id = atom.label_asym_id.c_str();
        cAtom.auth_atom_id = atom.auth_atom_id.c_str();
        cAtom.auth_comp_id = atom.auth_comp_id.c_str();
        cAtom.auth_asym_id = atom.auth_asym_id.c_str();
        cAtom.coords = atom.coords;
        cAtom.id = atom.id;
        cAtom.label_seq_id = atom.label_seq_id;
        cAtom.auth_seq_id = atom.auth_seq_id;
        cAtom.pdbx_PDB_ins_code = atom.pdbx_PDB_ins_code.c_str();
        cAtom.pdbx_PDB_model_num = atom.pdbx_PDB_model_num;
        cAtom.label_alt_id = atom.label_alt_id;
    }

    auto stepToCStep(const ClassifiedStep &step) -> LLKA_ClassifiedStep
    {
        return {
            .assignedNtC = step.assignedNtC,
            .assignedCANA = step.assignedCANA,
            .closestNtC = step.closestNtC,
            .closestCANA = step.closestCANA,
            .confalScore = step.confalScore,
            .euclideanDistanceNtCIdeal = step.euclideanDistanceNtCIdeal,
            .metrics = step.metrics,
            .differencesFromNtCAverages = step.differencesFromNtCAverages,
            .nuAngles_1 = step.nuAngles_1,
            .nuAngles_2 = step.nuAngles_2,
            .ribosePseudorotation_1 = step.ribosePseudorotation_1,
            .ribosePseudorotation_2 = step.ribosePseudorotation_2,
            .tau_1 = step.tau_1,
            .tau_2 = step.tau_2,
            .sugarPucker_1 = step.sugarPucker_1,
            .sugarPucker_2 = step.sugarPucker_2,
            .nuAngleDifferences_1 = step.nuAngleDifferences_1,
            .nuAngleDifferences_2 = step.nuAngleDifferences_2,
            .rmsdToClosestNtC = step.rmsdToClosestNtC,

            // WARNING: This requires that the step we are converting from stays alive long enough!
            .closestGoldenStep = step.closestGoldenStep.c_str(),

            .violations = step.violations,
            .violatingTorsionsAverage = step.violatingTorsionsAverage,
            .violatingTorsionsNearest = step.violatingTorsionsNearest,
        };
    }

    auto struToCStru(const Structure &stru, LLKA_Structure &cStru, std::vector<LLKA_Atom> &cAtoms) -> void
    {
        // TODO: Hide this behind an "unlikely" if size check?
        cAtoms.resize(stru.size());

        for (size_t idx = 0; idx < stru.size(); idx++) {
            const auto &atom = stru[idx];
            auto &cAtom = cAtoms[idx];
            atomToCAtom(atom, cAtom);
        }

        cStru.atoms = cAtoms.data();
        cStru.nAtoms = cAtoms.size();
    }

    auto struViewToCStru(const StructureView &view, LLKA_Structure &cStru, std::vector<LLKA_Atom> &cAtoms) -> void
    {
        // TODO: Hide this behind an "unlikely" if size check?
        cAtoms.resize(view.size());

        for (size_t idx = 0; idx < view.size(); idx++) {
            const auto &atom = *view[idx];
            auto &cAtom = cAtoms[idx];
            atomToCAtom(atom, cAtom);
        }

        cStru.atoms = cAtoms.data();
        cStru.nAtoms = cAtoms.size();
    }

    auto struToWrappedCStru(const Structure &stru) -> WrappedCStructure
    {
        WrappedCStructure wrap;
        struToCStru(stru, wrap.cStru, wrap.cAtoms);

        return wrap;
    }

    auto struViewToWrappedCStru(const StructureView &view) -> WrappedCStructure
    {
        WrappedCStructure wrap;
        struViewToCStru(view, wrap.cStru, wrap.cAtoms);

        return wrap;
    }

    auto strusToWrappedCStrus(const Structures &strus) -> WrappedCStructures {
        std::vector<LLKA_Structure> cStrus{strus.size()};

        for (size_t idx = 0; idx < strus.size(); idx++) {
            auto &s = strus[idx];
            auto &cS = cStrus[idx];

            auto cAtoms = std::unique_ptr<LLKA_Atom[]>(new LLKA_Atom[s.size()]);
            for (size_t jdx = 0; jdx < s.size(); jdx++)
                atomToCAtom(s[jdx], cAtoms[jdx]);

            cS.atoms = cAtoms.release();
            cS.nAtoms = s.size();
        }

        return WrappedCStructures{std::move(cStrus)};
    }

    auto cStruToStru(const LLKA_Structure &cStru) -> Structure
    {
        Structure stru{};
        cStruToStru(cStru, stru);

        return stru;
    }

    auto cStruToStru(const LLKA_Structure &cStru, Structure &stru) -> void
    {
        // TODO: Maybe just resize and overwrite instead of pushbacking?
        stru.clear();
        stru.reserve(cStru.nAtoms);

        for (size_t idx = 0; idx < cStru.nAtoms; idx++) {
            const auto &cAtom = cStru.atoms[idx];
            stru.emplace_back(
                makeAtom(
                    cAtom.type_symbol,
                    cAtom.label_atom_id,
                    cAtom.label_entity_id,
                    cAtom.label_comp_id,
                    cAtom.label_asym_id,
                    cAtom.auth_atom_id,
                    cAtom.auth_comp_id,
                    cAtom.auth_asym_id,
                    cAtom.coords,
                    cAtom.id,
                    cAtom.label_seq_id,
                    cAtom.auth_seq_id,
                    cAtom.pdbx_PDB_model_num,
                    cAtom.pdbx_PDB_ins_code,
                    cAtom.label_alt_id
                )
            );
        }
    }

    auto cStruViewToStruView(const LLKA_StructureView &cView, const WrappedCStructure &wStru, StructureView &view, const Structure &stru) -> void
    {
        /*
         * WARNING: Do not touch, do not touch, do not touch!!!
         */

        const LLKA_Atom *firstCAtomPtr = &wStru.cStru.atoms[0];
        const Atom *firstAtomPtr = stru.data();

        for (size_t idx = 0; idx < cView.nAtoms; idx++) {
            const LLKA_Atom *cAtomPtr = cView.atoms[idx];

            // This value represents the index of the viewed atom in the viewed structure
            // We take advantage of the fact that LLKA_Structure and C++ API Structure have
            // atoms ordered in the same way. With that, we can calculate the pointer into
            // the C++ API Structure vector,
            uintptr_t cDiff = cAtomPtr - firstCAtomPtr;
            const Atom *atomPtr = firstAtomPtr + cDiff;

            view.push_back(atomPtr);
        }
    }

    inline
    auto cCifDataToCifData(const LLKA_CifData *cCifData) -> CifData
    {
        CifData cifData;
        cifData.blocks.resize(cCifData->nBlocks);

        for (size_t blockIdx = 0; blockIdx < cCifData->nBlocks; blockIdx++) {
            const auto &cBlock = cCifData->blocks[blockIdx];
            auto &block = cifData.blocks[blockIdx];

            block.name = cBlock.name;

            auto cCat = cBlock.firstCategory;
            while (cCat != nullptr) {
                CifData::Category cat;

                cat.name = cCat->name;

                auto cItem = cCat->firstItem;
                while (cItem != nullptr) {
                    CifData::Item item;

                    item.keyword = cItem->keyword;

                    item.values.resize(cItem->nValues);
                    for (size_t valIdx = 0; valIdx < cItem->nValues; valIdx++) {
                        const auto cv = &cItem->values[valIdx];
                        auto &v = item.values[valIdx];

                        v.state = cv->state;

                        if (cv->state == LLKA_MINICIF_VALUE_SET)
                            v.text = cv->text;
                    }

                    cat.items.push_back(std::move(item));
                    cItem = LLKA_cifDataCategory_nextItem(cItem);
                }

                block.categories.push_back(std::move(cat));
                cCat = LLKA_cifDataBlock_nextCategory(cCat);
            }
        }

        return cifData;
    }

} // namespace helpers

} // namespace LLKA

#endif // __cplusplus
