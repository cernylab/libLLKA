/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include <llka_nucleotide.h>
#include "nucleotide.hpp"

LLKA_Structure LLKA_CC LLKA_extractNucleotide(const LLKA_Structure *stru, int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id)
{
    return LLKAInternal::extractNucleotide(stru, pdbx_PDB_model_num, label_asym_id, label_seq_id);
}

LLKA_StructureView LLKA_CC LLKA_extractNucleotideView(const LLKA_Structure *stru, int32_t pdbx_PDB_model_num, const char *label_asym_id, int32_t label_seq_id)
{
    return LLKAInternal::extractNucleotideView(stru, pdbx_PDB_model_num, label_asym_id, label_seq_id);
}

LLKA_RetCode LLKA_CC LLKA_extractRibose(const LLKA_Structure *stru, LLKA_Structure *riboseStru)
{
    if (stru->nAtoms < 1)
        return LLKA_E_MISSING_ATOMS;

    return LLKAInternal::extractRibose(stru, riboseStru);
}

LLKA_RetCode LLKA_CC LLKA_extractRiboseView(const LLKA_Structure *stru, LLKA_StructureView *riboseView)
{
    if (stru->nAtoms < 1)
        return LLKA_E_MISSING_ATOMS;

    return LLKAInternal::extractRibose(stru, riboseView);
}

LLKA_Bool LLKA_CC LLKA_isNucleotideCompound(const char *compId)
{
    return LLKAInternal::isNucleotideCompound(compId);
}

LLKA_SugarPucker LLKA_CC LLKA_nameToSugarPucker(const char *name)
{
    auto it = LLKAInternal::NAME_TO_SUGAR_PUCKER_MAPPING.find(name);
    if (it != LLKAInternal::NAME_TO_SUGAR_PUCKER_MAPPING.cend())
        return it->second;
    return LLKA_INVALID_SUGAR_PUCKER;
}

LLKA_RetCode LLKA_CC LLKA_riboseMetrics(const LLKA_Structure *stru, LLKA_RiboseMetrics *metrics)
{
    LLKA_StructureView riboseView;
    auto tRet = LLKAInternal::extractRibose(stru, &riboseView);
    if (tRet != LLKA_OK)
        return LLKA_E_INVALID_ARGUMENT;

    LLKAInternal::riboseMetrics(riboseView, *metrics);
    LLKA_destroyStructureView(&riboseView);

    return LLKA_OK;
}

LLKA_RetCode LLKA_CC LLKA_riboseMetricsView(const LLKA_StructureView *view, LLKA_RiboseMetrics *metrics)
{
    LLKA_StructureView riboseView;
    auto tRet = LLKAInternal::extractRibose(view, &riboseView);
    if (tRet != LLKA_OK)
        return LLKA_E_INVALID_ARGUMENT;

    LLKAInternal::riboseMetrics(riboseView, *metrics);
    LLKA_destroyStructureView(&riboseView);

    return LLKA_OK;
}

LLKA_RetCode LLKA_CC LLKA_sugarPucker(const LLKA_Structure *stru, LLKA_SugarPucker *pucker)
{
    LLKA_StructureView riboseView;
    auto tRet = LLKAInternal::extractRibose(stru, &riboseView);
    if (tRet != LLKA_OK)
        return LLKA_E_INVALID_ARGUMENT;

    LLKA_RiboseMetrics metrics;
    LLKAInternal::riboseMetrics(riboseView, metrics);

    *pucker = metrics.pucker;

    LLKA_destroyStructureView(&riboseView);

    return LLKA_OK;
}

LLKA_RetCode LLKA_CC LLKA_sugarPuckerView(const LLKA_StructureView *view, LLKA_SugarPucker *pucker)
{
    LLKA_StructureView riboseView;
    auto tRet = LLKAInternal::extractRibose(view, &riboseView);
    if (tRet != LLKA_OK)
        return LLKA_E_INVALID_ARGUMENT;

    LLKA_RiboseMetrics metrics;
    LLKAInternal::riboseMetrics(riboseView, metrics);

    *pucker = metrics.pucker;

    LLKA_destroyStructureView(&riboseView);

    return LLKA_OK;
}

const char * LLKA_CC LLKA_sugarPuckerToName(LLKA_SugarPucker pucker, LLKA_SugarPuckerNameBrevity brevity)
{
    if (pucker == LLKA_INVALID_SUGAR_PUCKER)
        return "";
    auto idx = static_cast<size_t>(pucker);

    switch (brevity) {
    case LLKA_SPN_VERY_TERSE:
        return LLKAInternal::SUGAR_PUCKER_NAMES_VERY_TERSE[idx];
    case LLKA_SPN_TERSE:
        return LLKAInternal::SUGAR_PUCKER_NAMES_TERSE[idx];
    case LLKA_SPN_FANCY:
        return LLKAInternal::SUGAR_PUCKER_NAMES_FANCY[idx];
    }

    return "";
}
