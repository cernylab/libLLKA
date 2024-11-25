/* vim: set sw=4 ts=4 sts=4 expandtab : */

#include "../include/dnatco.h"

#include <llka_util.h>

#include <algorithm>
#include <cmath>
#include <cstdio>


namespace DNATCOExtras {

Error::Error(LLKA_RetCode tRet) noexcept :
    m_tRet{tRet}
{
}

} // namespace DNATCOExtras

static
auto findSecondAtom(const LLKA::Structure &stru) -> const LLKA::Atom &
{
    const auto firstSeqId = stru.front().label_seq_id;
    for (const auto &atom : stru) {
        if (atom.label_seq_id != firstSeqId)
            return atom;
    }

    throw DNATCOExtras::Error{LLKA_E_BAD_DATA};
}

static
auto fmtNum(const char *format, auto number)
{
    char BUF[64];
    std::snprintf(BUF, 64, format, number);
    BUF[63] = '\0';

    return std::string{BUF};
}

static
auto getAltId(const LLKA::Structure &stru, int32_t seqId)
{
    for (const auto &atom : stru) {
        if (atom.label_seq_id != seqId)
            continue;

        if (atom.label_alt_id != LLKA_NO_ALTID)
            return atom.label_alt_id;
    }

    return LLKA_NO_ALTID;
}

static
auto makeItem(std::string keyword)
{
    LLKA::CifData::Item item;
    item.keyword = std::move(keyword);

    return item;
}

static
auto makeValue(std::string text, LLKA_CifDataValueState state)
{
    LLKA::CifData::Value v;
    v.text = std::move(text);
    v.state = state;

    return v;
}

static
auto makeNdbStructNtcOverallCategory()
{
    LLKA::CifData::Category cat;
    cat.name = "ndb_struct_ntc_overall";

    cat.items.push_back(makeItem("entry_id"));
    cat.items.push_back(makeItem("confal_score"));
    cat.items.push_back(makeItem("confal_percentile"));
    cat.items.push_back(makeItem("ntc_version"));
    cat.items.push_back(makeItem("cana_version"));
    cat.items.push_back(makeItem("num_steps"));
    cat.items.push_back(makeItem("num_classified"));
    cat.items.push_back(makeItem("num_unclassified"));
    cat.items.push_back(makeItem("num_unclassified_rmsd_close"));

    return cat;
}

static
auto makeNdbStructNtcStepCategory()
{
    LLKA::CifData::Category cat;
    cat.name = "ndb_struct_ntc_step";

    cat.items.push_back(makeItem("id"));
    cat.items.push_back(makeItem("name"));
    cat.items.push_back(makeItem("PDB_model_number"));
    cat.items.push_back(makeItem("label_entity_id_1"));
    cat.items.push_back(makeItem("label_asym_id_1"));
    cat.items.push_back(makeItem("label_seq_id_1"));
    cat.items.push_back(makeItem("label_comp_id_1"));
    cat.items.push_back(makeItem("label_alt_id_1"));
    cat.items.push_back(makeItem("label_entity_id_2"));
    cat.items.push_back(makeItem("label_asym_id_2"));
    cat.items.push_back(makeItem("label_seq_id_2"));
    cat.items.push_back(makeItem("label_comp_id_2"));
    cat.items.push_back(makeItem("label_alt_id_2"));
    cat.items.push_back(makeItem("auth_asym_id_1"));
    cat.items.push_back(makeItem("auth_seq_id_1"));
    cat.items.push_back(makeItem("auth_asym_id_2"));
    cat.items.push_back(makeItem("auth_seq_id_2"));
    cat.items.push_back(makeItem("PDB_ins_code_1"));
    cat.items.push_back(makeItem("PDB_ins_code_2"));

    return cat;
}

static
auto makeNdbStructNtcStepParametersCategory()
{
    LLKA::CifData::Category cat;
    cat.name = "ndb_struct_ntc_step_parameters";

    cat.items.push_back(makeItem("step_id"));
    cat.items.push_back(makeItem("tor_delta_1"));
    cat.items.push_back(makeItem("tor_epsilon_1"));
    cat.items.push_back(makeItem("tor_zeta_1"));
    cat.items.push_back(makeItem("tor_alpha_2"));
    cat.items.push_back(makeItem("tor_beta_2"));
    cat.items.push_back(makeItem("tor_gamma_2"));
    cat.items.push_back(makeItem("tor_delta_2"));
    cat.items.push_back(makeItem("tor_chi_1"));
    cat.items.push_back(makeItem("tor_chi_2"));
    cat.items.push_back(makeItem("dist_NN"));
    cat.items.push_back(makeItem("dist_CC"));
    cat.items.push_back(makeItem("tor_NCCN"));
    cat.items.push_back(makeItem("diff_tor_delta_1"));
    cat.items.push_back(makeItem("diff_tor_epsilon_1"));
    cat.items.push_back(makeItem("diff_tor_zeta_1"));
    cat.items.push_back(makeItem("diff_tor_alpha_2"));
    cat.items.push_back(makeItem("diff_tor_beta_2"));
    cat.items.push_back(makeItem("diff_tor_gamma_2"));
    cat.items.push_back(makeItem("diff_tor_delta_2"));
    cat.items.push_back(makeItem("diff_tor_chi_1"));
    cat.items.push_back(makeItem("diff_tor_chi_2"));
    cat.items.push_back(makeItem("diff_dist_NN"));
    cat.items.push_back(makeItem("diff_dist_CC"));
    cat.items.push_back(makeItem("diff_tor_NCCN"));
    cat.items.push_back(makeItem("confal_tor_delta_1"));
    cat.items.push_back(makeItem("confal_tor_epsilon_1"));
    cat.items.push_back(makeItem("confal_tor_zeta_1"));
    cat.items.push_back(makeItem("confal_tor_alpha_2"));
    cat.items.push_back(makeItem("confal_tor_beta_2"));
    cat.items.push_back(makeItem("confal_tor_gamma_2"));
    cat.items.push_back(makeItem("confal_tor_delta_2"));
    cat.items.push_back(makeItem("confal_tor_chi_1"));
    cat.items.push_back(makeItem("confal_tor_chi_2"));
    cat.items.push_back(makeItem("confal_dist_NN"));
    cat.items.push_back(makeItem("confal_dist_CC"));
    cat.items.push_back(makeItem("confal_tor_NCCN"));
    cat.items.push_back(makeItem("details"));

    return cat;
}

static
auto makeNdbStructNtCStepSummaryCategory()
{
    LLKA::CifData::Category cat;
    cat.name = "ndb_struct_ntc_step_summary";

    cat.items.push_back(makeItem("step_id"));
    cat.items.push_back(makeItem("assigned_CANA"));
    cat.items.push_back(makeItem("assigned_NtC"));
    cat.items.push_back(makeItem("confal_score"));
    cat.items.push_back(makeItem("euclidean_distance_NtC_ideal"));
    cat.items.push_back(makeItem("cartesian_rmsd_closest_NtC_representative"));
    cat.items.push_back(makeItem("closest_CANA"));
    cat.items.push_back(makeItem("closest_NtC"));
    cat.items.push_back(makeItem("closest_step_golden"));

    return cat;
}

static
auto makeNdbStructSugarStepParameters()
{
    LLKA::CifData::Category cat;
    cat.name = "ndb_struct_sugar_step_parameters";

    cat.items.push_back(makeItem("step_id"));
    cat.items.push_back(makeItem("P_1"));
    cat.items.push_back(makeItem("tau_1"));
    cat.items.push_back(makeItem("Pn_1"));
    cat.items.push_back(makeItem("P_2"));
    cat.items.push_back(makeItem("tau_2"));
    cat.items.push_back(makeItem("Pn_2"));
    cat.items.push_back(makeItem("nu_1_1"));
    cat.items.push_back(makeItem("nu_1_2"));
    cat.items.push_back(makeItem("nu_1_3"));
    cat.items.push_back(makeItem("nu_1_4"));
    cat.items.push_back(makeItem("nu_1_5"));
    cat.items.push_back(makeItem("nu_2_1"));
    cat.items.push_back(makeItem("nu_2_2"));
    cat.items.push_back(makeItem("nu_2_3"));
    cat.items.push_back(makeItem("nu_2_4"));
    cat.items.push_back(makeItem("nu_2_5"));
    cat.items.push_back(makeItem("diff_nu_1_1"));
    cat.items.push_back(makeItem("diff_nu_1_2"));
    cat.items.push_back(makeItem("diff_nu_1_3"));
    cat.items.push_back(makeItem("diff_nu_1_4"));
    cat.items.push_back(makeItem("diff_nu_1_5"));
    cat.items.push_back(makeItem("diff_nu_2_1"));
    cat.items.push_back(makeItem("diff_nu_2_2"));
    cat.items.push_back(makeItem("diff_nu_2_3"));
    cat.items.push_back(makeItem("diff_nu_2_4"));
    cat.items.push_back(makeItem("diff_nu_2_5"));

    return cat;
}


inline const std::vector<std::string> TORSION_NAMES = { "d1", "e1", "z1", "a2", "b2", "g2", "d2", "ch1", "ch2" };

static
auto makeDetailsString(const LLKA::ClassifiedStep &classifiedStep) {
    std::string detailsStr{};

    auto append = [&detailsStr](const std::string &s) {
        if (detailsStr.length() > 0)
            detailsStr += ";";
        detailsStr += s;
    };

    if (classifiedStep.violations & LLKA_CLASSIFICATION_E_AVERAGE_NEAREST_NEIGHBORS_TORSIONS_TOO_DIFFERENT) {
        for (size_t idx = 0; idx < 9; idx++) {
            if (classifiedStep.violatingTorsionsAverage & (1 << idx))
                append("cAn" + TORSION_NAMES[idx]);
        }
    }

    if (classifiedStep.violations & LLKA_CLASSIFICATION_E_NEAREST_NEIGHBOR_TORSIONS_TOO_DIFFERENT) {
        for (size_t idx = 0; idx < 9; idx++) {
            if (classifiedStep.violatingTorsionsNearest & (1 << idx))
                append("cNn" + TORSION_NAMES[idx]);
        }
    }

    if (classifiedStep.violations & LLKA_CLASSIFICATION_E_NN_TOO_LOW || classifiedStep.violations & LLKA_CLASSIFICATION_E_NN_TOO_HIGH)
        append("cNN");

    if (classifiedStep.violations & LLKA_CLASSIFICATION_E_CC_TOO_LOW || classifiedStep.violations & LLKA_CLASSIFICATION_E_CC_TOO_HIGH)
        append("cCC");

    if (classifiedStep.violations & LLKA_CLASSIFICATION_E_MU_TOO_LOW || classifiedStep.violations & LLKA_CLASSIFICATION_E_MU_TOO_HIGH)
        append("cmu");

    if (classifiedStep.violations & LLKA_CLASSIFICATION_E_MU_TOO_LOW || classifiedStep.violations & LLKA_CLASSIFICATION_E_MU_TOO_HIGH)
        append("cmu");

    if (classifiedStep.violations & LLKA_CLASSIFICATION_E_TOTAL_DISTANCE_TOO_HIGH)
        append("cMB");

    if (classifiedStep.violations & LLKA_CLASSIFICATION_E_FIRST_PSEUDOROTATION_TOO_DIFFERENT)
        append("cP");

    if (classifiedStep.violations & LLKA_CLASSIFICATION_E_SECOND_PSEUDOROTATION_TOO_DIFFERENT)
        append("cP1");

    return detailsStr;
}

static
auto toLowerCase(std::string str)
{
    //PERF: Have a special function to convert PDB ID to lower case
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

static
auto hasMultipleModels(const LLKA::CifData &cifData)
{
    const auto &cats = cifData.blocks[0].categories;
    const auto atomSite = std::find_if(
        cats.cbegin(),
        cats.cend(),
        [](const auto &cat) {
            auto lwr = toLowerCase(cat.name);
            return lwr == "atom_site";
        }
    );
    if (atomSite == cats.cend())
        return false; // This should never happen

    const auto modelNumItem = std::find_if(
        atomSite->items.cbegin(),
        atomSite->items.cend(),
        [](const auto &item) {
            auto lwr = toLowerCase(item.keyword);
            return lwr == "pdbx_pdb_model_num";
        }
    );
    if (modelNumItem == atomSite->items.cend())
        return false;

    const auto &values = modelNumItem->values;
    const auto num = values[0].state == LLKA_MINICIF_VALUE_SET ? std::stoi(values[0].text) : 0;
    for (size_t idx = 1; idx < modelNumItem->values.size(); idx++) {
        const auto _num = values[idx].state == LLKA_MINICIF_VALUE_SET ? std::stoi(values[idx].text) : 0;
        if (num != _num)
            return true;
    }
    return false;
}

static
auto makeDNATCOStepName(const std::string &entryId, bool hasMultipleModels, const LLKA::Atom &firstAtom, const LLKA::Atom &secondAtom, char firstAltId, char secondAltId) {
    const auto id = toLowerCase(entryId);
    const auto model = hasMultipleModels ? "-m" + std::to_string(firstAtom.pdbx_PDB_model_num) : "";
    const auto _firstAltId = firstAltId != LLKA_NO_ALTID ? "." + std::string{firstAltId} : "";
    const auto firstInsCode = firstAtom.pdbx_PDB_ins_code != LLKA_NO_INSCODE ? "." + std::string{firstAtom.pdbx_PDB_ins_code} : "";
    const auto _secondAltId = secondAltId != LLKA_NO_ALTID ? "." + std::string{secondAltId} : "";
    const auto secondInsCode = secondAtom.pdbx_PDB_ins_code != LLKA_NO_INSCODE ? "." + std::string{secondAtom.pdbx_PDB_ins_code} : "";

    return id + model + "_" + firstAtom.auth_asym_id + "_" + firstAtom.auth_comp_id + _firstAltId + "_" + std::to_string(firstAtom.auth_seq_id) + firstInsCode + "_" + secondAtom.auth_comp_id + _secondAltId + "_" + std::to_string(secondAtom.auth_seq_id) + secondInsCode;
}

namespace DNATCOExtras {

auto Error::retcode() const noexcept -> LLKA_RetCode
{
    return m_tRet;
}

auto Error::what() const noexcept -> const char *
{
    return LLKA_errorToString(m_tRet);
}

auto addDNATCOCategoriesToCif(LLKA::CifData cifData, const LLKA::AttemptedClassifiedSteps &attemptedSteps, const LLKA_AverageConfal &avgConfal, const LLKA::Structures &steps, const std::string &entryId) -> LLKA::CifData
{
    size_t totalStepsCount = 0;
    size_t assignedStepsCount = 0;
    size_t unassignedButCloseStepsCount = 0;

    const auto multipleModels = hasMultipleModels(cifData);

    auto overallCategory = makeNdbStructNtcOverallCategory();
    auto stepCategory = makeNdbStructNtcStepCategory();
    auto stepParametersCategory = makeNdbStructNtcStepParametersCategory();
    auto stepSummaryCategory = makeNdbStructNtCStepSummaryCategory();
    auto sugarCategory = makeNdbStructSugarStepParameters();

    if (attemptedSteps.size() != steps.size())
        throw Error{LLKA_E_MISMATCHING_SIZES};

    for (size_t idx = 0; idx < attemptedSteps.size(); idx++) {
        const auto &attemptedStep = attemptedSteps[idx];
        const auto &struStep = steps[idx];

        if (attemptedStep.status != LLKA_OK)
            continue;

        const auto &cs = attemptedStep.step;
        const auto &firstAtom = struStep.front();
        const auto &secondAtom = findSecondAtom(struStep);
        const auto firstAltId = getAltId(struStep, firstAtom.label_seq_id);
        const auto secondAltId = getAltId(struStep, secondAtom.label_seq_id);

        if (!cs.hasViolations())
            assignedStepsCount++;
        else if (cs.violations & LLKA_CLASSIFICATION_E_UNASSIGNED_BUT_CLOSE_ENOUGH)
            unassignedButCloseStepsCount++;

        const auto totalStepsCountStr = std::to_string(totalStepsCount + 1);

        // Fill out the summary category row
        stepSummaryCategory.items[0].values.push_back(makeValue(totalStepsCountStr, LLKA_MINICIF_VALUE_SET));
        stepSummaryCategory.items[1].values.push_back(makeValue(LLKA::CANAToName(cs.assignedCANA), LLKA_MINICIF_VALUE_SET));
        stepSummaryCategory.items[2].values.push_back(makeValue(LLKA::NtCToName(cs.assignedNtC), LLKA_MINICIF_VALUE_SET));
        stepSummaryCategory.items[3].values.push_back(makeValue(std::to_string(int(cs.confalScore.total + 0.5)), LLKA_MINICIF_VALUE_SET));
        stepSummaryCategory.items[4].values.push_back(makeValue(fmtNum("%.1f", cs.euclideanDistanceNtCIdeal), LLKA_MINICIF_VALUE_SET));
        stepSummaryCategory.items[5].values.push_back(makeValue(fmtNum("%.3f", cs.rmsdToClosestNtC), LLKA_MINICIF_VALUE_SET));
        stepSummaryCategory.items[6].values.push_back(makeValue(LLKA::CANAToName(cs.closestCANA), LLKA_MINICIF_VALUE_SET));
        stepSummaryCategory.items[7].values.push_back(makeValue(LLKA::NtCToName(cs.closestNtC), LLKA_MINICIF_VALUE_SET));
        stepSummaryCategory.items[8].values.push_back(makeValue(cs.closestGoldenStep, LLKA_MINICIF_VALUE_SET));

        // Fill out the NtC category row
        stepCategory.items[0].values.push_back(makeValue(totalStepsCountStr, LLKA_MINICIF_VALUE_SET));

        const auto stepName = makeDNATCOStepName(entryId, multipleModels, firstAtom, secondAtom, firstAltId, secondAltId);
        stepCategory.items[1].values.push_back(makeValue(stepName, LLKA_MINICIF_VALUE_SET));

        stepCategory.items[2].values.push_back(makeValue(std::to_string(firstAtom.pdbx_PDB_model_num), LLKA_MINICIF_VALUE_SET));
        stepCategory.items[3].values.push_back(makeValue(firstAtom.label_entity_id, LLKA_MINICIF_VALUE_SET));
        stepCategory.items[4].values.push_back(makeValue(firstAtom.label_asym_id, LLKA_MINICIF_VALUE_SET));
        stepCategory.items[5].values.push_back(makeValue(std::to_string(firstAtom.label_seq_id), LLKA_MINICIF_VALUE_SET));
        stepCategory.items[6].values.push_back(makeValue(firstAtom.label_comp_id, LLKA_MINICIF_VALUE_SET));
        if (firstAltId == LLKA_NO_ALTID)
            stepCategory.items[7].values.push_back(makeValue("", LLKA_MINICIF_VALUE_NONE));
        else
            stepCategory.items[7].values.push_back(makeValue(std::string{firstAltId}, LLKA_MINICIF_VALUE_SET));

        stepCategory.items[8].values.push_back(makeValue(secondAtom.label_entity_id, LLKA_MINICIF_VALUE_SET));
        stepCategory.items[9].values.push_back(makeValue(secondAtom.label_asym_id, LLKA_MINICIF_VALUE_SET));
        stepCategory.items[10].values.push_back(makeValue(std::to_string(secondAtom.label_seq_id), LLKA_MINICIF_VALUE_SET));
        stepCategory.items[11].values.push_back(makeValue(secondAtom.label_comp_id, LLKA_MINICIF_VALUE_SET));
        if (secondAltId == LLKA_NO_ALTID)
            stepCategory.items[12].values.push_back(makeValue("", LLKA_MINICIF_VALUE_NONE));
        else
            stepCategory.items[12].values.push_back(makeValue(std::string{secondAltId}, LLKA_MINICIF_VALUE_SET));

        stepCategory.items[13].values.push_back(makeValue(firstAtom.auth_asym_id, LLKA_MINICIF_VALUE_SET));
        stepCategory.items[14].values.push_back(makeValue(std::to_string(firstAtom.auth_seq_id), LLKA_MINICIF_VALUE_SET));
        stepCategory.items[15].values.push_back(makeValue(secondAtom.auth_asym_id, LLKA_MINICIF_VALUE_SET));
        stepCategory.items[16].values.push_back(makeValue(std::to_string(secondAtom.auth_seq_id), LLKA_MINICIF_VALUE_SET));

        if (firstAtom.pdbx_PDB_ins_code == LLKA_NO_INSCODE)
            stepCategory.items[17].values.push_back(makeValue("", LLKA_MINICIF_VALUE_NONE));
        else
            stepCategory.items[17].values.push_back(makeValue(firstAtom.pdbx_PDB_ins_code, LLKA_MINICIF_VALUE_SET));
        if (secondAtom.pdbx_PDB_ins_code == LLKA_NO_INSCODE)
            stepCategory.items[18].values.push_back(makeValue("", LLKA_MINICIF_VALUE_NONE));
        else
            stepCategory.items[18].values.push_back(makeValue(secondAtom.pdbx_PDB_ins_code, LLKA_MINICIF_VALUE_SET));

        // Fill step parameters row
        stepParametersCategory.items[0].values.push_back(makeValue(totalStepsCountStr, LLKA_MINICIF_VALUE_SET));

        stepParametersCategory.items[1].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.delta_1))), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[2].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.epsilon_1))), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[3].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.zeta_1))), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[4].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.alpha_2))), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[5].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.beta_2))), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[6].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.gamma_2))), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[7].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.delta_2))), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[8].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.chi_1))), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[9].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.chi_2))), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[10].values.push_back(makeValue(fmtNum("%.2f", cs.metrics.NN), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[11].values.push_back(makeValue(fmtNum("%.2f", cs.metrics.CC), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[12].values.push_back(makeValue(fmtNum("%.1f", LLKA_fullAngleFromDeg(LLKA_rad2deg(cs.metrics.mu))), LLKA_MINICIF_VALUE_SET));

        stepParametersCategory.items[13].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.delta_1)), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[14].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.epsilon_1)), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[15].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.zeta_1)), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[16].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.alpha_2)), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[17].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.beta_2)), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[18].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.gamma_2)), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[19].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.delta_2)), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[20].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.chi_1)), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[21].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.chi_2)), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[22].values.push_back(makeValue(fmtNum("%.2f", cs.differencesFromNtCAverages.NN), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[23].values.push_back(makeValue(fmtNum("%.2f", cs.differencesFromNtCAverages.CC), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[24].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.differencesFromNtCAverages.mu)), LLKA_MINICIF_VALUE_SET));

        stepParametersCategory.items[25].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.delta_1), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[26].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.epsilon_1), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[27].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.zeta_1), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[28].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.alpha_2), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[29].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.beta_2), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[30].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.gamma_2), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[31].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.delta_2), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[32].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.chi_1), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[33].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.chi_2), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[34].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.NN), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[35].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.CC), LLKA_MINICIF_VALUE_SET));
        stepParametersCategory.items[36].values.push_back(makeValue(fmtNum("%.1f", cs.confalScore.mu), LLKA_MINICIF_VALUE_SET));

        const auto details = makeDetailsString(cs);
        if (!details.empty())
            stepParametersCategory.items[37].values.push_back(makeValue(details, LLKA_MINICIF_VALUE_SET));
        else
            stepParametersCategory.items[37].values.push_back(makeValue("", LLKA_MINICIF_VALUE_NONE));

        // Fill out sugar category row
        sugarCategory.items[0].values.push_back(makeValue(totalStepsCountStr, LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[1].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.ribosePseudorotation_1))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[2].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.tau_1))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[3].values.push_back(makeValue(LLKA::sugarPuckerToName(cs.sugarPucker_1, LLKA_SPN_VERY_TERSE), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[4].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.ribosePseudorotation_2))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[5].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.tau_2))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[6].values.push_back(makeValue(LLKA::sugarPuckerToName(cs.sugarPucker_2, LLKA_SPN_VERY_TERSE), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[7].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_1.nu_0))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[8].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_1.nu_1))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[9].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_1.nu_2))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[10].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_1.nu_3))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[11].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_1.nu_4))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[12].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_2.nu_0))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[13].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_2.nu_1))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[14].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_2.nu_2))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[15].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_2.nu_3))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[16].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(LLKA_fullAngleFromRad(cs.nuAngles_2.nu_4))), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[17].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_1.nu_0)), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[18].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_1.nu_1)), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[19].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_1.nu_2)), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[20].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_1.nu_3)), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[21].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_1.nu_4)), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[22].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_2.nu_0)), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[23].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_2.nu_1)), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[24].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_2.nu_2)), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[25].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_2.nu_3)), LLKA_MINICIF_VALUE_SET));
        sugarCategory.items[26].values.push_back(makeValue(fmtNum("%.1f", LLKA_rad2deg(cs.nuAngleDifferences_2.nu_4)), LLKA_MINICIF_VALUE_SET));

        totalStepsCount++;
    }

    // Fill out the overall category
    overallCategory.items[0].values.push_back(makeValue(entryId, LLKA_MINICIF_VALUE_SET));
    overallCategory.items[1].values.push_back(makeValue(fmtNum("%.1f", avgConfal.score), LLKA_MINICIF_VALUE_SET));
    overallCategory.items[2].values.push_back(makeValue(fmtNum("%d", int(avgConfal.percentile)), LLKA_MINICIF_VALUE_SET));
    overallCategory.items[3].values.push_back(makeValue(LLKA_INTERNAL_NTC_VERSION, LLKA_MINICIF_VALUE_SET));
    overallCategory.items[4].values.push_back(makeValue(LLKA_INTERNAL_CANA_VERSION, LLKA_MINICIF_VALUE_SET));
    overallCategory.items[5].values.push_back(makeValue(std::to_string(totalStepsCount), LLKA_MINICIF_VALUE_SET));
    overallCategory.items[6].values.push_back(makeValue(std::to_string(assignedStepsCount), LLKA_MINICIF_VALUE_SET));
    overallCategory.items[7].values.push_back(makeValue(std::to_string(totalStepsCount - assignedStepsCount), LLKA_MINICIF_VALUE_SET));
    overallCategory.items[8].values.push_back(makeValue(std::to_string(unassignedButCloseStepsCount), LLKA_MINICIF_VALUE_SET));

    auto &block = cifData.blocks[0];
    block.categories.push_back(overallCategory);
    block.categories.push_back(stepCategory);
    block.categories.push_back(stepParametersCategory);
    block.categories.push_back(stepSummaryCategory);
    block.categories.push_back(sugarCategory);

    return cifData;
}

} // namespace DNATCOExtras
