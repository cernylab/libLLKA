/* vim: set sw=4 ts=4 sts=4 expandtab : */

#ifndef _LLKA_NTC_H
#define _LLKA_NTC_H

#include "llka_structure.h"

#define LLKA_INTERNAL_NTC_VERSION "3.5"
#define LLKA_INTERNAL_CANA_VERSION "2.3"

/*!
 * NtC step metrics that span across the first and second residue of the step
 */
typedef enum LLKA_CrossResidueMetric {
    LLKA_XR_DIST_CC,
    LLKA_XR_DIST_NN,
    LLKA_XR_TOR_MU
    ENUM_FORCE_INT32_SIZE(LLKA_CrossResidueMetric)
} LLKA_CrossResidueMetric;

/*!
 * Kinds of backbone torsion on a dinucleotide
 */
typedef enum LLKA_DinucleotideTorsion {
    LLKA_TOR_DELTA_1,
    LLKA_TOR_EPSILON_1,
    LLKA_TOR_ZETA_1,
    LLKA_TOR_ALPHA_2,
    LLKA_TOR_BETA_2,
    LLKA_TOR_GAMMA_2,
    LLKA_TOR_DELTA_2,
    LLKA_TOR_CHI_1,
    LLKA_TOR_CHI_2
    ENUM_FORCE_INT32_SIZE(LLKA_DinucleotideTorsion)
} LLKA_DinucleotideTorsion;

/*!
 * List of recognized CANA classes
 */
typedef enum LLKA_CANA {
    LLKA_AAA = 0, LLKA_AAw, LLKA_AAu,
    LLKA_AB, LLKA_BA,
    LLKA_BBB, LLKA_BBw, LLKA_B12, LLKA_BB2, LLKA_miB,
    LLKA_ICL,
    LLKA_OPN,
    LLKA_SYN,
    LLKA_ZZZ,
    LLKA_INVALID_CANA = 0x7FFFFFFF
} LLKA_CANA;
#define LLKA_LAST_CANA LLKA_ZZZ

/*!
 * List of recognized NtCs
 */
typedef enum LLKA_NtC {
    LLKA_AA00 = 0, LLKA_AA02, LLKA_AA03, LLKA_AA04, LLKA_AA08, LLKA_AA09, LLKA_AA01, LLKA_AA05, LLKA_AA06, LLKA_AA10, LLKA_AA11, LLKA_AA07, LLKA_AA12, LLKA_AA13,
    LLKA_AB01, LLKA_AB02, LLKA_AB03, LLKA_AB04, LLKA_AB05,
    LLKA_BA01, LLKA_BA05, LLKA_BA09, LLKA_BA08, LLKA_BA10, LLKA_BA13, LLKA_BA16, LLKA_BA17,
    LLKA_BB00, LLKA_BB01, LLKA_BB17, LLKA_BB02, LLKA_BB03, LLKA_BB11, LLKA_BB16, LLKA_BB04, LLKA_BB05, LLKA_BB07, LLKA_BB08, LLKA_BB10, LLKA_BB12, LLKA_BB13, LLKA_BB14, LLKA_BB15, LLKA_BB20,
    LLKA_IC01, LLKA_IC02, LLKA_IC03, LLKA_IC04, LLKA_IC05, LLKA_IC06, LLKA_IC07,
    LLKA_OP01, LLKA_OP02, LLKA_OP03, LLKA_OP04, LLKA_OP05, LLKA_OP06, LLKA_OP07, LLKA_OP08, LLKA_OP09, LLKA_OP10, LLKA_OP11, LLKA_OP12, LLKA_OP13, LLKA_OP14, LLKA_OP15, LLKA_OP16, LLKA_OP17, LLKA_OP18, LLKA_OP19, LLKA_OP20, LLKA_OP21, LLKA_OP22, LLKA_OP23, LLKA_OP24, LLKA_OP25, LLKA_OP26, LLKA_OP27, LLKA_OP28, LLKA_OP29, LLKA_OP30, LLKA_OP31, LLKA_OPS1, LLKA_OP1S,
    LLKA_AAS1, LLKA_AB1S, LLKA_AB2S,
    LLKA_BB1S, LLKA_BB2S, LLKA_BBS1,
    LLKA_ZZ01, LLKA_ZZ02, LLKA_ZZ1S, LLKA_ZZ2S, LLKA_ZZS1, LLKA_ZZS2,
    LLKA_INVALID_NTC = 0x7FFFFFFF
} LLKA_NtC;
#define LLKA_LAST_NTC LLKA_ZZS2

/*!
 * Identifies order of residue in step
 */
typedef enum LLKA_ResidueInStep {
    LLKA_FIRST_RESIDUE,
    LLKA_SECOND_RESIDUE
    ENUM_FORCE_INT32_SIZE(LLKA_ResidueInStep)
} LLKA_ResidueInStep;

/*!
 * Container for atom names.
 *
 * Used as a return type from functions that get atom names for torsions or cross-residue metrics
 * used to define NtC. Length of the name is limited to 7 characters.
 */
typedef struct LLKA_AtomNameQuad {
    char a[8];    /*!< First atom name. */
    char b[8];    /*!< Second atom name. */
    char c[8];    /*!< Third atom name. Set only for torsions. */
    char d[8];    /*!< Fourth atom name. Set only for torsions. */
} LLKA_AtomNameQuad;
LLKA_IS_POD(LLKA_AtomNameQuad)

/*!
 * Idenfifies backbone atom in step
 */
typedef struct LLKA_BackboneAtom {
    LLKA_ResidueInStep residue;    /*!< First or second residue */
    const char *name;              /*!< Name of the atom (C5', O3', ...) */
} LLKA_BackboneAtom;
LLKA_IS_POD(LLKA_BackboneAtom)

/*!
 * Properties common for a NtC step
 */
typedef struct LLKA_StepInfo {
    int32_t firstSeqId;              /*!< First residue sequence id */
    int32_t firstSeqIdAuth;          /*!< First residue author id */
    int32_t secondSeqId;             /*!< Second residue sequence id */
    int32_t secondSeqIdAuth;         /*!< Second residue author id */
    LLKA_BaseKind firstBaseKind;     /*!< Kind of the second resdue (purine, pyrimidine) */
    LLKA_BaseKind secondBaseKind;    /*!< Kind of the second resdue (purine, pyrimidine) */
} LLKA_StepInfo;
LLKA_IS_POD(LLKA_StepInfo)

/*!
 * Torions and distances used to classify a step as a concrete NtC class.
 * Torsions are in range -PI; +PI, distances are in Angstroms.
 */
typedef struct LLKA_StepMetrics {
    double delta_1;
    double epsilon_1;
    double zeta_1;
    double alpha_2;
    double beta_2;
    double gamma_2;
    double delta_2;
    double chi_1;
    double chi_2;
    double CC;
    double NN;
    double mu;
} LLKA_StepMetrics;
LLKA_IS_POD(LLKA_StepMetrics)

#define LLKA_INVALID_BKBN_ATOM_INDEX (size_t)-1

LLKA_BEGIN_API_FUNCTIONS

/*!
 * Returns index of a backbone atom in the stucture that contains backbone atoms of a single step.
 *
 * @param[in] bkbnAtom Atom whose index to return.
 * @param[in] backbone The backnone structure.
 *
 * @return Index of the backbone atom in the array of atoms of the structure
 */
LLKA_API size_t LLKA_CC LLKA_backboneAtomIndex(const LLKA_BackboneAtom *bkbnAtom, const LLKA_Structure *backbone);

/*!
 * Returns human-readable CANA name for a corresponding LLKA_CANA item.
 *
 * @param[in] ntc The CANA.
 *
 * @return String with the name
 */
LLKA_API const char * LLKA_CC LLKA_CANAToName(LLKA_CANA cana);

/*!
 * Calculates step metrics for a NtC step.
 *
 * @param[in] stru Structure to calculate the metrics for. This structure must be a valid single NtC step.
 * @param[out] metrics Calculated metrics.
 *
 * @retval LLKA_OK Metrics was calculated successfully
 * @retval LLKA_E_MISSING_ATOMS Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_DATA Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_SIZES Passed structure is not a step.
 * @retval LLKA_E_MULTIPLE_ALT_IDS Passed structure is not a valid step.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure is not a valid step.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_calculateStepMetrics(const LLKA_Structure *stru, LLKA_StepMetrics *metrics);

/*!
 * Calculates the difference between NtC step metrics of a given real step against a reference step
 *
 * @param[in] stru Structure to calculate the metrics for. This structure must be a valid single NtC step.
 * @param[in] ntc NtC ID of the reference step.
 * @param[out] metrics Calculated metrics difference.
 *             Note that the fields of the \p LLKA_StepMetrics struct contain the calculated
 *             differences in [-PI; +PI] range.
 *
 * @retval LLKA_OK Metrics was calculated successfully
 * @retval LLKA_E_MISSING_ATOMS Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_DATA Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_SIZES Passed structure is not a step.
 * @retval LLKA_E_MULTIPLE_ALT_IDS Passed structure is not a valid step.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure is not a valid step or an invalid NtC was passed as argument.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_calculateStepMetricsDifferenceAgainstReference(const LLKA_Structure *stru, LLKA_NtC ntc, LLKA_StepMetrics *metrics);

/*!
 * Returns atoms that are measured to calculate a particular NtC step metric.
 *
 * @param[in] metric The metric to calculate.
 * @param[in] stru Structure to calculate the metric for.
 * @param[out] metricStru Set of atoms needed to calcualte the metric. Atoms from \p stru are deep-copied to \p metricStru.
 *                        Atoms of the output structure are always in the same order regardless of their order in the input structure.
 *                        Atoms are deep-copied to the output structure.
 *
 * @retval LLKA_OK Atoms were extracted successfully.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure does not contain the required atoms.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_crossResidueMetric(LLKA_CrossResidueMetric metric, const LLKA_Structure *stru, LLKA_Structure *metricStru);

/*!
 * Returns names of atoms that make up a given cross-residue metric.
 *
 * This function correctly accounts for non-standard residues.
 *
 * @param[in] firstBase Three-letter code of the first base.
 * @param[in] secondBase Three-letter code of the second base.
 * @param[in] metrc Cross-residue for which to get the atom names.
 * @param[out] quad Quad with the atom names. If the cross-residue metric is a distance, only fields \p a and \p b are set, \p c and \p d are not modified.
 *
 * @retval LLKA_OK Success
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_crossResidueMetricAtomsFromBases(const char *firstBase, const char *secondBase, LLKA_CrossResidueMetric metric, LLKA_AtomNameQuad *quad);

/*!
 * Returns names of atoms that make up a given cross-residue metric.
 *
 * This function correctly accounts for non-standard residues
 *
 * @param[in] stru Structure to extract the information from. The structure must be a valid step.
 * @param[in] metrc Cross-residue for which to get the atom names.
 * @param[out] quad Quad with the atom names. If the cross-residue metric is a distance, only fields \p a and \p b are set, \p c and \p d are not modified.
 *
 * @retval LLKA_OK Success
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_crossResidueMetricAtomsFromStructure(const LLKA_Structure *stru, LLKA_CrossResidueMetric metric, LLKA_AtomNameQuad *quad);

/*!
 * Returns human-readable name of a cross-residue metric.
 *
 * @param[in] metric The metric to name.
 * @param[in] greek If true, the name will use greek letters if applicable.
 *
 * @return String with the name
 */
LLKA_API const char * LLKA_CC LLKA_crossResidueMetricName(LLKA_CrossResidueMetric metric, LLKA_Bool greek);

/*!
 * Returns atoms that are measured to calculate a particular dinucleotide torsion.
 *
 * @param[in] metric The torsion to calculate.
 * @param[in] stru Structure to calculate the torsion for.
 * @param[out] metricStru Set of atoms needed to calcualte the torsion. Atoms from \p stru are deep-copied to \p torsionStru.
 *                        Atoms of the output structure are always in the same order regardless of their order in the input structure.
 *                        Atoms are deep-copied to the output structure.
 *
 * @retval LLKA_OK Atoms were extracted successfully.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure does not contain the required atoms.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_dinucleotideTorsion(LLKA_DinucleotideTorsion torsion, const LLKA_Structure *stru, LLKA_Structure *torsionStru);

/*!
 * Returns names of atoms that make up a given torsion.
 *
 * This function correctly accounts for non-standard residues.
 *
 * @param[in] firstBase Three-letter code of the first base.
 * @param[in] secondBase Three-letter code of the second base.
 * @param[in] torsion Torsion for which to get the atom names.
 * @param[out] quad Quad with the atom names.
 *
 * @retval LLKA_OK Success
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_dinucleotideTorsionAtomsFromBases(const char *firstBase, const char *secondBase, LLKA_DinucleotideTorsion torsion, LLKA_AtomNameQuad *quad);

/*!
 * Returns names of atoms that make up a given torsion.
 *
 * This function correctly accounts for non-standard residues
 *
 * @param[in] stru Structure to extract the information from. The structure must be a valid step.
 * @param[in] torsion Torsion for which to get the atom names.
 * @param[out] quad Quad with the atom names.
 *
 * @retval LLKA_OK Success
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_dinucleotideTorsionAtomsFromStructure(const LLKA_Structure *stru, LLKA_DinucleotideTorsion torsion, LLKA_AtomNameQuad *quad);


/*!
 * Returns human-readable name of a torsion metric.
 *
 * @param[in] metric The torsion to name.
 * @param[in] greek If true, the name will use greek letters if applicable.
 *
 * @return Name string
 */
LLKA_API const char * LLKA_CC LLKA_dinucleotideTorsionName(LLKA_DinucleotideTorsion torsion, LLKA_Bool greek);

/*!
 * Extracts backbone from a step.
 *
 * @param[in] stru Structure to extract the backbone from. The structure must be a valid single NtC step.
 * @param[out] backbone The extracted backbone. Atoms of the extracted backbone are always
 *                      in the same order regardless of their order in the input structure.
 *                      Atoms are deep-copied to the output structure.
 *
 * @retval LLKA_OK Backbone extracted successfully.
 * @retval LLKA_E_MISSING_ATOMS Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_DATA Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_SIZES Passed structure is not a step.
 * @retval LLKA_E_MULTIPLE_ALT_IDS Passed structure is not a valid step.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure is not a valid step.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_extractBackbone(const LLKA_Structure *stru, LLKA_Structure *backbone);

/*!
 * Extracts backbone from a step as a view.
 *
 * @param[in] stru Structure to extract the backbone from. The structure must be a valid single NtC step.
 * @param[out] backbone The extracted backbone. Atoms of the extracted backbone are always
 *                      in the same order regardless of their order in the input structure.
 *                      Atoms are deep-copied to the output structure.
 *
 * @retval LLKA_OK Backbone extracted successfully.
 * @retval LLKA_E_MISSING_ATOMS Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_DATA Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_SIZES Passed structure is not a step.
 * @retval LLKA_E_MULTIPLE_ALT_IDS Passed structure is not a valid step.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure is not a valid step.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_extractBackboneView(const LLKA_Structure *stru, LLKA_StructureView *backbone);

/*!
 * Extracts extended backbone from a step. Extended backbone contains the backbone atoms and
 * C1' atoms from both residues and the N1 or N9 (depending on the kind of base) atoms for both residues.
 *
 * @param[in] stru Structure to extract the backbone from. The structure must be a valid single NtC step.
 * @param[out] backbone The extracted extended backbone. Atoms of the extracted extended backbone are always
 *                      in the same order regardless of their order in the input structure.
 *                        Atoms are deep-copied to the output structure.
 *
 * @retval LLKA_OK Extended backbone extracted successfully.
 * @retval LLKA_E_MISSING_ATOMS Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_DATA Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_SIZES Passed structure is not a step.
 * @retval LLKA_E_MULTIPLE_ALT_IDS Passed structure is not a valid step.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure is not a valid step.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_extractExtendedBackbone(const LLKA_Structure *stru, LLKA_Structure *backbone);

/*!
 * Extracts extended backbone from a step as a view. Extended backbone contains the backbone atoms and
 * C1' atoms from both residues and the N1 or N9 (depending on the kind of base) atoms for both residues.
 *
 * @param[in] stru Structure to extract the backbone from. The structure must be a valid single NtC step.
 * @param[out] backbone The extracted extended backbone. Atoms of the extracted extended backbone are always
 *                      in the same order regardless of their order in the input structure.
 *                        Atoms are deep-copied to the output structure.
 *
 * @retval LLKA_OK Extended backbone extracted successfully.
 * @retval LLKA_E_MISSING_ATOMS Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_DATA Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_SIZES Passed structure is not a step.
 * @retval LLKA_E_MULTIPLE_ALT_IDS Passed structure is not a valid step.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure is not a valid step.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_extractExtendedBackboneView(const LLKA_Structure *stru, LLKA_StructureView *backbone);

/*!
 * Extracts the atoms of the structure that are necessary to calculate NtC step metrics.
 *
 * @param[in] stru Structure to extract the atoms from. The structure must be a valid single NtC step.
 * @param[out] metricsStru Structure containing the extracted atoms. Atoms in the output structure are always in the same order
 *                         regardless of their order in the input structure. Atoms are deep-copied to the output structure.
 *
 * @retval LLKA_OK Backbone extracted successfully.
 * @retval LLKA_E_MISSING_ATOMS Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_DATA Passed structure is not a valid step.
 * @retval LLKA_E_MISMATCHING_SIZES Passed structure is not a step.
 * @retval LLKA_E_MULTIPLE_ALT_IDS Passed structure is not a valid step.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure is not a valid step.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_extractMetricsStructure(const LLKA_Structure *stru, LLKA_Structure *metricsStru);

/*!
 * Converts a CANA name from a string into LLKA_CANA item.
 *
 * @param[in] name Name of the CANA class (such as AAA, BB2, ...)
 *
 * @return LLKA_CANA item. LLKA_INVALID_CANA is returned when the name is not recognized.
 */
LLKA_API LLKA_CANA LLKA_CC LLKA_nameToCANA(const char *name);

/*!
 * Converts a NtC name from a string into LLKA_NtC item.
 *
 * @param[in] name Name of the NtC (such as AA00, ZZ1S, ...)
 *
 * @return LLKA_NtC item. LLKA_INVALID_NTC is returned when the name is not recognized.
 */
LLKA_API LLKA_NtC LLKA_CC LLKA_nameToNtC(const char *name);

/*!
 * Returns structure of a reference dinucletoide from a given NtC
 *
 * @param[in] NtC NtC whose reference dinucleotide to return
 *
 * @return LLKA_Structure object
 */
LLKA_API LLKA_Structure LLKA_CC LLKA_NtCStructure(LLKA_NtC ntc);

/*!
 * Returns human-readable NtC name for a corresponding LLKA_NtC item.
 *
 * @param[in] ntc The NtC.
 *
 * @return String with the name
 */
LLKA_API const char * LLKA_CC LLKA_NtCToName(LLKA_NtC ntc);

/*!
 * Checks if a structure represents a single NtC step.
 *
 * @param[in] stru Structure to check.
 * @param[out] info Common information for the step.
 *
 * @retval LLKA_OK Structure is a single NtC step.
 * @retval LLKA_E_MISSING_ATOMS Passed structure is not a NtC step because it is missing some atoms.
 * @retval LLKA_E_MISMATCHING_DATA Passed structure is not a NtC step because it contains wrong number of residues.
 * @retval LLKA_E_MISMATCHING_SIZES Passed structure is not a step because at least one of the residues has wrong number of atoms.
 * @retval LLKA_E_MULTIPLE_ALT_IDS Passed structure is not a valid step because it contains atoms from more that one alternate positions.
 * @retval LLKA_E_INVALID_ARGUMENT Passed structure is not a valid step because the kind of the base is unknown.
 */
LLKA_API LLKA_RetCode LLKA_CC LLKA_structureIsStep(const LLKA_Structure *stru, LLKA_StepInfo *info);

LLKA_END_API_FUNCTIONS

#endif /* _LLKA_NTC_H */
