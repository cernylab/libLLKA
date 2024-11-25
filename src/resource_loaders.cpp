/* vim: set sw=4 ts=4 sts=4 : */

#include <llka_resource_loaders.h>

#include "util/elementaries.h"
#include "util/csv.hpp"

namespace LLKAInternal {

inline constinit char CSV_DELIM{';'};
inline constinit char CSV_QUOTE{'"'};

namespace Csv {

template <>
struct Convert<LLKA_CANA> {
	template <typename Tag>
	static auto call(const std::string &s, const Tag &tag)
	{
		auto cana = LLKA_nameToCANA(s.c_str());
		if (cana == LLKA_INVALID_CANA)
			throw TokenParseError{tagStr(tag)};
		return cana;
	}
};

template <>
struct Convert<LLKA_NtC> {
	template <typename Tag>
	static auto call(const std::string &s, const Tag &tag)
	{
		auto ntc = LLKA_nameToNtC(s.c_str());
		if (ntc == LLKA_INVALID_NTC)
			throw TokenParseError{tagStr(tag)};
		return ntc;
	}
};

template <>
struct Convert<LLKA_SugarPucker> {
	template <typename Tag>
	static auto call(const std::string &s, const Tag &tag)
	{
		auto pucker = LLKA_nameToSugarPucker(s.c_str());
		if (pucker == LLKA_INVALID_SUGAR_PUCKER)
			throw TokenParseError{tagStr(tag)};
		return pucker;
	}
};

} // namespace Csv

namespace Setters {

template <LLKA_ClassificationMetric LLKA_ClassificationCluster::* MetricPtr, double LLKA_ClassificationMetric:: *ValuePtr>
struct ClassificationMetric {
	using Type = double;
	static auto set(LLKA_ClassificationCluster &cluster, double v) { cluster.*MetricPtr.*ValuePtr = v; }
};

template <typename Container, typename T, T Container::* ClsPtr>
struct Direct {
	using Type = T;
	static auto set(Container &c, T &&v) { c.*ClsPtr = std::move(v); }
};
template <typename Container, const char * Container::* ClsPtr>
struct Direct<Container, const char *, ClsPtr> {
	using Type = const char *;
	static auto set(Container &c, const char *cc) { c.*ClsPtr = cc; }
};

template <typename Container, LLKA_StepMetrics Container::* SMClsPtr, double LLKA_StepMetrics::* MetricClsPtr>
struct Metric {
	using Type = double;
	static auto set(Container &c, double torsion) { c.*SMClsPtr.*MetricClsPtr = torsion; }
};

template <typename Container, LLKA_NuAngles Container::* NAClsPtr, double LLKA_NuAngles::* NuAngleClsPtr>
struct NuAngle {
	using Type = double;
	static auto set(Container &c, double angle) { c.*NAClsPtr.*NuAngleClsPtr = angle; }
};

template <typename Container, LLKA_NuAnglesMetrics Container::* CNAClsPtr, LLKA_ClassificationMetric LLKA_NuAnglesMetrics::* NuAngleClsPtr, double LLKA_ClassificationMetric::* MClsPtr>
struct NuAngleMetrics {
	using Type = double;
	static auto set(Container &c, double angle) { c.*CNAClsPtr.*NuAngleClsPtr.*MClsPtr = angle; }
};

} // Setters

// CSV schemas

using ClusterNuAnglesSchema = Csv::Schema<
	Csv::SchemaKind::Named,
	LLKA_ClusterNuAngles,
	Csv::NamedField<"clusterNumber", Setters::Direct<LLKA_ClusterNuAngles, int32_t, &LLKA_ClusterNuAngles::clusterNumber>>,
	// First nucleotide
	Csv::NamedField<"nu0_1_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_0, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu0_1_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_0, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"nu1_1_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_1, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu1_1_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_1, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"nu2_1_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_2, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu2_1_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_2, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"nu3_1_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_3, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu3_1_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_3, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu4_1_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_4, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu4_1_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAnglesMetrics::nu_4, &LLKA_ClassificationMetric::deviation>>,
	// Second nucletotide
	Csv::NamedField<"nu0_2_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_0, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu0_2_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_0, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"nu1_2_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_1, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu1_2_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_1, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"nu2_2_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_2, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu2_2_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_2, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"nu3_2_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_3, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu3_2_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_3, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"nu4_2_meanValue", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_4, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"nu4_2_deviation", Setters::NuAngleMetrics<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAnglesMetrics::nu_4, &LLKA_ClassificationMetric::deviation>>
>;

using ConfalPercentilesSchema = Csv::Schema<
	Csv::SchemaKind::Numbered,
	LLKA_ConfalPercentile,
	Csv::NumberedField<0, Setters::Direct<LLKA_ConfalPercentile, double, &LLKA_ConfalPercentile::value>>
>;

/*
Example of a schema with numbered fields

using ClusterNuAnglesSchema = Csv::Schema<
	Csv::SchemaKind::Numbered,
	LLKA_ClusterNuAngles,
	Csv::NumberedField<0, Setters::Direct<LLKA_ClusterNuAngles, int32_t, &LLKA_ClusterNuAngles::clusterNumber>>,
	Csv::NumberedField<4, Setters::NuAngle<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::firstNucleotide,  &LLKA_NuAngles::nu_3>>,
	Csv::NumberedField<10, Setters::NuAngle<LLKA_ClusterNuAngles, &LLKA_ClusterNuAngles::secondNucleotide, &LLKA_NuAngles::nu_4>>
>;
*/

using ClustersSchema = Csv::Schema<
	Csv::SchemaKind::Named,
	LLKA_ClassificationCluster,
	Csv::NamedField<"clusterNumber",          Setters::Direct<LLKA_ClassificationCluster, int32_t, &LLKA_ClassificationCluster::number>>,
	Csv::NamedField<"ntc",                    Setters::Direct<LLKA_ClassificationCluster, LLKA_NtC, &LLKA_ClassificationCluster::NtC>>,
	Csv::NamedField<"cana",                   Setters::Direct<LLKA_ClassificationCluster, LLKA_CANA, &LLKA_ClassificationCluster::CANA>>,
	Csv::NamedField<"delta_1_deviation",      Setters::ClassificationMetric<&LLKA_ClassificationCluster::delta_1, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"delta_1_meanValue",      Setters::ClassificationMetric<&LLKA_ClassificationCluster::delta_1, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"epsilon_1_deviation",    Setters::ClassificationMetric<&LLKA_ClassificationCluster::epsilon_1, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"epsilon_1_meanValue",    Setters::ClassificationMetric<&LLKA_ClassificationCluster::epsilon_1, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"zeta_1_deviation",       Setters::ClassificationMetric<&LLKA_ClassificationCluster::zeta_1, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"zeta_1_meanValue",       Setters::ClassificationMetric<&LLKA_ClassificationCluster::zeta_1, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"alpha_2_deviation",      Setters::ClassificationMetric<&LLKA_ClassificationCluster::alpha_2, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"alpha_2_meanValue",      Setters::ClassificationMetric<&LLKA_ClassificationCluster::alpha_2, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"beta_2_deviation",       Setters::ClassificationMetric<&LLKA_ClassificationCluster::beta_2, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"beta_2_meanValue",       Setters::ClassificationMetric<&LLKA_ClassificationCluster::beta_2, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"gamma_2_deviation",      Setters::ClassificationMetric<&LLKA_ClassificationCluster::gamma_2, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"gamma_2_meanValue",      Setters::ClassificationMetric<&LLKA_ClassificationCluster::gamma_2, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"delta_2_deviation",      Setters::ClassificationMetric<&LLKA_ClassificationCluster::delta_2, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"delta_2_meanValue",      Setters::ClassificationMetric<&LLKA_ClassificationCluster::delta_2, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"chi_1_deviation",        Setters::ClassificationMetric<&LLKA_ClassificationCluster::chi_1, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"chi_1_meanValue",        Setters::ClassificationMetric<&LLKA_ClassificationCluster::chi_1, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"chi_2_deviation",        Setters::ClassificationMetric<&LLKA_ClassificationCluster::chi_2, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"chi_2_meanValue",        Setters::ClassificationMetric<&LLKA_ClassificationCluster::chi_2, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"CC_deviation",           Setters::ClassificationMetric<&LLKA_ClassificationCluster::CC, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"CC_meanValue",           Setters::ClassificationMetric<&LLKA_ClassificationCluster::CC, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"NN_deviation",           Setters::ClassificationMetric<&LLKA_ClassificationCluster::NN, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"NN_meanValue",           Setters::ClassificationMetric<&LLKA_ClassificationCluster::NN, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"mu_deviation",           Setters::ClassificationMetric<&LLKA_ClassificationCluster::mu, &LLKA_ClassificationMetric::deviation>>,
	Csv::NamedField<"mu_meanValue",           Setters::ClassificationMetric<&LLKA_ClassificationCluster::mu, &LLKA_ClassificationMetric::meanValue>>,
	Csv::NamedField<"ribosePseudorotation_1", Setters::Direct<LLKA_ClassificationCluster, double, &LLKA_ClassificationCluster::ribosePseudorotation_1>>,
	Csv::NamedField<"ribosePseudorotation_2", Setters::Direct<LLKA_ClassificationCluster, double, &LLKA_ClassificationCluster::ribosePseudorotation_2>>
>;


using ConfalsSchema = Csv::Schema<
	Csv::SchemaKind::Named,
	LLKA_Confal,
	Csv::NamedField<"clusterNumber", Setters::Direct<LLKA_Confal, int32_t, &LLKA_Confal::clusterNumber>>,
	Csv::NamedField<"delta_1",       Setters::Direct<LLKA_Confal, double, &LLKA_Confal::delta_1>>,
	Csv::NamedField<"epsilon_1",     Setters::Direct<LLKA_Confal, double, &LLKA_Confal::epsilon_1>>,
	Csv::NamedField<"zeta_1",        Setters::Direct<LLKA_Confal, double, &LLKA_Confal::zeta_1>>,
	Csv::NamedField<"alpha_2",       Setters::Direct<LLKA_Confal, double, &LLKA_Confal::alpha_2>>,
	Csv::NamedField<"beta_2",        Setters::Direct<LLKA_Confal, double, &LLKA_Confal::beta_2>>,
	Csv::NamedField<"gamma_2",       Setters::Direct<LLKA_Confal, double, &LLKA_Confal::gamma_2>>,
	Csv::NamedField<"delta_2",       Setters::Direct<LLKA_Confal, double, &LLKA_Confal::delta_2>>,
	Csv::NamedField<"chi_1",         Setters::Direct<LLKA_Confal, double, &LLKA_Confal::chi_1>>,
	Csv::NamedField<"chi_2",         Setters::Direct<LLKA_Confal, double, &LLKA_Confal::chi_2>>,
	Csv::NamedField<"CC",            Setters::Direct<LLKA_Confal, double, &LLKA_Confal::CC>>,
	Csv::NamedField<"NN",            Setters::Direct<LLKA_Confal, double, &LLKA_Confal::NN>>,
	Csv::NamedField<"mu",            Setters::Direct<LLKA_Confal, double, &LLKA_Confal::mu>>,
	Csv::NamedField<"nu0_1",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusFirst, &LLKA_NuAngles::nu_0>>,
	Csv::NamedField<"nu1_1",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusFirst, &LLKA_NuAngles::nu_1>>,
	Csv::NamedField<"nu2_1",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusFirst, &LLKA_NuAngles::nu_2>>,
	Csv::NamedField<"nu3_1",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusFirst, &LLKA_NuAngles::nu_3>>,
	Csv::NamedField<"nu4_1",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusFirst, &LLKA_NuAngles::nu_4>>,
	Csv::NamedField<"nu0_2",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusSecond, &LLKA_NuAngles::nu_0>>,
	Csv::NamedField<"nu1_2",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusSecond, &LLKA_NuAngles::nu_1>>,
	Csv::NamedField<"nu2_2",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusSecond, &LLKA_NuAngles::nu_2>>,
	Csv::NamedField<"nu3_2",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusSecond, &LLKA_NuAngles::nu_3>>,
	Csv::NamedField<"nu4_2",         Setters::NuAngle<LLKA_Confal, &LLKA_Confal::nusSecond, &LLKA_NuAngles::nu_4>>
>;

template <typename T, T LLKA_GoldenStep::* ClsPtr>
using GSSDirect = Setters::Direct<LLKA_GoldenStep, T, ClsPtr>;
template <double LLKA_StepMetrics::* TorsionClsPtr>
using GSSMetric = Setters::Metric<LLKA_GoldenStep, &LLKA_GoldenStep::metrics, TorsionClsPtr>;
template <double LLKA_NuAngles::* NAClsPtr>
using GSSNuAngle1 = Setters::NuAngle<LLKA_GoldenStep, &LLKA_GoldenStep::nuAngles_1, NAClsPtr>;
template <double LLKA_NuAngles::* NAClsPtr>
using GSSNuAngle2 = Setters::NuAngle<LLKA_GoldenStep, &LLKA_GoldenStep::nuAngles_2, NAClsPtr>;

using GoldenStepsSchema = Csv::Schema<
	Csv::SchemaKind::Named,
	LLKA_GoldenStep,
	Csv::NamedField<"clusterNumber", GSSDirect<int32_t, &LLKA_GoldenStep::clusterNumber>>,
	Csv::NamedField<"delta_1",       GSSMetric<&LLKA_StepMetrics::delta_1>>,
	Csv::NamedField<"epsilon_1",     GSSMetric<&LLKA_StepMetrics::epsilon_1>>,
	Csv::NamedField<"zeta_1",        GSSMetric<&LLKA_StepMetrics::zeta_1>>,
	Csv::NamedField<"alpha_2",       GSSMetric<&LLKA_StepMetrics::alpha_2>>,
	Csv::NamedField<"beta_2",        GSSMetric<&LLKA_StepMetrics::beta_2>>,
	Csv::NamedField<"gamma_2",       GSSMetric<&LLKA_StepMetrics::gamma_2>>,
	Csv::NamedField<"delta_2",       GSSMetric<&LLKA_StepMetrics::delta_2>>,
	Csv::NamedField<"chi_1",         GSSMetric<&LLKA_StepMetrics::chi_1>>,
	Csv::NamedField<"chi_2",         GSSMetric<&LLKA_StepMetrics::chi_2>>,
	Csv::NamedField<"CC",            GSSMetric<&LLKA_StepMetrics::CC>>,
	Csv::NamedField<"NN",            GSSMetric<&LLKA_StepMetrics::NN>>,
	Csv::NamedField<"mu",            GSSMetric<&LLKA_StepMetrics::mu>>,
	Csv::NamedField<"name",          GSSDirect<const char *, &LLKA_GoldenStep::name>>,
	Csv::NamedField<"pucker_1",      GSSDirect<LLKA_SugarPucker, &LLKA_GoldenStep::pucker_1>>,
	Csv::NamedField<"pucker_2",      GSSDirect<LLKA_SugarPucker, &LLKA_GoldenStep::pucker_2>>,
	Csv::NamedField<"nu0_1",         GSSNuAngle1<&LLKA_NuAngles::nu_0>>,
	Csv::NamedField<"nu1_1",         GSSNuAngle1<&LLKA_NuAngles::nu_1>>,
	Csv::NamedField<"nu2_1",         GSSNuAngle1<&LLKA_NuAngles::nu_2>>,
	Csv::NamedField<"nu3_1",         GSSNuAngle1<&LLKA_NuAngles::nu_3>>,
	Csv::NamedField<"nu4_1",         GSSNuAngle1<&LLKA_NuAngles::nu_4>>,
	Csv::NamedField<"nu0_2",         GSSNuAngle2<&LLKA_NuAngles::nu_0>>,
	Csv::NamedField<"nu1_2",         GSSNuAngle2<&LLKA_NuAngles::nu_1>>,
	Csv::NamedField<"nu2_2",         GSSNuAngle2<&LLKA_NuAngles::nu_2>>,
	Csv::NamedField<"nu3_2",         GSSNuAngle2<&LLKA_NuAngles::nu_3>>,
	Csv::NamedField<"nu4_2",         GSSNuAngle2<&LLKA_NuAngles::nu_4>>
>;

template <typename Schema, typename Input>
static
auto read(const Input &input, typename Schema::Line * &data, size_t &count)
{
	try {
		auto [ lines, warnings ] = LLKAInternal::Csv::parse<Schema>(input, LLKAInternal::CSV_DELIM, LLKAInternal::CSV_QUOTE);

		if (!warnings.empty())
			return LLKA_E_BAD_DATA;

		if (lines.empty())
			return LLKA_E_BAD_DATA;

		data = new typename Schema::Line[lines.size()];
		count = lines.size();

		for (size_t idx = 0; idx < lines.size(); idx++)
			data[idx] = lines[idx];

		return LLKA_OK;
	} catch (const LLKA_RetCode tRet) {
		return tRet;
	}
}

} // namespace LLKA

void LLKA_CC LLKA_destroyResource(LLKA_Resource *resource)
{
	if (resource->count == 0)
		return;

	switch (resource->type) {
	case LLKA_RES_AVERAGE_NU_ANGLES:
		delete [] resource->data.clusterNuAngles;
		break;
	case LLKA_RES_CLUSTERS:
		delete [] resource->data.clusters;
		break;
	case LLKA_RES_CONFALS:
		delete [] resource->data.confals;
		break;
	case LLKA_RES_GOLDEN_STEPS:
		for (size_t idx = 0; idx < resource->count; idx++)
			LLKAInternal::destroyString(resource->data.goldenSteps[idx].name);
		delete [] resource->data.goldenSteps;
		break;
	case LLKA_RES_CONFAL_PERCENTILES:
		delete [] resource->data.confalPercentiles;
		break;
	}
}

struct LoadClusterNuAngles {
	using Schema = LLKAInternal::ClusterNuAnglesSchema;
	static constexpr LLKA_ClusterNuAngles * LLKA_UResource::* Data = &LLKA_UResource::clusterNuAngles;
};
struct LoadClusters {
	using Schema = LLKAInternal::ClustersSchema;
	static constexpr LLKA_ClassificationCluster * LLKA_UResource::* Data = &LLKA_UResource::clusters;
};
struct LoadConfals {
	using Schema = LLKAInternal::ConfalsSchema;
	static constexpr LLKA_Confal * LLKA_UResource::* Data = &LLKA_UResource::confals;
};
struct LoadConfalPercentiles {
	using Schema = LLKAInternal::ConfalPercentilesSchema;
	static constexpr LLKA_ConfalPercentile * LLKA_UResource::* Data = &LLKA_UResource::confalPercentiles;
};
struct LoadGoldenSteps {
	using Schema = LLKAInternal::GoldenStepsSchema;
	static constexpr LLKA_GoldenStep * LLKA_UResource::* Data = &LLKA_UResource::goldenSteps;
};

enum class ResourceAs {
	FILE,
	TEXT
};

template <ResourceAs As, typename Loader, typename ParamStr>
static
LLKA_RetCode loadResourceByType(const ParamStr param, LLKA_Resource *resource)
{
	if constexpr (As == ResourceAs::FILE)
		return LLKAInternal::read<typename Loader::Schema, std::filesystem::path>(param, resource->data.*Loader::Data, resource->count);
	else if constexpr (As == ResourceAs::TEXT)
		return LLKAInternal::read<typename Loader::Schema, std::string>(param, resource->data.*Loader::Data, resource->count);
	else
		return LLKA_E_INVALID_ARGUMENT;
}

template <ResourceAs As, typename ParamStr>
static
LLKA_RetCode loadResource(const ParamStr param, LLKA_Resource *resource)
{
	if (resource->type == LLKA_RES_AVERAGE_NU_ANGLES)
		return loadResourceByType<As, LoadClusterNuAngles, ParamStr>(param, resource);
	else if (resource->type == LLKA_RES_CLUSTERS) {
		auto tRet = loadResourceByType<As, LoadClusters, ParamStr>(param, resource);
		if (tRet != LLKA_OK)
			return tRet;

		// minValue and maxValue will be initialized later, set then to NaNs to indicate uninitialized data
		for (size_t idx = 0; idx < resource->count; idx++) {
			auto &cluster = resource->data.clusters[idx];
			for (
				auto &clsPtr :
				{ &LLKA_ClassificationCluster::delta_1, &LLKA_ClassificationCluster::epsilon_1, &LLKA_ClassificationCluster::zeta_1,
				  &LLKA_ClassificationCluster::alpha_2, &LLKA_ClassificationCluster::beta_2, &LLKA_ClassificationCluster::gamma_2,
				  &LLKA_ClassificationCluster::delta_2, &LLKA_ClassificationCluster::chi_1, &LLKA_ClassificationCluster::chi_2,
				  &LLKA_ClassificationCluster::CC, &LLKA_ClassificationCluster::NN, &LLKA_ClassificationCluster::mu }
			) {
				auto &min = (cluster.*clsPtr).minValue;
				auto &max = (cluster.*clsPtr).maxValue;

				min = LLKAInternal::NaN<std::remove_reference_t<decltype(min)>>();
				max = LLKAInternal::NaN<std::remove_reference_t<decltype(max)>>();
			}
		}

		return LLKA_OK;
	} else if (resource->type == LLKA_RES_CONFALS)
		return loadResourceByType<As, LoadConfals, ParamStr>(param, resource);
	else if (resource->type == LLKA_RES_GOLDEN_STEPS)
		return loadResourceByType<As, LoadGoldenSteps, ParamStr>(param, resource);
	else if (resource->type == LLKA_RES_CONFAL_PERCENTILES)
		return loadResourceByType<As, LoadConfalPercentiles, ParamStr>(param, resource);

	return LLKA_E_INVALID_ARGUMENT;
}

#ifndef LLKA_FILESYSTEM_ACCESS_DISABLED
LLKA_RetCode LLKA_CC LLKA_loadResourceFile(const LLKA_PathChar *path, LLKA_Resource *resource)
{
	return loadResource<ResourceAs::FILE>(path, resource);
}
#endif /* LLKA_FILESYSTEM_ACCESS_DISABLED */

LLKA_RetCode LLKA_CC LLKA_loadResourceText(const char *text, LLKA_Resource *resource)
{
	return loadResource<ResourceAs::TEXT>(text, resource);
}
