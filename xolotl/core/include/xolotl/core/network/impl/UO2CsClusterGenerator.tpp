#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
KOKKOS_INLINE_FUNCTION
bool
UO2CsClusterGenerator::refine(const Region& region, BoolArray& result) const
{
	result[0] = true;
	result[1] = true;
	result[2] = true;

	// I is never grouped
	if (region[Species::I].begin() > 0) {
		return true;
	}

	// V is never grouped
	if (region[Species::V].end() > 1 && region[Species::V].begin() < 3 &&
		region[Species::Cs].begin() == 0 && region[Species::I].begin() == 0) {
		return true;
	}

	// Cs is never grouped
	if (region[Species::Cs].end() > 1 && region[Species::Cs].begin() < 2 &&
		region[Species::V].begin() == 0 && region[Species::I].begin() == 0) {
		return true;
	}

	// CsV
	if (region[Species::Cs].begin() < _groupingMin &&
		region[Species::V].begin() < _groupingMin) {
		return true;
	}

	// Edges
	if (region[Species::Cs].end() > _maxCs + 1 &&
		region[Species::V].end() > _maxV + 1) {
		return true;
	}

	// Middle
	if (region[Species::Cs].length() <
		util::max((double)(_groupingWidthCs + 1),
			region[Species::Cs].begin() * 1.0e-1)) {
		result[0] = false;
	}
	if (region[Species::V].length() <
		util::max((double)(_groupingWidthV + 1),
			region[Species::V].begin() * 1.0e-1)) {
		result[1] = false;
	}

	// Edges
	if (region[Species::Cs].begin() == 0) {
		result[0] = true;
	}
	if (region[Species::V].begin() == 0) {
		result[1] = true;
	}
	if (region[Species::Cs].end() > _maxCs + 1) {
		result[0] = true;
	}
	if (region[Species::V].end() > _maxV + 1) {
		result[1] = true;
	}

	if (!result[0] && !result[1]) {
		return false;
	}

	return true;
}

KOKKOS_INLINE_FUNCTION
bool
UO2CsClusterGenerator::select(const Region& region) const
{
	// Remove 0
	if (region[Species::Cs].end() == 1 && region[Species::V].end() == 1 &&
		region[Species::I].end() == 1) {
		return false;
	}

	// Interstitials
	if (region[Species::I].begin() > 0 &&
		(region[Species::Cs].begin() > 0 || region[Species::V].begin() > 0)) {
		return false;
	}

	// Cesium
	if (region[Species::Cs].begin() > 1 && region[Species::V].end() == 1 &&
		region[Species::I].end() == 1) {
		return false;
	}
	if (region[Species::Cs].begin() > _maxCs) {
		return false;
	}

	// Vacancy
	if (region[Species::V].begin() > 2 && region[Species::Cs].end() == 1 &&
		region[Species::I].end() == 1) {
		return false;
	}
	if (region[Species::V].begin() > _maxV) {
		return false;
	}

	// Cs_1V
	if (region[Species::Cs].begin() == 1 && region[Species::V].end() > 9 &&
		region[Species::I].end() == 1) {
		return false;
	}

	return true;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
UO2CsClusterGenerator::getFormationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	// Not used?
	return 0.0;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
UO2CsClusterGenerator::getMigrationEnergy(
	const Cluster<PlsmContext>& cluster) const noexcept
{
	return util::infinity<double>;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
UO2CsClusterGenerator::getDiffusionFactor(
	const Cluster<PlsmContext>& cluster, double latticeParameter) const noexcept
{
	return 0.0;
}

template <typename PlsmContext>
KOKKOS_INLINE_FUNCTION
double
UO2CsClusterGenerator::getReactionRadius(const Cluster<PlsmContext>& cluster,
	double latticeParameter, double interstitialBias,
	double impurityRadius) const noexcept
{
	const auto& reg = cluster.getRegion();
	double radius = 0.0;
	double FourPi = 4.0 * ::xolotl::core::pi;
	double omega =
		0.25 * latticeParameter * latticeParameter * latticeParameter;
	if (reg.isSimplex()) {
		Composition comp(reg.getOrigin());
		if (comp.isOnAxis(Species::I)) {
			radius = latticeParameter / 2.0;
		}
		else if (comp.isOnAxis(Species::Cs)) {
			radius = impurityRadius;
		}
		else if (comp.isOnAxis(Species::V)) {
			radius = latticeParameter * sqrt(2.0) / 2.0;
		}
		else {
			radius = latticeParameter * sqrt(2.0) / 2.0 +
				cbrt((3.0 * omega * (double)comp[Species::V]) / FourPi) -
				cbrt((3.0 * omega) / FourPi);
		}
	}
	else {
		// Loop on the V range
		for (auto j : makeIntervalRange(reg[Species::V])) {
			radius += latticeParameter * sqrt(2.0) / 2.0 +
				cbrt((3.0 * omega * (double)j) / FourPi) -
				cbrt((3.0 * omega) / FourPi);
		}
		// Average the radius
		radius /= reg[Species::V].length();
	}

	return radius;
}
} // namespace network
} // namespace core
} // namespace xolotl
