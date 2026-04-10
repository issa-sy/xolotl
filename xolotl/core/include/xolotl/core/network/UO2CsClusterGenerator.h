#pragma once

#include <plsm/refine/Detector.h>

namespace xolotl
{
namespace core
{
namespace network
{
class UO2CsClusterGenerator : public plsm::refine::Detector<UO2CsClusterGenerator>
{
public:
	using Superclass = plsm::refine::Detector<UO2CsClusterGenerator>;
	using NetworkType = UO2CsReactionNetwork;
	using Species = UO2CsSpecies;

	template <typename PlsmContext>
	using Cluster = Cluster<NetworkType, PlsmContext>;

	using Region = typename NetworkType::Region;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;
	using BoolArray = plsm::refine::BoolVec<Region>;

	UO2CsClusterGenerator(const options::IOptions& options);

	UO2CsClusterGenerator(
		const options::IOptions& options, std::size_t refineDepth);

	KOKKOS_INLINE_FUNCTION
	bool
	refine(const Region& region, BoolArray& result) const;

	KOKKOS_INLINE_FUNCTION
	bool
	select(const Region& region) const;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getFormationEnergy(const Cluster<PlsmContext>& cluster) const noexcept;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getMigrationEnergy(const Cluster<PlsmContext>& cluster) const noexcept;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getDissociationEnergy(const Cluster<PlsmContext>& cluster) const noexcept;
	double
	getDiffusionFactor(const Cluster<PlsmContext>& cluster,
		double latticeParameter) const noexcept;

	template <typename PlsmContext>
	KOKKOS_INLINE_FUNCTION
	double
	getReactionRadius(const Cluster<PlsmContext>& cluster,
		double latticeParameter, double interstitialBias,
		double impurityRadius) const noexcept;

private:
	AmountType _maxCs;
	AmountType _maxV{0};
	AmountType _groupingMin;
	AmountType _groupingWidthCs;
	AmountType _groupingWidthV;
	double _csDiffusivity;
	bool _csDiffusive;
	double _density;
};
} // namespace network
} // namespace core
} // namespace xolotl
