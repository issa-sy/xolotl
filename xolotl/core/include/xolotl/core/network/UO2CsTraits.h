#pragma once

#include <tuple>

#include <Kokkos_UnorderedMap.hpp>

#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/detail/ClusterData.h>

namespace xolotl
{
namespace core
{
namespace network
{
class UO2CsProductionReaction;
class UO2CsDissociationReaction;
class UO2CsReSolutionReaction;
class UO2CsNucleationReaction;
class UO2CsSinkReaction;
class UO2CsReactionNetwork;
class UO2CsClusterGenerator;
namespace detail
{
class UO2CsClusterUpdater;
}

enum class UO2CsSpecies
{
	Cs,
	V,
	I
};

inline const std::string&
toLabelString(UO2CsSpecies species)
{
	static const std::string labelArray[] = {"Cs", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(UO2CsSpecies species)
{
	static const std::string nameArray[] = {"Cesium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<UO2CsSpecies> : std::integral_constant<std::size_t, 3>
{
};

template <>
struct NumberOfInterstitialSpecies<UO2CsSpecies> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<UO2CsSpecies> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<UO2CsSpecies, 3>
{
	using Sequence = EnumSequence<UO2CsSpecies, 3>;
	static constexpr auto first = Sequence(UO2CsSpecies::Cs);
	static constexpr auto last = Sequence(UO2CsSpecies::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<UO2CsSpecies>
	mapToMomentId(EnumSequence<UO2CsSpecies, 3> value)
	{
		if (value == UO2CsSpecies::I)
			return 1;
		return value();
	}
};

template <>
struct ReactionNetworkTraits<UO2CsReactionNetwork>
{
	using Species = UO2CsSpecies;

	static constexpr std::size_t numSpecies = 3;

	// using ReactionType = UO2CsReaction;
	using ProductionReactionType = UO2CsProductionReaction;
	using DissociationReactionType = UO2CsDissociationReaction;
	using ReSolutionReactionType = UO2CsReSolutionReaction;
	using NucleationReactionType = UO2CsNucleationReaction;
	using SinkReactionType = UO2CsSinkReaction;

	using ReactionTypeList =
		std::tuple<ProductionReactionType, DissociationReactionType,
			ReSolutionReactionType, NucleationReactionType, SinkReactionType>;

	using ClusterGenerator = UO2CsClusterGenerator;
	using ClusterUpdater = detail::UO2CsClusterUpdater;
};

namespace detail
{
template <typename PlsmContext>
struct ClusterDataExtra<UO2CsReactionNetwork, PlsmContext>
{
	using NetworkType = UO2CsReactionNetwork;

	template <typename TData>
	using View = ViewType<TData, PlsmContext>;
	using MapType = Kokkos::UnorderedMap<int, int, PlsmContext>;

	using IndexType = detail::ReactionNetworkIndexType;

	ClusterDataExtra() = default;

	template <typename PC>
	KOKKOS_INLINE_FUNCTION
	ClusterDataExtra(const ClusterDataExtra<NetworkType, PC>& data) :
		constantRates(data.constantRates),
		fileClusterMap(data.fileClusterMap)
	{
	}

	template <typename PC>
	void
	deepCopy(const ClusterDataExtra<NetworkType, PC>& data)
	{
		if (!data.constantRates.is_allocated()) {
			return;
		}

		if (!constantRates.is_allocated()) {
			constantRates = create_mirror_view(data.constantRates);
		}

		if (!fileClusterMap.is_allocated()) {
			fileClusterMap = MapType(data.fileClusterMap.capacity());
		}

		deep_copy(constantRates, data.constantRates);
		deep_copy(fileClusterMap, data.fileClusterMap);
		fileClusterNumber = data.fileClusterNumber;
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret = 0;

		ret += constantRates.required_allocation_size(constantRates.extent(0),
			constantRates.extent(1), constantRates.extent(2));
		ret += sizeof(fileClusterMap);
		ret += sizeof(fileClusterNumber);

		return ret;
	}

	void
	initialize(IndexType size)
	{
		constantRates = View<double***>("Constant Rates", size, size + 1, 2);
		fileClusterMap = MapType(size);
	}

	View<double***> constantRates;
	MapType fileClusterMap;
	IdType fileClusterNumber{0};
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
