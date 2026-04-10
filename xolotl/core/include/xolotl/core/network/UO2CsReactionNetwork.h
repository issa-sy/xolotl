#pragma once

#include <xolotl/core/network/UO2CsReaction.h>
#include <xolotl/core/network/UO2CsTraits.h>
#include <xolotl/core/network/ReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class UO2CsReactionGenerator;

class UO2CsClusterUpdater;
} // namespace detail

class UO2CsReactionNetwork : public ReactionNetwork<UO2CsReactionNetwork>
{
	friend class ReactionNetwork<UO2CsReactionNetwork>;

public:
	using Superclass = ReactionNetwork<UO2CsReactionNetwork>;
	using Subpaving = typename Superclass::Subpaving;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using IndexType = typename Superclass::IndexType;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;

	using Superclass::Superclass;

	IndexType
	checkLargestClusterId();

	void
	initializeExtraClusterData(const options::IOptions& options);

private:
	double
	checkLatticeParameter(double latticeParameter);

	double
	computeAtomicVolume(double latticeParameter)
	{
		// 4 atoms per cell
		return 0.25 * latticeParameter * latticeParameter * latticeParameter;
	}

	double
	checkImpurityRadius(double impurityRadius);

	detail::UO2CsReactionGenerator
	getReactionGenerator() const noexcept;

	void
	readClusters(const std::string filename = "reaction.dat");

	void
	readReactions(
		double temperature, const std::string filename = "reaction.dat");

	double
	getDegeneracy(IndexType v)
	{
		switch (v) {
		case 2:
			return 6.0;
		case 3:
			return 8.0;
		case 4:
			return 2.0;
		case 5:
			return 6.0;
		case 7:
			return 8.0;
		case 8:
			return 8.0;
		default:
			return 1.0;
		}
	}
};

namespace detail
{
class UO2CsReactionGenerator :
	public ReactionGenerator<UO2CsReactionNetwork, UO2CsReactionGenerator>
{
	friend class ReactionGeneratorBase<UO2CsReactionNetwork, UO2CsReactionGenerator>;

public:
	using NetworkType = UO2CsReactionNetwork;
	using Subpaving = typename NetworkType::Subpaving;
	using Superclass =
		ReactionGenerator<UO2CsReactionNetwork, UO2CsReactionGenerator>;

	using Superclass::Superclass;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	operator()(IndexType i, IndexType j, TTag tag) const;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addSinks(IndexType i, TTag tag) const;

private:
	ReactionCollection<NetworkType>
	getReactionCollection() const;
};

class UO2CsClusterUpdater
{
public:
	using NetworkType = UO2CsReactionNetwork;
	using ClusterData = typename NetworkType::ClusterData;
	using IndexType = typename NetworkType::IndexType;

	KOKKOS_INLIUO2Cs_FUNCTION
	void
	updateDiffusionCoefficient(const ClusterData& data, IndexType clusterId,
		IndexType gridIndex) const;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/UO2CsClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/UO2CsReactionNetwork.tpp>
#endif
