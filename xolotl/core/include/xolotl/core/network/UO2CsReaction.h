#pragma once

#include <xolotl/core/network/UO2CsTraits.h>
#include <xolotl/core/network/NucleationReaction.h>
#include <xolotl/core/network/ReSolutionReaction.h>
#include <xolotl/core/network/SinkReaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
class UO2CsReactionNetwork;

class UO2CsProductionReaction :
	public ProductionReaction<UO2CsReactionNetwork, UO2CsProductionReaction>
{
	friend class Reaction<UO2CsReactionNetwork, UO2CsProductionReaction>;

public:
	using Superclass =
		ProductionReaction<UO2CsReactionNetwork, UO2CsProductionReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	UO2CsProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1,
		IndexType cluster2 = invalidIndex, IndexType cluster3 = invalidIndex) :
		ProductionReaction(reactionData, clusterData, reactionId, cluster0,
			cluster1, cluster2, cluster3)
	{
		if (this->_clusterData->extraData.fileClusterMap.exists(
				_reactants[0]) and
			this->_clusterData->extraData.fileClusterMap.exists(
				_reactants[1])) {
			auto& clusterMap = this->_clusterData->extraData.fileClusterMap;
			this->_deltaG0 = reactionData.reactionEnergies(
				clusterMap.value_at(clusterMap.find(_reactants[0])),
				clusterMap.value_at(clusterMap.find(_reactants[1])));
		}
	}

	KOKKOS_INLINE_FUNCTION
	UO2CsProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		UO2CsProductionReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2,
			clusterSet.cluster3)
	{
	}

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy();

private:
	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex, double time = 0.0);

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);
};

class UO2CsDissociationReaction :
	public DissociationReaction<UO2CsReactionNetwork, UO2CsDissociationReaction>
{
	friend class Reaction<UO2CsReactionNetwork, UO2CsDissociationReaction>;

public:
	using Superclass =
		DissociationReaction<UO2CsReactionNetwork, UO2CsDissociationReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	UO2CsDissociationReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1, IndexType cluster2) :
		DissociationReaction(
			reactionData, clusterData, reactionId, cluster0, cluster1, cluster2)
	{
		if (this->_clusterData->extraData.fileClusterMap.exists(
				_products[0]) and
			this->_clusterData->extraData.fileClusterMap.exists(_products[1])) {
			auto& clusterMap = this->_clusterData->extraData.fileClusterMap;
			this->_deltaG0 = reactionData.reactionEnergies(
				clusterMap.value_at(clusterMap.find(_products[0])),
				clusterMap.value_at(clusterMap.find(_products[1])));
		}
	}

	KOKKOS_INLINE_FUNCTION
	UO2CsDissociationReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		UO2CsDissociationReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2)
	{
	}

	KOKKOS_INLINE_FUNCTION
	double
	getRateForProduction(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeBindingEnergy();

private:
	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex, double time = 0.0);

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);
};

class UO2CsReSolutionReaction :
	public ReSolutionReaction<UO2CsReactionNetwork, UO2CsReSolutionReaction>
{
public:
	using Superclass =
		ReSolutionReaction<UO2CsReactionNetwork, UO2CsReSolutionReaction>;

	using Superclass::Superclass;
};

class UO2CsNucleationReaction :
	public NucleationReaction<UO2CsReactionNetwork, UO2CsNucleationReaction>
{
public:
	using Superclass =
		NucleationReaction<UO2CsReactionNetwork, UO2CsNucleationReaction>;

	using Superclass::Superclass;
};

class UO2CsSinkReaction : public SinkReaction<UO2CsReactionNetwork, UO2CsSinkReaction>
{
	friend class Reaction<UO2CsReactionNetwork, UO2CsSinkReaction>;

public:
	using Superclass = SinkReaction<UO2CsReactionNetwork, UO2CsSinkReaction>;

	using Superclass::Superclass;

	KOKKOS_INLINE_FUNCTION
	UO2CsSinkReaction(ReactionDataRef reactionData, const ClusterData& clusterData,
		IndexType reactionId, IndexType cluster0) :
		SinkReaction(reactionData, clusterData, reactionId, cluster0)
	{
		if (this->_clusterData->extraData.fileClusterMap.exists(_reactant)) {
			auto& clusterMap = this->_clusterData->extraData.fileClusterMap;
			this->_deltaG0 = reactionData.reactionEnergies(
				clusterMap.value_at(clusterMap.find(_reactant)),
				this->_clusterData->extraData.fileClusterNumber);
		}
	}

	KOKKOS_INLINE_FUNCTION
	UO2CsSinkReaction(ReactionDataRef reactionData, const ClusterData& clusterData,
		IndexType reactionId, const detail::ClusterSet& clusterSet) :
		UO2CsSinkReaction(
			reactionData, clusterData, reactionId, clusterSet.cluster0)
	{
	}

	KOKKOS_INLINE_FUNCTION
	double
	getSinkBias();

	KOKKOS_INLINE_FUNCTION
	double
	getSinkStrength();

private:
	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);
};
} // namespace network
} // namespace core
} // namespace xolotl
