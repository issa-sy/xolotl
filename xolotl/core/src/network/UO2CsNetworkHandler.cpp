#include <xolotl/core/network/UO2CsNetworkHandler.h>
#include <xolotl/core/network/UO2CsReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
auto uo2csNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		UO2CsNetworkHandler>({"UO2Cs"});
}

auto uo2csNetworkGenerator = [](const options::IOptions& options) {
	using NetworkType = core::network::UO2CsReactionNetwork;

	// Get the boundaries from the options
	NetworkType::AmountType maxV = options.getMaxV();
	NetworkType::AmountType maxI = options.getMaxI();
	NetworkType::AmountType maxCs = options.getMaxImpurity();
	NetworkType::AmountType groupingWidthCs = options.getGroupingWidthA();
	NetworkType::AmountType groupingWidthV = options.getGroupingWidthB();
	// Take care of the case with no grouping
	if (options.getGroupingMin() > maxV) {
		groupingWidthCs = maxCs + 1;
		groupingWidthV = maxV + 1;
	}
	else {
		// Adapt maxCs and maxV
		int i = 0;
		while (maxCs + 1 > pow(groupingWidthCs, i)) {
			++i;
		}
		maxCs = pow(groupingWidthCs, i) - 1;
		i = 0;
		while (maxV + 1 > pow(groupingWidthV, i)) {
			++i;
		}
		maxV = pow(groupingWidthV, i) - 1;
	}

	std::vector<NetworkType::AmountType> maxSpeciesAmounts = {
		maxCs, maxV, maxI};
	std::vector<NetworkType::SubdivisionRatio> subdivRatios = {
		{groupingWidthCs, groupingWidthV, maxI + 1}};
	auto network = std::make_shared<NetworkType>(
		maxSpeciesAmounts, subdivRatios, 1, options);

	return network;
};

UO2CsNetworkHandler::UO2CsNetworkHandler(const options::IOptions& options) :
	NetworkHandler(options, uo2csNetworkGenerator)
{
}
} // namespace network
} // namespace core
} // namespace xolotl
