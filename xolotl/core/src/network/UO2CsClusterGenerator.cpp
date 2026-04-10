#include <xolotl/core/network/UO2CsReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace network
{
UO2CsClusterGenerator::UO2CsClusterGenerator(const options::IOptions& options) :
	_maxCs(options.getMaxImpurity()),
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthCs(options.getGroupingWidthA()),
	_groupingWidthV(options.getGroupingWidthB()),
	_csDiffusivity(options.getCesiumDiffusivity()),
	_csDiffusive(_csDiffusivity > 0.0),
	_density(options.getDensity())
{
}

UO2CsClusterGenerator::UO2CsClusterGenerator(
	const options::IOptions& options, std::size_t refineDepth) :
	Superclass(refineDepth),
	_maxCs(options.getMaxImpurity()),
	_maxV(options.getMaxV()),
	_groupingMin(options.getGroupingMin()),
	_groupingWidthCs(options.getGroupingWidthA()),
	_groupingWidthV(options.getGroupingWidthB()),
	_csDiffusivity(options.getCesiumDiffusivity()),
	_csDiffusive(_csDiffusivity > 0.0),
	_density(options.getDensity())
{
}
} // namespace network
} // namespace core
} // namespace xolotl
