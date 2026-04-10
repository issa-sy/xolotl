#pragma once

#include <xolotl/core/network/NetworkHandler.h>
#include <xolotl/factory/network/NetworkHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace network
{
class UO2CsNetworkHandler : public NetworkHandler
{
public:
	UO2CsNetworkHandler(const options::IOptions& options);
};
} // namespace network
} // namespace core
} // namespace xolotl
