#include <xolotl/core/material/UO2CsMaterialHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
namespace detail
{
auto UO2CsMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		UO2CsMaterialHandler>({"UO2Cs"});
}
} // namespace material
} // namespace core
} // namespace xolotl