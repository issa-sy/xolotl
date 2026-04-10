#pragma once

// ============================================================================
// À placer dans : xolotl/core/include/xolotl/core/material/UO2CsMaterialHandler.h
// ============================================================================

#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/core/flux/UO2CsFitFluxHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class UO2CsMaterialHandler : public MaterialHandler
{
public:
    UO2CsMaterialHandler(const options::IOptions& options) :
        MaterialHandler(options,
            MaterialSubHandlerGenerator<core::flux::UO2CsFitFluxHandler>{})
    {
    }
};
} // namespace material
} // namespace core
} // namespace xolotl
