/* @HEADER@ */
/* @HEADER@ */

#include "SundanceSpatiallyConstantExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

SpatiallyConstantExpr::SpatiallyConstantExpr(const double& value)
	: EvaluatableExpr(), value_(value)
{}


