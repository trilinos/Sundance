/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUnaryExpr.hpp"
#include "SundanceExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;


UnaryExpr::UnaryExpr(const RefCountPtr<ScalarExpr>& arg)
	: ExprWithChildren(tuple(arg))
{}


