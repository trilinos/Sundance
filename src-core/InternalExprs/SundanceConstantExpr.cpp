/* @HEADER@ */
/* @HEADER@ */

#include "SundanceConstantExpr.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

ConstantExpr::ConstantExpr(const double& value)
	: SpatiallyConstantExpr(value)
{}



ostream& ConstantExpr::toText(ostream& os, bool /* paren */) const 
{
	os << value();
	return os;
}

ostream& ConstantExpr::toLatex(ostream& os, bool /* paren */) const 
{
	os << value();
	return os;
}


XMLObject ConstantExpr::toXML() const 
{
	XMLObject rtn("Constant");
	rtn.addAttribute("value", Teuchos::toString(value()));
	return rtn;
}


