/* @HEADER@ */
/* @HEADER@ */

#include "SundanceParameter.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

Parameter::Parameter(const double& value, const string& name)
	: FuncElementBase(name), 
    SpatiallyConstantExpr(value)
{}


XMLObject Parameter::toXML() const 
{
	XMLObject rtn("Parameter");
	rtn.addAttribute("name", name());
	rtn.addAttribute("value", Teuchos::toString(value()));
	return rtn;
}



