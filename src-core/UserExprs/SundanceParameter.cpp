/* @HEADER@ */
/* @HEADER@ */

#include "SundanceParameter.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

Parameter::Parameter(const string& name)
	: FuncElementBase(name), 
    SpatiallyConstantExpr(1.0)
{}


XMLObject Parameter::toXML() const 
{
	XMLObject rtn("Parameter");
	rtn.addAttribute("name", name());
	rtn.addAttribute("value", Teuchos::toString(value()));
	return rtn;
}



