/* @HEADER@ */
/* @HEADER@ */

#include "SundanceUnknownFuncElement.hpp"
#include "SundanceUnknownFunctionBase.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::FrameworkInterface;
using namespace Teuchos;

UnknownFuncElement::UnknownFuncElement(const UnknownFunctionBase* master,
                                       const string& name,
                                       int myIndex)
	: SymbolicFuncElement(name, myIndex), master_(master)
{}


XMLObject UnknownFuncElement::toXML() const 
{
	XMLObject rtn("UnknownFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}

