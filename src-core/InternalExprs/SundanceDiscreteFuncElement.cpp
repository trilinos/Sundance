/* @HEADER@ */
/* @HEADER@ */

#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

DiscreteFuncElement::DiscreteFuncElement(DiscreteFunctionStub* master, 
                                         const string& name,
                                         int myIndex)
	: LeafExpr(), 
    FuncElementBase(name),
    master_(master),
    myIndex_(myIndex)
{}


XMLObject DiscreteFuncElement::toXML() const 
{
	XMLObject rtn("DiscreteFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}

