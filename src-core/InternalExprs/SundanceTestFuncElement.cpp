/* @HEADER@ */
/* @HEADER@ */

#include "SundanceTestFuncElement.hpp"
#include "SundanceTestFunctionStub.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;

TestFuncElement::TestFuncElement(const TestFunctionStub* master,
                                 const string& name,
                                 const string& suffix,
                                 int myIndex)
	: SymbolicFuncElement(name, suffix, myIndex), master_(master)
{}


XMLObject TestFuncElement::toXML() const 
{
	XMLObject rtn("TestFuncElement");
	rtn.addAttribute("name", name());
	return rtn;
}


