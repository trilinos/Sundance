/* @HEADER@ */
/* @HEADER@ */

#include "SundanceCoordExpr.hpp"
#include "SundanceCoordDeriv.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceSparsityPattern.hpp"

using namespace SundanceCore;
using namespace SundanceUtils;

using namespace SundanceCore::Internal;
using namespace Teuchos;

CoordExpr::CoordExpr(int dir, const string& name)
  : FuncElementBase(coordName(dir, name)), 
    LeafExpr(), dir_(dir)
{;}

XMLObject CoordExpr::toXML() const 
{
  XMLObject rtn("CoordExpr");
  rtn.addAttribute("dir", Teuchos::toString(dir_));
  rtn.addAttribute("name", name());
  return rtn;
}

string CoordExpr::coordName(int dir, const string& name)
{
  if (name.length() > 0) return name;
  switch(dir)
    {
    case 0:
      return "x";
    case 1:
      return "y";
    case 2:
      return "z";
    default:
      TEST_FOR_EXCEPTION(true, RuntimeError,
                         "CoordExpr::coordName direction out of range [0,2]");
      return "error";
    }
}

bool CoordExpr::hasNonzeroDeriv(const MultipleDeriv& d) const
{
  if (d.order()==0) return true;

  if (d.order()==1)
    {
      const CoordDeriv* sd = (*(d.begin())).coordDeriv();
      return (sd != 0 && sd->dir()==dir_);
    }
  return false;
}




