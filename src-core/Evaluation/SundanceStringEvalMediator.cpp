/* @HEADER@ */
/* @HEADER@ */


#include "SundanceStringEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceExceptions.hpp"

using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace SundanceCore;
using namespace SundanceUtils;
using namespace TSFExtended;
using namespace Teuchos;
using namespace std;



StringEvalMediator::StringEvalMediator()
  : AbstractEvalMediator() 
{}

void StringEvalMediator::evalCoordExpr(const CoordExpr* expr,
                                       RefCountPtr<EvalVector>& vec) const
{
  SUNDANCE_OUT(verbosity() > VerbSilent, "evaluating coord expr " << expr->toXML().toString());
  
  vec->setString(expr->name());
}

void StringEvalMediator
::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
                          const Array<MultiIndex>& mi,
                          Array<RefCountPtr<EvalVector> >& vec) const 
{
  static Array<string> coordNames;

  if (coordNames.size() != 3)
    {
      coordNames.resize(3);
      coordNames[0] = "x";
      coordNames[1] = "y";
      coordNames[2] = "z";
    }

  string funcName = expr->name();
  
  for (unsigned int i=0; i<mi.size(); i++)
    {
      if (mi[i].order()==0)
        {
          vec[i]->setString(funcName);
        }
      else
        {
          int dir = mi[i].firstOrderDirection();
          string deriv = "D[" + funcName + ", " + coordNames[dir] + "]";
          vec[i]->setString(deriv);
        }
    }
}
