/* @HEADER@ */
/* @HEADER@ */

#include "SundanceExprFieldWrapper.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteFuncElement.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;


ExprFieldWrapper::ExprFieldWrapper(const Expr& expr)
  : expr_(expr),
    vector_(),
    discreteSpace_(),
    map_(),
    indices_()
{
  if (expr.size()==1)
    {
      const DiscreteFunction* df 
        = dynamic_cast<const DiscreteFunction*>(expr[0].ptr().get());
      if (df != 0)
        {
          vector_ = df->vector();
          discreteSpace_ = df->discreteSpace();
          map_ = df->map();
          indices_ = tuple(0);
        }
      const DiscreteFuncElement* dfe 
        = dynamic_cast<const DiscreteFuncElement*>(expr[0].ptr().get());
      if (dfe != 0)
        {
          const DiscreteFunction* f
            = dynamic_cast<const DiscreteFunction*>(dfe->master());
          TEST_FOR_EXCEPTION(f == 0, RuntimeError,
                             "ExprFieldWrapper ctor argument "
                             << expr << " is not a discrete function");
          vector_ = f->vector();
          discreteSpace_ = f->discreteSpace();
          map_ = f->map();
          indices_ = tuple(dfe->myIndex());
        }

      TEST_FOR_EXCEPTION(df == 0 && dfe == 0, RuntimeError,
                         "ExprFieldWrapper ctor argument is not a discrete function");
    }
  else
    {
      TEST_FOR_EXCEPTION(expr.size() != 1, RuntimeError,
                         "non-scalar expr given to ExprFieldWrapper ctor");
    }
}


double ExprFieldWrapper::getData(int cellDim, int cellID, int elem) const
{
  Array<int> dofs;
  map_->getDOFsForCell(cellDim, cellID, indices_[elem], dofs);
  TEST_FOR_EXCEPTION(dofs.size() > 1, RuntimeError,
                     "too many DOFs found in ExprFieldWrapper::getData()");

  return vector_.getElement(dofs[0]);
}
    
