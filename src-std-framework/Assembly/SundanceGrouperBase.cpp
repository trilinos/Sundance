/* @HEADER@ */
/* @HEADER@ */

#include "SundanceGrouperBase.hpp"
#include "SundanceOut.hpp"
#include "SundanceTabs.hpp"
#include "SundanceFuncWithBasis.hpp"
#include "SundanceFunctionalDeriv.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFuncElement.hpp"
#include "SundanceUnknownFunction.hpp"
#include "SundanceTestFunction.hpp"

using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceUtils;
using namespace Teuchos;
using namespace TSFExtended;



void GrouperBase::extractWeakForm(const EquationSet& eqn,
                                  const MultipleDeriv& functionalDeriv,
                                  BasisFamily& varBasis, 
                                  BasisFamily& unkBasis,
                                  MultiIndex& miVar, MultiIndex& miUnk,
                                  int& varID, int& unkID, 
                                  bool& isOneForm) const
{
  Tabs tab;

  MultipleDeriv::const_iterator iter;

  isOneForm = false;  

  if (functionalDeriv.size()==0) return;
  TEST_FOR_EXCEPTION(functionalDeriv.size() > 2, InternalError,
                     "WeakFormBatch::extractWeakForm detected a functional "
                     "derivative of order > 2: " 
                     << functionalDeriv.toString());

  bool foundUnk = false;
  bool foundVar = false;

  for (iter = functionalDeriv.begin(); iter != functionalDeriv.end(); iter++)
    {
      const Deriv& d = *iter;
      
      TEST_FOR_EXCEPTION(!d.isFunctionalDeriv(), InternalError,
                         "WeakFormBatch::extractWeakForm "
                         "detected a non-functional derivative: "
                         << functionalDeriv.toString());
      
      const FunctionalDeriv* f = d.funcDeriv();
      
      const SymbolicFuncElement* s 
        = dynamic_cast<const SymbolicFuncElement*>(f->func());

      TEST_FOR_EXCEPTION(s==0, InternalError, 
                         "WeakFormBatch::extractWeakForm failed to cast "
                         "function to SymbolicFuncElement");
      

      int funcID = f->funcID();
      int myIndex = s->myIndex();

      if (eqn.hasUnkID(funcID))
        {
          const UnknownFuncElement* u
            = dynamic_cast<const UnknownFuncElement*>(s);
          TEST_FOR_EXCEPTION(u==0, InternalError, 
                             "WeakFormBatch::extractWeakForm could not cast "
                             "unknown function to UnknownFuncElement");
          foundUnk = true;
          unkID = eqn.reducedUnkID(funcID);

          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found unkID=" << unkID);

          const FuncWithBasis* fb 
            = dynamic_cast<const  FuncWithBasis*>(u->master());
          unkBasis = fb->basis()[myIndex];
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found unkBasis=" << unkBasis);

          miUnk = f->multiIndex();
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found unk multi index=" << miUnk.toString());
        }
      else
        {
          foundVar = true;
          varID = eqn.reducedVarID(funcID);

          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found varID=" << varID);

          const UnknownFuncElement* u
            = dynamic_cast<const UnknownFuncElement*>(s);

          const TestFuncElement* t
            = dynamic_cast<const TestFuncElement*>(s);

          TEST_FOR_EXCEPTION(u==0 && t==0, InternalError, 
                             "WeakFormBatch::extractWeakForm could not cast "
                             "variational function to either an "
                             "UnknownFuncElement or a TestFuncElement");

          if (t != 0) 
            {
              const FuncWithBasis* fb 
                = dynamic_cast<const  FuncWithBasis*>(t->master());
              varBasis = fb->basis()[myIndex];
            }
          else
            {
              const FuncWithBasis* fb 
                = dynamic_cast<const  FuncWithBasis*>(u->master());
              varBasis = fb->basis()[myIndex];
            }
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found varBasis=" << varBasis);

          miVar = f->multiIndex();
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found var multi index=" << miVar.toString());
        }
    }

  if (!foundUnk) isOneForm = true;
}
