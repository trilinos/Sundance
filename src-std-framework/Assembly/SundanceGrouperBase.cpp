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
                                  BasisFamily& testBasis, 
                                  BasisFamily& unkBasis,
                                  MultiIndex& miTest, MultiIndex& miUnk,
                                  int& testID, int& unkID, 
                                  bool& isOneForm) const
{
  Tabs tab;

  MultipleDeriv::const_iterator iter;

  isOneForm = false;  

  if (functionalDeriv.size()==0) return;
  TEST_FOR_EXCEPTION(functionalDeriv.size()>2, InternalError,
                     "WeakFormBatch::extractWeakForm detected a functional "
                     "derivative of order > 2: " 
                     << functionalDeriv.toString());

  bool foundTest = false;
  bool foundUnk = false;

  for (iter = functionalDeriv.begin(); iter != functionalDeriv.end(); iter++)
    {
      const Deriv& d = *iter;
      
      TEST_FOR_EXCEPTION(!d.isFunctionalDeriv(), InternalError,
                         "WeakFormBatch::extractWeakForm "
                         "detected a non-functional derivative: "
                         << functionalDeriv.toString());
      
      const FunctionalDeriv* f = d.funcDeriv();
      
      const UnknownFuncElement* u 
        = dynamic_cast<const UnknownFuncElement*>(f->func());
      if (u != 0)
        {
          foundUnk = true;
          unkID = f->funcID();
          unkID = eqn.reducedUnkID(unkID);
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found unkID=" << unkID);
          const UnknownFunction* uf 
            = dynamic_cast<const  UnknownFunction*>(u->master());
          TEST_FOR_EXCEPTION(uf==0, InternalError, 
                             "WeakFormBatch::extractWeakForm failed to cast "
                             "UnknownFunctionStub to UnknownFunction");
          const FuncWithBasis* fb 
            = dynamic_cast<const  FuncWithBasis*>(uf);
          TEST_FOR_EXCEPTION(fb==0, InternalError, 
                             "WeakFormBatch::extractWeakForm failed to cast "
                             "UnknownFunction to funcWithBasis");
          int myIndex = u->myIndex();
          unkBasis = fb->basis()[myIndex];
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found unkBasis=" << unkBasis);
          miUnk = f->multiIndex();
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found unk multi index=" << miUnk.toString());
          continue;
        }
      
      const TestFuncElement* t 
        = dynamic_cast<const TestFuncElement*>(f->func());
      if (t != 0)
        {
          foundTest = true;
          testID = f->funcID();
          testID = eqn.reducedTestID(testID);
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found testID=" << testID);
          
          const TestFunction* tf 
            = dynamic_cast<const  TestFunction*>(t->master());
          TEST_FOR_EXCEPTION(tf==0, InternalError, 
                             "WeakFormBatch::extractWeakForm failed to cast "
                             "TestFunctionStub to TestFunction");
          const FuncWithBasis* fb 
            = dynamic_cast<const  FuncWithBasis*>(tf);
          TEST_FOR_EXCEPTION(fb==0, InternalError, 
                             "WeakFormBatch::extractWeakForm failed to cast "
                             "TestFunction to funcWithBasis");

          int myIndex = t->myIndex();
          testBasis = fb->basis()[myIndex];
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found testBasis=" << testBasis);
          miTest = f->multiIndex();
          SUNDANCE_OUT(verbosity() > VerbMedium, 
                       tab << "found test multi index=" << miTest.toString());
          continue;
        }
    }

  if (!foundUnk) isOneForm = true;
}
