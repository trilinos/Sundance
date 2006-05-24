#include "SundanceExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceStringEvalMediator.hpp"

using namespace SundanceUtils;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;
using SundanceCore::List;

using SundanceCore::Internal::UnknownFuncElement;

static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

void doit(const Expr& e, 
          const Expr& vars,
          const Expr& varEval,
          const Expr& unks,
          const Expr& unkEval,
          const Expr& fixed,
          const Expr& fixedEval, 
          const EvalContext& region)
{
  EvalManager mgr;
  mgr.setRegion(region);

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  Array<Set<MultiSet<int> > > funcDerivs(3);
  Array<Set<MultiIndex> > spatialDerivs(3);

  Expr vf = vars.flatten();
  Expr uf = unks.flatten();
  Expr ff = fixed.flatten();
  Expr v0 = varEval.flatten();
  Expr u0 = unkEval.flatten();
  Expr f0 = fixedEval.flatten();

  funcDerivs[0].put(MultiSet<int>());

  for (unsigned int i=0; i<ff.size(); i++)
    {
      const SymbolicFuncElement* fPtr
        = dynamic_cast<const SymbolicFuncElement*>(ff[i].ptr().get());
      RefCountPtr<DiscreteFuncElement> f0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(f0[i].ptr());
      RefCountPtr<ZeroExpr> f0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(f0[i].ptr());
      if (f0Ptr.get() != 0)
        {
          fPtr->substituteFunction(f0Ptr);
        }
      else if (f0ZeroPtr.get() != 0)
        {
          fPtr->substituteZero();
        }
      else
        {
          TEST_FOR_EXCEPT(true);
        }
    }      
  

  for (unsigned int i=0; i<vf.size(); i++)
    {
      const SymbolicFuncElement* vPtr
        = dynamic_cast<const SymbolicFuncElement*>(vf[i].ptr().get());
      RefCountPtr<DiscreteFuncElement> v0Ptr
        = rcp_dynamic_cast<DiscreteFuncElement>(v0[i].ptr());
      RefCountPtr<ZeroExpr> v0ZeroPtr
        = rcp_dynamic_cast<ZeroExpr>(v0[i].ptr());
      if (v0Ptr.get() != 0)
        {
          vPtr->substituteFunction(v0Ptr);
        }
      else if (v0ZeroPtr.get() != 0)
        {
          vPtr->substituteZero();
        }
      else
        {
          TEST_FOR_EXCEPT(true);
        }
      
      int vid = vPtr->funcID();
      
      funcDerivs[1].put(makeMultiSet<int>(vid));

      for (unsigned int j=0; j<uf.size(); j++)
        {
          const SymbolicFuncElement* uPtr
            = dynamic_cast<const SymbolicFuncElement*>(uf[j].ptr().get());
          RefCountPtr<DiscreteFuncElement> u0Ptr
            = rcp_dynamic_cast<DiscreteFuncElement>(u0[j].ptr());
          RefCountPtr<ZeroExpr> u0ZeroPtr
            = rcp_dynamic_cast<ZeroExpr>(u0[j].ptr());
          if (u0Ptr.get() != 0)
            {
              uPtr->substituteFunction(u0Ptr);
            }
          else if (u0ZeroPtr.get() != 0)
            {
              uPtr->substituteZero();
            }
          else
            {
              TEST_FOR_EXCEPT(true);
            }
          int uid = uPtr->funcID();
          funcDerivs[2].put(makeMultiSet<int>(vid, uid));
        }
    }

  Array<Set<MultipleDeriv> > RInput = ev->computeInputR(region, funcDerivs, spatialDerivs);

  ev->determineR(region, RInput);

  
  

  for (unsigned int i=0; i<RInput.size(); i++)
    {
      Tabs tabs0;
      cout << tabs0 << "---------------------------------------------------------------------"
           << endl;
      cout << tabs0 << "order = " << i << endl;
      cout << tabs0 << "constants " << endl;
      const Set<MultipleDeriv>& C = ev->findC(i, region);
      for (Set<MultipleDeriv>::const_iterator it=C.begin(); it != C.end(); it++)
        {
          Tabs tab;
          cout << tab << *it << endl;
        }
      cout << tabs0 << "variables " << endl;
      const Set<MultipleDeriv>& V = ev->findV(i, region);
      for (Set<MultipleDeriv>::const_iterator it=V.begin(); it != V.end(); it++)
        {
          Tabs tab;
          cout << tab << *it << endl;
        }
      cout << tabs0 << "all " << endl;
      const Set<MultipleDeriv>& R = ev->findR(i, region);
      for (Set<MultipleDeriv>::const_iterator it=R.begin(); it != R.end(); it++)
        {
          Tabs tab;
          cout << tab << *it << endl;
        }
    }

  ev->displayNonzeros(cout, region);
}



int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());

      int maxDiffOrder = 2;

      verbosity<SymbolicTransformation>() = VerbSilent;
      verbosity<Evaluator>() = VerbSilent;
      verbosity<EvalVector>() = VerbSilent;
      verbosity<EvaluatableExpr>() = VerbExtreme;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionStub("u");
			Expr w = new UnknownFunctionStub("w");
			Expr v = new TestFunctionStub("v");
			Expr alpha = new UnknownFunctionStub("alpha");

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr u0 = new DiscreteFunctionStub("u0");
      Expr w0 = new DiscreteFunctionStub("w0");
      Expr alpha0 = new DiscreteFunctionStub("alpha0");
      Expr zero = new ZeroExpr();

      Array<Expr> tests;

      

      tests.append(v*dx*(u*x));


      for (int i=0; i<tests.length(); i++)
        {
          RegionQuadCombo rqc(rcp(new CellFilterStub()), 
                              rcp(new QuadratureFamilyStub(1)));
          EvalContext context(rqc, maxDiffOrder, EvalContext::nextID());
          doit(tests[i], 
                   SundanceCore::List(v),
                   SundanceCore::List(zero),
                   SundanceCore::List(u, w),
                   SundanceCore::List(u0, w0),
                   SundanceCore::List(alpha),
                   SundanceCore::List(alpha0),
                   context);
        }


    }
	catch(exception& e)
		{
			Out::println(e.what());
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
