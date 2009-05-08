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
#include "SundanceProductTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
#include "SundanceIntegral.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceStringEvalMediator.hpp"

using namespace SundanceUtils;
using namespace SundanceCore;
using namespace Teuchos;
using namespace TSFExtended;
using SundanceCore::List;


static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

static Time& doitTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("doit"); 
  return *rtn;
}



void doit(const Expr& e, 
          const Expr& tests,
          const Expr& unks,
          const Expr& u0, 
          const EvalContext& region)
{
  TimeMonitor t0(doitTimer());
  EvalManager mgr;
  mgr.setRegion(region);

  static RefCountPtr<AbstractEvalMediator> mediator 
    = rcp(new StringEvalMediator());

  mgr.setMediator(mediator);
  mgr.setVerbosity(5);

  Expr params;
  Expr fixed;

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  DerivSet d = SymbPreprocessor::setupFwdProblem(e[0], 
                                                 tests,
                                                 unks,
                                                 u0,
                                                 params,
                                                 params,
                                                 fixed, fixed,
                                                 fixed, fixed,
                                                 region);

  Tabs tab;
  Out::os() << tab << *ev->sparsitySuperset(region) << endl;

  Array<double> constantResults;
  Array<RefCountPtr<EvalVector> > vectorResults;

  Out::os() << tab << "evaluating..." << endl;
  ev->evaluate(mgr, constantResults, vectorResults);

  Out::os() << tab << "done evaluating..." << endl;
  ev->sparsitySuperset(region)->print(Out::os(), vectorResults, constantResults);

}



void testExpr(const Expr& e,  
              const Expr& tests,
              const Expr& unks,
              const Expr& u0, 
              const EvalContext& region)
{
  Out::os() << endl 
       << "------------------------------------------------------------- " << endl;
  Out::os()  << "-------- testing " << e.toString() << " -------- " << endl;
  Out::os() << endl 
       << "------------------------------------------------------------- " << endl;

  try
    {
      doit(e, tests, unks, u0, region);
    }
  catch(exception& ex)
    {
      Out::os() << "EXCEPTION DETECTED!" << endl;
      Out::os() << ex.what() << endl;
      exit(1);
    }
}

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      Expr::showAllParens() = true;
//      ProductTransformation::optimizeFunctionDiffOps()=true;

      EvalVector::shadowOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

			Expr u = new UnknownFunctionStub("u");
			Expr v = new TestFunctionStub("v");

      Handle<CellFilterStub> interior = rcp(new CellFilterStub());
      Handle<QuadratureFamilyStub> quad = rcp(new QuadratureFamilyStub(1));
      
      Expr u0 = new DiscreteFunctionStub("u0");
      Expr w0 = new DiscreteFunctionStub("w0");
      Expr zero = new ZeroExpr();

      Array<Expr> tests;



      int maxDiffOrder=2;
//      tests.append(v*u);
//      tests.append((dx*v)*(dx*u));
      tests.append((dx*v));
//      tests.append(v*(dx*(u - u0)));


      for (int i=0; i<tests.length(); i++)
        {
          RegionQuadCombo rqc(rcp(new CellFilterStub()), 
                              rcp(new QuadratureFamilyStub(1)));
          EvalContext context(rqc, maxDiffOrder, EvalContext::nextID());
          context.setSetupVerbosity(5);
          testExpr(tests[i], 
                   SundanceCore::List(v),
                   SundanceCore::List(u),
                   SundanceCore::List(zero),
                   context);
        }

      

      TimeMonitor::summarize();
    }
	catch(exception& e)
		{
			Out::println(e.what());
		}


  
}
