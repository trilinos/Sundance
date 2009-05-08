#include "SundanceExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceUnknownParameter.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
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
using namespace SundanceCore;
using namespace Teuchos;
using namespace TSFExtended;
using SundanceCore::List;

using SundanceCore::UnknownFuncElement;

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
          const Expr& unkParams,
          const Expr& paramVals, 
          const EvalContext& region)
{
  TimeMonitor t0(doitTimer());
  EvalManager mgr;
  mgr.setRegion(region);

  static RefCountPtr<AbstractEvalMediator> mediator 
    = rcp(new StringEvalMediator());

  mgr.setMediator(mediator);

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  Expr fixed;
  Expr fixed0;

  DerivSet d = SymbPreprocessor::setupSensitivities(e[0], 
                                                    tests,
                                                    unks,
                                                    u0,
                                                    unkParams,
                                                    paramVals,
                                                    fixed,
                                                    fixed0,
                                                    fixed0,
                                                    fixed0,
                                                    region);

  Tabs tab;
  cerr << tab << *ev->sparsitySuperset(region) << endl;
  //  ev->showSparsity(cerr, region);

  // RefCountPtr<EvalVectorArray> results;

  Array<double> constantResults;
  Array<RefCountPtr<EvalVector> > vectorResults;

  ev->evaluate(mgr, constantResults, vectorResults);

  ev->sparsitySuperset(region)->print(cerr, vectorResults, constantResults);

  
  // results->print(cerr, ev->sparsitySuperset(region).get());
}



void testExpr(const Expr& e,  
              const Expr& tests,
              const Expr& unks,
              const Expr& u0, 
          const Expr& unkParams,
          const Expr& paramVals, 
              const EvalContext& region)
{
  cerr << endl 
       << "------------------------------------------------------------- " << endl;
  cerr  << "-------- testing " << e.toString() << " -------- " << endl;
  cerr << endl 
       << "------------------------------------------------------------- " << endl;

  try
    {
      doit(e, tests, unks, u0, unkParams, paramVals, region);
    }
  catch(exception& ex)
    {
      cerr << "EXCEPTION DETECTED!" << endl;
      cerr << ex.what() << endl;
      // cerr << "repeating with increased verbosity..." << endl;
      //       cerr << "-------- testing " << e.toString() << " -------- " << endl;
      //       Evaluator::verbosity() = 2;
      //       EvalVector::verbosity() = 2;
      //       EvaluatableExpr::verbosity() = 2;
      //       Expr::showAllParens() = true;
      //       doit(e, region);
      exit(1);
    }
}

int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);

      TimeMonitor t(totalTimer());

      int maxDiffOrder = 2;

      verbosity<SymbolicTransformation>() = VerbSilent;
      verbosity<Evaluator>() = VerbExtreme;
      verbosity<EvalVector>() = VerbExtreme;
      verbosity<EvaluatableExpr>() = VerbExtreme;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionStub("u");
			Expr alpha = new UnknownParameter("alpha");
			Expr alpha0 = new Parameter(3.14, "alpha0");
			Expr beta = new UnknownParameter("beta");
			Expr beta0 = new Parameter(2.72, "beta0");
			Expr v = new TestFunctionStub("v");

      cerr << "u=" << u << endl;
      cerr << "v=" << v << endl;
      cerr << "alpha=" << alpha << endl;

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr u0 = new DiscreteFunctionStub("u0");
      Expr zero = new ZeroExpr();

      Array<Expr> tests;

      

      tests.append(v*dx*(u*u - x*u) + dx*(v*alpha*u) + v*sin(alpha*u) + 2.0*v + v*exp(u*beta));


      for (int i=0; i<tests.length(); i++)
        {
          RegionQuadCombo rqc(rcp(new CellFilterStub()), 
                              rcp(new QuadratureFamilyStub(1)));
          EvalContext context(rqc, maxDiffOrder, EvalContext::nextID());
          testExpr(tests[i], 
                   SundanceCore::List(v),
                   SundanceCore::List(u),
                   SundanceCore::List(u0),
                   SundanceCore::List(alpha, beta),
                   SundanceCore::List(alpha0, beta0),
                   context);
        }


    }
	catch(exception& e)
		{
			Out::println(e.what());
		}
  TimeMonitor::summarize();

  
}
