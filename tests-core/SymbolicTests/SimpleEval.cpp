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
              const EvalContext& region)
{
  cerr << endl 
       << "------------------------------------------------------------- " << endl;
  cerr  << "-------- testing " << e.toString() << " -------- " << endl;
  cerr << endl 
       << "------------------------------------------------------------- " << endl;

  try
    {
      doit(e, tests, unks, u0, region);
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

      verbosity<SymbolicTransformation>() = VerbExtreme;
      verbosity<Evaluator>() = VerbSilent;
      verbosity<EvalVector>() = VerbSilent;
      verbosity<EvaluatableExpr>() = VerbSilent;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

			Expr u = new UnknownFunctionStub("u");
			Expr v = new TestFunctionStub("v");

			Expr ux = new UnknownFunctionStub("ux");
			Expr vx = new TestFunctionStub("vx");

			Expr uy = new UnknownFunctionStub("uy");
			Expr vy = new TestFunctionStub("vy");

			Expr p = new UnknownFunctionStub("p");
			Expr q = new TestFunctionStub("q");double pi = 4.0*atan(1.0);
      Expr sx = sin(pi*x);
      Expr cx = cos(pi*x);
      Expr sy = sin(pi*y);
      Expr cy = cos(pi*y);
      Expr psiExact = pow(pi, -3.0) * sx*sy;
      Expr uExact = pow(pi, -2.0)*List(-sx*cy, cx*sy);
      Expr fy = 4.0*cx*sy;

      Handle<CellFilterStub> interior = rcp(new CellFilterStub());
      Handle<QuadratureFamilyStub> quad = rcp(new QuadratureFamilyStub(1));

      Expr grad = List(dx, dy);
      Expr eqn = Integral(interior, (grad*vx)*(grad*ux)  
                          + (grad*vy)*(grad*uy)  - p*(dx*vx+dy*vy)
                          + q*(dx*ux+dy*uy) - vy*fy,
                          quad);

			Expr w = new UnknownFunctionStub("w");
			Expr s = new TestFunctionStub("s");

      cerr << "u=" << u << endl;
      cerr << "v=" << v << endl;


      Expr u0 = new DiscreteFunctionStub("u0");
      Expr w0 = new DiscreteFunctionStub("w0");
      Expr zero = new ZeroExpr();

      Array<Expr> tests;

      Expr z = Complex(x,y);
      Expr I = Complex(0.0, 1.0);

      cout << "z = " << x + I*y << endl;

      cout << "u + y + v + s + x + w = " << u + y + v + s + x + w << endl;

      cout << "|z|^2 = " << (dx*(x + I*y)) * (dx*(x - I*y)) << endl;



#ifdef BLAh
      

      tests.append(v*(dx*(u - u0)));


      for (int i=0; i<tests.length(); i++)
        {
          RegionQuadCombo rqc(rcp(new CellFilterStub()), 
                              rcp(new QuadratureFamilyStub(1)));
          EvalContext context(rqc, maxDiffOrder, EvalContext::nextID());
          testExpr(tests[i], 
                   SundanceCore::List(v, s),
                   SundanceCore::List(u, w),
                   SundanceCore::List(zero, zero),
                   context);
        }

      Expr uu0;
      {
        Expr uu = new UnknownFunctionStub("uu", 2);
        uu0 = uu[0];
      }


      cerr << "uu0 = " << uu0 << endl;
#endif
      TimeMonitor::summarize();
    }
	catch(exception& e)
		{
			Out::println(e.what());
		}


  
}
