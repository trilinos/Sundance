#include "SundanceExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
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
#include "SundanceBruteForceEvaluator.hpp"
#include "SundanceEvalVectorArray.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceSymbPreprocessor.hpp"


using namespace SundanceUtils;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}
static Time& ioTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("result output"); 
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

  RefCountPtr<EvaluatorFactory> factory 
    = rcp(new BruteForceEvaluatorFactory());

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  DerivSet d = SymbPreprocessor::setupExpr(e[0], 
                                           tests,
                                           unks,
                                           u0,
                                           region, factory.get(), 2);

  RefCountPtr<EvalVectorArray> results;

  ev->evaluate(mgr, results);

  results->print(cerr, d);
}

void testExpr(const Expr& e,  
              const Expr& tests,
              const Expr& unks,
              const Expr& u0, 
              const EvalContext& region)
{
  cerr << endl 
       << "-------- testing " << e.toString() << " -------- " << endl;

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

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());
      verbosity<SymbolicTransformation>() = VerbSilent;
      verbosity<Evaluator>() = VerbSilent;
      verbosity<EvalVector>() = VerbSilent;
      verbosity<EvaluatableExpr>() = VerbSilent;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionStub("u");
			Expr w = new UnknownFunctionStub("w");
			Expr v = new TestFunctionStub("v");
			Expr s = new TestFunctionStub("s");

      cerr << "u=" << u << endl;
      cerr << "v=" << v << endl;

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr u0 = new DiscreteFunctionStub("u0");
      Expr w0 = new DiscreteFunctionStub("w0");

      Array<Expr> tests;

      tests.append(v);

      tests.append(s + v);

      tests.append(u*v);

      tests.append(u*v + s);

      tests.append(s + u*v);

      tests.append(s + v+u);

      tests.append(v + v*u*u);

      tests.append(v*u*u + v);

      tests.append(v*w + v*u*u);

      tests.append(v*u*u + v*w);

      tests.append(v*(u+w));

      tests.append((u+w)*v);

      tests.append((v+s)*(u+w));

      tests.append(dx*v);

      tests.append(dx*v + dx*s);

      tests.append((dx*u)*(dx*v));

      tests.append(u*(dx*v));

      tests.append((dx*v)*u);

      tests.append((dx*u)*v);

      tests.append(v*(dx*u));

      tests.append(v*u*dx*u + v*w*dy*u);


      for (int i=0; i<tests.length(); i++)
        {
          RegionQuadCombo rqc(rcp(new CellFilterStub()), 
                              rcp(new QuadratureFamilyStub(0)));
          EvalContext context(rqc, EvalContext::nextID());
          testExpr(tests[i], 
                   SundanceCore::List(v, s),
                   SundanceCore::List(u, w),
                   SundanceCore::List(u0, w0),
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
