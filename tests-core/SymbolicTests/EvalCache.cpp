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
#include "SundanceSymbPreprocessor.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
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

      Expr v = new TestFunctionStub("v");

			Expr a = new UnknownFunctionStub("a");
			Expr b = new UnknownFunctionStub("b");
			Expr c = new UnknownFunctionStub("c");

			Expr a0 = new DiscreteFunctionStub("a0");
			Expr b0 = new DiscreteFunctionStub("b0");
			Expr c0 = new DiscreteFunctionStub("c0");

      Expr dx = new Derivative(0);

			Expr d = a + b;
      Expr e = c*d;
      Expr expr = v*(a + dx*e + d*d + e*d);




      

      RegionQuadCombo region(rcp(new CellFilterStub()),
                        rcp(new QuadratureFamilyStub(0)));
      EvalContext context(region, 0);
      EvalManager mgr;
      mgr.setRegion(context);

      RefCountPtr<EvaluatorFactory> factory 
        = rcp(new BruteForceEvaluatorFactory());
      
      DerivSet ds = SymbPreprocessor::setupExpr(expr,
                                                SundanceCore::List(v),
                                                SundanceCore::List(a,b,c),
                                                SundanceCore::List(a0,b0,c0),
                                                context,
                                                factory.get(), 2);
      
      RefCountPtr<EvalVectorArray> results;


      const EvaluatableExpr* ev 
        = dynamic_cast<const EvaluatableExpr*>(expr[0].ptr().get());

      cerr << "---- first pass -------" << endl;
      ev->evaluate(mgr, results);
      results->print(cerr, ds);

      
      cerr << "---- second pass -------" << endl;
      ev->flushResultCache();
      ev->evaluate(mgr, results);
      results->print(cerr, ds);


    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
