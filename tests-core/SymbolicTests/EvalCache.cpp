#include "SundanceExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionBase.hpp"
#include "SundanceTestFunctionBase.hpp"
#include "SundanceDiscreteFunctionBase.hpp"
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
#include "SundanceEvalRegion.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceBruteForceEvaluator.hpp"
#include "SundanceEvalVectorArray.hpp"
#include "SundanceSymbPreprocessor.hpp"


using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;

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
      SymbolicTransformation::verbosity() = 0;
      Evaluator::verbosity() = 0;
      EvalVector::verbosity() = 0;
      EvaluatableExpr::verbosity() = 0;
      Expr::showAllParens() = true;

      Expr v = new TestFunctionBase("v");

			Expr a = new UnknownFunctionBase("a");
			Expr b = new UnknownFunctionBase("b");
			Expr c = new UnknownFunctionBase("c");

			Expr a0 = new DiscreteFunctionBase("a0");
			Expr b0 = new DiscreteFunctionBase("b0");
			Expr c0 = new DiscreteFunctionBase("c0");

      Expr dx = new Derivative(0);

			Expr d = a + b;
      Expr e = c*d;
      Expr expr = v*(a + dx*e + d*d + e*d);




      

      EvalRegion region("region");
      EvalManager mgr;
      mgr.setRegion(region);

      RefCountPtr<EvaluatorFactory> factory 
        = rcp(new BruteForceEvaluatorFactory());
      
      DerivSet ds = SymbPreprocessor::setupExpr(expr,
                                                SundanceCore::List(v),
                                                SundanceCore::List(a,b,c),
                                                SundanceCore::List(a0,b0,c0),
                                                region,
                                                factory.get());
      
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
