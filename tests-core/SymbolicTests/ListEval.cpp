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
#include "SundanceEquationSet.hpp"


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
      EquationSet::classVerbosity() = VerbHigh;
      Expr::showAllParens() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionBase("u", 2);
			Expr v = new TestFunctionBase("v", 2);
			Expr p = new UnknownFunctionBase("p");
			Expr q = new TestFunctionBase("q");


      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      Expr grad = SundanceCore::List(dx, dy);

      Expr u0 = new DiscreteFunctionBase("u0", 2);
      Expr p0 = new DiscreteFunctionBase("p0");

      Expr navierStokes = v*(u*grad)*u + p*(grad*v)
        + (grad*v[0])*(grad*u[0]) + (grad*v[1])*(grad*u[1]) 
        + q*(grad*u);

      Expr bc = v[0]*(u[0]-1.0) + v[1]*(u[1]-2.0);

      EvalRegion internalRegion(rcp(new CellFilterBase()), 
                                rcp(new QuadratureFamilyBase(0)));

      EvalRegion bcRegion(rcp(new CellFilterBase()), 
                                rcp(new QuadratureFamilyBase(0)));

      EvalManager mgr;
      mgr.setRegion(internalRegion);

      RefCountPtr<EvaluatorFactory> factory 
        = rcp(new BruteForceEvaluatorFactory());
      
      DerivSet dIn = SymbPreprocessor::setupExpr(navierStokes, 
                                               List(v,q),
                                               List(u,p),
                                               List(u0,p0),
                                               internalRegion,
                                               factory.get());
      
      DerivSet dBC = SymbPreprocessor::setupExpr(bc,
                                                 List(v,q),
                                                 List(u,p),
                                                 List(u0,p0),
                                                 bcRegion,
                                                 factory.get());
      RefCountPtr<EvalVectorArray> results;


      const EvaluatableExpr* ev 
        = dynamic_cast<const EvaluatableExpr*>(navierStokes[0].ptr().get());
      ev->evaluate(mgr, results);
      
      cerr << "---- internal -------" << endl;
      results->print(cerr, dIn);

      mgr.setRegion(bcRegion);

      const EvaluatableExpr* ev2 
        = dynamic_cast<const EvaluatableExpr*>(bc[0].ptr().get());
      
      ev2->evaluate(mgr, results);
      
      cerr << "---- bc -------" << endl;
      results->print(cerr, dBC);

      Expr weak = Integral(new CellFilterBase(), navierStokes);
      Expr weakBC = EssentialBC(new CellFilterBase(), bc);
      EquationSet eqns(weak, weakBC, List(v,q), List(u,p), List(u0,p0),
                       rcp(new BruteForceEvaluatorFactory()));
      
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
