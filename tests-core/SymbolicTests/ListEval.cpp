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
#include "SundanceEquationSet.hpp"

using SundanceCore::List;
using namespace SundanceCore;
using namespace SundanceUtils;
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

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionStub("u", 2);
			Expr v = new TestFunctionStub("v", 2);
			Expr p = new UnknownFunctionStub("p");
			Expr q = new TestFunctionStub("q");


      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      Expr grad = SundanceCore::List(dx, dy);

      Expr u0 = new DiscreteFunctionStub("u0", 2);
      Expr p0 = new DiscreteFunctionStub("p0");

      Expr navierStokes = v*(u*grad)*u + p*(grad*v)
        + (grad*v[0])*(grad*u[0]) + (grad*v[1])*(grad*u[1]) 
        + q*(grad*u);

      Expr bc = v[0]*(u[0]-1.0) + v[1]*(u[1]-2.0);

      RegionQuadCombo internalRegion(rcp(new CellFilterStub()), 
                                rcp(new QuadratureFamilyStub(0)));
      EvalContext internalContext(internalRegion, 0);

      RegionQuadCombo bcRegion(rcp(new CellFilterStub()), 
                                rcp(new QuadratureFamilyStub(0)));
      EvalContext bcContext(bcRegion, 0);

      EvalManager mgr;
      mgr.setRegion(internalContext);

      RefCountPtr<EvaluatorFactory> factory 
        = rcp(new BruteForceEvaluatorFactory());
      
      DerivSet dIn = SymbPreprocessor::setupExpr(navierStokes, 
                                               List(v,q),
                                               List(u,p),
                                               List(u0,p0),
                                               internalContext,
                                               factory.get(), 2);
      
      DerivSet dBC = SymbPreprocessor::setupExpr(bc,
                                                 List(v,q),
                                                 List(u,p),
                                                 List(u0,p0),
                                                 bcContext,
                                                 factory.get(), 2);
      RefCountPtr<EvalVectorArray> results;


      const EvaluatableExpr* ev 
        = dynamic_cast<const EvaluatableExpr*>(navierStokes[0].ptr().get());
      ev->evaluate(mgr, results);
      
      cerr << "---- internal -------" << endl;
      results->print(cerr, dIn);

      mgr.setRegion(bcContext);

      const EvaluatableExpr* ev2 
        = dynamic_cast<const EvaluatableExpr*>(bc[0].ptr().get());
      
      ev2->evaluate(mgr, results);
      
      cerr << "---- bc -------" << endl;
      results->print(cerr, dBC);

      Expr weak = Integral(new CellFilterStub(), navierStokes);
      Expr weakBC = EssentialBC(new CellFilterStub(), bc);
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
