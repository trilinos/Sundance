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
#include "SundanceIntegral.hpp"
#include "SundanceEssentialBC.hpp"

using namespace SundanceUtils;
using namespace SundanceCore;
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

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionStub("u");
			Expr v = new TestFunctionStub("v");

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      Expr grad = List(dx, dy);

      Handle<CellFilterStub> interior = new CellFilterStub();
      Handle<CellFilterStub> left = new CellFilterStub();
      Handle<CellFilterStub> right = new CellFilterStub();
      Handle<CellFilterStub> top = new CellFilterStub();

      Handle<QuadratureFamilyStub> quad8 = new QuadratureFamilyStub(8);
      
      Expr eqn = Integral(interior, (grad*u)*(grad*v))
        + Integral(interior, y*v, quad8)
        + Integral(left, v*x*x);

      Expr bc = EssentialBC(right, v*u) + EssentialBC(top, v*(u-x));

      cerr << "eqn = " << endl << eqn.toString() << endl;
      cerr << "bc = " << endl << bc.toString() << endl;

      
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
