#include "SundanceExpr.hpp"
#include "SundanceStdMathOps.hpp"
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
#include "SundanceEvalVector.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceStringEvalMediator.hpp"
#include "SundanceEvaluationTester.hpp"

using namespace SundanceUtils;
using namespace SundanceTesting;
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




int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);
      Tabs tabs;
      TimeMonitor timer(totalTimer());

      //       verbosity<SymbolicTransformation>() = VerbSilent;
      //       verbosity<EvaluationTester>() = VerbExtreme;
      //       verbosity<Evaluator>() = VerbExtreme;
      //       verbosity<EvalVector>() = VerbExtreme;
      //       verbosity<EvaluatableExpr>() = VerbExtreme;
      //       verbosity<AbstractEvalMediator>() = VerbExtreme;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = false;

      Time stopwatch("test");

      for (int n=1; n<20; n++)
        {
          double t0 = stopwatch.wallTime();
          stopwatch.start();
          int depth = n;

          ADField U(ADBasis(1), sqrt(2.0));
          Expr u = new TestUnknownFunction(U, "u");
          Expr x = new CoordExpr(0);

          Expr expr = u + x;

          for (int i=0; i<depth; i++)
            {
              expr = expr * expr;
            }

          EvaluationTester tester(expr);
          Array<double> df;
          Array<Array<double> > df2;
          tester.evaluate(df, df2);
          stopwatch.stop();
          double t1 = stopwatch.wallTime();
          cerr << n << "   "  << tester.numNodes() << "    " << t1-t0 << endl;

        }
    }
	catch(exception& e)
		{
			Out::println(e.what());
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
