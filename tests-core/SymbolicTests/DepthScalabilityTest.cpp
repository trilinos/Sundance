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
#include "Teuchos_GlobalMPISession.hpp"
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




int main(int argc, char** argv)
{
  
  try
		{
      GlobalMPISession session(&argc, &argv);
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

      for (int n=1; n<10; n++)
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

      for (int n=2; n<40; n+=4)
        {
          double t0 = stopwatch.wallTime();
          stopwatch.start();
          int depth = n;

          ADField U(ADBasis(1), sqrt(2.0));
          Array<Expr> unks(n);
          Array<ADField> Unks(n);
          for (int i=0; i<n; i++) 
            {
              Unks[i] = ADField(ADBasis(1), sqrt(double(1+i)));
              unks[i] = new TestUnknownFunction(Unks[i], "u" + Teuchos::toString(i));
            }
          Expr x = new CoordExpr(0);

          Expr expr = x;
          for (int i=0; i<depth; i++)
            {
              for (int j=0; j<depth; j++)
                {
                  expr = expr + unks[i]*unks[j];
                }
            }

          EvaluationTester tester(expr);
          Array<double> df;
          Array<Array<double> > df2;
          tester.evaluate(df, df2);
          stopwatch.stop();
          double t1 = stopwatch.wallTime();
          cerr << n << "   "  << tester.numNodes() << "    " << t1-t0 << endl;

        }
      TimeMonitor::summarize();
    }
	catch(exception& e)
		{
			Out::println(e.what());
		}


  
}
