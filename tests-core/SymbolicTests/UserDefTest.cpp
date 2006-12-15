#include "SundanceExpr.hpp"
#include "SundanceUserDefOp.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceListExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceProductTransformation.hpp"
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
using SundanceCore::List;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}


#define LOUD()                                          \
  {                                                     \
    verbosity<EvaluationTester>() = VerbExtreme;        \
    verbosity<Evaluator>() = VerbExtreme;               \
    verbosity<SparsitySuperset>() = VerbSilent;         \
    verbosity<EvalVector>() = VerbSilent;               \
    verbosity<EvaluatableExpr>() = VerbExtreme;         \
    verbosity<AbstractEvalMediator>() = VerbExtreme;    \
  }

#define QUIET()                                         \
  {                                                     \
    verbosity<EvaluationTester>() = VerbSilent;         \
    verbosity<Evaluator>() = VerbSilent;                \
    verbosity<SparsitySuperset>() = VerbSilent;         \
    verbosity<EvalVector>() = VerbSilent;               \
    verbosity<EvaluatableExpr>() = VerbSilent;          \
    verbosity<AbstractEvalMediator>() = VerbSilent;     \
  }


class MyFunc : public UserDefFunctor
{
public:
  MyFunc() : UserDefFunctor("MyFunc"){;}
  virtual ~MyFunc(){;}
  double eval1(const Array<double>& vars, double* df) const ;
  double eval0(const Array<double>& vars) const ;
  int numArgs() const {return 3;}
};


double MyFunc::eval1(const Array<double>& vars, double* df) const
{
  double rtn = vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
  df[0] = sin(vars[1]) + 2.0*vars[2];
  df[1] = vars[0]*cos(vars[1]);
  df[2] = 2.0*vars[0];
  return rtn;
}

double MyFunc::eval0(const Array<double>& vars) const
{
  return vars[0]*sin(vars[1]) + 2.0*vars[2]*vars[0];
}


int main(int argc, char** argv)
{
  try
		{
      GlobalMPISession session(&argc, &argv);
      Tabs tabs;
      TimeMonitor timer(totalTimer());

      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

      ProductTransformation::optimizeFunctionDiffOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);

      ADField U(ADBasis(1), sqrt(2.0));
      ADField W(ADBasis(2), sqrt(3.0));

      ADCoord X(0);
      ADCoord Y(1);
      ADCoord Z(2);

      ADDerivative Dx(0);
      ADDerivative Dy(1);
      ADDerivative Dz(2);

      ADReal C_old = sin(X)*sin(Y);

			Expr u = new TestUnknownFunction(U, "u");
			Expr w = new TestUnknownFunction(W, "w");

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr p = 2.0*(dx*u) + w;
      Expr q = 2.0*u*w + x*w;
      Expr f = new UserDefOp(List(p, q, x), rcp(new MyFunc()));

      //      verbosity<Evaluator>() = VerbExtreme;


      


      Array<double> df1(2);
      Array<double> df2(2);

      cout << "evaluating symb expr, eval1()" << endl;
      EvaluationTester fTester2(p*sin(q) + 2.0*p*x, 1);
      double f2 = fTester2.evaluate(df2);
      cout << "value (symb)    = " << f2 << endl;
      cout << "deriv (symb)    = " << df2 << endl;

      cout << "evaluating user def expr, eval1()" << endl;

      EvaluationTester fTester1(f, 1);
      double f1 = fTester1.evaluate(df1);

      cout << "value (functor) = " << f1 << endl;
      cout << "deriv (functor) = " << df1 << endl;


      cout << "evaluating symb expr, eval0()" << endl;
      EvaluationTester fTester3(p*sin(q) + 2.0*p*x, 0);
      double f3 = fTester3.evaluate();
      cout << "value (symb)    = " << f3 << endl;


      cout << "evaluating user def expr, eval0()" << endl;
      EvaluationTester fTester4(f, 0);
      double f4 = fTester4.evaluate();
      cout << "value (functor)    = " << f4 << endl;


      

      double tol = 1.0e-10;

      bool isOK = ::fabs(f2 - f1) < tol;
      for (unsigned int i=0; i<df1.size(); i++)
        {
          isOK = isOK && ::fabs(df2[i]-df1[i]) < tol;
        }


      if (isOK)
        {
          cerr << "all tests PASSED!" << endl;
        }
      else
        {
          cerr << "overall test FAILED!" << endl;
        }
      TimeMonitor::summarize();
    }
	catch(exception& e)
		{
      cerr << "overall test FAILED!" << endl;
      cerr << "detected exception: " << e.what() << endl;
		}

  
}
