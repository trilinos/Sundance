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


#define TESTER(expr, adExpr)\
{\
Tabs tabs1;\
cerr << tabs1 << endl << tabs1\
<< "------------- Testing " << (expr).toString() << " -----------"\
<< endl << tabs1 << endl;\
EvaluationTester tester((expr));\
bool thisTestIsOK = true;\
double f = tester.fdEvaluate(fdStep, tol1, tol2, thisTestIsOK);\
if (!thisTestIsOK)\
{\
  isOK = false;\
}\
double adf = (adExpr).value();\
cerr << tabs1 << "expr value = " << f << " check=" << adf \
<< " |f-check|=" << fabs(f-adf) << endl;\
double fError = fabs(f-adf);\
if (fError > tol1)\
{\
thisTestIsOK=false;\
cerr << "value computation FAILED" << endl;\
isOK = false;\
}\
if (!thisTestIsOK)\
{\
failures.append(#expr);\
}\
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
      //   verbosity<Evaluator>() = VerbExtreme;
      verbosity<EvalVector>() = VerbSilent;
      // verbosity<EvaluatableExpr>() = VerbExtreme;
      verbosity<AbstractEvalMediator>() = VerbSilent;
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

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

			Expr u = new TestUnknownFunction(U, "u");
			Expr w = new TestUnknownFunction(W, "w");

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      double tol1 = 1.0e-6;
      double tol2 = 1.0e-6;
      double fdStep = 1.0e-3;
      bool isOK = true;
      Array<string> failures;


      TESTER(u, U);

      TESTER(-u, -U);

      /* ----------- cases of sum expressions ------------------- */

      /* tests const-const and vec-vec sums */
      TESTER( u + u, U + U );

      /* tests const-const and vec-vec subtractions */
      TESTER( u - u, U - U );

      /* tests vec-vec and const-0 sums */
      TESTER( u + x, U + X );

      /* tests vec-vec and const-0 subtractions */
      TESTER( u - x, U - X );

      /* tests vec-vec and 0-const sums */
      TESTER( x + u, X + U );

      /* tests vec-vec and 0-const subtractions */
      TESTER( x - u, X - U );

      /* tests vec-vec and const-vec sums */
      TESTER( u + x*u, U + X*U );

      /* tests vec-vec and const-vec subtractions */
      TESTER( u - x*u, U - X*U );

      /* tests vec-vec and vec-const sums */
      TESTER( x*u + u, X*U + U );

      /* tests vec-vec and vec-const subtractions */
      TESTER( x*u - u, X*U - U );

      /* ----------- cases of product expressions ------------------- */

      /* */
      TESTER( u*u, U*U );

      /* */
      TESTER( u*u*u, U*U*U );

      /* */
      TESTER( u*w, U*W );

      /* */
      TESTER( u*u*w, U*U*W );

      /* */
      TESTER( w*u*w, W*U*W );

      /* */
      TESTER( w*u*u, W*U*U );

      /* */
      TESTER( (u+w)*u, (U+W)*U );

      /* */
      TESTER( u*(u+w), U*(U+W) );

      /* */
      TESTER( (u+w*u)*u, (U+W*U)*U );

      /* */
      TESTER( u*(u+w*u), U*(U+W*U) );

      /* */
      TESTER( (u+w)*(u+w), (U+W)*(U+W) );

      /* */
      TESTER( (u+x)*(u+w), (U+X)*(U+W) );

      /* */
      TESTER( (u+w*u)*(u+w), (U+W*U)*(U+W) );

      /* */
      TESTER( (u+w)*(u+w*u), (U+W)*(U+W*U) );


      /* */
      TESTER( (x*u)*(y*u), (X*U)*(Y*U) );

      /* */
      TESTER( (x*u)*(y*w), (X*U)*(Y*W) );


      /* */
      TESTER( (2.0*u)*(y*u), (2.0*U)*(Y*U) );

      /* */
      TESTER( (x*u)*(2.0*u), (X*U)*(2.0*U) );

      /* */
      TESTER( (2.0*u)*(4.0*u), (2.0*U)*(4.0*U) );

      /* */
      TESTER( (2.0*u*u)*(y*u), (2.0*U*U)*(Y*U) );

      /* */
      TESTER( (x*u)*(2.0*u*u), (X*U)*(2.0*U*U) );

      /* */
      TESTER( 2.0*(y*u), 2.0*(Y*U) );

      /* */
      TESTER( (x*u)*2.0, (X*U)*2.0 );

      /* */
      TESTER( u*(y*u), U*(Y*U) );

      /* */
      TESTER( (x*u)*u, (X*U)*U );

      

      /* -------------- tests of diff ops ----------------------- */


     

      /* */
      TESTER((dx*u), (Dx*U));

      /* */
      TESTER((dx*x), (Dx*X));

      /* */
      TESTER((dx*y), (Dx*Y));

      
      TESTER((dx*(u+w)), (Dx*(U+W)));

      TESTER((dx*(u-w)), (Dx*(U-W)));

      TESTER((dx*(x+w)), (Dx*(X+W)));

      TESTER((dx*(x-w)), (Dx*(X-W)));

      TESTER((dx*(x*w+w)), (Dx*(X*W+W)));

      TESTER((dx*(x*w-w)), (Dx*(X*W-W)));

      TESTER((dx*(w + x*w)), (Dx*(W + X*W)));


      TESTER((dx*(w - x*w)), (Dx*(W - X*W)));

      TESTER((dx*(x*w)), (Dx*(X*W)));

      TESTER((dx*(y*w)), (Dx*(Y*W)));

      TESTER((dx*(u*w)), (Dx*(U*W)));

      Expr g = x*x + y*y;

      TESTER((u*(dx*(g) + dy*(g))), (U*(Dx*(X*X + Y*Y) + Dy*(X*X + Y*Y))));
       


      TESTER((u*g + w*g), (U*(X*X + Y*Y) + W*(X*X + Y*Y)));



      
      TESTER((dx*(2.0*u+4.0*w)), (Dx*(2.0*U+4.0*W)));

      TESTER((dx*(2.0*u-4.0*w)), (Dx*(2.0*U-4.0*W)));

      TESTER((dx*(2.0*x+4.0*w)), (Dx*(2.0*X+4.0*W)));

      TESTER((dx*(2.0*x-4.0*w)), (Dx*(2.0*X-4.0*W)));



      TESTER((dx*(x*w+u)), (Dx*(X*W+U)));

      TESTER((dx*(w + x*w)), (Dx*(W + X*W)));

      TESTER((dx*(x*w+u*y)), (Dx*(X*W+U*Y))); 

      TESTER(((dx*u)*(dx*w)), ((Dx*U)*(Dx*W)));

      TESTER((x*(dx*u)*(dx*w)), (X*(Dx*U)*(Dx*W))); 

      TESTER((u*(dx*u)*(dx*w)), (U*(Dx*U)*(Dx*W)));

      TESTER((x*(dx*u)*(dx*w)), (X*(Dx*U)*(Dx*W)));

      TESTER((u*(dx*u)*(dx*w)), (U*(Dx*U)*(Dx*W)));

      TESTER((x*(dx*u)*(dx*u)), (X*(Dx*U)*(Dx*U)));

      TESTER((u*(dx*u)*(dx*u)), (U*(Dx*U)*(Dx*U)));

      TESTER((u*(dx*u)), (U*(Dx*U)));

      TESTER((u*(dx*(u*x))), (U*(Dx*(U*X))));

      /* Unary operators */
      TESTER(sin(x), sin(X));

      TESTER(w*sin(u), W*sin(U));

      TESTER(sin(u*cos(x)), sin(U*cos(X)));

      TESTER(sin(u*x), sin(U*X));

      TESTER(sin(u*x+0.5*w*u), sin(U*X+0.5*W*U));

      TESTER(sin(cos(u*x)+0.5*w*u), sin(cos(U*X)+0.5*W*U));


      if (isOK)
        {
          cerr << "all tests PASSED!" << endl;
        }
      else
        {
          cerr << "test FAILED!" << endl;
          cerr << endl << "failed exprs: " << endl
               << endl;
          for (unsigned int i=0; i<failures.size(); i++)
            {
              cerr << failures[i] << endl;
            }
          cerr << endl;
        }
    }
	catch(exception& e)
		{
			Out::println(e.what());
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
