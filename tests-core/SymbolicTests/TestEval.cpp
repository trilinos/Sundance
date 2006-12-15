#include "SundanceExpr.hpp"
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

#define TESTER(expr, adExpr)                                            \
  {                                                                     \
    Tabs tabs1;                                                         \
    cerr << tabs1 << endl << tabs1                                      \
         << "------------- Testing " << #expr << " -----------"        \
         << endl << tabs1 << endl;                                      \
    bool thisTestIsOK = true;                                           \
    try                                                                 \
      {                                                                 \
      EvaluationTester tester((expr));                                  \
      double f = tester.fdEvaluate(fdStep, tol1, tol2, thisTestIsOK);   \
      if (!thisTestIsOK)                                                \
        {                                                               \
          isOK = false;                                                 \
        }                                                               \
      double adf = (adExpr).value();                                    \
      cerr << tabs1 << "expr value = " << f << " check=" << adf         \
           << " |f-check|=" << fabs(f-adf) << endl;                     \
      double fError = fabs(f-adf);                                      \
      if (fError > tol1)                                                \
        {                                                               \
          thisTestIsOK=false;                                           \
          cerr << "value computation FAILED" << endl;                   \
          isOK = false;                                                 \
        }                                                               \
      }                                                                 \
    catch(std::exception& ex)                                           \
      {                                                                 \
        thisTestIsOK = false;                                           \
        isOK=false;                                                     \
        cerr << "exception: " << ex.what() << endl ;                    \
      }                                                                 \
    if (!thisTestIsOK)                                                  \
      {                                                                 \
        failures.append(#expr);                                         \
        cerr << "test " << (expr).toString() << " FAILED" << endl << endl;\
      }                                                                 \
    else                                                                \
      {                                                                 \
        cerr << "test " << (expr).toString() << " PASSED" << endl << endl; \
      }\
  }


#define XTESTER(X,Y)                            \
  LOUD();                                       \
  TESTER(X,Y);                                   \
  QUIET();                                      \
  goto finish;





int main(int argc, char** argv)
{
  try
		{
      GlobalMPISession session(&argc, &argv);
      Tabs tabs;
      TimeMonitor timer(totalTimer());

      //verbosity<SymbolicTransformation>() = VerbExtreme;
      //#define BLAHBLAH
#ifdef BLAHBLAH
      verbosity<EvaluationTester>() = VerbExtreme;
      verbosity<Evaluator>() = VerbExtreme;
      verbosity<SparsitySuperset>() = VerbSilent;
      verbosity<EvalVector>() = VerbSilent;
      verbosity<EvaluatableExpr>() = VerbExtreme;
      verbosity<AbstractEvalMediator>() = VerbExtreme;
#endif
      Expr::showAllParens() = true;

      EvalVector::shadowOps() = true;

      ProductTransformation::optimizeFunctionDiffOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr dz = new Derivative(2);

      ADField U(ADBasis(1), sqrt(2.0));
      ADField V(ADBasis(2), sqrt(2.5));
      ADField W(ADBasis(2), sqrt(3.0));

      ADCoord X(0);
      ADCoord Y(1);
      ADCoord Z(2);

      ADDerivative Dx(0);
      ADDerivative Dy(1);
      ADDerivative Dz(2);

      ADReal C_old = sin(X)*sin(Y);

			Expr u = new TestUnknownFunction(U, "u");
			Expr v = new TestUnknownFunction(V, "v");
			Expr w = new TestUnknownFunction(W, "w");

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr g = x*x + y*y;
      Expr f = x*x;
      Expr h = x+y;
      Expr c_old = sin(x)*sin(y);
      double dt = 0.01;

      Expr grad = List(dx, dy);

      double tol1 = 1.0e-5;
      double tol2 = 1.0e-5;
      double fdStep = 1.0e-5;
      bool isOK = true;
      Array<string> failures;



      TESTER(u, U);


      TESTER(-u, -U);

      
      /* ----------- tests of symbolic simplifications -------------*/
      TESTER( u - u, U - U );

      TESTER( u + u, 2.0*U );



      /* ----------- distinct cases of sum expressions ------------------- */

      /* tests const-const and vec-vec sums */
      TESTER( u + w, U + W );

      /* tests const-const and vec-vec subtractions */
      TESTER( u - w, U - W );

      /* tests vec-vec and const-0 sums */
      TESTER( u + x, U + X );

      /* tests vec-vec and const-0 subtractions */
      TESTER( u - x, U - X );

      /* tests vec-vec and const-0 subtractions */
      TESTER( u - (w + v), U - (W + V) );

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
      TESTER( w*(u*u-2.0), W*(U*U-2.0) );

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


      //#endif     

      /* */
      TESTER((dx*u), (Dx*U));

      /* */
      TESTER((dx*x), (Dx*X));

      TESTER((dx*(x*x)), (Dx*(X*X)));

      TESTER((dx*(y*y)), (Dx*(Y*Y)));

      TESTER((dx*(x*x + y*y)), (Dx*(X*X + Y*Y)));

      TESTER((dx*(x*w)), (Dx*(X*W)));

      /* */
      TESTER((dx*y), (Dx*Y));

      TESTER((dx*u+dx*w), (Dx*U+Dx*W));
      
      TESTER((dx*(u+w)), (Dx*(U+W)));

      TESTER((dx*(u-w)), (Dx*(U-W)));

      TESTER((dx*(x+w)), (Dx*(X+W)));

      TESTER((dx*(x-w)), (Dx*(X-W)));


      TESTER((dx*(x*w+w)), (Dx*(X*W+W)));
      //#ifdef BLARF
      TESTER((dx*(x*w-w)), (Dx*(X*W-W)));

      TESTER((dx*(w + x*w)), (Dx*(W + X*W)));


      TESTER((dx*(w - x*w)), (Dx*(W - X*W)));

      //#endif


      //#ifdef BLARF
      TESTER((dx*(y*w)), (Dx*(Y*W)));

      TESTER((dx*(u*w)), (Dx*(U*W)));

      TESTER((dx*(3.0*u*w)), (Dx*(3.0*U*W)));

      TESTER((u*(dx*u)), (U*(Dx*U)));

      TESTER((dx*(u*u)), (Dx*(U*U)));


      TESTER((dx*(h)), (Dx*(X+Y)));
      TESTER((dx*(x*h)), (Dx*(X*(X+Y))));

      TESTER((dx*(h) + dy*(h)), (Dx*(X+Y) + Dy*(X+Y)));

      TESTER((dx*(x*h) + dy*(y*h)), (Dx*(X*(X+Y)) + Dy*(Y*(X+Y))));

      TESTER((dx*(y*h) + dy*(x*h)), (Dx*(Y*(X+Y)) + Dy*(X*(X+Y))));

      TESTER((u*(dx*(f) + dy*(f))), (U*(Dx*(X*X) + Dy*(X*X))));


      TESTER((u*(dx*(x*x+y*y) + dy*(x*x+y*y))), (U*(Dx*(X*X + Y*Y) + Dy*(X*X + Y*Y))));

      TESTER((u*(dx*(g) + dy*(g))), (U*(Dx*(X*X + Y*Y) + Dy*(X*X + Y*Y))));
       


      TESTER((u*g + w*g), (U*(X*X + Y*Y) + W*(X*X + Y*Y)));



      
      TESTER((dx*(2.0*u+4.0*w)), (Dx*(2.0*U+4.0*W)));

      TESTER((dx*(2.0*u-4.0*w)), (Dx*(2.0*U-4.0*W)));

      TESTER((dx*(2.0*x+4.0*w)), (Dx*(2.0*X+4.0*W)));

      TESTER((dx*(2.0*x-4.0*w)), (Dx*(2.0*X-4.0*W)));



      TESTER((dx*(x*w+u)), (Dx*(X*W+U)));

      TESTER((dx*(w + x*w)), (Dx*(W + X*W)));

      TESTER((dx*(x*w+u*y)), (Dx*(X*W+U*Y))); 

      TESTER(((dx*u)+(dx*w) + w), ((Dx*U)+(Dx*W) + W));

      TESTER(((dx*u)*(dx*w) + 2.0*w), ((Dx*U)*(Dx*W) + 2.0*W));

      TESTER(((dx*u)*(dx*w) + (dy*u)*(dy*w)), ((Dx*U)*(Dx*W) + (Dy*U)*(Dy*W)));


      TESTER((x*(dx*u)*(dx*w)), (X*(Dx*U)*(Dx*W))); 

      TESTER((u*(dx*u)*(dx*w)), (U*(Dx*U)*(Dx*W)));

      TESTER((x*(dx*u)*(dx*w)), (X*(Dx*U)*(Dx*W)));

      TESTER((u*(dx*u)*(dx*w)), (U*(Dx*U)*(Dx*W)));

      TESTER((x*(dx*u)*(dx*u)), (X*(Dx*U)*(Dx*U)));

      TESTER((u*(dx*u)*(dx*u)), (U*(Dx*U)*(Dx*U)));

      TESTER((u*(dx*u)), (U*(Dx*U)));

      TESTER((u*(dx*(u*x))), (U*(Dx*(U*X))));

      /* Unary operators */
      TESTER(sin(u), sin(U));

      TESTER(sin(0.5*u), sin(0.5*U));

      TESTER(sin(u+w), sin(U+W));

      TESTER(sin(u*w), sin(U*W));

      TESTER(sin(0.5+u), sin(0.5+U));

      TESTER(sin(x), sin(X));

      TESTER(sin(u*x), sin(U*X));

      TESTER(sin(0.5*u*x), sin(0.5*U*X));

      TESTER(w*sin(u), W*sin(U));

      TESTER(w*sin(u-x), W*sin(U-X));

      TESTER(sin(dx*u), sin(Dx*U));

      TESTER(w*sin(dx*u), W*sin(Dx*U));


      TESTER(sin(u*cos(x)), sin(U*cos(X)));

      TESTER(sin(0.5*(w+u)), sin(0.5*(W+U)));

      TESTER(sin(0.5*w*u), sin(0.5*W*U));

      TESTER(sin(u*x+0.5*w*u), sin(U*X+0.5*W*U));

      TESTER(sin(cos(u*x)+0.5*w*u), sin(cos(U*X)+0.5*W*U));


      TESTER(sin(u)/u, sin(U)/U);

      TESTER(sin(cos(u)), sin((cos(U))));

      TESTER(sin(cos(x)), sin((cos(X))));

      TESTER((dx*sin(x)), (Dx*(sin(X))));

      TESTER((dx*sin(u)), (Dx*(sin(U))));

      TESTER((dx*exp(u)), (Dx*(exp(U))));

      TESTER((dx*exp(u+w)), (Dx*(exp(U+W))));

      TESTER((dx*exp(u*u)), (Dx*(exp(U*U))));

      TESTER((dx*exp(u*w)), (Dx*(exp(U*W))));

      TESTER((dx*exp(2.0*u)), (Dx*(exp(2.0*U))));

      TESTER((dx*exp(-u)), (Dx*(exp(-U))));

      TESTER((dx*exp(-1.0*u)), (Dx*(exp(-1.0*U))));

      TESTER((dx*(cos(x)*sin(x))), (Dx*(cos(X)*sin(X))));

      TESTER((dx*(cos(x)*sin(y))), (Dx*(cos(X)*sin(Y))));

      TESTER((dx*sin(cos(u))), (Dx*(sin(cos(U)))));

      TESTER((dx*sin(2.0*cos(u))), (Dx*(sin(2.0*cos(U)))));

      TESTER(w*(dx*sin(u)), W*(Dx*(sin(U))));

      TESTER(dx*(pow(u, 2.0)), Dx*(pow(U, 2.0)));

      TESTER(dx*(pow(u, 4.0)), Dx*(pow(U, 4.0)));

      TESTER(log(u), log(U));

      TESTER(dx*(log(u)), Dx*(log(U)));

      TESTER(dx*(exp(log(u))), Dx*(exp(log(U))));

      TESTER(dx*(exp(-log(u))), Dx*(exp(-log(U))));

      TESTER(dx*(exp(2.0*log(u))), Dx*(exp(2.0*log(U))));

      TESTER(w*((u-c_old)/dt) + (grad*w)*(grad*(u + c_old)/2.0),
             W*((U-C_old)/dt) + (Dx*W)*(Dx*(U + C_old)/2.0) + (Dy*W)*(Dy*(U + C_old)/2.0));
             


    finish:
      if (isOK)
        {
          cerr << "all tests PASSED!" << endl;
        }
      else
        {
          cerr << "overall test FAILED!" << endl;
          cerr << endl << "failed exprs: " << endl
               << endl;
          for (unsigned int i=0; i<failures.size(); i++)
            {
              cerr << failures[i] << endl;
            }
          cerr << endl;
        }
      TimeMonitor::summarize();
    }
	catch(exception& e)
		{
      cerr << "overall test FAILED!" << endl;
      cerr << "detected exception: " << e.what() << endl;
		}
}
