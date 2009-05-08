#include "SundanceExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace SundanceUtils;
using namespace SundanceCore;
using namespace Teuchos;

bool checkStringForms(const Expr& e1, const Expr& e2)
{
  string s1 = e1[0].toString();
  string s2 = e2[0].toString();

  cerr << "1: " << s1 << endl;
  cerr << "2: " << s2 << endl;

  return s1==s2;
}

int main(int argc, char** argv)
{

  try
		{
      /* initialize MPI */
      GlobalMPISession session(&argc, &argv);

      
      //      verbosity<SymbolicTransformation>() = VerbExtreme;

      Expr::showAllParens() = false;
      EvalVector::shadowOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionStub("u");
			Expr v = new UnknownFunctionStub("v");
			Expr w = new UnknownFunctionStub("w");

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr alpha = new Parameter(1.0, "alpha");
      Expr beta = new Parameter(3.14, "beta");
      Expr zero = new ZeroExpr();

      bool ok = true;

      ok = ok && checkStringForms( u + zero, u );
      
      ok = ok && checkStringForms( zero + u, u );

      ok = ok && checkStringForms( u - zero, u );
      
      ok = ok && checkStringForms( -zero + u, u );

      ok = ok && checkStringForms( -x + v, -(x-v));

      ok = ok && checkStringForms( u + alpha, alpha + u );

      ok = ok && checkStringForms( u - alpha, -alpha + u );

      ok = ok && checkStringForms( u+alpha + beta, beta + alpha + u);

      ok = ok && checkStringForms( (u-alpha) + beta, (beta + -alpha) + u);

      ok = ok && checkStringForms( (u+alpha) - beta, (-beta + alpha) + u);

      ok = ok && checkStringForms( (u-alpha) - beta, (-beta + -alpha) + u);

      ok = ok && checkStringForms( (u+alpha) + v, alpha + (u+v));

      ok = ok && checkStringForms( alpha + (beta + u), (alpha + beta) + u );

      ok = ok && checkStringForms( u + (alpha + v), alpha + (u + v) );

      ok = ok && checkStringForms( u * alpha, alpha * u );

      ok = ok && checkStringForms( u * 0.0, 0.0);

      ok = ok && checkStringForms( 0.0 * u, 0.0);

      ok = ok && checkStringForms( u * 1.0, u);

      ok = ok && checkStringForms( 1.0 * u, u);

      ok = ok && checkStringForms( -1.0 * u, -u);

      ok = ok && checkStringForms( -(-u), u);

      ok = ok && checkStringForms( u * 1.0 + u*0.0, u);

      ok = ok && checkStringForms( 1.0 * u + u*0.0, u);

      ok = ok && checkStringForms( u * 1.0 + 0.0*u, u);

      ok = ok && checkStringForms( 1.0 * u + 0.0*u, u);

      ok = ok && checkStringForms( u * (alpha * v), alpha * (u * v) );

      ok = ok && checkStringForms( alpha * (beta * u), (alpha * beta) * u);

      ok = ok && checkStringForms( (alpha * u) * beta, (beta * alpha) * u);

      ok = ok && checkStringForms( (alpha * u) * v, alpha * (u * v) );

      ok = ok && checkStringForms( (-u)*(-v), u*v);

      ok = ok && checkStringForms( u*(-v), -u*v);

      ok = ok && checkStringForms( (-u)*v, -u*v);

      ok = ok && checkStringForms( x + (-u)*v, x-u*v);

      ok = ok && checkStringForms( -x + (-u)*v, -(x+u*v));

      ok = ok && checkStringForms( -x + v, v - x);

      ok = ok && checkStringForms( u*(-1.0*v), -u*v);

      ok = ok && checkStringForms( (-1.0*u)*v, -u*v);

      ok = ok && checkStringForms( dx * (alpha * u), alpha * (dx * u));

      ok = ok && checkStringForms( dx * (alpha + u), dx * u);

      ok = ok && checkStringForms( dx * alpha, zero );

      ok = ok && checkStringForms( (dx + dy)*u, dx*u + dy*u );

      ok = ok && checkStringForms( (u*dx)*v, u*(dx*v) );




      if (ok) cerr << "All are OK" << endl;
      else cerr << "failures detected!" << endl;

      TimeMonitor::summarize();
    }
	catch(exception& e)
		{
			Out::println(e.what());
		}


  return 0;
}
