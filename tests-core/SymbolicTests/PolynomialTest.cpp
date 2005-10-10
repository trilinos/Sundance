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
#include "SundanceFunctionalPolynomial.hpp"

#include "SundanceSymbolicTransformation.hpp"
#include "SundanceProductTransformation.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceEquationSet.hpp"

using SundanceCore::List;
using namespace SundanceCore;
using namespace SundanceUtils;
using namespace SundanceCore::Internal;
using namespace Teuchos;
using namespace TSFExtended;

RefCountPtr<ScalarExpr> expr2scalar(const Expr& e)
{
  return rcp_dynamic_cast<ScalarExpr>(e[0].ptr());
}

RefCountPtr<FunctionalPolynomial> expr2poly(const Expr& e)
{
  return FunctionalPolynomial::toPoly(expr2scalar(e));
}

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      Expr::showAllParens() = true; 
      ProductTransformation::optimizeFunctionDiffOps() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

      Expr u = new UnknownFunctionStub("u");
      Expr v = new UnknownFunctionStub("v");
      Expr w = new UnknownFunctionStub("w");
      
      RefCountPtr<FunctionalPolynomial> p 
        = rcp(new FunctionalPolynomial(expr2scalar(u)));
      
      RefCountPtr<FunctionalPolynomial> q 
        = rcp(new FunctionalPolynomial(expr2scalar(w)));
      
      RefCountPtr<FunctionalPolynomial> r
        = rcp(new FunctionalPolynomial(expr2scalar(w)));

      p = p->multiplyPoly(expr2poly(v).get());
      p = p->multiplyPoly(expr2poly(v).get());
      p = p->multiplyPoly(expr2poly(u).get());

      q = q->multiplyPoly(expr2poly(u).get());
      q = q->multiplyPoly(expr2poly(v).get());

      p = p->addPoly(q.get(), 1);

      p = p->addPoly(r.get(), 1);

          
      cerr << "p = " << p->toXML() << endl;
      p->toText(cerr, false) << endl;
      
    }
	catch(exception& e)
		{
      cerr << e.what() << endl;
		}
  TimeMonitor::summarize();
  MPISession::finalize();
}
