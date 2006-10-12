#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "Sundance.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"

int main(int argc, void** argv)
{
  try
    {
      Sundance::init(&argc, &argv);
      
      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      Expr::showAllParens() = true;

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      int ndim = 2;
      int order = 2;

      SpectralBasis SB(ndim, order);
     
      Array<Expr> Ex1(6);
      Ex1[0] = 1.0;
      Ex1[1] = 1.0;
      Ex1[2] = 0.0;
      Ex1[3] = 1.0;
      Ex1[4] = 0.0;
      Ex1[5] = 0.0;

      Expr SE1 = new SpectralExpr(SB, Ex1);


      Array<Expr> Ex2(6);
      Ex2[0] = -1.0;
      Ex2[1] = 0.0;
      Ex2[2] = -1.0;
      Ex2[3] = 0.0;
      Ex2[4] = -1.0;
      Ex2[5] = 0.0;

      Expr SE2 = new SpectralExpr(SB, Ex2);

      Expr G = x*x;

      Expr Sum  = SE1-SE2;

      const SpectralExpr* se = dynamic_cast<const SpectralExpr*>(Sum.ptr().get());
      SpectralBasis basis = se->getSpectralBasis();

      for(int i=0; i< basis.nterms(); i++)
	cout << se->getCoeff(i) << endl;


      Expr J = Sum + SE1; 
      cout << J << endl; 

    }
  
  catch(exception& e)
    {
      Sundance::handleException(e);
    }
  Sundance::finalize();
}

