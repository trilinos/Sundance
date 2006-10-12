#include "SundanceExpr.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceSpectralBasis.hpp"
#include "SundanceSpectralExpr.hpp"
#include "SundanceUnknownFunctionStub.hpp"

int main(int argc, void** argv)
{
  try
    {
      MPISession::init(&argc, &argv);
      

      Expr::showAllParens() = true;

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr z = new CoordExpr(2);

      int ndim = 2;
      int order = 2;

      SpectralBasis SB(ndim, order);

      Expr u = new UnknownFunctionStub("u", SB);
      Expr w = new UnknownFunctionStub("w", SB);
     
      Array<Expr> Ex1(6);
      Ex1[0] = 1.0;
      Ex1[1] = x;
      Ex1[2] = 0.0;
      Ex1[3] = x*y;
      Ex1[4] = x+y;
      Ex1[5] = y;

      Expr SE1 = new SpectralExpr(SB, Ex1);


      Array<Expr> Ex2(6);
      Ex2[0] = -3.0*x;
      Ex2[1] = 0.0;
      Ex2[2] = -y;
      Ex2[3] = x-y;
      Ex2[4] = -4.0*y + 2.0*x;
      Ex2[5] = 0.0;

      Expr SE2 = new SpectralExpr(SB, Ex2);

      Expr G = x*x;

      Expr Sum  = w * u * z ;

      const SpectralExpr* se = dynamic_cast<const SpectralExpr*>(Sum.ptr().get());
      SpectralBasis basis = se->getSpectralBasis();

      for(int i=0; i< basis.nterms(); i++)
        cout << se->getCoeff(i) << endl;


      Expr J = Sum + SE1; 
      cout << J << endl; 

    }
  
  catch(exception& e)
    {
      cerr << e.what() << endl;
    }
  MPISession::finalize();
}

