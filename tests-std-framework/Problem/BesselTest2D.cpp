#include "Sundance.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceUserDefOp.hpp"

/** 
 * 
 */

class BesselJFunc : public UserDefFunctor
{
public: 
  /** */
  BesselJFunc(int n) : UserDefFunctor("J_" + Teuchos::toString(n)), n_(n) {;}

  /** */
  virtual double eval0(const Array<double>& vars) const ;

private:
  int n_;
};

Expr BesselJ(int n, const Expr& r)
{
  return  new UserDefOp(r, rcp(new BesselJFunc(n)));
}


/** 
 * To plug in a user-defined function, you write a functor object,
 * deriving from UserDefFunctor, and implement your function in 
 * the eval() method of the functor. Any constant parameters (i.e., 
 * NOT expr objects) to your function, you can pass in as ctor 
 * arguments. 
 *
 * A user-defined functor will be accessed through a UserDefOp
 * object, whose ctor takes Sundance expression specifying the 
 * arguments of the function. This expression can be list-valued
 * in the event that your function is of more than one variable.
 *
 * The eval() method evaluates the function at a single point
 * in space. The argument to eval() is an array of doubles, the elements
 * of which are the values of the user-level argument Expr at
 * the current evaluation point. 
 */
class DrumFuncFunctor : public UserDefFunctor
{
public:
  /** */
  DrumFuncFunctor(const double& A, int n) 
    : UserDefFunctor("Drum_" + Teuchos::toString(n)), A_(A), n_(n) {;}

  /** */
  virtual double eval0(const Array<double>& vars) const ;

private:
  double A_;
  int n_;
};


/* Definition of eval() function. This is where the "meat" is. */

double DrumFuncFunctor::eval0(const Array<double>& vars) const
{
  double x = vars[0];
  double y = vars[1];

  double theta = atan2(y,x);
  double r = sqrt(x*x + y*y);

  return A_*cos(n_*theta)*jn(n_, r);
}



/* A little wrapper to be able to build your user-defined function
* in pretty syntax. */
Expr DrumFunc(int n, double A, const Expr& x, const Expr& y)
{
  return new UserDefOp(List(x,y), rcp(new DrumFuncFunctor(A, n)));
}





int main(int argc, void** argv)
{
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasisSimplicialMesh, and will
       * be built using a PartitionedLineMesher. */
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 10.0, 64, 1,
                                                         0.0, 10.0, 64, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();


      DiscreteSpace discSpace(mesh, new Lagrange(2), vecType);

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);
      Expr r = sqrt(x*x + y*y);

      int n = 2;
      Expr Jn = BesselJ(n, r);

      Expr Dn = DrumFunc(n, 1.0, x, y);

      Expr JnDisc = L2Projector(discSpace, Jn).project();
      Expr DnDisc = L2Projector(discSpace, Dn).project();

       /* Write the field in VTK format */
      FieldWriter w = new VTKWriter("Bessel2d");
      w.addMesh(mesh);
      w.addField("J_n", new ExprFieldWrapper(JnDisc));
      w.addField("Drum_n", new ExprFieldWrapper(DnDisc));
      w.write();

      double tol = 1.0e-12;
      Sundance::passFailTest(0, tol);
    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}



/* Definition of BesselJFunc::eval() */

double BesselJFunc::eval0(const Array<double>& vars) const
{
  return jn(n_, vars[0]);
}
