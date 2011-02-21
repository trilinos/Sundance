#include "Sundance.hpp"


// written by Michael Ulbrich, TU Muenchen, Dec 30, 2010

int main(int argc, char** argv)
{
  
  try
  {
    Sundance::init(&argc, &argv);

    /* We will do our linear algebra using Epetra */
    VectorType<double> vecType = new EpetraVectorType();

    MeshType meshType = new BasicSimplicialMeshType();
    int nx = 1;

    MeshSource mesher = new PartitionedLineMesher(0.0, 1.0, nx, meshType);
    Mesh mesh = mesher.getMesh();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter(); 

    /* Create unknown function */
    BasisFamily L1=new Lagrange(1);
    Expr u = new UnknownFunction(L1, "u");
    Expr w = new TestFunction(L1, "w");

    Expr dx = new Derivative(0);
      
    Expr x = new CoordExpr(0);

    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad = new GaussianQuadrature(4);

    DiscreteSpace discSpaceL1(mesh, L1, vecType);
    Expr ud = sqrt(2.0)*x; //proj.project();
    
    WatchFlag watch("watch");
    watch.setParam("evaluation", 5);
    watch.setParam("evaluator setup", 5);
    watch.setParam("discrete function evaluation", 1);
    watch.setParam("integration setup", 0);
    watch.setParam("symbolic preprocessing", 0);
    watch.setParam("integration", 0);

    Expr eqn1 = Integral(interior, w*(dx*(u+ud)), quad, watch);
    Expr eqn2 = Integral(interior, w*(dx*u)+w*(dx*ud), 
      quad, watch);
    Expr bc;

    Out::root() << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" 
                << " @@@@@@@@@@@@@@@ creating BAD operator @@@@@@@@@@@@@\n" 
                << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" 
                << endl;
    
    LinearProblem prob1(mesh, eqn1, bc, w, u, vecType);
    LinearOperator<double> A1 = prob1.getOperator();

    
    Out::root() << "BAD operator = " << endl << A1 << endl << endl << endl;


    Out::root() << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" 
                << " @@@@@@@@@@@@@@@ creating GOOD operator @@@@@@@@@@@@\n" 
                << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
                << endl;
    LinearProblem prob2(mesh, eqn2, bc, w, u, vecType);
    LinearOperator<double> A2 = prob2.getOperator();
    Out::root() << "GOOD operator = " << endl << A2 << endl;

    Vector<double> b = A2.domain().createMember();
    b.randomize();

    Vector<double> r = 2.0*(A2*b) - A1*b;

    Out::root() << "difference in operator application = " 
                << r.norm2() << endl;
    

  }
	catch(std::exception& e)
  {
    std::cerr << e.what() << std::endl;
  }
  Sundance::finalize();
//  return Sundance::testStatus(); 
  return(0);
}
