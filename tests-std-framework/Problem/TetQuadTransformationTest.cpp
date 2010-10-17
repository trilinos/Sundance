
#include "Sundance.hpp"

using Sundance::List;




#if defined(HAVE_SUNDANCE_EXODUS) && defined(Trilinos_DATA_DIR)

int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);

    VectorType<double> vecType = new EpetraVectorType();

    MeshType meshType = new BasicSimplicialMeshType();
      
    string meshFile = "tetPrism";
    MeshSource meshSrc
      =  new ExodusMeshReader(meshFile, meshType);
    Mesh mesh = meshSrc.getMesh();

    CellFilter interior = new MaximalCellFilter();

    CellFilter faces = new DimensionalCellFilter(2);
    CellFilter top = faces.labeledSubset(1);
    CellFilter bottom = faces.labeledSubset(2);



    BasisFamily basis = new Lagrange(2);
    Expr T = new UnknownFunction(basis, "T");
    Expr vT = new TestFunction(basis, "vT");

    Expr z = new CoordExpr(2);

    QuadratureFamily quad4 = new GaussianQuadrature(4);

    WatchFlag watch1("interior");
    WatchFlag watch2("surface");
    //    watchMe.setParam("integration setup", 5);
    watch1.setParam("fill", 5);
    watch2.setParam("fill", 5);
    watch1.setParam("assembly loop", 5);
    watch2.setParam("assembly loop", 5);

    watch2.setParam("integration", 5);
    watch2.setParam("integral transformation", 5);

    watch1.deactivate();

    Expr exact = z;

    Expr eqn1 = Integral(interior, vT*(T-z), quad4, watch1);
    Expr bc1 = EssentialBC(bottom, vT*(T-z), quad4, watch2);

    Expr eqn2 = Integral(interior, vT*(T-z), quad4, watch1)
      + Integral(bottom, vT*(T-z), quad4, watch2);
    Expr bc2;

    Expr eqn3 = Integral(interior, vT*(T-z), quad4, watch1);
    Expr bc3;

    Array<LinearProblem> prob(3);
    Out::os() << "=============== Problem 1 ctor" << endl;
    prob[0] = LinearProblem(mesh, eqn1, bc1, vT, T, vecType);
    Out::os() << "=============== Problem 2 ctor" << endl;
    prob[1] = LinearProblem(mesh, eqn2, bc2, vT, T, vecType);
    Out::os() << "=============== Problem 3 ctor" << endl;
    prob[2] = LinearProblem(mesh, eqn3, bc3, vT, T, vecType);

    FieldWriter w1 = new VerboseFieldWriter("dump");
    w1.addMesh(mesh);
    w1.write();

    Array<LinearOperator<double> > A(3);
    Array<Vector<double> > b(3);

    LinearSolver<double> solver 
      = LinearSolverBuilder::createSolver("amesos.xml");

    Array<Expr> soln(3);
    
    for (int i=0; i<3; i++)
    {
      Out::os() << "=============== Problem #" << i+1 << " =========="
                << endl;

      string matFile = "mtx-" + Teuchos::toString(i+1) + ".dat";
      string vecFile = "vec-" + Teuchos::toString(i+1) + ".dat";
      string mapFile = "map-" + Teuchos::toString(i+1) + ".dat";
      ofstream maps(mapFile.c_str());
      ofstream vecs(vecFile.c_str());
      ofstream mats(matFile.c_str());
      prob[i].rowMap(0)->print(maps);

      A[i] = prob[i].getOperator();
      b[i] = prob[i].getSingleRHS();
      mats << "matrix = " << endl << A[i] << endl;
      vecs << "vector = " << endl << b[i] << endl;

      soln[i] = prob[i].solve(solver);
      Out::os() << "=============== Done Problem #" << i+1 << " =========="
                << endl;


    }

    Expr err0 = Integral(interior, 
      pow(soln[0]-exact, 2.0),
      new GaussianQuadrature(4));

    Expr err1 = Integral(interior, 
      pow(soln[1]-exact, 2.0),
      new GaussianQuadrature(4)); 

    Expr err2 = Integral(interior, 
      pow(soln[2]-exact, 2.0),
      new GaussianQuadrature(4));

    double error0Sq = evaluateIntegral(mesh, err0);
    double error1Sq = evaluateIntegral(mesh, err1);
    double error2Sq = evaluateIntegral(mesh, err2);

    std::cerr << "prob 0 error norm = " << sqrt(error0Sq) << std::endl << std::endl;
    std::cerr << "prob 1 error norm = " << sqrt(error1Sq) << std::endl << std::endl;
    std::cerr << "prob 2 error norm = " << sqrt(error2Sq) << std::endl << std::endl;

    FieldWriter w = new VTKWriter("tetPrism");
    w.addMesh(mesh);
    w.addField("bc", new ExprFieldWrapper(soln[0]));
    w.addField("penalty", new ExprFieldWrapper(soln[1]));
    w.addField("interior", new ExprFieldWrapper(soln[2]));
    w.write();
    
    double tol = 1.0e-12;
    double errorSq = error0Sq + error1Sq + error2Sq;
    Sundance::passFailTest(sqrt(errorSq), tol);
  }
  catch(exception& e) 
  {
    Sundance::handleException(e);
  }
  Sundance::finalize(); return Sundance::testStatus(); 

}


#else

int main(int argc, char** argv)
{
  Sundance::init(&argc, &argv);
  std::cout << "dummy TetQuadTransformationTest PASSED. Enable exodus to run the actual test" <<
 std::endl;
  Sundance::finalize();
  return 0;
}

#endif
