/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "Sundance.hpp"
#include "TSFExplicitlyTransposedLTIProblemFactory.hpp"


class ThermalLTIProblemFactory 
  : public ExplicitlyTransposedLTIProblemFactory<double>
{
public:
  /** */
  ThermalLTIProblemFactory(int nSteps, double deltaT);

  /** */
  const DiscreteSpace& solutionSpace() const 
    {return *(timestepProb_.solnSpace()[0]);}

  /** */
  const Mesh& mesh() const {return mesh_;}
private:
  Mesh mesh_;
  LinearProblem timestepProb_;
  LinearProblem massProb_;

  /** */
  void setup(double deltaT);

  
};

int main(int argc, char** argv)
{
  try
  {
    Sundance::init(&argc, &argv);

    int nSteps = 10;
    const double pi = 4.0*atan(1.0);
    double dt = 0.005/(2.0*pi*pi);

    RefCountPtr<ThermalLTIProblemFactory> pf 
      = rcp(new ThermalLTIProblemFactory(nSteps, dt));

    const DiscreteSpace& discSpace = pf->solutionSpace();

    /* Set up a function for the initial state */

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);
    Expr f = sin(pi*x)*sin(pi*y) + 0.5*sin(2.0*pi*x)*sin(3.0*pi*y)
      + 0.3*sin(3.0*pi*x)*sin(4.0*pi*y);

    /* Discretize the initial state */
    L2Projector proj(discSpace, f);
    Expr u0 = proj.project();

    /* Extract the state vector from discrete function form of 
     * the state variable */
    Vector<double> X0 = DiscreteFunction::discFunc(u0)->getVector();

    /* Get the forward timestepping operator */
    LinearOperator<double> bigAInv = pf->getBigAInv();

    /* Create the zero-padded operand for the timestepping operator */
    LinearOperator<double> bigF = pf->getBigF();
    Vector<double> X = bigF.range().createMember();
    X.zero();
    X.setBlock(0, X0);

    /* Do the timestepping. The whole process is encapsulated by 
     * the application of the timestepping operator bigAInv. */
    Vector<double> U = bigAInv * X;

    /* Write the time history */
    bool writeIt = false;
    if (writeIt)
    {
      for (int i=0; i<U.space().numBlocks(); i++)
      {
        Expr uOut = new DiscreteFunction(discSpace, U.getBlock(i));
        FieldWriter w = new VTKWriter("ThermalLTI-" + Teuchos::toString(i));
        w.addMesh(pf->mesh());
        w.addField("soln", new ExprFieldWrapper(uOut[0]));
        w.write();
      }
    }
    

    /* ------ Now apply the Hessian, just to show how it's done ------- */

    /* Get the Hessian */
    LinearOperator<double> H = pf->getH();

    /* Multiply X0 by F^T in order to put the first block of X0 into
     * a single-block block vector. This is needed for space compatibility. */
    Vector<double> x0 = pf->getBigF().domain().createMember();
    Thyra::randomize(-1.0, 1.0, x0.ptr().get());

    /* Do the Hessian application */
    /* Note: this doesn't work yet. I need to get transpose solves working
     * with the aztec adapters. */

    cerr << "applying Hessian" << endl;
    Vector<double> Hx = H * x0;
    
  }
	catch(exception& e)
  {
    Sundance::handleException(e);
  }
  Sundance::finalize();
}



ThermalLTIProblemFactory::ThermalLTIProblemFactory(int nSteps, double deltaT) 
  : ExplicitlyTransposedLTIProblemFactory<double>(nSteps),
    mesh_(),
    timestepProb_(),
    massProb_()
{
  setup(deltaT);
}



void ThermalLTIProblemFactory::setup(double deltaT)
{
  int nx = 64;
  int ny = 64;
  int np = MPIComm::world().getNProc();

  /* We will do our linear algebra using Epetra */
  VectorType<double> vecType = new EpetraVectorType();

  /* Create a mesh. It will be of type BasisSimplicialMesh, and will
   * be built using a PartitionedRectangleMesher. */
  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, nx, np,
    0.0, 1.0, ny, 1,
    meshType);
  mesh_ = mesher.getMesh();

  /* Create a cell filter that will identify the maximal cells
   * in the interior of the domain */
  CellFilter interior = new MaximalCellFilter();
  CellFilter bdry = new BoundaryCellFilter();

  /* Create unknown and test functions, discretized using first-order
   * Lagrange interpolants */
  BasisFamily basis = new Lagrange(1);
  Expr u = new UnknownFunction(basis, "u");
  Expr v = new TestFunction(basis, "v");

  /* Create differential operator and coordinate functions */
  Expr dx = new Derivative(0);
  Expr dy = new Derivative(1);
  Expr grad = List(dx, dy);

  /* We need a quadrature rule for doing the integrations */
  QuadratureFamily quad2 = new GaussianQuadrature(2);

  /* Set up the backward Euler timestep operator (M - deltaT*K) */
  Expr timestepEqn 
    = Integral(interior, u*v + deltaT*(grad*v)*(grad*u), quad2);
  Expr timestepBC = EssentialBC(bdry, u*v, quad2);
      
  timestepProb_ = LinearProblem(mesh_, timestepEqn, timestepBC, 
    v, u, vecType);

  /* Set up a problem to for the BC-corrected mass matrix */
  Expr massEqn 
    = Integral(interior, u*v, quad2);
  Expr massBC = EssentialBC(bdry, u*v, quad2);

  massProb_ =  LinearProblem(mesh_, massEqn, massBC, 
    v, u, vecType);


      
  string solverFile = "aztec.xml";
  ParameterXMLFileReader reader(searchForFile("SolverParameters/" + solverFile));
  ParameterList solverParams = reader.getParameters();
  LinearSolver<double> solver 
    = LinearSolverBuilder::createSolver(solverParams);

  /* BE is (M - deltaT * K), the matrix that is solved upon each Backwards
   * Euler timestep */
  LinearOperator<double> BE = timestepProb_.getOperator();
  /* Lacking a transpose solve in AztecOO, we form the transpose of BE
   * explicitly */
  LinearOperator<double> BEt = formExplicitTranspose(BE);

  /* M is the mass matrix */
  LinearOperator<double> M = massProb_.getOperator();

  /* Form A = BE*M */
  LinearOperator<double> A = BE.inverse(solver) * M;
  /* Form A^T = M^T * BE^{-T}, where we have used the explicit
   * transpose of BE.  */
  LinearOperator<double> At = M.transpose() * BEt.inverse(solver);

  LinearOperator<double> C = identityOperator(BE.range());

  this->init(A, At, C);
}

