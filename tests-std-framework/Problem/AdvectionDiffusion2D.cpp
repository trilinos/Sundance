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

//
// Send a TSFExtended::LinearOperator<double> to a file!
//

#include "EpetraExt_OperatorOut.h"
#include "TSFEpetraMatrix.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace {

void printLinearOperatorType(
  const TSFExtended::LinearOperator<double> &A,
  std::ostream &out
  )
{
  out
    << "\nConcrete type of underlying LinearOperator is: "
    << typeid(*A.ptr()).name() << "\n";
}

Teuchos::RefCountPtr<const Epetra_Operator>
get_Epetra_Operator( const TSFExtended::LinearOperator<double> &A )
{
  // Dynamic cast to an interface where you can extract an Epetra_Operator
  // view of the underlying linear operator.  Note that
  // TSFExtended::EpetraMatrix derives from this interface when
  // HAVE_THYRA_EPETRA is defined.  The actually matrix that
  // TSFExtended::EpetraMatrix stores is an Epetra_CrsMatrix object so we
  // actually loose information when we call this function.  However, a simple
  // dynamic cast will get this back again.
  Teuchos::RefCountPtr<const Thyra::EpetraLinearOpBase>
    teA = Teuchos::rcp_dynamic_cast<const Thyra::EpetraLinearOpBase>(
      A.ptr(), true
      );
  Teuchos::RefCountPtr<const Epetra_Operator> epetraOp;
  Thyra::ETransp epetraOpTransp;
  Thyra::EApplyEpetraOpAs epetraOpApplyAs;
  Thyra::EAdjointEpetraOp epetraOpAdjointSupport;
  teA->getEpetraOpView(
    &epetraOp, &epetraOpTransp, &epetraOpApplyAs, &epetraOpAdjointSupport );
  // Note, above the last argument 
  TEST_FOR_EXCEPT(is_null(epetraOp));
  TEST_FOR_EXCEPT(epetraOpTransp != Thyra::NOTRANS);
  TEST_FOR_EXCEPT(epetraOpApplyAs != Thyra::EPETRA_OP_APPLY_APPLY);
  TEST_FOR_EXCEPT(epetraOpAdjointSupport != Thyra::EPETRA_OP_ADJOINT_SUPPORTED);
  // If we get here
  return epetraOp;
}

void writeLinearOperatorMatrixMarketFile(
  const std::string &fileName,
  const TSFExtended::LinearOperator<double> &A,
  const std::string &matrixName = "",
  const std::string &matrixDescription = "",
  bool writeHeader = 0
  )
{
  printLinearOperatorType(A,std::cout);
  TEST_FOR_EXCEPT(0==fileName.length());
  EpetraExt::OperatorToMatrixMarketFile(
    fileName.c_str(),
    *Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(
      get_Epetra_Operator(A)
      ),
    matrixName.length() ? matrixName.c_str() : 0,
    matrixDescription.length() ? matrixDescription.c_str() : 0,
    writeHeader
    );
}

} // namespace

/** 
 * Solves the advection-diffusion equation in 2D, with a velocity
 * field computed from a potential flow model.
 */

CELL_PREDICATE(LeftPointTest, {return fabs(x[0]) < 1.0e-10;})
CELL_PREDICATE(BottomPointTest, {return fabs(x[1]) < 1.0e-10;})
CELL_PREDICATE(RightPointTest, {return fabs(x[0]-1.0) < 1.0e-10;})
CELL_PREDICATE(TopPointTest, {return fabs(x[1]-1.0) < 1.0e-10;})



int main(int argc, char** argv)
{
  
  try
		{
      Sundance::init(&argc, &argv);
      int np = MPIComm::world().getNProc();

      /* We will do our linear algebra using Epetra */
      VectorType<double> vecType = new EpetraVectorType();

      /* Create a mesh. It will be of type BasicSimplicialMesh, and will
       * be built using a PartitionedRectangleMesher. */
      int n = 1;
      MeshType meshType = new BasicSimplicialMeshType();
      MeshSource mesher = new PartitionedRectangleMesher(0.0, 1.0, n, np,
                                                         0.0, 1.0, n, 1,
                                                         meshType);
      Mesh mesh = mesher.getMesh();

      /* Create a cell filter that will identify the maximal cells
       * in the interior of the domain */
      CellFilter interior = new MaximalCellFilter();
      CellFilter edges = new DimensionalCellFilter(1);

      CellFilter left = edges.subset(new LeftPointTest());
      CellFilter right = edges.subset(new RightPointTest());
      CellFilter top = edges.subset(new TopPointTest());
      CellFilter bottom = edges.subset(new BottomPointTest());


      
      /* Create unknown and test functions, discretized using first-order
       * Lagrange interpolants */
      int order = 2;
      Expr u = new UnknownFunction(new Lagrange(order), "u");
      Expr v = new TestFunction(new Lagrange(order), "v");

      /* Create differential operator and coordinate functions */
      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);
      Expr grad = List(dx, dy);
      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      /* We need a quadrature rule for doing the integrations */
      QuadratureFamily quad2 = new GaussianQuadrature(2);
      QuadratureFamily quad4 = new GaussianQuadrature(4);

      /* Define the weak form for the potential flow equation */
      Expr flowEqn = Integral(interior, (grad*v)*(grad*u), quad2);

      /* Define the Dirichlet BC */
      Expr flowBC = EssentialBC(bottom, v*(u-0.5*x*x), quad4)
        + EssentialBC(top, v*(u - 0.5*(x*x - 1.0)), quad4)
        + EssentialBC(left, v*(u + 0.5*y*y), quad4)
        + EssentialBC(right, v*(u - 0.5*(1.0-y*y)), quad4);

      /* We can now set up the linear problem! */
      LinearProblem flowProb(mesh, flowEqn, flowBC, v, u, vecType);

      ParameterXMLFileReader reader(searchForFile("SolverParameters/bicgstab.xml"));
      ParameterList solverParams = reader.getParameters();
      cerr << "params = " << solverParams << endl;
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);

      /* solve the problem */
      Expr u0 = flowProb.solve(solver);

      /* Write the operator to a matrix market file */
      writeLinearOperatorMatrixMarketFile(
        "AdvectionDiffusion2D.mtx",flowProb.getOperator(),
        "AdvectionDiffusion2D","Just to show how to do this"
        ,true
        );

      /* Now set up and solve the advection-diffusion equation for r */
      Expr r = new UnknownFunction(new Lagrange(order), "u");
      Expr s = new TestFunction(new Lagrange(order), "v");

      Expr velocity = grad*u0;
      //      Expr velocity = grad*(0.5*(x*x - y*y));
      Expr adEqn = Integral(interior, (grad*s)*(grad*r), quad2)
        + Integral(interior, s*velocity*(grad*r), quad4);
        

      Expr adBC = EssentialBC(bottom, s*r, quad4)
        + EssentialBC(top, s*(r-x), quad4)
        + EssentialBC(left, s*r, quad4)
        + EssentialBC(right, s*(r-y), quad4);

      LinearProblem adProb(mesh, adEqn, adBC, s, r, vecType);
      Expr r0 = adProb.solve(solver);

      FieldWriter w = new VTKWriter("AD-2D");
      w.addMesh(mesh);
      w.addField("potential", new ExprFieldWrapper(u0[0]));
      w.addField("concentration", new ExprFieldWrapper(r0[0]));
      w.write();

      Expr exactPotential = 0.5*(x*x - y*y);
      Expr exactConcentration = x*y;

      Expr uErr = exactPotential - u0;
      Expr uErrExpr = Integral(interior, 
                              uErr*uErr,
                              quad4);

      Expr rErr = exactConcentration-r0;
      Expr rErrExpr = Integral(interior, 
                               rErr*rErr, 
                               quad4);

      FunctionalEvaluator uInt(mesh, uErrExpr);
      FunctionalEvaluator rInt(mesh, rErrExpr);

      double uErrorSq = uInt.evaluate();
      cerr << "potential error norm = " << sqrt(uErrorSq) << endl << endl;

      double rErrorSq = rInt.evaluate();
      cerr << "concentration error norm = " << sqrt(rErrorSq) << endl << endl;

      Sundance::passFailTest(uErrorSq + rErrorSq, 1.0e-11);

    }
	catch(exception& e)
		{
      Sundance::handleException(e);
		}
  Sundance::finalize();
}
