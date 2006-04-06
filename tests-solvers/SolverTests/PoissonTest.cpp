//@HEADER
// ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "Teuchos_MPISession.hpp"
#include "TSFVector.hpp"
#include "TSFLinearCombination.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "TSFLinearSolver.hpp"
#include "TSFAztecSolver.hpp"
#include "TSFMatrixLaplacian1D.hpp"

using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;


int main(int argc, void *argv[]) 
{
  typedef Teuchos::ScalarTraits<double> ST;

  try
    {
      int verbosity = 1;

      MPISession::init(&argc, &argv);

      MPIComm::world().synchronize();

      VectorType<double> type = new EpetraVectorType();

      /* create the range space  */
      int nLocalRows = 10;

      MatrixLaplacian1D builder(nLocalRows, type);

      LinearOperator<double> A = builder.getOp();

      Vector<double> x = A.domain().createMember();
      Thyra::randomize(-ST::one(),+ST::one(),x.ptr().get());

      Vector<double> y = A*x;
      Vector<double> ans = A.range().createMember();

      cerr << "x=" << x << endl;
      cerr << "y=" << y << endl;

      ParameterList params;
      params.set("Method", "GMRES");
      params.set("Precond", "ML");
      params.set("ML Levels", 2);
      params.set("Max Iterations", 100);
      params.set("Tolerance", 1.0e-12);

      LinearSolver<double> solver = new AztecSolver(params);

      SolverState<double> state = solver.solve(A, y, ans);
      
      cerr << state << endl;

      cerr << "solver is " << solver << endl;

      double err = (x-ans).norm2();
      cerr << "error norm = " << err << endl;

      double tol = 1.0e-10;
      if (err > tol)
        {
          cerr << "Poisson solve test FAILED" << endl;
        }
      else
        {
          cerr << "Poisson solve test PASSED" << endl;
        }
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();
}
