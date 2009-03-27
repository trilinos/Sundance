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

#include "Teuchos_GlobalMPISession.hpp"
#include "TSFPoissonBoltzmannOp.hpp"
#include "TSFEpetraVectorType.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "SundancePathUtils.hpp"
#include "TSFNOXSolver.H"

using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;




int main(int argc, char *argv[]) 
{
  try
    {
      GlobalMPISession session(&argc, &argv);


      MPIComm::world().synchronize();

      /* create the nonlinear operator */
      VectorType<double> type = new EpetraVectorType();
      int nProc = MPIComm::world().getNProc();
      int nLocalRows = 128/nProc;
      PoissonBoltzmannOp* prob = new PoissonBoltzmannOp(nLocalRows, type);
      NonlinearOperator<double> F = prob;

      /* create the nox solver */

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(SundanceUtils::searchForFile("SolverParameters/nox.xml"));
#else
      ParameterXMLFileReader reader("nox.xml");
#endif

      ParameterList noxParams = reader.getParameters();

      cerr << "solver params = " << noxParams << endl;

      NOXSolver solver(noxParams, F);

      NOX::StatusTest::StatusType stat = solver.solve();
      TEST_FOR_EXCEPTION(stat != NOX::StatusTest::Converged,
        runtime_error, "solve failed");

      Vector<double> soln = solver.getSolution();

      cerr << "solution = " << endl << soln << endl;

      Vector<double> exact = prob->exactSoln();

      cerr << "exact solution = " << endl << exact << endl;

//bvbw reddish port hack
      double temp_val = nLocalRows*nProc;
      double err = (exact-soln).norm2()/sqrt(temp_val);
      cerr << "error norm = " << err << endl;
      

      double tol = 1.0e-6;
      if (err > tol)
        {
          cerr << "NOX Poisson-Boltzmann test FAILED" << endl;
          return 1;
        }
      else
        {
          cerr << "NOX Poisson-Boltzmann test PASSED" << endl;
          return 0;
        }
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
      return -1;
    }
}

