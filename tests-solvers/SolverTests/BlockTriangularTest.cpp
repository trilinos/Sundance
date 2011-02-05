//@HEADER
// ***********************************************************************
// 
//           Playa: Trilinos Solver Framework Extended
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
#include "PlayaVectorDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaEpetraVectorType.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaBlockTriangularSolverDecl.hpp"
#include "PlayaLinearCombinationDecl.hpp"
#include "PlayaBlockVectorSpaceDecl.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaMatrixLaplacian1D.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "SundancePathUtils.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaSimpleBlockOpDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaSimpleBlockOpImpl.hpp"
#include "PlayaBlockVectorSpaceImpl.hpp"
#include "PlayaBlockTriangularSolverImpl.hpp"

#endif


using namespace Teuchos;
using namespace Playa;
using namespace PlayaOps;


int main(int argc, char *argv[]) 
{
  typedef Teuchos::ScalarTraits<double> ST;

  try
    {
      GlobalMPISession session(&argc, &argv);
      

      MPIComm::world().synchronize();

      VectorType<double> type = new EpetraVectorType();

      /* create the range space  */
      int nLocalRows = 10;

      MatrixLaplacian1D builder(nLocalRows, type);

      LinearOperator<double> A = builder.getOp();

      int nBlocks = 3;
      Array<Vector<double> > x(nBlocks);
      Array<VectorSpace<double> > space(nBlocks);
      for (int i=0; i<nBlocks; i++)
        {
          space[i] = A.domain();
          x[i] = A.domain().createMember();
          Thyra::randomize(-ST::one(),+ST::one(),x[i].ptr().ptr());          
        }

      VectorSpace<double> blockSpace = blockSpace(space);

      LinearOperator<double> bigA = makeBlockOperator(blockSpace, blockSpace);
      Vector<double> bigRHS = blockSpace.createMember();
      Vector<double> bigX = blockSpace.createMember();
      
      for (int i=0; i<nBlocks; i++)
        {
          bigX.setBlock(i, x[i]);
          for (int j=i; j<nBlocks; j++)
            {
              MatrixLaplacian1D builder(nLocalRows, type);
              LinearOperator<double> Aij = builder.getOp();
              bigA.setBlock(i,j,Aij);
            }
        }
      bigA.endBlockFill();
      
      bigRHS = bigA * bigX;
      Vector<double> bigSoln = blockSpace.createMember();

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(Sundance::searchForFile("SolverParameters/poissonParams.xml"));
#else
      ParameterXMLFileReader reader("poissonParams.xml");
#endif

      ParameterList solverParams = reader.getParameters();
      LinearSolver<double> solver 
        = LinearSolverBuilder::createSolver(solverParams);
      LinearSolver<double> blockSolver 
        = new BlockTriangularSolver<double>(solver);
      
      SolverState<double> state = blockSolver.solve(bigA, bigRHS, bigSoln);
      
      std::cerr << state << std::endl;

      double err = (bigSoln - bigX).norm2();
      std::cerr << "error norm = " << err << std::endl;

      double tol = 1.0e-8;
      if (err > tol)
        {
          std::cerr << "Poisson solve test FAILED" << std::endl;
          return 1;
        }
      else
        {
          std::cerr << "Poisson solve test PASSED" << std::endl;
          return 0;
        }
    }
  catch(std::exception& e)
    {
      std::cerr << "Caught exception: " << e.what() << std::endl;
      return -1;
    }
  return 0;
}

