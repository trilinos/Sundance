//@HEADER@

//@HEADER@

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
#include "PlayaDefaultBlockVectorSpaceDecl.hpp"
#include "PlayaLinearCombinationDecl.hpp"
#include "PlayaBlockVectorSpaceDecl.hpp"
#include "PlayaAztecSolver.hpp"
#include "PlayaMatrixLaplacian1D.hpp"
#include "PlayaLinearSolverBuilder.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaSimpleBlockOpDecl.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaLinearSolverImpl.hpp"
#include "PlayaVectorImpl.hpp"
#include "PlayaSimpleBlockOpImpl.hpp"
#include "PlayaDefaultBlockVectorImpl.hpp"
#include "PlayaBlockTriangularSolverImpl.hpp"

#endif


using namespace Teuchos;
using namespace Playa;
using namespace PlayaExprTemplates;


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
          x[i].randomize();
        }

      VectorSpace<double> blkSpace = blockSpace(space);

      LinearOperator<double> bigA = makeBlockOperator(blkSpace, blkSpace);
      Vector<double> bigRHS = blkSpace.createMember();
      Vector<double> bigX = blkSpace.createMember();
      
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
      Vector<double> bigSoln = blkSpace.createMember();

#ifdef HAVE_CONFIG_H
      ParameterXMLFileReader reader(Playa::searchForFile("SolverParameters/poissonParams.xml"));
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

