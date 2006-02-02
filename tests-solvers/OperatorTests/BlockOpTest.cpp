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


#include <cstdlib>
#include "Teuchos_MPISession.hpp"
#include "TSFVector.hpp"
#include "TSFLinearCombination.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFLoadableMatrix.hpp"
#include "TSFVectorType.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFEpetraVectorType.hpp"
#include "TSFEpetraVectorSpace.hpp"
//#include "TSFCoreEpetraVectorSpace.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPIComm.hpp"
//#include "TSFLinearSolver.hpp"
#include "TSFBICGSTABSolver.hpp"
#include "TSFProductVectorSpace.hpp"
#include "TSFTransposeOperator.hpp"
//#include "TSFInverseOperator.hpp"
#include "TSFLinearOperator.hpp"
#include "TSFEpetraMatrix.hpp"
//#include "Thyra_LinearOpBase.hpp"
#include "TSFIdentityOperator.hpp"
#include "TSFZeroOperator.hpp"
#include "TSFBlockOperator.hpp"
#include "TSFDiagonalOperator.hpp"
#include "TSFScaledOperator.hpp"
#include "TSFSumOperator.hpp"
#include "TSFComposedOperator.hpp"
#include "TSFMatrixLaplacian1D.hpp"
#include "TSFRandomSparseMatrix.hpp"
#include "TSFRandomBlockMatrix.hpp"
#include "TSFCompoundTester.hpp"


using namespace Teuchos;
using namespace TSFExtended;
using namespace TSFExtendedOps;
using Thyra::TestSpecifier;

int main(int argc, void *argv[]) 
{
  try
    {
      MPISession::init(&argc, &argv);
      MPIComm::world().synchronize();

      cerr << "go!" << endl;
      VectorType<double> type = new EpetraVectorType();

      Array<int> domainBlockSizes = tuple(2,3,4);
      Array<int> rangeBlockSizes = tuple(2,2);

      Array<VectorSpace<double> > domainBlocks(domainBlockSizes.size());
      Array<VectorSpace<double> > rangeBlocks(rangeBlockSizes.size());

      for (unsigned int i=0; i<domainBlocks.size(); i++)
        {
          domainBlocks[i] = type.createEvenlyPartitionedSpace(domainBlockSizes[i]);
        }

      for (unsigned int i=0; i<rangeBlocks.size(); i++)
        {
          rangeBlocks[i] = type.createEvenlyPartitionedSpace(rangeBlockSizes[i]);
        }
      
      VectorSpace<double> domain = productSpace(domainBlocks);
      VectorSpace<double> range = productSpace(rangeBlocks);

      double blockDensity = 0.75;
      double onProcDensity = 0.5;
      double offProcDensity = 0.1;
      
      RandomBlockMatrix<double> builder(domain, range, 
                                        blockDensity,
                                        onProcDensity,
                                        offProcDensity,
                                        type);

      LinearOperator<double> A = builder.getOp();

      cerr << "A num block rows = " << A.numBlockRows() << endl;
      cerr << "A num block cols = " << A.numBlockCols() << endl;

      Vector<double> x = domain.createMember();
      cerr << "randomizing trial vector" << endl;
      Thyra::randomize(-1.0, 1.0, x.ptr().get());

      Array<Vector<double> > xBlock(domain.numBlocks());
      for (unsigned int i=0; i<xBlock.size(); i++)
        {
          xBlock[i] = x.getBlock(i);
        }

      cerr << "------------------------------------------------------------" << endl;
      cerr << "computing A*x..." << endl;
      Vector<double> y0 = A * x;
      for (int i=0; i<y0.space().numBlocks(); i++)
        {
          cerr << "y0[" << i << "] = " << endl << y0.getBlock(i) << endl;
        }
      

      Vector<double> y1 = range.createMember();
      cerr << "------------------------------------------------------------" << endl;
      cerr << "computing A*x block-by-block..." << endl;
      Array<Vector<double> > yBlock(range.numBlocks());
      for (unsigned int i=0; i<yBlock.size(); i++)
        {
          yBlock[i] = range.getBlock(i).createMember();
          yBlock[i].zero();
          for (unsigned int j=0; j<xBlock.size(); j++)
            {
              LinearOperator<double> Aij = A.getBlock(i,j);
              if (Aij.ptr().get() != 0)
                {
                  cerr << "A(" << i << ", " << j << ") = " << endl 
                       << Aij << endl;
                }
              else
                {
                  cerr << "A(" << i << ", " << j << ") = 0 " << endl;
                }
              cerr << "x[" << j << "] = " << endl << xBlock[j] << endl;
              if (Aij.ptr().get()==0) continue;
              yBlock[i] = yBlock[i] + Aij * xBlock[j];
            }
          y1.setBlock(i, yBlock[i]);
        }

      for (int i=0; i<y1.space().numBlocks(); i++)
        {
          cerr << "y1[" << i << "] = " << endl << y1.getBlock(i) << endl;
        }
      double err = (y1 - y0).norm2();
      cerr << "error = " << err << endl;

      double tol = 1.0e-13;
      if (err < tol)
        {
          cerr << "block op test PASSED" << endl;
        }
      else
        {
          cerr << "block op test FAILED" << endl;
        }
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();
}



