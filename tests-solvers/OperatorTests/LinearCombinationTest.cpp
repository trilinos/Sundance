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
#include "TSFDiagonalOperator.hpp"
#include "TSFScaledOperator.hpp"
#include "TSFSumOperator.hpp"
#include "TSFComposedOperator.hpp"
#include "TSFMatrixLaplacian1D.hpp"
#include "TSFRandomSparseMatrix.hpp"
#include "TSFLinearCombinationTester.hpp"

using namespace Teuchos;
using namespace std;
using namespace TSFExtended;
using namespace TSFExtendedOps;
using Thyra::TestSpecifier;

int main(int argc, void *argv[]) 
{
  try
    {
      MPISession::init(&argc, &argv);

      VectorType<double> type = new EpetraVectorType();

      int nLocalRows = 1;
      
      double onProcDensity = 0.5;
      double offProcDensity = 0.1;

      LinearCombinationTester<double> tester(nLocalRows,
                                             onProcDensity,
                                             offProcDensity,
                                             type,
                                             TestSpecifier<double>(true, 1.0e-13, 1.0e-10));


      tester.runAllTests();
    }
  catch(std::exception& e)
    {
      cerr << "Caught exception: " << e.what() << endl;
    }
  MPISession::finalize();
}



