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
#include "SundanceTriangleQuadrature.hpp"
#include "SundanceTetQuadrature.hpp"

using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
using namespace SundanceStdFwk;
using namespace SundanceStdFwk::Internal;
using namespace Teuchos;
using namespace TSFExtended;


int main(int argc, void** argv)
{
	try
		{
			MPISession::init(&argc, &argv);

      Array<int> validTetOrders = tuple(1, 2, 4, 6);
      Array<int> validTriOrders = tuple(1,2,3,4,5,6);

      cerr << "------------- testing triangle rules -------------------"  << endl;
      for (int i=0; i<validTriOrders.size(); i++)
				{
          int p = validTriOrders[i];
					bool pass = TriangleQuadrature::test(p);
					if (pass) cerr << "order " << p << " PASSED" << endl;
					else cerr << "order " << p << " FAILED" << endl;
				}
      cerr << "------------- testing tet rules -------------------"  << endl;
          
      for (int i=0; i<validTetOrders.size(); i++)
				{
          int p = validTetOrders[i];
					bool pass = TetQuadrature::test(p);
					if (pass) cerr << "order " << p << " PASSED" << endl;
					else cerr << "order " << p << " FAILED" << endl;
				}

		}
	catch(exception& e)
		{
      cerr << "Detected exception: " << e.what() << endl;
      cerr << "Quadrature test FAILED" << endl;
		}


  MPISession::finalize();
}
