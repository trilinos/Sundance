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


int main(int argc, char** argv)
{
  int stat = 0 ;
	try
		{
			GlobalMPISession session(&argc, &argv);

      Array<int> validTetOrders = tuple(1, 2, 4, 6);
      int maxorder = 15;
      Array<int> validTriOrders(maxorder);
      for (int i=0; i<maxorder; i++) validTriOrders[i] = i+1;

      Array<int> triFailures;
      Array<int> tetFailures;

      cerr << "------------- testing triangle rules -------------------"  << endl;
      for (unsigned int i=0; i<validTriOrders.size(); i++)
				{
          int p = validTriOrders[i];
					bool pass = TriangleQuadrature::test(p);
					if (pass) cerr << "order " << p << " PASSED" << endl;
					else 
            {
              cerr << "order " << p << " FAILED" << endl;
              triFailures.append(p);
            }
				}
      cerr << "------------- testing tet rules -------------------"  << endl;
          
      for (unsigned int i=0; i<validTetOrders.size(); i++)
				{
          int p = validTetOrders[i];
					bool pass = TetQuadrature::test(p);
					if (pass) cerr << "order " << p << " PASSED" << endl;
					else 
            {
              cerr << "order " << p << " FAILED" << endl;
              tetFailures.append(p);
            }
				}

      if (tetFailures.size()>0) 
        {
          cout << "failures detected for tets: orders " << tetFailures << endl;
          cout << "tet tests FAILED" << endl;
          stat = -1;
        }
      else
        {
          cout << "tet tests PASSED" << endl;
        }

      if (triFailures.size()>0) 
        {
          cout << "failures detected for tris: orders " << triFailures << endl;
          cout << "tri tests FAILED" << endl;
          stat = -1;
        }
      else
        {
          cout << "tri tests PASSED" << endl;
        }


		}
	catch(exception& e)
		{
      cerr << "Detected exception: " << e.what() << endl;
      cerr << "Quadrature test FAILED" << endl;
      stat = -1;
		}

  return stat;
  
}
