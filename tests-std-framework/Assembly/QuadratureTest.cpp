/* @HEADER@ */
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
		}


  MPISession::finalize();
}
