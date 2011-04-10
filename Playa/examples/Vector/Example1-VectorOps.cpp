/* @HEADER@ */
//   
/* @HEADER@ */

#include "PlayaEpetraVectorType.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaOut.hpp"
#include "PlayaLinearCombinationImpl.hpp"
#include "PlayaVectorOpsImpl.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_MPIComm.hpp"


/* ------------------------------------------------------------------------
 *
 * This example shows the construction and use of simple vectors
 * 
 * ------------------------------------------------------------------------ */

using namespace Playa;

int main(int argc, char *argv[])
{
  int rtn = 0;

  try
  {
    /* Initialize MPI */
    GlobalMPISession session(&argc, &argv);


    /* The VectorType object will be used when we create vector spaces, 
     * specifying what type of low-level linear algebra implementation
     * will be used. */
    VectorType<double> vecType = new EpetraVectorType();

    /* Construct a vector space  */
    int n = 100;
    VectorSpace<double> vs = vecType.createEvenlyPartitionedSpace(MPIComm::world(), n);


    /* Make some vectors */
    Vector<double> x = vs.createMember();
    x.randomize();

    Vector<double> y = vs.createMember();
    y.randomize();

    Vector<double> z = vs.createMember();
    z.randomize();

    Vector<double> u = vs.createMember();
    u.randomize();

    Vector<double> v = vs.createMember();
    v.randomize();

    
    double err = 1.0e10;
    double tol = 1.0e-10;
    int fails = 0;

    /* Make sure we've not created a zero-norm random vector */
    Out::root() << "||x||=" << norm2(x) << endl;


    /* Test the deep copying of a vector */
    Vector<double> xCopy = x.copy();
    err = norm2(x-xCopy);
    Out::root() << "||copy error|| = " << err << endl;
    if (!(err < tol)) fails++;

    /* Do the computation w/o making a temporary */
    err = norm2Dist(x, xCopy); 
    Out::root() << "||copy error|| = " << err << endl;
    if (!(err < tol)) fails++;

    /* Test scalar multiplication */
    x *= 2.0;
 
    err = ::fabs( (x-xCopy).norm2() - xCopy.norm2() );
    Out::root() << "|norm(2x-x) - norm(x)| = " << err << endl;
    if (!(err < tol)) fails++;

    /* Test scalar division */
    x /= 2.0;
 
    err = norm2Dist(x, xCopy);
    Out::root() << "norm(2x/2-x)= " << err << endl;
    if (!(err < tol)) fails++;

    
     
    if (rtn == 0)
    {
      Out::root() << "test PASSED" << endl;
    }
    else
    {
      Out::root() << "test FAILED" << endl;
    }
  }
  catch(std::exception& e)
  {
    std::cerr << "Caught exception: " << e.what() << endl;
    rtn = -1;
  }
  return rtn;
}
