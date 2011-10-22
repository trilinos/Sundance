/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaGlobalAnd.hpp"
#include "PlayaMPIComm.hpp"
#include "Teuchos_Assert.hpp"
#include "PlayaOut.hpp"
#ifdef HAVE_MPI
#include "mpi.h"
#endif


namespace Playa
{
bool globalAnd(bool localVal)
{
  int out = localVal;

  
#ifdef HAVE_MPI
  int in = localVal;
  int ierr = ::MPI_Allreduce((void*) &in, (void*) &out, 1, MPI_INT,
    MPI_LAND, MPI_COMM_WORLD);
  TEUCHOS_TEST_FOR_EXCEPTION(ierr != 0, std::runtime_error,
    "MPI_Allreduce returned error code=" << ierr);
#endif

  return out;
}
}
