#include "SundanceCollectiveExceptionCheck.hpp"
#include "PlayaMPIComm.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace Playa;

namespace Sundance
{
  void reportFailure(const MPIComm& comm)
  {
    int myBad = 1;
    int anyBad = 0;
    comm.allReduce((void*) &myBad, (void*) &anyBad, 1, MPIDataType::intType(),
                     MPIOp::sumOp());
  }

  bool checkForFailures(const MPIComm& comm)
  {
    int myBad = 0;
    int anyBad = 0;
    comm.allReduce((void*) &myBad, (void*) &anyBad, 1, MPIDataType::intType(),
                     MPIOp::sumOp());
    return anyBad > 0;
  }


}



	





