// @HEADER
// @HEADER

#include "PlayaErrorPolling.hpp"
#include "PlayaMPIComm.hpp"

namespace Playa
{
  void ErrorPolling::reportFailure(const MPIComm& comm)
  {
    if (isActive())
      {
        int myBad = 1;
        int anyBad = 0;
        comm.allReduce((void*) &myBad, (void*) &anyBad, 1, MPIDataType::intType(),
                       MPIOp::sumOp());
      }
  }

  bool ErrorPolling::pollForFailures(const MPIComm& comm)
  {
    /* bypass if inactive */
    if (!isActive()) return true;

    int myBad = 0;
    int anyBad = 0;
    try
      {
        comm.allReduce((void*) &myBad, (void*) &anyBad, 1, MPIDataType::intType(),
                       MPIOp::sumOp());
      }
    catch(const std::exception&)
      {
        return true;
      }
    return anyBad > 0;
  }
}



	





