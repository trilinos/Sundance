/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_RAND_H
#define PLAYA_RAND_H

#include "PlayaDefs.hpp"

namespace Teuchos
{
class MPIComm;
}

namespace Playa
{
using Teuchos::MPIComm;

/** */
class Rand
{
public:
  /** */
  static void setLocalSeed(const MPIComm& comm, int seed);

  /** */
  static double val();
};

}

#endif
