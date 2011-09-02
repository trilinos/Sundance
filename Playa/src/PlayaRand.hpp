/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_RAND_H
#define PLAYA_RAND_H

#include "PlayaDefs.hpp"


namespace Playa
{
class MPIComm;

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
