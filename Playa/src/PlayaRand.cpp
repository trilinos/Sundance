/* @HEADER@ */
//   
/* @HEADER@ */


#include "PlayaRand.hpp"
#include "Teuchos_MPIComm.hpp"
#include "Teuchos_TestForException.hpp"
#include "PlayaOut.hpp"

namespace Playa
{
using Teuchos::MPIComm;

void Rand::setLocalSeed(const MPIComm& comm, int seed)
{
  int rank = comm.getRank();
  int lseed = seed;
  for (int i=0; i<rank; i++) lseed = (lseed * 371761) % 5476181;

#ifdef _WIN32
  srand((long) lseed);
#else
  srand48(lseed);
#endif
}

double Rand::val()
{
#ifdef _WIN32
  return ((double) rand())/RAND_MAX;
#else
  return drand48();
#endif
}


}
