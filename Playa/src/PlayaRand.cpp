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
  srand48(lseed);
}

double Rand::val()
{
  return drand48();
}


}
