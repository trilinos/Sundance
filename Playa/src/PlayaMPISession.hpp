// @HEADER

// @HEADER

#ifndef PLAYA_MPISESSION_H
#define PLAYA_MPISESSION_H

/*! \file Playa_MPISession.hpp
    \brief A MPI utilities class, providing methods for initializing,
	finalizing, and querying the global MPI session
*/
#include "PlayaDefs.hpp"
#include <stack>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

namespace Playa
{
  /**
   * \brief This class provides methods for initializing, finalizing, 
   * and querying the global MPI session. 
   */
  class MPISession
    {
    public:
      //! Initializer, calls MPI_Init() if necessary
      static void init(int* argc, void*** argv);

      //! Returns the process rank relative to MPI_COMM_WORLD
      static int getRank() {return rank_;}

      //! Returns the number of processors in MPI_COMM_WORLD 
      static int getNProc() {return nProc_;}

      //! Finalizer, calls MPI_Finalize() if necessary
      static void finalize();

      /** Set to true if a message should be written by each processor
       * at startup. */
      static bool& showStartupMessage() {static bool rtn=false; return rtn;}

    private:
      static int rank_;
      static int nProc_;
    };
}

#endif
