/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "Sundance.hpp"

int main(int argc, char** argv)
{
  try
  {
    /* Initialization */
    Sundance::init(&argc, &argv);

    /* ---- BEGIN CODE BODY --- */

    /* The main simulation code goes here. In this example, all we do
     * is to print some information about the processor ranks. */

    MPIComm comm = MPIComm::world();

    /* Print a header from the root processor only. Although this executes on
     * all processors, anything written to the output stream Out::root() 
     * is ignored on all non-root processors (rank != 0). 
     * After writing, synchronize to keep this message from getting jumbled
     * together with the subsequent messages. 
     */
    Out::root() << "Example: getting started" << endl;
    comm.synchronize();

    /* Every processor now speaks up and identifies itself */
    int myRank = comm.getRank();
    int nProc = comm.getNProc();
    Out::os() << "Processor " << myRank 
              << " of " << nProc << " checking in" << endl;

    /* ---- END CODE BODY --- */

    /* Test success or failure. Most examples you'll see will do this 
     * as part of the Trilinos regression testing system. 
     * If you write a simulation code that won't become part of Trilinos,
     * you often can bypass this step.
     *
     * Here the test is a trival one: every processor's rank must be
     * smaller than the total number of processors. If this fails,
     * your MPI installation is probably broken!
     * */
    Sundance::passFailTest(myRank < nProc);
  }
	catch(std::exception& e) /* exception handling */
  {
    cerr << "exception!" << endl;
    Sundance::handleException(e);
  }
  /* Finalization */
  Sundance::finalize();

  return Sundance::testStatus();
}
