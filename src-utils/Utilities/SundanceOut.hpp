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

#ifndef SUNDANCE_OUT_H
#define SUNDANCE_OUT_H

#include "SundanceDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include "TSFObjectWithVerbosity.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
  /**
   *
   */
  class Out
    {
    public:
      
      static void println(const string& str) {cerr << str << endl;}
    private:
    };

}

#define SUNDANCE_OUT(test, msg) \
{ \
  if (test) \
    { \
      TeuchosOStringStream omsg; \
      omsg << msg; \
      SundanceUtils::Out::println(string(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg))); \
    } \
}


#define SUNDANCE_VERB_EXTREME(msg) SUNDANCE_OUT(verbosity() > VerbHigh, msg)
#define SUNDANCE_VERB_HIGH(msg) SUNDANCE_OUT(verbosity() > VerbMedium, msg)
#define SUNDANCE_VERB_MEDIUM(msg) SUNDANCE_OUT(verbosity() > VerbLow, msg)
#define SUNDANCE_VERB_LOW(msg) SUNDANCE_OUT(verbosity() > VerbSilent, msg)




#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
