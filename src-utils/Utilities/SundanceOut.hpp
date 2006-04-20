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
#include "Teuchos_RefCountPtr.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
  using namespace Teuchos;
  /**
   *
   */
  class Out
    {
    public:
      
      static void println(const string& str) 
      {
        if (hasLogFile()) *logFile() << str << endl;
        if (!suppressStdout()) cout << str << endl;
      }

      static void setLogFile(const string& filename)
      {
        logFile() = rcp(new ofstream(filename.c_str()));
        hasLogFile() = true;
      }

      static bool& suppressStdout() {static bool rtn=false; return rtn;}

    private:
      static bool& hasLogFile() {static bool rtn=false; return rtn;}
      static RefCountPtr<ostream>& logFile() {static RefCountPtr<ostream> rtn; return rtn;}
      
    };

}

#define SUNDANCE_OUT(test, msg) \
{ \
  if (test) \
    { \
      TeuchosOStringStream omsg; \
      omsg << msg; \
      Out::println(omsg.str());                 \
    } \
}


#define SUNDANCE_VERB_EXTREME(msg) SUNDANCE_OUT(this->verbosity() > VerbHigh, msg)
#define SUNDANCE_VERB_HIGH(msg) SUNDANCE_OUT(this->verbosity() > VerbMedium, msg)
#define SUNDANCE_VERB_MEDIUM(msg) SUNDANCE_OUT(this->verbosity() > VerbLow, msg)
#define SUNDANCE_VERB_LOW(msg) SUNDANCE_OUT(this->verbosity() > VerbSilent, msg)

#define SUNDANCE_HEADER_LINE "\n------------------------------------------------------------------\n"



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
