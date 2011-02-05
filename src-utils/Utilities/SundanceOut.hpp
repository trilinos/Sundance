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
#include "SundanceObjectWithVerbosity.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_MPIComm.hpp"

using Playa::Out;
using Playa::Tabs;

#define SUNDANCE_OUT(test, msg) PLAYA_OUT(test, msg)


#define SUNDANCE_VERB_EXTREME(msg) PLAYA_MSG4(this->verb(), msg)
#define SUNDANCE_VERB_HIGH(msg) PLAYA_MSG3(this->verb(), msg)
#define SUNDANCE_VERB_MEDIUM(msg) PLAYA_MSG2(this->verb(), msg)
#define SUNDANCE_VERB_LOW(msg) PLAYA_MSG1(this->verb(), msg)

#define SUNDANCE_HEADER_LINE "\n------------------------------------------------------------------\n"

#define SUNDANCE_MSG(context, level, msg) PLAYA_OUT(this->verbLevel(context) >= level, msg)

#define SUNDANCE_LEVEL1(context, msg) PLAYA_MSG(context, 1, msg)

#define SUNDANCE_LEVEL2(context, msg) PLAYA_MSG(context, 2, msg)

#define SUNDANCE_LEVEL3(context, msg) PLAYA_MSG(context, 3, msg)

#define SUNDANCE_LEVEL4(context, msg) PLAYA_MSG(context, 4, msg)

#define SUNDANCE_LEVEL5(context, msg) PLAYA_MSG(context, 5, msg)


#define SUNDANCE_MSG1(level, msg) PLAYA_OUT(level >= 1, msg)

#define SUNDANCE_MSG2(level, msg) PLAYA_OUT(level >= 2, msg)

#define SUNDANCE_MSG3(level, msg) PLAYA_OUT(level >= 3, msg)

#define SUNDANCE_MSG4(level, msg) PLAYA_OUT(level >= 4, msg)

#define SUNDANCE_MSG5(level, msg) PLAYA_OUT(level >= 5, msg)

#define SUNDANCE_BANNER1(level, tab, msg) \
  PLAYA_MSG1(level, tab \
    << "===================================================================");\
  PLAYA_MSG1(level, tab << std::endl << tab \
    << "  " << msg); \
  PLAYA_MSG1(level, tab << std::endl << tab\
    << "===================================================================");


#define SUNDANCE_BANNER2(level, tab, msg) \
  SUNDANCE_MSG2(level, tab \
    << "-------------------------------------------------------------------");\
  SUNDANCE_MSG2(level, tab << msg); \
  SUNDANCE_MSG2(level, tab\
    << "-------------------------------------------------------------------");



#define SUNDANCE_BANNER3(level, tab, msg) \
  SUNDANCE_MSG3(level, tab \
    << "-------------------------------------------------------------------");\
  SUNDANCE_MSG3(level, tab << std::endl << tab \
    << msg); \
  SUNDANCE_MSG3(level, tab << std::endl << tab\
    << "-------------------------------------------------------------------");

#define SUNDANCE_TRACE(e) \
{ \
  TeuchosOStringStream omsg; \
        omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
        throw std::runtime_error(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define SUNDANCE_TRACE_MSG(e, msg)                      \
{ \
  TeuchosOStringStream omsg; \
        omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
  omsg << msg << std::endl; \
  throw std::runtime_error(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}


#define SUNDANCE_ERROR(msg) \
{ \
  TeuchosOStringStream omsg; \
        omsg << __FILE__ << ":" << __LINE__ << ": " \
       << ": " << msg; \
  const std::string &omsgstr = omsg.str(); \
  TestForException_break(omsgstr); \
  throw std::runtime_error(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}


#endif
