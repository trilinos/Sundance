/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_EXCEPTIONS_H
#define PLAYA_EXCEPTIONS_H

#include "PlayaDefs.hpp"
#include "PlayaDebug.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include <stdexcept>


//bvbw for backard compatibility reasons
//     I could not get this to work with ifdefs hence the hack

#define PLAYA_ERROR7(msg) \
{ \
  TestForException_break(); \
  std::ostringstream omsg; \
	omsg << __FILE__ << ":" << __LINE__ << ": " \
       << ": " << msg; \
  throw Playa::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define PLAYA_ERROR(msg) \
{ \
  std::ostringstream omsg; \
	omsg << __FILE__ << ":" << __LINE__ << ": " \
       << ": " << msg; \
  const std::string &omsgstr = omsg.str(); \
  TestForException_break(omsgstr); \
  throw Playa::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}


#define PLAYA_TRACE(e) \
{ \
  std::ostringstream omsg; \
	omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
  throw Playa::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define PLAYA_TRACE_MSG(e, msg)                      \
{ \
  std::ostringstream omsg; \
	omsg << e.what() << std::endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << std::endl ; \
  omsg << msg << std::endl; \
  throw Playa::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define PLAYA_BOUNDSCHECK(i, low, high, location) \
{ \
  TEST_FOR_EXCEPTION( i < low || i > high, Playa::RuntimeError, \
    "Bounds violation in " << location << ": " \
    << #i << "is out of range [" \
    << #low << ", " << #high << "]") \
}

#define PLAYA_CHECK_ARRAY_SIZE_MATCH(a1, a2) \
  {\
    TEST_FOR_EXCEPTION(a1.size() != a2.size(), Playa::RuntimeError, \
y      "Mismatched array sizes: size(" << #a1 << ")=" << a1.size() \
      << " size(" << #a2 << ")=" << a2.size() << ". Expected equal sizes");\
  }



namespace Playa
{
  /**
   * InternalError is thrown when an "impossible" case is detected
   * in Playa's internals. Occurance of an InternalError indicates
   * either a bug in Playa or runtime memory corruption that is
   * scrambling an object's virtual tables.
   */
  class InternalError : public std::logic_error
    {
    public:
      /** */
      InternalError(const std::string& msg);
    };

  /**
   * RuntimeError is an exception that occurs as a result of invalid
   * user-level code.
   */
  class RuntimeError : public std::runtime_error
    {
    public:
      /** */
      RuntimeError(const std::string& msg);
    };
  
}

                  

#endif
