/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EXCEPTIONS_H
#define SUNDANCE_EXCEPTIONS_H

#include "SundanceDefs.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include <stdexcept>


#ifndef DOXYGEN_DEVELOPER_ONLY

#define SUNDANCE_ERROR(msg) \
{ \
  TestForException_break(); \
  TeuchosOStringStream omsg; \
	omsg << __FILE__ << ":" << __LINE__ << ": " \
       << ": " << msg; \
  throw SundanceUtils::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define SUNDANCE_TRACE(e) \
{ \
  TeuchosOStringStream omsg; \
	omsg << e.what() << endl \
  << "caught in " << __FILE__ << ":" << __LINE__ << endl ; \
  throw SundanceUtils::RuntimeError(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
}

#define SUNDANCE_BOUNDSCHECK(i, low, high, msg) \
{ \
  TEST_FOR_EXCEPTION( i < low || i > high, SundanceUtils::RuntimeError, \
                     "Bounds violation: " << #i << "is out of range [" \
                      << #low << ", " << #high << "]") \
}


namespace SundanceUtils
{
  /**
   * InternalError is thrown when an "impossible" case is detected
   * in Sundance's internals. Occurance of an InternalError indicates
   * either a bug in Sundance or runtime memory corruption that is
   * scrambling an object's virtual tables.
   */
  class InternalError : public std::logic_error
    {
    public:
      /** */
      InternalError(const string& msg);
    };

  /**
   * RuntimeError is an exception that occurs as a result of invalid
   * user-level code.
   */
  class RuntimeError : public std::runtime_error
    {
    public:
      /** */
      RuntimeError(const string& msg);
    };

  /**
   * BadSymbolicsError is thrown when a mathematically nonsensical
   * expression is detected. Unfortunately, it is possible to form expressions
   * that are legal C++ but illegal mathematics. For example, one can
   * happily divide a differential operator by another expression,
   * \code
   * Expr f = new UnknownFunction();
   * Expr dx = new Derivative(0);
   * Expr e = dx/f;
   * \endcode
   */
  class BadSymbolicsError : public RuntimeError
    {
    public:
      /** */
      BadSymbolicsError(const string& msg);
    };
}

                  
#endif  /* DOXYGEN_DEVELOPER_ONLY */  


#endif
