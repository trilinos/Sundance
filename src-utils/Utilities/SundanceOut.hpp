/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_OUT_H
#define SUNDANCE_OUT_H

#include "SundanceDefs.hpp"
#include "Teuchos_TestForException.hpp"


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
      Out::println(string(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg))); \
    } \
}



#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
