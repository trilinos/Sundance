/* @HEADER@ */
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
