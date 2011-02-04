/* @HEADER@ */
//   
/* @HEADER@ */

#ifndef PLAYA_OUT_HPP
#define PLAYA_OUT_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_MPIComm.hpp"

namespace Teuchos
{
template <class T> class Array;
}

namespace Playa
{
class Tabs;
using namespace Teuchos;

/**
 * Class Out provides standardized access to fancy streams for writing
 * diagnostic output. 
 */
class Out
{
public:
      
  /** Print a line followed by termination */
  static void println(const std::string& str) 
    {
      if (!suppressStdout()) os() << str << std::endl;
    }

  /** Provide a fancy ostream wrapper to cout */
  static FancyOStream& os()
    {
      static RCP<std::ostream> os = rcp(&std::cout, false);
      static RCP<FancyOStream> rtn = fancyOStream(os);
      static bool first = true;
      if (first)
      {
        rtn->setShowProcRank(true);
        first = false;
      }
      return *rtn;
    }

  /** Provide a fancy ostream wrapper to cout, ignoring output from
   * non-root processors*/
  static FancyOStream& root()
    {
      static bool isRoot = MPIComm::world().getRank()==0;
      static RCP<FancyOStream> blackHole
        = rcp(new FancyOStream(rcp(new oblackholestream())));

      if (isRoot)
      {
        return os();
      }
      else
      {
        return *blackHole;
      }
    }

  /** */
  static bool& suppressStdout() {static bool rtn=false; return rtn;}
};


/** Write an array formatted to show a specified number of columns */
void writeTable(std::ostream& os, const Tabs& tab,
  const Array<double>& a, int cols);
/** Write an array formatted to show a specified number of columns */
void writeTable(std::ostream& os, const Tabs& tab, 
  const Array<int>& a, int cols);

}

#define PLAYA_OUT(test, msg) \
  { \
    if (test) \
    { \
      TeuchosOStringStream omsg; \
      omsg << msg; \
      Out::println(omsg.str());                 \
    } \
  }


#define PLAYA_VERB_EXTREME(msg) PLAYA_MSG4(this->verb(), msg)
#define PLAYA_VERB_HIGH(msg) PLAYA_MSG3(this->verb(), msg)
#define PLAYA_VERB_MEDIUM(msg) PLAYA_MSG2(this->verb(), msg)
#define PLAYA_VERB_LOW(msg) PLAYA_MSG1(this->verb(), msg)

#define PLAYA_HEADER_LINE "\n------------------------------------------------------------------\n"

#define PLAYA_MSG(context, level, msg) PLAYA_OUT(this->verbLevel(context) >= level, msg)

#define PLAYA_LEVEL1(context, msg) PLAYA_MSG(context, 1, msg)

#define PLAYA_LEVEL2(context, msg) PLAYA_MSG(context, 2, msg)

#define PLAYA_LEVEL3(context, msg) PLAYA_MSG(context, 3, msg)

#define PLAYA_LEVEL4(context, msg) PLAYA_MSG(context, 4, msg)

#define PLAYA_LEVEL5(context, msg) PLAYA_MSG(context, 5, msg)


#define PLAYA_MSG1(level, msg) PLAYA_OUT(level >= 1, msg)

#define PLAYA_MSG2(level, msg) PLAYA_OUT(level >= 2, msg)

#define PLAYA_MSG3(level, msg) PLAYA_OUT(level >= 3, msg)

#define PLAYA_MSG4(level, msg) PLAYA_OUT(level >= 4, msg)

#define PLAYA_MSG5(level, msg) PLAYA_OUT(level >= 5, msg)

#define PLAYA_BANNER1(level, tab, msg) \
  PLAYA_MSG1(level, tab \
    << "===================================================================");\
  PLAYA_MSG1(level, tab << std::endl << tab \
    << "  " << msg); \
  PLAYA_MSG1(level, tab << std::endl << tab\
    << "===================================================================");


#define PLAYA_BANNER2(level, tab, msg) \
  PLAYA_MSG2(level, tab \
    << "-------------------------------------------------------------------");\
  PLAYA_MSG2(level, tab << msg); \
  PLAYA_MSG2(level, tab\
    << "-------------------------------------------------------------------");



#define PLAYA_BANNER3(level, tab, msg) \
  PLAYA_MSG3(level, tab \
    << "-------------------------------------------------------------------");\
  PLAYA_MSG3(level, tab << std::endl << tab \
    << msg); \
  PLAYA_MSG3(level, tab << std::endl << tab\
    << "-------------------------------------------------------------------");

#endif
