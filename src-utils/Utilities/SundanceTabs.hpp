/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_TABS_H
#define SUNDANCE_TABS_H

#include "SundanceDefs.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceUtils
{
  /**
   *
   */
  class Tabs
    {
    public:
      /** */
      Tabs();

      /** */
      ~Tabs();

      /** */
      void print(ostream& os) const ;

    private:
      static int& tabLevel() {static int rtn = 0; return rtn;}
    };
}

namespace std
{
  /** \relates Tabs */
  inline ostream& operator<<(ostream& os, const SundanceUtils::Tabs& t)
    {
      t.print(os);
      return os;
    }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
