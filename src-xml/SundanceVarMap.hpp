#ifndef SUNDANCE_VARMAP_H
#define SUNDANCE_VARMAP_H

#include "Sundance.hpp"

namespace SundanceXML
{
  template <class T> class VarMap : public SundanceUtils::Map<string, T>
  {
  public:
    /** */
    VarMap() : SundanceUtils::Map<string, T>() {;}

    /** */
    const T& lookup(const string& varName) const 
      {
        TEST_FOR_EXCEPTION(!containsKey(varName), InternalError,
                           "var name " << varName 
                           << " not found in map ");
      }
    
  };
}

#endif
