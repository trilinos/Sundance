/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FIELDBASE_H
#define SUNDANCE_FIELDBASE_H


#ifndef DOXYGEN_DEVELOPER_ONLY


#include "SundanceDefs.hpp"

namespace SundanceStdMesh
{
  namespace Internal
  {
    /**
     */
    class FieldBase 
    {
    public:
      /** */
      FieldBase(){;}

      /** virtual dtor */
      virtual ~FieldBase(){;}
    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
