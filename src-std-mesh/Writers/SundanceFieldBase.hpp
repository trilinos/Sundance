/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_FIELDBASE_H
#define SUNDANCE_FIELDBASE_H


#ifndef DOXYGEN_DEVELOPER_ONLY


#include "SundanceDefs.hpp"
#include "TSFHandleable.hpp"

namespace SundanceStdMesh
{
  namespace Internal
  {
    /**
     *
     */
    class FieldBase : public TSFExtended::Handleable<FieldBase>
    {
    public:
      /** */
      FieldBase(){;}

      /** virtual dtor */
      virtual ~FieldBase(){;}

      /** */
      virtual int numElems() const {return 1;}

      /** */
      virtual double getData(int cellDim, int cellID, int elem) const = 0 ;

    };
  }
}


#endif /* DOXYGEN_DEVELOPER_ONLY */

#endif
