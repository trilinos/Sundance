/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_NULLCELLFILTER_H
#define SUNDANCE_NULLCELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellFilterBase.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * NullCellFilter is used as a placeholder cell filter
     * in those equations defined
     * independently of geometry, i.e., equations involving only
     * global parameters.
     **/
    class NullCellFilter : public CellFilterBase 
    {
    public:
      /** */
      NullCellFilter();

      /** */
      virtual ~NullCellFilter(){;}

      /** */
      virtual XMLObject toXML() const ;

      /** */
      virtual string typeName() const {return "NullCellFilter";}

      /** */
      virtual bool lessThan(const CellFilterBase* other) const ;

      /** Return the dimension of the cells that will be identified
       * by this filter when acting on the given mesh */
      virtual int dimension(const Mesh& mesh) const ;

      /* */
      GET_RCP(CellFilterBase);
    
    protected:

      /** */
      virtual CellSet internalGetCells(const Mesh& mesh) const ;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
