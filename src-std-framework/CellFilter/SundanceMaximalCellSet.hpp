/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_MAXIMALCELLSET_H
#define SUNDANCE_MAXIMALCELLSET_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellSetBase.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  namespace Internal {}
  using namespace Internal;
  using namespace Teuchos;

  /** */
  class MaximalCellSet : public CellSetBase 
  {
  public:
    /** */
    MaximalCellSet();

    /** */
    virtual ~MaximalCellSet(){;}

    /** */
    virtual XMLObject toXML() const ;

    /** */
    virtual string typeName() const {return "MaximalCellSet";}

    /** */
    virtual bool lessThan(const CellSetBase* other) const ;

    /** */
    virtual RefCountPtr<CellSetBase> getRcp() {return rcp(this);}
  };

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
