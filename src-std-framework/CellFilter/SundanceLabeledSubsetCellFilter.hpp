/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_LABELEDSUBSETCELLFILTER_H
#define SUNDANCE_LABELEDSUBSETCELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellFilterBase.hpp"

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
  class LabeledSubsetCellFilter : public CellFilterBase 
  {
  public:
    /** */
    LabeledSubsetCellFilter(const CellFilter& superset,
                            const string& label);

    /** */
    virtual ~LabeledSubsetCellFilter(){;}

    /** */
    virtual XMLObject toXML() const ;

    /** */
    virtual string typeName() const {return "LabeledSubsetCellFilter";}

    /** */
    virtual bool lessThan(const CellFilterStub* other) const ;

    /** */
    virtual RefCountPtr<CellFilterBase> getRcp() {return rcp(this);}
    

  protected:
    /** */
    virtual CellSet internalGetCells(const Mesh& mesh) const ;

    /** */
    string label_;

    /** */
    CellFilter superset_;
  };

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
