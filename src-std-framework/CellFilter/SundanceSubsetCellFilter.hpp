/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_SUBSETCELLFILTER_H
#define SUNDANCE_SUBSETCELLFILTER_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellPredicate.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  using namespace TSFExtended;
  namespace Internal 
  {
    using namespace Teuchos;

    /** */
    class SubsetCellFilter : public CellFilterBase 
    {
    public:
      /** */
      SubsetCellFilter(const CellFilter& superset,
                       const CellPredicate& predicate);

      /** */
      virtual ~SubsetCellFilter(){;}

      /** Return the dimension of the cells that will be identified
       * by this filter when acting on the given mesh */
      virtual int dimension(const Mesh& mesh) const 
      {return superset_.dimension(mesh);}

      /** */
      virtual XMLObject toXML() const ;

      /** */
      virtual string typeName() const {return "SubsetCellFilter";}

      /** */
      virtual bool lessThan(const CellFilterBase* other) const ;

      /* */
      GET_RCP(CellFilterBase);
    

    protected:
      /** */
      virtual CellSet internalGetCells(const Mesh& mesh) const ;

    private:
      /** */
      CellFilter superset_;

      /** */
      CellPredicate predicate_;
    };
  }
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
