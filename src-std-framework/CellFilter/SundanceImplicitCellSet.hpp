/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_IMPLICITCELLSET_H
#define SUNDANCE_IMPLICITCELLSET_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellSetBase.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceCellType.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMesh.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * ImplicitCellSet is a cell set subtype where the set of cell LIDs
     * is never stored. Iteration is done by simply advancing the 
     * LID by one. 
     * 
     * @see CellFilter, CellSetBase, CellIterator 
     **/
    class ImplicitCellSet : public CellSetBase
    {
    public:

      /** Construct with a mesh */
      ImplicitCellSet(const Mesh& mesh, int cellDim,
                      const CellType& cellType);

      /** Returns an iterator pointing to the first element
       * in the set. */
      virtual CellIterator begin() const ;

      /** Returns a past-the-end iterator */
      virtual CellIterator end() const ;

      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(ostream& os) const ;
      //@}

      /* Handleable interface */
      GET_RCP(CellSetBase);


    private:
      
    };
  }

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
