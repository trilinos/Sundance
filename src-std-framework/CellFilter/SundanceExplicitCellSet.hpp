/* @HEADER@ */
/* @HEADER@ */

#ifndef SUNDANCE_EXPLICITCELLSET_H
#define SUNDANCE_EXPLICITCELLSET_H

#ifndef DOXYGEN_DEVELOPER_ONLY

#include "SundanceDefs.hpp"
#include "SundanceCellSetBase.hpp"

namespace SundanceStdFwk
{
 using namespace SundanceUtils;
using namespace SundanceStdMesh;
using namespace SundanceStdMesh::Internal;
  namespace Internal
  {
    using namespace Teuchos;

    /** 
     * ExplicitCellSet is a cell set subtype where the cell LIDs
     * are stored explicitly in an STL set. 
     * 
     * @see CellFilter, CellSet, CellSetBase, CellIterator 
     **/
    class ExplicitCellSet : public CellSetBase
    {
    public:

      /** Construct with a mesh, initializing to an empty set */
      ExplicitCellSet(const Mesh& mesh, int cellDim,
                      const CellType& cellType);

      /** Returns an iterator pointing to the first element
       * in the set. */
      virtual CellIterator begin() const ;

      /** Returns a past-the-end iterator */
      virtual CellIterator end() const ;

      /** Returns a modifiable reference to the set of cells */
      Set<int>& cells() {return cells_;}

      /** \name Printable interface */
      //@{
      /** Print to a stream */
      virtual void print(ostream& os) const ;
      //@}

      /* Handleable interface */
      GET_RCP(CellSetBase);

    private:

      /** The set of cell LIDs */
      Set<int> cells_;

      
    };
  }

}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
